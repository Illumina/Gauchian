#!/usr/bin/env python3
#
# Gauchian: GBA variant caller
# Copyright 2021 Illumina, Inc.
# All rights reserved.
#
# Author: Xiao Chen <xchen2@illumina.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

import os
import argparse
import json
import logging
import datetime
from collections import namedtuple
from .caller.gba_caller import GbaCaller
from .depth_calling.utilities import (
    parse_gmm_file,
    parse_region_file,
)
from .depth_calling.snp_count import get_snp_position

MAD_THRESHOLD = 0.11
resource_info = namedtuple(
    "resource_info",
    "genome gmm_parameter region_dic snp_db var_db var_homo_db haplotype_db var_list",
)


def load_parameters():
    """Return parameters."""
    parser = argparse.ArgumentParser(
        description="WGS-based targeted variant caller for GBA."
    )
    parser.add_argument(
        "-m", "--manifest",
        help="Manifest listing absolute paths to input BAM/CRAM files",
        required=True,
    )
    parser.add_argument(
        "-g", "--genome", help="Reference genome, select from 19, 37, or 38", required=True
    )
    parser.add_argument("-o", "--outDir", help="Output directory", required=True)
    parser.add_argument("-p", "--prefix", help="Prefix to output file", required=True)
    parser.add_argument(
        "-t", "--threads",
        help="Number of threads to use. Default is 1",
        type=int,
        default=1,
        required=False,
    )
    parser.add_argument(
        "--reference",
        help="Optional path to reference fasta file for CRAM",
        required=False,
    )

    args = parser.parse_args()
    if args.genome not in ["19", "37", "38"]:
        raise Exception("Genome not recognized. Select from 19, 37, or 38")

    return args


def prepare_data_files(datadir, gene, genome):
    """Check data files and prepare resource data for GBA calling"""
    region_file = os.path.join(datadir, "%s_region_%s.bed" % (gene, genome))
    snp_file = os.path.join(datadir, "%s_SNP_%s.txt" % (gene, genome))
    gmm_file = os.path.join(datadir, "%s_gmm.txt" % gene)
    variant_file = os.path.join(
        datadir, "%s_target_variant_%s.txt" % (gene, genome)
    )
    variant_homology_file = os.path.join(
        datadir, "%s_target_variant_homology_region_%s.txt" % (gene, genome)
    )
    haplotype_file = os.path.join(datadir, "%s_haplotype_%s.txt" % (gene, genome))

    for required_file in [
        region_file,
        snp_file,
        variant_file,
        variant_homology_file,
        haplotype_file,
        gmm_file,
    ]:
        if os.path.exists(required_file) == 0:
            raise Exception("File %s not found." % required_file)

    region_dic = parse_region_file(region_file)

    snp_db = None
    if os.path.exists(snp_file):
        snp_db = get_snp_position(snp_file)

    var_list = []
    var_db = None
    if os.path.exists(variant_file):
        var_db = get_snp_position(variant_file)
        with open(variant_file) as f:
            for line in f:
                if line[0] != "#":
                    var_name = line.split()[-1]
                    var_list.append(var_name)

    var_homo_db = None
    if os.path.exists(variant_homology_file):
        var_homo_db = get_snp_position(variant_homology_file)
        with open(variant_homology_file) as f:
            for line in f:
                if line[0] != "#":
                    var_name = line.split()[-1]
                    var_list.append(var_name)

    HAPLOTYPE_VAR = set()
    haplotype_db = {}
    if os.path.exists(haplotype_file):
        with open(haplotype_file) as f:
            for line in f:
                at = line.split()
                HAPLOTYPE_VAR.add(at[-1])
        for variant in HAPLOTYPE_VAR:
            haplotype_db.setdefault(
                variant,
                get_snp_position(haplotype_file, variant)
                )

    gmm_parameter = None
    if os.path.exists(gmm_file):
        gmm_parameter = parse_gmm_file(gmm_file)

    call_parameters = resource_info(
        genome,
        gmm_parameter,
        region_dic,
        snp_db,
        var_db,
        var_homo_db,
        haplotype_db,
        var_list,
    )
    return call_parameters


def write_to_tsv(final_output, out_tsv):
    """Prepare tsv output"""
    header = [
        "Sample",
        "is_biallelic_GBAP1-like_variant_exon9-11",
        "is_carrier_GBAP1-like_variant_exon9-11",
        "total_CN",
        "deletion_breakpoint_in_GBA_gene",
        "GBAP1-like_variant_exon9-11",
        "other_variants"
        ]
    with open(out_tsv, "w") as tsv_output:
        tsv_output.write("\t".join(header) + "\n")
        for sample_id in final_output:
            final_call = final_output[sample_id]
            if final_call["haplotypes"] is not None:
                p1like_variants = None
                if final_call["GBAP1_like_variant_exon9_11"] != "":
                    p1like_variants = final_call["GBAP1_like_variant_exon9_11"]
                other_variants = None
                if final_call["other_variants"] != "":
                    other_variants = final_call["other_variants"]
                output_per_sample = [
                    sample_id,
                    final_call["haplotypes"]["is_biallelic"],
                    final_call["haplotypes"]["is_carrier"],
                    final_call["haplotypes"]["total_cn"],
                    final_call["haplotypes"]["deletion_bp_in_gene"],
                    p1like_variants,
                    other_variants
                ]
                tsv_output.write("\t".join([str(a) for a in output_per_sample]) + "\n")
            else:
                tsv_output.write("\t".join([str(a) for a in [sample_id]+[None]*6]) + "\n")


def run():
    parameters = load_parameters()
    manifest = parameters.manifest
    outdir = parameters.outDir
    genome = parameters.genome
    prefix = parameters.prefix
    threads = parameters.threads
    reference_fasta = parameters.reference
    logging.basicConfig(level=logging.DEBUG)

    datadir = os.path.join(os.path.dirname(__file__), "data")
    if os.path.exists(outdir) == 0:
        os.makedirs(outdir)
    out_json = os.path.join(outdir, prefix + ".json")
    out_tsv = os.path.join(outdir, prefix + ".tsv")
    call_parameters = prepare_data_files(datadir, "GBA", genome)

    final_output = {}
    with open(manifest) as read_manifest:
        for line in read_manifest:
            bam_name = line.strip()
            sample_id = os.path.splitext(os.path.basename(bam_name))[0]
            if os.path.exists(bam_name) == 0:
                logging.warning(
                    "Input alignmet file for sample %s does not exist.", sample_id
                )
            else:
                logging.info(
                    "Processing sample %s at %s", sample_id, datetime.datetime.now()
                )
                gba_caller = GbaCaller()
                gba_caller.set_par(
                    bam_name,
                    call_parameters,
                    threads,
                    reference_fasta,
                )
                gba_call = gba_caller.call()._asdict()
                # Use normalized coverage MAD across stable regions
                # as a sample QC measure.
                if gba_call["Coverage_MAD"] > MAD_THRESHOLD:
                    logging.warning(
                        "Sample %s has uneven coverage. CN calls may be unreliable.",
                        sample_id,
                    )
                final_output.setdefault(sample_id, gba_call)

    # Write to json
    logging.info("Writing to json at %s", datetime.datetime.now())
    with open(out_json, "w") as json_output:
        json.dump(final_output, json_output)

    # Write to tsv
    logging.info('Writing to tsv at %s', datetime.datetime.now())
    write_to_tsv(final_output, out_tsv)


if __name__ == "__main__":
    run()

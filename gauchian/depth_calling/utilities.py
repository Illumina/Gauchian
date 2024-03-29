#!/usr/bin/env python3
#
# Gauchian: GBA variant caller
# Copyright 2021 Illumina, Inc.
# All rights reserved.
#
# Author: Xiao Chen <xchen2@illumina.com>
#
# This program is licensed under the terms of the Polyform strict license
#
# ***As far as the law allows, the software comes as is, without
# any warranty or condition, and the licensor will not be liable
# to you for any damages arising out of these terms or the use
# or nature of the software, under any kind of legal claim.***
#
# You should have received a copy of the PolyForm Strict License 1.0.0
# along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
#
#


import pysam


def parse_region_file(region_file):
    """Return the set of regions for counting from a bed file."""
    region_dic = {}
    with open(region_file) as read_region:
        for line in read_region:
            nchr, region_start, region_end, region_name, region_type, region_gc = (
                line.strip().split()
            )
            region_start = int(region_start)
            region_end = int(region_end)
            region = (nchr, region_start, region_end, region_name)
            region_dic.setdefault(region_type, []).append((region, region_gc))
    return region_dic


def parse_gmm_file(gmm_file):
    """Return the gmm parameters stored in input file."""
    dpar_tmp = {}
    with open(gmm_file) as read_gmm:
        for line in read_gmm:
            split_line = line.strip().split()
            dpar_tmp.setdefault(split_line[0], {})
            list_value = [a.split(":")[-1] for a in split_line[2:]]
            dpar_tmp[split_line[0]].setdefault(split_line[1], list_value)
    return dpar_tmp


def open_alignment_file(alignment_file, reference_fasta=None):
    if alignment_file.endswith("cram"):
        return pysam.AlignmentFile(
            alignment_file, "rc", reference_filename=reference_fasta
        )
    return pysam.AlignmentFile(alignment_file, "rb")

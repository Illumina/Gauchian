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
import sys
import numpy as np
from collections import namedtuple
dir_name = os.path.dirname(os.path.dirname(__file__))
sys.path.append(dir_name)

from depth_calling.utilities import open_alignment_file
from depth_calling.gene_caller import GeneCaller
from .phasing_gba import PhasingGba


class GbaCaller(GeneCaller):
    sample_call = namedtuple(
        "sample_call",
        "Coverage_MAD Median_depth deletion_CN deletion_CN_raw \
        GBAP1_like_variant_exon9_11 other_variants variant_raw_count \
        snp_call snp_raw haplotypes",
    )

    def __init__(self):
        super().__init__()

    def improve_total_cn(self):
        """
        Improve GBA total cn when cn is high
        """
        # there could be noise in high CN calls. pick the most likely total CN
        # so that sites at the beginning of GBA (should always be CN(GBA)=2)
        # are called close to CN2
        self.total_cn = None
        if self.gmm_call.cn is not None and self.gmm_call.cn <= 4:
            self.total_cn = self.gmm_call.cn + 2
        elif self.gmm_call.depth_value > 4.25:
            total_cn_low = int(np.floor(self.gmm_call.depth_value + 2))
            total_cn_high = int(np.ceil(self.gmm_call.depth_value + 2))
            gene_sites_low = [a*total_cn_low for a in self.gene1_fraction[-20:]]
            gene_sites_high = [a*total_cn_high for a in self.gene1_fraction[-20:]]
            # substract the expected CN, which is 2
            gene_sites_low_diff = abs(np.median([a - 2 for a in gene_sites_low]))
            gene_sites_high_diff = abs(np.median([a - 2 for a in gene_sites_high]))
            if gene_sites_low_diff < 0.1 and gene_sites_high_diff > 0.2:
                self.total_cn = total_cn_low
            elif gene_sites_high_diff < 0.1 and gene_sites_low_diff > 0.2:
                self.total_cn = total_cn_high
        if self.total_cn is None and self.gmm_call.cn is not None:
            self.total_cn = self.gmm_call.cn + 2

    def call_gba_recombinants(self, debug=False):
        """
        Call variants in homology region by analyzing haplotypes
        """
        phase_gba = PhasingGba()
        bam_handle = open_alignment_file(self.bamfile, self.reference)
        phase_gba.set_par(
            bam_handle,
            self.call_parameters.haplotype_db["Exon9-11"],
            self.total_cn,
            self.reference
        )
        bam_handle.close()
        return phase_gba.call_variants_from_haplotype(debug=debug)

    def call(self):
        """Perform variant calling"""
        # 1. read counting, normalization
        self.normalize_depth()
        if self.normalized_depth.normalized["GBAdel"] is None:
            return self.sample_call(
                self.normalized_depth.mad,
                self.normalized_depth.mediandepth,
                None, None, None, None, None, None, None, None)
        # 2. call copy number
        self.call_cn("GBAdel")
        # 3. Get read counts and SNP ratios
        self.count_reads_diff_sites_and_variants()
        self.improve_total_cn()
        if self.total_cn is None:
            return self.sample_call(
                self.normalized_depth.mad,
                self.normalized_depth.mediandepth,
                self.gmm_call.cn, self.gmm_call.depth_value,
                None, None, self.raw_count, None, None, None
            )
        # 4. Call variants (simple)
        self.call_cn_diff_sites_and_variants()
        # site 48 onwards gets into GBA gene.
        # Mark here so we can easily visualize breakpoint location
        self.raw_gene1_cn_out = self.raw_gene1_cn
        self.raw_gene1_cn_out.insert(48, 'gene')
        # 5. Call variants (use the phasing approach in homology region Exon9-11)
        recombinant_variants = []
        variants_called_from_haplotype = self.call_gba_recombinants()
        if (
            variants_called_from_haplotype.is_biallelic or
            variants_called_from_haplotype.is_carrier
        ):
            recombinant_variants.append(
                ",".join(variants_called_from_haplotype.variants_allele1)
                )
            recombinant_variants.append(
                ",".join(variants_called_from_haplotype.variants_allele2)
                )
        # 6. Return final call
        self.cn_call_snp_out = self.cn_call_snp
        self.cn_call_snp_out.insert(48, 'gene')
        sample_gba_call = self.sample_call(
            self.normalized_depth.mad,
            self.normalized_depth.mediandepth,
            self.gmm_call.cn,
            self.gmm_call.depth_value,
            "/".join(recombinant_variants),
            ",".join(self.total_callset),
            self.raw_count,
            ",".join([str(a) for a in self.cn_call_snp_out]),
            ",".join([str(a) for a in self.raw_gene1_cn_out]),
            variants_called_from_haplotype._asdict()
        )
        return sample_gba_call

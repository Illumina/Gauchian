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

from collections import OrderedDict
from .snp_count import (
    get_supporting_reads,
    get_fraction,
    get_supporting_reads_single_region,
)
from .gmm import Gmm
from .bin_count import get_normed_depth
from .copy_number_call import (
    call_cn_snp,
    call_cn_var,
    call_cn_var_homo,
    get_called_variants,
)
from .utilities import open_alignment_file


class GeneCaller:
    """
    A class for calling variant in paralogs
    """
    def __init__(self):
        self.bamfile = None
        self.call_parameters = None
        self.reference = None
        self.threads = None
        self.normalized_depth = None
        self.gmm_call = None
        self.total_cn = None
        self.raw_count = OrderedDict()
        self.gene1_read_count = []
        self.gene2_read_count = []
        self.gene1_fraction = []
        self.raw_gene1_cn = []
        self.var_homo_alt = []
        self.var_homo_ref = []
        self.var_alt = []
        self.var_ref = []
        self.var_list = []
        self.cn_call_snp = []
        self.total_callset = []

    def set_par(self, bamfile, call_parameters, threads=None, reference=None):
        """Establish parameter values"""
        self.bamfile = bamfile
        self.call_parameters = call_parameters
        self.threads = threads
        self.reference = reference
        self.var_list = self.call_parameters.var_list

    def normalize_depth(self, gc_correct=True):
        """Count reads and perform depth normalization"""
        self.normalized_depth = get_normed_depth(
            self.bamfile,
            self.call_parameters.region_dic,
            self.threads,
            self.reference,
            gc_correct=gc_correct
        )

    def call_cn(self, region_name):
        """Call CN for regions of interest"""
        gmm = Gmm()
        gmm_parameter = self.call_parameters.gmm_parameter
        gmm.set_gmm_par(gmm_parameter, region_name)
        self.gmm_call = gmm.gmm_call(self.normalized_depth.normalized[region_name])

    def count_reads_diff_sites_and_variants(self):
        """Get read counts at differentiating sites and variant sites"""
        bam_handle = open_alignment_file(self.bamfile, self.reference)
        snp_db = self.call_parameters.snp_db
        self.gene1_read_count, self.gene2_read_count = get_supporting_reads(
            bam_handle, snp_db.dsnp1, snp_db.dsnp2, snp_db.nchr, snp_db.dindex
        )
        self.gene1_fraction = get_fraction(self.gene1_read_count, self.gene2_read_count)
        var_db = self.call_parameters.var_db
        self.var_alt, self.var_ref, var_alt_forward, var_alt_reverse = get_supporting_reads_single_region(
            bam_handle, var_db.dsnp1, var_db.nchr, var_db.dindex
        )
        var_homo_db = self.call_parameters.var_homo_db
        self.var_homo_alt, self.var_homo_ref = get_supporting_reads(
            bam_handle,
            var_homo_db.dsnp1,
            var_homo_db.dsnp2,
            var_homo_db.nchr,
            var_homo_db.dindex,
        )
        bam_handle.close()
        self.raw_count = OrderedDict()
        for i in range(len(self.var_list)):
            if i < len(self.var_alt):
                self.raw_count.setdefault(self.var_list[i], "%i,%i" % (self.var_alt[i], self.var_ref[i]))
            else:
                self.raw_count.setdefault(
                    self.var_list[i],
                    "%i,%i"
                    % (self.var_homo_alt[i - len(self.var_alt)], self.var_homo_ref[i - len(self.var_alt)]),
                )

    def call_cn_diff_sites_and_variants(self):
        """Call CN at differentiating sites and call variants"""
        # differentiating sites
        self.cn_call_snp = call_cn_snp(self.total_cn, self.gene1_read_count, self.gene2_read_count)
        self.raw_gene1_cn = [round(self.total_cn * a, 3) for a in self.gene1_fraction]
        # homology region
        cn_call_var_homo = call_cn_var_homo(
            self.total_cn, self.var_homo_alt, self.var_homo_ref, max_cn=2
            )
        # non-homology region
        cn_call_var = call_cn_var(2, self.var_alt, self.var_ref)
        # update variant names
        self.total_callset = get_called_variants(self.var_list, cn_call_var)
        called_var_homo = get_called_variants(self.var_list, cn_call_var_homo, len(cn_call_var))
        self.total_callset += called_var_homo


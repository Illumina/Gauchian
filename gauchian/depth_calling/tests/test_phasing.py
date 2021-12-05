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

import sys
import os
import pytest

from ..phasing import Phasing

test_data_dir = os.path.join(os.path.dirname(__file__), "test_data")


class TestPhasing(object):
    def test_get_variants(self):
        phase = Phasing()
        phase.variant_sites = [1, 3]
        phase.variant_names = ['var1', 'var2']
        assert phase.get_variants('12111') == ['var1']
        assert phase.get_variants('11121') == ['var2']
        assert phase.get_variants('12121') == ['var1', 'var2']

    def test_assess_assembled_haps(self):
        phase = Phasing()
        phase.variant_sites = [1]
        haps = ['12x', '111', '121']
        fully_assembled, var_sites_haps = phase.assess_assembled_haps(haps)
        assert fully_assembled is False
        assert var_sites_haps == ['1', '2']

    def test_carrier_from_assembled_haps(self):
        phase = Phasing()
        phase.total_cn = 3
        phase.variant_sites = [1]
        haps = ['122', '111', '121']
        fully_assembled, var_sites_haps = phase.assess_assembled_haps(haps)
        assert fully_assembled is True
        assert phase.carrier_from_assembled_haps(var_sites_haps) is True

        haps = ['122', '111', '112']
        fully_assembled, var_sites_haps = phase.assess_assembled_haps(haps)
        assert fully_assembled is True
        assert phase.carrier_from_assembled_haps(var_sites_haps) is False

    def test_filter_haps_get_var_sites(self):
        phase = Phasing()
        phase.variant_sites = [1]
        phase.num_sites = 3
        phase.depth_calls = [2, 1, 2]
        haps = ['222', '111', '121']
        var_index_in_haps = phase.filter_haps_get_var_sites(haps)
        assert var_index_in_haps == {"121": [1]}

        phase.num_sites = 4
        phase.depth_calls = [2, 1, 1, 1]
        haps = ['2222', '1111', '1222']
        var_index_in_haps = phase.filter_haps_get_var_sites(haps)
        assert var_index_in_haps == {"1222": [1]}

        phase.num_sites = 4
        phase.depth_calls = [3, 2, 2, 2]
        haps = ['2222', '1111', '1222']
        var_index_in_haps = phase.filter_haps_get_var_sites(haps)
        assert var_index_in_haps == {}

    def test_get_base_cns(self):
        phase = Phasing()
        phase.depth_calls = [2, 1, 2, 2, 1]
        phase.get_base_cn(0, len(phase.depth_calls)-1)
        assert phase.depth_base_cn == 2

        phase = Phasing()
        phase.depth_calls = [3, 1, 3, 2, 1]
        phase.get_base_cn(0, 2)
        assert phase.depth_base_cn is None

    def test_get_candidate_gene_hap(self):
        phase = Phasing()
        haps = ['11111', '11211', '22222']
        phase.get_candidate_gene_hap(haps)
        assert phase.candidate_gene_hap == ['11111', '11211']

        phase = Phasing()
        haps = ['111111', '112221', '222221']
        phase.get_candidate_gene_hap(haps)
        assert phase.candidate_gene_hap == ['111111']

    def test_compare_haplotype_against_depth(self):
        phase = Phasing()
        phase.depth_calls = [2, 1, 2, 2]
        assert phase.compare_haplotype_against_depth('1211') == 1
        assert phase.compare_haplotype_against_depth('1221') == 0.75

    def test_assess_compound_het(self):
        phase = Phasing()
        phase.variant_sites = list(range(5))
        haps = ['11111', '12111', '22222']
        assert phase.assess_compound_het(haps) is False

        haps = ['21111', '12111', '22222']
        assert phase.assess_compound_het(haps) is True

    def test_assess_homozygous(self):
        phase = Phasing()
        phase.variant_sites = list(range(5))
        phase.variant_names = ['']*5
        phase.depth_calls = [0, 1, 1, 1, 2]
        assert phase.assess_homozygous() is True

        phase.depth_calls = [1, 1, 1, 1, 2]
        assert phase.assess_homozygous() is False

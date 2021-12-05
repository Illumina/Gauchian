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


import pytest
from ..phasing_gba import PhasingGba


class TestPhasingGBA(object):
    def test_get_variants_gba(self):
        phase = PhasingGba()
        hap = "1112111111"
        variants = phase.get_variants_gba(hap)
        assert variants == ["L483P"]

        hap = "1222111111"
        variants = phase.get_variants_gba(hap)
        assert variants == ["RecNciI"]

        hap = "2222222111"
        variants = phase.get_variants_gba(hap)
        assert variants == ["RecTL"]

    def test_check_deletion_bp_in_gene(self):
        phase = PhasingGba()
        phase.total_cn = 4
        haps = {
            "1112111111": 0,
            "1111111111": 0,
        }
        phase.check_deletion_bp_in_gene(haps)
        assert phase.deletion_bp_in_gene is None

        phase.total_cn = 3
        haps = {
            "2222111111": 0,
            "1111111111": 0,
        }
        phase.check_deletion_bp_in_gene(haps)
        assert phase.deletion_bp_in_gene is True

    def test_assemble_haplotypes(self):
        phase = PhasingGba()
        # one hybrid site
        phase.haplotype_per_read = {
            'read0': '1x2xxxxxxx',
            'read1': 'x121xxxxxx',
            'read2': 'xx211xxxxx',
            'read3': '111xxxxxxx',
            'read4': 'xx11xxxxxx',
            'read5': 'xx11x1xxxx',
            'read6': 'xxx111xxxx',
            'read7': 'xxxx1111xx',
            'read8': 'xxxxxx111x',
            'read9': 'xxxxx1x11x',
            'read10': 'xxxxxxx111',
        }
        assembled_haplotypes = phase.assemble_haplotypes()
        assert sorted(assembled_haplotypes.full_haplotypes) == ['1111111111', '1121111111']

        # two hybrid sites, phased
        phase.haplotype_per_read = {
            'read0':  '1x2xxxxxxx',
            'read1':  'x121xxxxxx',
            'read2':  'xx211xxxxx',
            'read3':  '111xxxxxxx',
            'read4':  'xx11xxxxxx',
            'read5':  'xx11x1xxxx',
            'read6':  'xxx112xxxx',
            'read7':  'xxxx1211xx',
            'read8':  'xxxxxx111x',
            'read9':  'xxxxx2x11x',
            'read10': 'xxxxxxx111',
            'read11': 'xxxxxx11x1',
            'read12': 'xxxx121xxx',
            'read13': 'xxx1111xxx',
            'read14': 'xx111xxxxx',
        }
        assembled_haplotypes = phase.assemble_haplotypes()
        assert sorted(assembled_haplotypes.hap_5p) == ['1111111xxx', '1121121xxx']
        assert sorted(assembled_haplotypes.hap_3p) == ['xx11111111', 'xxx1121111']
        assert sorted(assembled_haplotypes.full_haplotypes) == ['1111111111', '1121121111']

        # two hybrid sites, generate all possible combinations
        phase.haplotype_per_read = {
            'read0':  '1x2xxxxxxx',
            'read1':  'x121xxxxxx',
            'read2':  'xx211xxxxx',
            'read3':  '111xxxxxxx',
            'read4':  'xx11xxxxxx',
            'read5':  'xx111xxxxx',
            'read6':  'xxx112xxxx',
            'read7':  'xxxx1211xx',
            'read8':  'xxxxxx111x',
            'read9':  'xxxxx2x11x',
            'read10': 'xxxxxxx111',
            'read11': 'xxxxxx11x1',
            'read12': 'xxxx121xxx',
            'read13': 'xxx1111xxx',
            'read14': 'xx111xxxxx',
        }
        assembled_haplotypes = phase.assemble_haplotypes()
        assert sorted(assembled_haplotypes.hap_5p) == ['11111xxxxx', '11211xxxxx']
        assert sorted(assembled_haplotypes.hap_3p) == ['xxx1111111', 'xxx1121111']
        assert sorted(assembled_haplotypes.full_haplotypes) == ['1111111111', '1111121111', '1121111111', '1121121111']

        # two hybrid sites, farther away, generate all possible combinations
        phase.haplotype_per_read = {
            'read0': '1x2xxxxxxx',
            'read1': 'x121xxxxxx',
            'read2': 'xx211xxxxx',
            'read3': '111xxxxxxx',
            'read4': 'xx11xxxxxx',
            'read5': 'xx11x1xxxx',
            'read6': 'xxx111xxxx',
            'read7': 'xxxx1111xx',
            'read8': 'xxxxxx111x',
            'read9': 'xxxxx1x11x',
            'read10': 'xxxxxxx112',
            'read11': 'xxxxxx11x1',
            'read12': 'xxxx111xxx',
            'read13': 'xxx1111xxx',
            'read14': 'xxxx111xxx',
            'read15': 'xx111xxxxx',
            'read16': 'xxxx111xxx',
            'read17': 'xxx1111xxx',
        }
        assembled_haplotypes = phase.assemble_haplotypes()
        assert sorted(assembled_haplotypes.hap_5p) == ['1111111xxx', '1121111xxx']
        assert sorted(assembled_haplotypes.hap_3p) == ['xxx1111111', 'xxx1111112']
        assert sorted(assembled_haplotypes.full_haplotypes) == ['1111111111', '1111111112', '1121111111', '1121111112']

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

from ..haplotype import (
    haplotype_groups,
    get_base1_base2,
    get_haplotype_per_read,
    extract_hap,
    extend_hap_5p,
    extend_hap_3p,
    get_block_read_count,
    get_possible_hap_cns,
    call_hap_cn,
    get_read_count_threshold,
    join_subblocks,
    filter_hap,
    compare_two_haps,
    get_seq_consensus,
)
from ..snp_count import get_snp_position

test_data_dir = os.path.join(os.path.dirname(__file__), "test_data")


class TestHaplotype(object):
    def test_get_bases(self):
        snp_file = os.path.join(test_data_dir, "SMN_SNP_37.txt")
        base_db = get_snp_position(snp_file)
        target_positions = [11, 12, 13]
        base1, base2 = get_base1_base2(base_db, target_positions)
        assert base1 == ["G", "C", "A"]
        assert base2 == ["A", "T", "G"]

    def test_get_haplotype_per_read(self):
        snp_file = os.path.join(test_data_dir, "SMN_SNP_37.txt")
        base_db = get_snp_position(snp_file)
        target_positions = [11, 12]
        base1, base2 = get_base1_base2(base_db, target_positions)
        assert base1 == ["G", "C"]
        assert base2 == ["A", "T"]
        dread = {
            "read1": {0: "G", 1: "C"},
            "read2": {0: "A", 1: "T"},
            "read3": {0: "G", 1: "T"},
            "read4": {0: "G", 1: "T"},
        }
        dhap = get_haplotype_per_read(dread, base1, base2)
        hap_count = extract_hap(dhap, [0, 1])
        assert hap_count["11"] == [1]
        assert hap_count["22"] == [1]
        assert hap_count["12"] == [1, 1]

        hap_count = extract_hap(dhap, [1])
        assert hap_count["1"] == [1]
        assert hap_count["2"] == [1, 1, 1]

    def test_get_block_read_count(self):
        haps = {
            "111": 10,
            "222": 20
        }
        block_count = get_block_read_count(haps, 0)
        assert block_count == [10, 20]

        haps = {
            "11": 10,
            "12": 20
        }
        block_count = get_block_read_count(haps, 0)
        assert block_count == [10, 20]

        haps = {
            "12": 10,
            "22": 20
        }
        block_count = get_block_read_count(haps, 0)
        assert block_count is None

        haps = {
            "111": 10,
            "121": 15,
            "222": 20
        }
        block_count = get_block_read_count(haps, 0)
        assert block_count == [10, 15, 20]

        haps = {
            "111": 10,
            "122": 15,
            "121": 20
        }
        block_count = get_block_read_count(haps, 0)
        assert block_count == [10, 15, 20]

        haps = {
            "111": 10,
            "122": 15,
            "121": 20,
            "222": 25
        }
        block_count = get_block_read_count(haps, 0)
        assert block_count == [10, 15, 20, 25]

        haps = {
            "211": 10,
            "222": 15,
        }
        block_count = get_block_read_count(haps, 0)
        assert block_count == [0]

    def test_get_possible_hap_cns(self):
        possible_hap_cns = get_possible_hap_cns(4, 3)
        assert sorted(list(possible_hap_cns)) == [
            (1, 1, 2), (1, 2, 1), (2, 1, 1)
            ]
        possible_hap_cns = get_possible_hap_cns(5, 4)
        assert sorted(list(possible_hap_cns)) == [
            (1, 1, 1, 2), (1, 1, 2, 1), (1, 2, 1, 1), (2, 1, 1, 1)
            ]

    def test_call_hap_cn(self):
        cn_call = call_hap_cn([10, 10], 2)
        assert cn_call > 0.9

        cn_call = call_hap_cn([10, 10, 10], 2)
        assert cn_call == 0

        cn_call = call_hap_cn([10, 10, 10], 3)
        assert cn_call > 0.9

        cn_call = call_hap_cn([10, 2, 10], 3)
        assert cn_call > 0.9

        cn_call = call_hap_cn([10, 10, 20], 4)
        assert cn_call > 0.9

        cn_call = call_hap_cn([20, 10, 10], 4)
        assert cn_call < 0.4

        cn_call = call_hap_cn([10, 20, 10], 4)
        assert cn_call > 0.9

        cn_call = call_hap_cn([10, 20, 20], 5)
        assert cn_call > 0.9

        cn_call = call_hap_cn([20, 10, 20], 5)
        assert cn_call < 0.4

        cn_call = call_hap_cn([10, 10, 10, 10], 4)
        assert cn_call > 0.9

        cn_call = call_hap_cn([10, 10, 20, 10], 5)
        assert cn_call > 0.9

        cn_call = call_hap_cn([10, 20, 10, 10], 5)
        assert cn_call > 0.9

        cn_call = call_hap_cn([20, 10, 10, 10], 5)
        assert cn_call < 0.4

        cn_call = call_hap_cn([10, 10, 30, 10], 6)
        assert cn_call > 0.9

        cn_call = call_hap_cn([10, 20, 20, 10], 6)
        assert cn_call > 0.9

        cn_call = call_hap_cn([10, 20, 30, 10], 7)
        assert cn_call > 0.9

        cn_call = call_hap_cn([10, 20, 40, 10], 8)
        assert cn_call > 0.9

        cn_call = call_hap_cn([20, 10, 40, 10], 8)
        assert cn_call < 0.4

        cn_call = call_hap_cn([10, 10, 20, 10], 10)
        assert cn_call == 0

    def test_get_seq_consensus(self):
        con = get_seq_consensus("121xx", "x2111")
        assert con == "12111"

        con = get_seq_consensus("12xxx", "xx111")
        assert con is False

        con = get_seq_consensus("12xx2", "xx111")
        assert con is False

    def test_join_subblocks(self):
        dread = {
            'read1': "12xxxx",
            'read2': "x211xx",
            'read3': "xxx121",
            'read4': "x2xx21"
        }
        lhap1 = ['1211xx']
        lhap2 = ['xx1121']
        merged, dcount = join_subblocks(dread, lhap1, lhap2, 2, 4)
        assert merged == ['121121']
        assert dcount == {'121121': ['read4']}

    def test_extend_hap_5p(self):
        matching_haplotype_groups = haplotype_groups(
            {
                'x11x': ['21xx', '211x', '21xx']
            }, None
        )
        new_hap_list = extend_hap_5p(matching_haplotype_groups)
        assert new_hap_list == ['211x']

        matching_haplotype_groups = haplotype_groups(
            {
                'x11x': ['21xx', '111x', '21xx']
            }, None
        )
        new_hap_list = extend_hap_5p(matching_haplotype_groups)
        assert sorted(new_hap_list) == ['111x', '211x']
        new_hap_list = extend_hap_5p(matching_haplotype_groups, min_read=2)
        assert sorted(new_hap_list) == ['211x']

        matching_haplotype_groups = haplotype_groups(
            {
                'x11x': ['21xx', '111x', '21xx'],
                'x21x': ['x21x', '121x']
            }, None
        )
        # 'x21x' shouldn't be dropped due to #reads < min_read
        new_hap_list = extend_hap_5p(matching_haplotype_groups, min_read=2)
        assert sorted(new_hap_list) == ['211x', 'x21x']

    def test_extend_hap_3p(self):
        matching_haplotype_groups = haplotype_groups(
            {
                'x11x': ['x112', '211x', '21x2']
            }, None
        )
        new_hap_list = extend_hap_3p(matching_haplotype_groups)
        assert new_hap_list == ['x112']

        matching_haplotype_groups = haplotype_groups(
            {
                'x11x': ['21xx', '111x', '21xx']
            }, None
        )
        new_hap_list = extend_hap_3p(matching_haplotype_groups)
        assert new_hap_list == ['x11x']

        matching_haplotype_groups = haplotype_groups(
            {
                'x11x': ['21xx', '1111', '21xx']
            }, None
        )
        new_hap_list = extend_hap_3p(matching_haplotype_groups)
        assert new_hap_list == ['x111']
        # 'x11x' shouldn't be dropped due to #reads < min_read
        new_hap_list = extend_hap_3p(matching_haplotype_groups, min_read=2)
        assert new_hap_list == ['x11x']

        matching_haplotype_groups = haplotype_groups(
            {
                'x11x': ['21x2', '1111', '21x2', '1112', 'x112']
            }, None
        )
        # min_read in this case is 1, so x111
        new_hap_list = extend_hap_3p(matching_haplotype_groups)
        assert sorted(new_hap_list) == ['x111', 'x112']
        matching_haplotype_groups = haplotype_groups(
            {
                'x11x': ['21x2', '1111', '21x2', '1112', 'x112', 'x112']
            }, None
        )
        # min_read in this case is 2, so no x111
        new_hap_list = extend_hap_3p(matching_haplotype_groups)
        assert new_hap_list == ['x112']

    def test_compare_two_haps(self):
        hap1 = 'xx1121xx'
        hap2 = 'x1112x2x'
        match, mismatch, extend = compare_two_haps(hap1, hap2)
        assert match == 3
        assert mismatch == 0
        assert extend == 2

        hap1 = 'xx1221xx'
        hap2 = 'x1112x2x'
        match, mismatch, extend = compare_two_haps(hap1, hap2)
        assert match == 2
        assert mismatch == 1
        assert extend == 2

    def test_get_read_count_threshold(self):
        bases = ['1', '1', '1', '2', '2']
        assert get_read_count_threshold(bases) == 1

        bases = ['1', '1', '1', '1', '1', '2']
        assert get_read_count_threshold(bases) == 2

    def test_filter_hap(self):
        haplotype_read_count = {
            '11': [1, 1, 1],
            '22': [1, 1]
        }
        haplotype_read_count_filtered = filter_hap(haplotype_read_count)
        assert haplotype_read_count_filtered == {'11': 3, '22': 2}
        haplotype_read_count_filtered = filter_hap(haplotype_read_count, min_read=3)
        assert haplotype_read_count_filtered == {'11': 3}

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


import sys
import os
import pytest

from ..copy_number_call import (
    call_reg1_cn,
    process_raw_call_gc,
    process_raw_call_denovo,
    get_called_variants
)


class TestCallCN(object):
    def test_call_reg1_cn(self):
        call = call_reg1_cn(4, 30, 30)
        assert call == [2]
        call = call_reg1_cn(4, 10, 30)
        assert call == [1]
        call = call_reg1_cn(4, 30, 10)
        assert call == [3]
        call = call_reg1_cn(3, 30, 30)
        assert call == [2, 0.689, 1, 0.311]

    def test_process_raw_call_gc(self):
        lcn = [[1], [1, 0.8, 2, 0.2]]
        filtered_call = process_raw_call_gc(lcn, 0.7)
        assert filtered_call == [1, 1]
        filtered_call = process_raw_call_gc(lcn, 0.9)
        assert filtered_call == [1, None]
        filtered_call = process_raw_call_gc(lcn, 0.9, keep_none=False)
        assert filtered_call == [1]

    def test_process_raw_call_denovo(self):
        lcn = [[1], [0, 0.8, 1, 0.2], [2, 0.6, 1, 0.4]]
        filtered_call = process_raw_call_denovo(lcn, 0.9, 0.7)
        assert filtered_call == [1, None, 1]
        filtered_call = process_raw_call_denovo(lcn, 0.9, 0.55)
        assert filtered_call == [1, None, 2]

        lcn = [[1], [1, 0.65, 0, 0.35], [2, 0.6, 1, 0.4]]
        filtered_call = process_raw_call_denovo(lcn, 0.9, 0.62, [1, 1, 1])
        assert filtered_call == [1, 1, 1]
        filtered_call = process_raw_call_denovo(lcn, 0.9, 0.7, [1, 1, 1])
        assert filtered_call == [1, 0, 1]

    def test_get_called_variants(self):
        var_list = ["var1", "var2", "var3", "var4"]
        cn_called = [0, None, 1, 2]
        var_called = get_called_variants(var_list, cn_called)
        assert var_called == ["var3", "var4", "var4"]

        var_list = ["var1", "var2", "var3", "var4", "var5"]
        cn_called = [0, None, 1, 2]
        var_called = get_called_variants(var_list, cn_called, 1)
        assert var_called == ["var4", "var5", "var5"]

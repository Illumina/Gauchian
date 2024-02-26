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

from ..utilities import parse_region_file

test_data_dir = os.path.join(os.path.dirname(__file__), "test_data")


class TestUtilities(object):
    def test_parse_reigon_file(self):
        region_file = os.path.join(test_data_dir, "SMN_region_19_short.bed")
        region_dic = parse_region_file(region_file)
        assert len(region_dic["norm"]) == 500
        assert len(region_dic["exon16"]) == 2
        assert len(region_dic["exon78"]) == 2

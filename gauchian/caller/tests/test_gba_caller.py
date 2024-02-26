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

import pytest
from collections import namedtuple
from ..gba_caller import GbaCaller


class TestGbaCaller(object):
    cn_call = namedtuple("cn_call", "cn depth_value")
    def test_improve_total_cn(self):
        caller = GbaCaller()
        caller.total_cn = None
        caller.gmm_call = self.cn_call(2, 1.0)
        caller.improve_total_cn()
        assert caller.total_cn == 4

        caller.gmm_call = self.cn_call(None, 4.5)
        caller.gene1_fraction = [0.33]*35
        caller.improve_total_cn()
        assert caller.total_cn == 6

        caller.gene1_fraction = [0.28]*35
        caller.improve_total_cn()
        assert caller.total_cn == 7

        caller.gmm_call = self.cn_call(5, 4.5)
        caller.gene1_fraction = [0.33]*35
        caller.improve_total_cn()
        assert caller.total_cn == 6

        caller.gmm_call = self.cn_call(5, 4.5)
        caller.gene1_fraction = [0.3]*35
        caller.improve_total_cn()
        assert caller.total_cn == 7

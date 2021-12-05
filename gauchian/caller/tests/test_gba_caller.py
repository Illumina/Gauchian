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

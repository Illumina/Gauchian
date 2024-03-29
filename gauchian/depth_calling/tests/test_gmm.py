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
import numpy as np


from ..gmm import Gmm
from ..utilities import parse_gmm_file

test_data_dir = os.path.join(os.path.dirname(__file__), "test_data")


class TestGMM(object):
    def test_gmm_parameter(self):
        gmm_file = os.path.join(test_data_dir, "SMN_gmm.txt")
        dpar_tmp = parse_gmm_file(gmm_file)

        test_gmm = Gmm()
        test_gmm.set_gmm_par(dpar_tmp, "exon1-6")
        assert test_gmm.value_shift == 0.994
        assert len(test_gmm.mu_state) == 11
        assert len(test_gmm.prior_state) == 11
        assert test_gmm.mu_state[0:4] == [0, 0.5, 1, 1.495]
        assert round(test_gmm.mu_state[4], 4) == 1.99
        assert round(test_gmm.mu_state[5], 4) == 2.485
        # priors sum up to 1
        assert round(sum(test_gmm.prior_state), 2) == 1
        assert test_gmm.prior_state[2] == 0.026
        assert round(test_gmm.prior_state[0], 4) == 0.001
        assert round(test_gmm.prior_state[8], 4) == 0.0003
        assert test_gmm.sigma_state[0:4] == [
            0.032,
            0.051 / np.sqrt(2),
            0.051,
            0.051 * np.sqrt(1.5),
        ]

    def test_gmmcall(self):
        gmm_file = os.path.join(test_data_dir, "SMN_gmm.txt")
        dpar_tmp = parse_gmm_file(gmm_file)
        test_gmm = Gmm()
        test_gmm.set_gmm_par(dpar_tmp, "exon1-6")
        cncall = test_gmm.gmm_call(2.1)
        assert cncall[0] == 2
        cncall = test_gmm.gmm_call(6.48)
        assert cncall[0] is None
        # when depth value > default_GMM_nstate (11), a no-call will be made
        cncall = test_gmm.gmm_call(12.05)
        assert cncall[0] is None

        test_gmm = Gmm()
        test_gmm.set_gmm_par(dpar_tmp, "exon7-8")
        cncall = test_gmm.gmm_call(0.18)
        assert cncall[0] == 0
        cncall = test_gmm.gmm_call(0.95)
        assert cncall[0] == 1
        cncall = test_gmm.gmm_call(1.635)
        assert cncall[0] == 2
        cncall = test_gmm.gmm_call(2.26)
        assert cncall[0] == 2
        cncall = test_gmm.gmm_call(2.364)
        assert cncall[0] is None
        cncall = test_gmm.gmm_call(2.53)
        assert cncall[0] == 3
        cncall = test_gmm.gmm_call(3.35)
        assert cncall[0] == 3
        cncall = test_gmm.gmm_call(3.39)
        assert cncall[0] is None
        cncall = test_gmm.gmm_call(3.515)
        assert cncall[0] == 4
        cncall = test_gmm.gmm_call(4.391)
        assert cncall[0] == 4
        cncall = test_gmm.gmm_call(4.4)
        assert cncall[0] is None
        cncall = test_gmm.gmm_call(4.6)
        assert cncall[0] == 5
        cncall = test_gmm.gmm_call(5.38)
        assert cncall[0] is None

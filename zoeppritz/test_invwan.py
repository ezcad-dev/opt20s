# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.

import unittest
import numpy as np
from zoeppritz.utils import elapar_hs2delta
from zoeppritz.modaki import inc2ave_angle
from zoeppritz.modwan import wang1999
from zoeppritz.invwan import wan1itr


class Test(unittest.TestCase):
    def test_1(self):
        # Two half spaces elastic model
        # vp1, vp2 = 3.0, 2.0
        # vs1, vs2 = 1.5, 1.0
        # ro1, ro2 = 2.3, 2.0
        vp1, vp2 = 3.0, 3.3
        vs1, vs2 = 1.5, 1.7
        ro1, ro2 = 2.3, 2.4

        # Change parameterization
        ro_rd, vp_rd, vs_rd, vs_vp_ratio = \
            elapar_hs2delta(vp1, vs1, ro1, vp2, vs2, ro2)

        # Define angles
        angles = np.arange(0, 60, 6)

        # Calculate the reflection amplitude or b in Ax=b
        ave_angles = inc2ave_angle(angles, vp_rd)
        rpp = wang1999(vs_vp_ratio, ro_rd, vp_rd, vs_rd, ave_angles)

        print("Target model:", ro_rd, vp_rd, vs_rd)
        ro_rd_ini = ro_rd * 0.9
        vp_rd_ini = vp_rd * 0.9
        vs_rd_ini = vs_rd * 1.1
        x_ini = (ro_rd_ini, vp_rd_ini, vs_rd_ini)
        print("Initial model:", x_ini)

        for i in range(9):
            x_new = wan1itr(angles, rpp, x_ini, vs_vp_ratio)
            print("Updated", x_new)
            x_ini = x_new

        err = ro_rd - x_new[0]
        self.assertLessEqual(err, 0.001)

    def test_2(self):
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()

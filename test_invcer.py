# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.

import unittest
import numpy as np
from utils import elapar_hs2ratio
from modcer import rpp_cer1977, rps_cer1977
from invcer import cer1itr


class Test(unittest.TestCase):
    def test_1(self):
        # Two half spaces elastic model
        # vp1, vp2 = 3.0, 2.0
        # vs1, vs2 = 1.5, 1.0
        # ro1, ro2 = 2.3, 2.0
        vp1, vp2 = 4.0, 2.0
        vs1, vs2 = 2.0, 1.0
        ro1, ro2 = 2.4, 2.0

        # Change parameterization
        r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)

        # Define angles
        angles = np.arange(1, 60, 6)

        # Calculate the reflection amplitude or b in Ax=b
        m = len(angles)
        rpp = np.zeros(m)
        rps = np.zeros(m)
        for i in range(m):
            angle = angles[i]
            amp, pha = rpp_cer1977(r1, r2, r3, r4, angle)
            rpp[i] = amp
            amp, pha = rps_cer1977(r1, r2, r3, r4, angle, amp_type='abs')
            rps[i] = amp

        print("Target model:", r1, r2, r3, r4)
        r1_ini = 2.4 / 4.0
        r2_ini = 2.2 / 4.0
        r3_ini = 1.3 / 4.0
        r4_ini = 1.6 / 2.4
        x_ini = (r1_ini, r2_ini, r3_ini, r4_ini)
        print("Initial model:", x_ini)

        for i in range(5):
            x_new = cer1itr(angles, rpp, x_ini, rps=rps)
            print("Updated", x_new)
            x_ini = x_new

        err = r1 - x_new[0]
        self.assertLessEqual(err, 0.001)


if __name__ == '__main__':
    unittest.main()

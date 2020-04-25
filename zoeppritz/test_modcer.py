# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.

import unittest
# import numpy as np
from zoeppritz.utils import elapar_hs2ratio
from zoeppritz.modcer import rpp_cer1977, rps_cer1977


class Test(unittest.TestCase):
    def test_1(self):
        # Two half spaces elastic model
        vp1, vp2 = 3.0, 2.0
        vs1, vs2 = 1.5, 1.0
        ro1, ro2 = 2.3, 2.0
        # vp1, vp2 = 2.0, 4.0
        # vs1, vs2 = 0.88, 1.54
        # ro1, ro2 = 2.0, 2.3

        # Change parameterization
        r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)

        # Define angles
        # angles = np.arange(0, 90, 1)
        # for angle in angles:
        #     amp, pha = rpp_cer1977(r1, r2, r3, r4, angle, amp_type='abs')
        #     amp, pha = rpp_cer1977(r1, r2, r3, r4, angle, amp_type='real')
        #     print("ang,amp,pha =", angle, amp, pha)

        angle = 0
        amp, pha = rpp_cer1977(r1, r2, r3, r4, angle, amp_type='real')
        amp_truth = -0.266055
        err = amp_truth - amp
        self.assertLessEqual(err, 0.001)

    def test_2(self):
        # Two half spaces elastic model
        vp1, vp2 = 5.72, 2.87
        vs1, vs2 = 2.93, 1.61
        ro1, ro2 = 2.86, 2.14

        # Change parameterization
        r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)

        # Define angles
        # angles = np.arange(0, 90, 1)
        # for angle in angles:
        #     amp, pha = rps_cer1977(r1, r2, r3, r4, angle, amp_type='abs')
        #     amp, pha = rps_cer1977(r1, r2, r3, r4, angle, amp_type='real')
        #     print("ang,amp,pha =", angle, amp, pha)

        angle = 10
        amp, pha = rps_cer1977(r1, r2, r3, r4, angle, amp_type='real')
        amp_truth = -0.14823324
        err = amp_truth - amp
        self.assertLessEqual(err, 0.001)


if __name__ == '__main__':
    unittest.main()

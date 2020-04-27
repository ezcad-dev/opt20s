# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.

import unittest
import numpy as np
from zoeppritz.utils import elapar_hs2ratio
from zoeppritz.modcer import rpp_cer1977, rps_cer1977
from zoeppritz.invcer import cer1itr


class Test(unittest.TestCase):
    def test_pp_clean(self):
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

        # print("Target model:", r1, r2, r3, r4)
        r1_ini = 2.4 / 4.0
        r2_ini = 2.2 / 4.0
        r3_ini = 1.3 / 4.0
        r4_ini = 1.6 / 2.4
        x_ini = (r1_ini, r2_ini, r3_ini, r4_ini)
        # print("Initial model:", x_ini)

        for i in range(5):
            x_new = cer1itr(angles, rpp, x_ini, rps=rps)
            # print("Updated", x_new)
            x_ini = x_new

        self.assertLessEqual(r1 - x_new[0], 0.001)
        self.assertLessEqual(r2 - x_new[1], 0.001)
        self.assertLessEqual(r3 - x_new[2], 0.001)
        self.assertLessEqual(r4 - x_new[3], 0.001)

    def test_pp_noise(self):
        # Two half spaces elastic model
        vp1, vp2 = 4.0, 2.0
        vs1, vs2 = 2.0, 1.0
        ro1, ro2 = 2.4, 2.0

        # Change parameterization
        r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)

        # Define angles
        angles = np.arange(1, 60, 1)

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

        # Add noise to data
        # ramp = np.max(rpp) - np.min(rpp)
        # mu, sigma = 0, 0.05 * ramp  # mean and standard deviation
        # noise = np.random.normal(mu, sigma, m)

        # hard code noise to make test result consistent
        noise = np.array([
            -3.40745756e-03, 4.41326891e-03, -1.84969329e-02, -8.57665267e-03,
            -4.64722728e-03, 1.81164323e-02, 2.35764041e-03, 8.66820650e-03,
            2.61862025e-03, 5.60835320e-03, -1.32386200e-02, -1.13868325e-02,
            -3.85411636e-03, 8.30156732e-04, 6.08364262e-03, -1.13107829e-02,
            7.51568819e-03, 8.32391400e-03, 7.18915187e-03, 2.48970883e-03,
            1.42114394e-02, 2.45652884e-04, -4.69414374e-03, 4.60964000e-03,
            1.43935631e-02, -5.88788401e-03, 3.13041871e-03, -6.68177919e-04,
            -6.20489672e-03, -1.68069368e-04, -1.78392131e-02, 8.38724551e-04,
            1.30622636e-03, -9.83497743e-03, -1.17627106e-02, -1.62056738e-02,
            4.62611536e-03, 1.48628494e-02, -1.24973356e-02, -1.01725440e-02,
            7.38562227e-03, 9.21933387e-03, -6.69923701e-03, 6.42089408e-03,
            -4.77129595e-03, 2.33900064e-03, 3.29402557e-05, 9.54770479e-04,
            -1.49280387e-02, -6.65381602e-03, -1.58004300e-02, -7.08064272e-03,
            5.65539007e-04, -2.76684435e-03, -5.60120257e-03, 8.84405490e-03,
            -3.24883460e-03, 5.64724034e-03, -9.45532624e-03,
        ])
        rpp_noisy = rpp + noise

        # ramp = np.max(rps) - np.min(rps)
        # mu, sigma = 0, 0.05 * ramp  # mean and standard deviation
        # noise = np.random.normal(mu, sigma, m)

        # hard code noise to make test result consistent
        noise = np.array([
            -0.0309984, 0.00092359, -0.00770345, -0.03662312, 0.00336188,
            0.00583431, -0.02101242, -0.0248055, -0.00333648, 0.02492424,
            -0.00099495, 0.00944948, -0.00325943, 0.01934984, -0.00704765,
            0.01490579, 0.00779604, 0.02183828, -0.00405295, -0.01820525,
            -0.00446887, 0.01793082, 0.03251096, 0.0026122, 0.01377384,
            -0.01452418, 0.02901279, -0.00881719, 0.02308159, 0.01260138,
            -0.00522267, 0.00769085, 0.02171298, -0.01478435, 0.01349567,
            -0.00778548, -0.01922285, -0.01798599, -0.02126122, -0.00327526,
            0.01550364, 0.00130878, 0.00680895, 0.02670106, -0.05456112,
            0.02081972, 0.02333233, 0.03656901, 0.01069452, -0.01197574,
            0.02639394, 0.01850353, 0.0232636, -0.00037154, -0.01148699,
            0.03056004, 0.006255, 0.01079065, 0.02806546,
        ])
        rps_noisy = rps + noise

        # print("Target model:", r1, r2, r3, r4)
        r1_ini = 2.4 / 4.0
        r2_ini = 2.2 / 4.0
        r3_ini = 1.3 / 4.0
        r4_ini = 1.6 / 2.4
        x_ini = (r1_ini, r2_ini, r3_ini, r4_ini)
        # print("Initial model:", x_ini)

        for i in range(10):
            # x_new = cer1itr(angles, rpp_noisy, x_ini, rps=None)  # fails
            x_new = cer1itr(angles, rpp_noisy, x_ini, rps=rps_noisy)
            # print("Updated", i, x_new)
            x_ini = x_new

        self.assertLessEqual(0.54935233 - x_new[0], 0.001)
        self.assertLessEqual(0.48925867 - x_new[1], 0.001)
        self.assertLessEqual(0.25932287 - x_new[2], 0.001)
        self.assertLessEqual(0.76031868 - x_new[3], 0.001)


if __name__ == '__main__':
    unittest.main()

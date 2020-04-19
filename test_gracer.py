# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.

import unittest
from utils import elapar_hs2ratio
from modcer import rpp_cer1977
from gracer import pdr1, pdr2, pdr3, pdr4


class Test(unittest.TestCase):
    def test_1(self):
        # Two half spaces elastic model
        vp1, vp2 = 3.0, 2.0
        vs1, vs2 = 1.5, 1.0
        ro1, ro2 = 2.3, 2.0

        vp2_a = 1.2
        vp2_b = 2.5
        dx = 0.1
        nx = int((vp2_b - vp2_a) / dx)
        angle = 20.

        for i in range(nx):
            vp2 = vp2_a + (i + 1) * dx
            r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)
            fr_ana = pdr1(r1, r2, r3, r4, angle)

            delta = 0.001
            a1, _ = rpp_cer1977(r1, r2, r3, r4, angle)
            a2, _ = rpp_cer1977(r1 + delta, r2, r3, r4, angle)
            fr_num = (a2 - a1) / delta
            res = fr_ana - fr_num

            print("pdr1 =", vp2, r1, fr_ana, fr_num, abs(res))

    def test_2(self):
        # Two half spaces elastic model
        vp1, vp2 = 3.0, 2.0
        vs1, vs2 = 1.5, 1.0
        ro1, ro2 = 2.3, 2.0

        vs1_a = 1.0
        vs1_b = 2.0
        dx = 0.1
        nx = int((vs1_b - vs1_a) / dx)
        angle = 20.

        for i in range(nx):
            vs1 = vs1_a + (i + 1) * dx
            r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)
            fr_ana = pdr2(r1, r2, r3, r4, angle)

            delta = 0.001
            a1, _ = rpp_cer1977(r1, r2, r3, r4, angle)
            a2, _ = rpp_cer1977(r1, r2 + delta, r3, r4, angle)
            fr_num = (a2 - a1) / delta
            res = fr_ana - fr_num

            print("pdr2 =", vs1, r2, fr_ana, fr_num, abs(res))

    def test_3(self):
        # Two half spaces elastic model
        vp1, vp2 = 3.0, 2.0
        vs1, vs2 = 1.5, 1.0
        ro1, ro2 = 2.3, 2.0

        vs2_a = 0.8
        vs2_b = 1.4
        dx = 0.05
        nx = int((vs2_b - vs2_a) / dx)
        angle = 20.

        for i in range(nx):
            vs2 = vs2_a + (i + 1) * dx
            r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)
            fr_ana = pdr3(r1, r2, r3, r4, angle)

            delta = 0.001
            a1, _ = rpp_cer1977(r1, r2, r3, r4, angle)
            a2, _ = rpp_cer1977(r1, r2, r3 + delta, r4, angle)
            fr_num = (a2 - a1) / delta
            res = fr_ana - fr_num

            print("pdr3 =", vs2, r3, fr_ana, fr_num, abs(res))

    def test_4(self):
        # Two half spaces elastic model
        vp1, vp2 = 3.0, 2.0
        vs1, vs2 = 1.5, 1.0
        ro1, ro2 = 2.3, 2.0

        ro2_a = 1.5
        ro2_b = 2.5
        dx = 0.1
        nx = int((ro2_b - ro2_a) / dx)
        angle = 20.

        for i in range(nx):
            ro2 = ro2_a + (i + 1) * dx
            r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)
            fr_ana = pdr4(r1, r2, r3, r4, angle)

            delta = 0.001
            a1, _ = rpp_cer1977(r1, r2, r3, r4, angle)
            a2, _ = rpp_cer1977(r1, r2, r3, r4 + delta, angle)
            fr_num = (a2 - a1) / delta
            res = fr_ana - fr_num

            print("pdr4 =", ro2, r4, fr_ana, fr_num, abs(res))


if __name__ == '__main__':
    unittest.main()

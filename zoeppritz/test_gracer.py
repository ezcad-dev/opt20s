# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.

import unittest
from zoeppritz.utils import elapar_hs2ratio
from zoeppritz.gracer import gradient


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

        mode, rid = 'PP', 1
        for i in range(nx):
            vp2 = vp2_a + (i + 1) * dx
            r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)
            self.cmp2m(r1, r2, r3, r4, angle, mode, rid)

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

        mode, rid = 'PP', 2
        for i in range(nx):
            vs1 = vs1_a + (i + 1) * dx
            r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)
            self.cmp2m(r1, r2, r3, r4, angle, mode, rid)

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

        mode, rid = 'PP', 3
        for i in range(nx):
            vs2 = vs2_a + (i + 1) * dx
            r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)
            self.cmp2m(r1, r2, r3, r4, angle, mode, rid)

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

        mode, rid = 'PP', 4
        for i in range(nx):
            ro2 = ro2_a + (i + 1) * dx
            r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)
            self.cmp2m(r1, r2, r3, r4, angle, mode, rid)

    def test_5(self):
        # Two half spaces elastic model
        vp1, vp2 = 5.72, 2.87
        vs1, vs2 = 2.93, 1.61
        ro1, ro2 = 2.86, 2.14

        vp2_a = 1.7
        vp2_b = 5.5
        dx = 0.5
        nx = int((vp2_b - vp2_a) / dx)
        angle = 30.

        mode, rid = 'PS', 1
        for i in range(nx):
            vp2 = vp2_a + (i + 1) * dx
            r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)
            self.cmp2m(r1, r2, r3, r4, angle, mode, rid)

    def test_6(self):
        # Two half spaces elastic model
        vp1, vp2 = 5.72, 2.87
        vs1, vs2 = 2.93, 1.61
        ro1, ro2 = 2.86, 2.14

        vs1_a = 1.8
        vs1_b = 3.5
        dx = 0.4
        nx = int((vs1_b - vs1_a) / dx)
        angle = 30.

        mode, rid = 'PS', 2
        for i in range(nx):
            vs1 = vs1_a + (i + 1) * dx
            r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)
            self.cmp2m(r1, r2, r3, r4, angle, mode, rid)

    def test_7(self):
        # Two half spaces elastic model
        vp1, vp2 = 5.72, 2.87
        vs1, vs2 = 2.93, 1.61
        ro1, ro2 = 2.86, 2.14

        vs2_a = 1.7
        vs2_b = 2.7
        dx = 0.2
        nx = int((vs2_b - vs2_a) / dx)
        angle = 30.

        mode, rid = 'PS', 3
        for i in range(nx):
            vs2 = vs2_a + (i + 1) * dx
            r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)
            self.cmp2m(r1, r2, r3, r4, angle, mode, rid)

    def test_8(self):
        # Two half spaces elastic model
        vp1, vp2 = 5.72, 2.87
        vs1, vs2 = 2.93, 1.61
        ro1, ro2 = 2.86, 2.14

        ro2_a = 0.9
        ro2_b = 2.8
        dx = 0.5
        nx = int((ro2_b - ro2_a) / dx)
        angle = 30.

        mode, rid = 'PS', 4
        for i in range(nx):
            ro2 = ro2_a + (i + 1) * dx
            r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)
            self.cmp2m(r1, r2, r3, r4, angle, mode, rid)

    @staticmethod
    def cmp2m(r1, r2, r3, r4, angle, mode, rid):
        method = 'analytic'
        fr_ana = gradient(r1, r2, r3, r4, angle, mode, rid, method=method)
        method = 'numeric'
        fr_num = gradient(r1, r2, r3, r4, angle, mode, rid, method=method)
        res = fr_ana - fr_num
        print(r1, r2, r3, r4, angle, mode, rid, fr_ana, fr_num, abs(res))


if __name__ == '__main__':
    unittest.main()

# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.
"""
Input space partition, input domain modeling (IDM)
"""

import unittest
from zoeppritz.modcer import rpp_cer1977, rps_cer1977


class Test(unittest.TestCase):
    def test_idm_interface(self):
        # Interface-based IDM
        # r1 float, blocks: <0, =0, (0,1), =1, >1
        r1s = [-1, 0, 0.5, 1, 2]
        r2s = [-1, 0, 0.5, 1, 2]
        r3s = [-1, 0, 0.5, 1, 2]
        r4s = [-1, 0, 0.5, 1, 2]
        angles = [-1, 0, 100]  # blocks: <0, =0, >0
        # angles = [-1, 0, 10]  # blocks: <0, =0, >0
        amp_types = ['real', 'abs', 'else']

        # Total number of test cases is 5**4 * 4 * 3 = 7500
        # To illustrate, we fix these parameters, iterate r1 and angle
        r2, r3, r4, amp_type = 0.5, 0.5, 1.0, 'real'

        r1 = -1  # illegal value
        for angle in angles:
            with self.assertRaises(ValueError):
                a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
            with self.assertRaises(ValueError):
                a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        r1 = 0  # illegal value
        for angle in angles:
            with self.assertRaises(ValueError):
                a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
            with self.assertRaises(ValueError):
                a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        r1 = 0.5  # illegal r3 = 0.5 > 0.707 * r1
        for angle in angles:
            with self.assertRaises(ValueError):
                a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
            with self.assertRaises(ValueError):
                a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        r1 = 1
        angle = -1  # illegal value
        with self.assertRaises(ValueError):
            a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        with self.assertRaises(ValueError):
            a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        # No elastic contrast
        r1 = 1
        angle = 0
        a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        self.assertEqual(a, 0)
        self.assertEqual(p, 0)
        a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        self.assertEqual(a, 0)
        self.assertEqual(p, 0)

        r1 = 1
        angle = 100  # illegal value
        with self.assertRaises(ValueError):
            a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        with self.assertRaises(ValueError):
            a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        r1 = 2
        angle = -1  # illegal value
        with self.assertRaises(ValueError):
            a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        with self.assertRaises(ValueError):
            a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        r1 = 2
        angle = 0
        a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        self.assertLessEqual(0.33333333 - a, 0.001)
        self.assertEqual(p, 0)
        a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        self.assertEqual(a, 0)
        self.assertEqual(p, 0)

        r1 = 2
        angle = 100  # illegal value
        with self.assertRaises(ValueError):
            a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        with self.assertRaises(ValueError):
            a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

    def test_idm_functionality(self):
        # Functionality-based IDM

        # r1 physical non negative, blocks: (0,1), >=1
        r1s = [0.9, 1.1]
        angles = [-1, 0, 30, 100]  # blocks: <0, =0, (0,90), >=90

        # base choice
        r1, r2, r3, r4, angle, amp_type = 0.9, 0.5, 0.5, 0.9, 30, 'real'

        r1, angle = 0.9, -1
        with self.assertRaises(ValueError):
            a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        with self.assertRaises(ValueError):
            a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        r1, angle = 0.9, 0
        a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        self.assertLessEqual(-0.104972376 - a, 0.001)
        self.assertEqual(p, 180)
        a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        self.assertEqual(a, 0)
        self.assertEqual(p, 0)

        r1, angle = 0.9, 30
        a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        self.assertLessEqual(-0.10705785 - a, 0.001)
        self.assertEqual(p, 180)
        a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        self.assertLessEqual(-0.04646549 - a, 0.001)
        self.assertEqual(p, 180)

        r1, angle = 0.9, 100
        with self.assertRaises(ValueError):
            a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        with self.assertRaises(ValueError):
            a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        r1, angle = 1.1, -1
        with self.assertRaises(ValueError):
            a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        with self.assertRaises(ValueError):
            a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        r1, angle = 1.1, 0
        a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        self.assertLessEqual(-0.0050251256 - a, 0.001)
        self.assertEqual(p, 180)
        a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        self.assertEqual(a, 0)
        self.assertEqual(p, 0)

        r1, angle = 1.1, 30
        a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        self.assertLessEqual(0.026191689 - a, 0.001)
        self.assertEqual(p, 0)
        a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        self.assertLessEqual(-0.04659759 - a, 0.001)
        self.assertEqual(p, 180)

        r1, angle = 1.1, 100
        with self.assertRaises(ValueError):
            a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        with self.assertRaises(ValueError):
            a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        # Add tests to increase coverage

        r1, r2, r3, r4, angle, amp_type = -1, 0.5, 0.5, 0.9, 30, 'real'
        with self.assertRaises(ValueError):
            a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        r1, r2, r3, r4, angle, amp_type = 0.9, 0.8, 0.5, 0.9, 30, 'real'
        with self.assertRaises(ValueError):
            a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        r1, r2, r3, r4, angle, amp_type = 0.9, 0.5, 0.8, 0.9, 30, 'real'
        with self.assertRaises(ValueError):
            a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        r1, r2, r3, r4, angle, amp_type = 0.9, 0.5, 0.5, -0.9, 30, 'real'
        with self.assertRaises(ValueError):
            a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        r1, r2, r3, r4, angle, amp_type = 0.9, 0.5, 0.5, 0.9, 0, 'abs'
        a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        r1, r2, r3, r4, angle, amp_type = 0.9, 0.5, 0.5, 0.9, 0, 'else'
        with self.assertRaises(ValueError):
            a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        r1, r2, r3, r4, angle, amp_type = 0.9, 0.5, 0.5, 0.9, 30, 'abs'
        a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        r1, r2, r3, r4, angle, amp_type = 0.9, 0.5, 0.5, 0.9, 30, 'else'
        with self.assertRaises(ValueError):
            a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        with self.assertRaises(ValueError):
            a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        # post-critical angle
        r1, r2, r3, r4, angle, amp_type = 2, 0.5, 0.9, 1.2, 45, 'real'
        a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        # Add tests of geophysical interest

        # No contrast: r1=1, r2=r3, r4=1
        r1, r2, r3, r4, angle, amp_type = 1, 0.5, 0.5, 1, 0, 'real'
        a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        self.assertEqual(a, 0)
        self.assertEqual(p, 0)
        a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        self.assertEqual(a, 0)
        self.assertEqual(p, 0)
        r1, r2, r3, r4, angle, amp_type = 1, 0.5, 0.5, 1, 30, 'real'
        a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        self.assertEqual(a, 0)
        self.assertEqual(p, 0)
        a, p = rps_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        self.assertEqual(a, 0)
        self.assertEqual(p, 0)

        # Only density contrast: r1=1, r2=r3, r4!=1

        # NOT handle acoustic media
        # r1, r2, r3, r4, angle, amp_type = 1, 0, 0, 1.1, 30, 'real'
        # a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)

        # It is known that acoustic RC has NO AVA on density-only contrast.
        # Is there AVA with elastic density contrast? Seems there is.
        r1, r2, r3, r4, angle, amp_type = 1, 0.5, 0.5, 1.1, 0, 'real'
        a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        # print('test1', a, p)
        r1, r2, r3, r4, angle, amp_type = 1, 0.5, 0.5, 1.1, 30, 'real'
        a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        # print('test2', a, p)
        r1, r2, r3, r4, angle, amp_type = 1, 0.5, 0.5, 1.1, 50, 'real'
        a, p = rpp_cer1977(r1, r2, r3, r4, angle, amp_type=amp_type)
        # print('test3', a, p)


if __name__ == '__main__':
    unittest.main()

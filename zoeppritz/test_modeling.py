# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.

import unittest
from zoeppritz.modeling import modeling


class Test(unittest.TestCase):
    def test_1(self):
        # Two half spaces elastic model
        vp1, vp2 = 3.0, 2.0
        vs1, vs2 = 1.5, 1.0
        ro1, ro2 = 2.3, 2.0

        model = vp1, vs1, ro1, vp2, vs2, ro2
        angles = '1,2,3'
        equation, reflection = 'linear', 'PP'
        r = modeling(model, angles, equation, reflection)
        self.assertEqual(r.shape, (3, 3))

        angles = '1,2,3'
        equation, reflection = 'linear', 'PS'
        with self.assertRaises(NotImplementedError):
            r = modeling(model, angles, equation, reflection)

        angles = '1,2,3'
        equation, reflection = 'quadratic', 'PP'
        r = modeling(model, angles, equation, reflection)
        self.assertEqual(r.shape, (3, 3))

        angles = '1,2,3'
        equation, reflection = 'quadratic', 'PS'
        with self.assertRaises(NotImplementedError):
            r = modeling(model, angles, equation, reflection)

        angles = '1,2,3'
        equation, reflection = 'zoeppritz', 'PP'
        r = modeling(model, angles, equation, reflection)
        self.assertEqual(r.shape, (3, 3))

        angles = '1,2,3'
        equation, reflection = 'zoeppritz', 'PS'
        r = modeling(model, angles, equation, reflection)
        self.assertEqual(r.shape, (3, 3))

        angles = '1,2,3'
        equation, reflection = 'zoeppritz', 'PT'
        with self.assertRaises(NotImplementedError):
            r = modeling(model, angles, equation, reflection)

        angles = '1,2,3'
        equation, reflection = 'exponential', 'PP'
        with self.assertRaises(NotImplementedError):
            r = modeling(model, angles, equation, reflection)

        angles = '0-60(10)'
        equation, reflection = 'zoeppritz', 'PP'
        r = modeling(model, angles, equation, reflection)
        self.assertEqual(r.shape, (6, 3))


if __name__ == '__main__':
    unittest.main()

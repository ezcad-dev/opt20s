# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.

import unittest
from gopoint.tool.zoep_modeling import zoep_modeling


class Test(unittest.TestCase):
    def test_1(self):
        # Two half spaces elastic model
        vp1, vp2 = 3.0, 2.0
        vs1, vs2 = 1.5, 1.0
        ro1, ro2 = 2.3, 2.0

        model = vp1, vs1, ro1, vp2, vs2, ro2
        angles = '1,2,3'
        equation, reflection = 'zoeppritz', 'PP'
        complexity = 'amplitude'
        object_name = 'testo'
        d = zoep_modeling(model, angles, equation, reflection, complexity,
                          object_name)
        self.assertEqual(d.name, object_name)
        self.assertEqual(d.vertexes['xyz'].shape, (3, 3))

        complexity = 'phase'
        angles = '0-60(10)'
        object_name = 'testp'
        d = zoep_modeling(model, angles, equation, reflection, complexity,
                          object_name)
        self.assertEqual(d.name, object_name)
        self.assertEqual(d.vertexes['xyz'].shape, (6, 3))

        complexity = 'else'
        with self.assertRaises(ValueError):
            d = zoep_modeling(model, angles, equation, reflection, complexity,
                              object_name)


if __name__ == '__main__':
    unittest.main()

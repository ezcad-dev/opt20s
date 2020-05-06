# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.
"""
This module creates points.
"""

import numpy as np
from ..new.new import from_data_array
from zoeppritz.modeling import modeling


def zoep_modeling(model, inc_angles, equation, reflection, complexity,
                  object_name):
    arr = modeling(model, inc_angles, equation, reflection)
    # m, n = ar.shape
    # z = np.zeros([m, 1])
    # arr = np.append(ar, z, axis=1)
    if complexity == 'amplitude':
        pass
    elif complexity == 'phase':
        angles, amplitude, phase = arr[:, 0], arr[:, 1], arr[:, 2]
        arr = np.vstack((angles, phase, amplitude)).T
    else:
        raise ValueError("Unknown value")
    props = ('X', 'Y', 'Z')
    dob = from_data_array(arr, props=props, object_name=object_name)
    return dob

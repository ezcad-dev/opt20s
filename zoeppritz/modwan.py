# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.
"""
Rpp modeling and inversion with quadratic approximation, Wang 1999.
"""

from math import pi
import numpy as np
from zoeppritz.modaki import aki1980_coe


def wang1999(vs_vp_ratio, ro_rd, vp_rd, vs_rd, average_angles,
             amp_type='real'):
    """
    Calculate PP reflection amplitude using Wang (1999)
    quadratic approximation (equation 10)

    For arguments, refer function aki1980()
    """
    A = aki1980_coe(vs_vp_ratio, average_angles)
    x = np.array([ro_rd, vp_rd, vs_rd])
    R = np.matmul(A, x)

    angles = average_angles / 180. * pi
    quad_coef = vs_vp_ratio ** 3 * np.cos(angles) * np.sin(angles) ** 2
    quad_term = (ro_rd + 2 * vs_rd)**2
    R += quad_coef * quad_term
    if amp_type is 'real':
        return R
    elif amp_type is 'abs':
        return np.abs(R)
    else:
        raise ValueError("Unknown amplitude type")

# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.
"""
Linearized inversion
"""

import numpy as np
from modcer import rpp_cer1977
from gracer import gradient


def main():
    pass


def cer1itr(angles, rpp, x_ini, scale=1, constraints={}):
    """
    One iteration of linearized inversion.

    Parameters
    ----------
    angles : array
        incident angles in degrees.
    rpp : array
        Rpp amplitude at the angles, also the b in Ax=b.
    x_ini : tuple
        Initial or starting model of this iteration.
    scale : float
        Scale to the model update
    constraints : dict
        constraints e.g. {'r2': 0.5}

    Returns
    -------
    x_new : tuple
        Updated model
    """
    x_ini_copy = [x for x in x_ini]
    if 'r1' in constraints:
        x_ini_copy[0] = constraints['r1']
    if 'r2' in constraints:
        x_ini_copy[1] = constraints['r2']
    if 'r3' in constraints:
        x_ini_copy[2] = constraints['r3']
    if 'r4' in constraints:
        x_ini_copy[3] = constraints['r4']
    r1_ini, r2_ini, r3_ini, r4_ini = x_ini_copy

    m = len(angles)
    rpp_ini = np.zeros(m)
    A = np.zeros((m, 4))
    for i in range(m):
        angle = angles[i]
        amp, pha = rpp_cer1977(r1_ini, r2_ini, r3_ini, r4_ini, angle)
        rpp_ini[i] = amp

        # Calculate the Jacobian matrix A in Ax=b
        fr1 = gradient(r1_ini, r2_ini, r3_ini, r4_ini, angle, 'PP', 1)
        fr2 = gradient(r1_ini, r2_ini, r3_ini, r4_ini, angle, 'PP', 2)
        fr3 = gradient(r1_ini, r2_ini, r3_ini, r4_ini, angle, 'PP', 3)
        fr4 = gradient(r1_ini, r2_ini, r3_ini, r4_ini, angle, 'PP', 4)
        A[i] = [fr1, fr2, fr3, fr4]

    # A *= -1  # needed when we take abs of negative rpp
    b_dif = rpp - rpp_ini
    lstsq = np.linalg.lstsq(A, b_dif, rcond=None)
    x_dif = lstsq[0]
    if 'r1' in constraints:
        x_dif[0] = 0
    if 'r2' in constraints:
        x_dif[1] = 0
    if 'r3' in constraints:
        x_dif[2] = 0
    if 'r4' in constraints:
        x_dif[3] = 0
    x_new = x_ini_copy + x_dif * scale
    return x_new


if __name__ == '__main__':
    main()

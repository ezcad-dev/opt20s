# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.
"""
Linearized inversion
"""

import numpy as np
from modcer import rpp_cer1977, rps_cer1977
from gracer import gradient


def main():
    pass


def cer1itr(angles, rpp, x_ini, rps=None, fm='numeric', scale=1,
    constraints={}):
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
    rps : array
        Rps amplitude at the angles, append b in Ax=b.
    fm : str
        the method to calculate Frechet derivatives, numeric or analytic
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
    n = 4
    rpp_ini = np.zeros(m)
    A1 = np.zeros((m, n))
    for i in range(m):
        angle = angles[i]
        amp, pha = rpp_cer1977(r1_ini, r2_ini, r3_ini, r4_ini, angle)
        rpp_ini[i] = amp
        # Calculate the Jacobian matrix A in Ax=b
        for j in range(n):
            fre = gradient(r1_ini, r2_ini, r3_ini, r4_ini, angle,
                'PP', j+1, method=fm)
            A1[i, j] = fre
    # A *= -1  # needed when we take abs of negative rpp
    b1 = rpp - rpp_ini

    if rps is None:
        A = A1
        b_dif = b1
    else:
        rps_ini = np.zeros(m)
        A2 = np.zeros((m, n))
        for i in range(m):
            angle = angles[i]
            amp, pha = rps_cer1977(r1_ini, r2_ini, r3_ini, r4_ini, angle,
                amp_type='abs')
            rps_ini[i] = amp
            # Calculate the Jacobian matrix A in Ax=b
            for j in range(n):
                fre = gradient(r1_ini, r2_ini, r3_ini, r4_ini, angle,
                    'PS', j+1, method=fm)
                A2[i, j] = fre
        b2 = rps - rps_ini
        A = np.concatenate((A1, A2), axis=0)
        b_dif = np.concatenate((b1, b2), axis=0)

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

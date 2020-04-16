# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.
"""
Rpp modeling and inversion with quadratic approximation, Wang 1999.
"""

from math import pi
import numpy as np
from utils import elapar_hs2delta
from aki1980 import aki1980_coe, inc2ave_angle


def main():
    wang1999_inv()


def wang1999_inv():
    # Two half spaces elastic model
    # vp1, vp2 = 3.0, 2.0
    # vs1, vs2 = 1.5, 1.0
    # ro1, ro2 = 2.3, 2.0
    vp1, vp2 = 3.0, 3.3
    vs1, vs2 = 1.5, 1.7
    ro1, ro2 = 2.3, 2.4

    # Change parameterization
    ro_rd, vp_rd, vs_rd, vs_vp_ratio = \
        elapar_hs2delta(vp1, vs1, ro1, vp2, vs2, ro2)

    # Define angles
    angles = np.arange(0, 60, 6)

    # Calculate the reflection amplitude or b in Ax=b
    ave_angles = inc2ave_angle(angles, vp_rd)
    rpp = wang1999(vs_vp_ratio, ro_rd, vp_rd, vs_rd, ave_angles)

    print("Target model:", ro_rd, vp_rd, vs_rd)
    ro_rd_ini = ro_rd * 0.9
    vp_rd_ini = vp_rd * 0.9
    vs_rd_ini = vs_rd * 1.1
    x_ini = (ro_rd_ini, vp_rd_ini, vs_rd_ini)
    print("Initial model:", x_ini)

    for i in range(9):
        x_new = inv1itr(angles, rpp, x_ini, vs_vp_ratio)
        print("Updated", x_new)
        x_ini = x_new


def inv1itr(angles, rpp, x_ini, vs_vp_ratio=0.5):
    ro_rd_ini, vp_rd_ini, vs_rd_ini = x_ini

    ave_angles = inc2ave_angle(angles, vp_rd_ini)
    rpp_ini = wang1999(vs_vp_ratio, ro_rd_ini, vp_rd_ini, vs_rd_ini, ave_angles)

    # Calculate the Jacobian matrix A in Ax=b
    A = wang1999_jac(vs_vp_ratio, ro_rd_ini, vs_rd_ini, ave_angles)

    b_dif = rpp - rpp_ini
    lstsq = np.linalg.lstsq(A, b_dif, rcond=None)
    x_dif = lstsq[0]
    x_new = x_ini + x_dif
    return x_new


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


def wang1999_jac(vs_vp_ratio, ro_rd, vs_rd, average_angles):
    """
    Calculate Jacobian matrix of partial derivatives using Wang (1999)
    quadratic approximation (equation 10)

    Parameters
    ----------
    vs_vp_ratio : float
        Vs over Vp ratio of background model
    ro_rd : float
        relative difference of density.
        Say background (average) density is bd = (d2 + d1) / 2.
        ro_rd = (d2 - d1) / bd
    vs_rd : float
        relative difference of Vs.
    average_angles : array
        average of incident and transmission angles.
        The unit is degree. Length is m.

    Returns
    -------
    A : array
        The Jacobian matrix, partial derivatives, shape (m, 3).
        The three columns are for density, Vp, and Vs, respectively.
    """
    angles = average_angles / 180. * pi
    quad_coef = vs_vp_ratio ** 3 * np.cos(angles) * np.sin(angles) ** 2
    quad_ro_pd = 2 * (ro_rd + 2 * vs_rd)
    quad_vs_pd = 2 * quad_ro_pd

    A = aki1980_coe(vs_vp_ratio, average_angles)
    A[:, 0] += quad_coef * quad_ro_pd
    A[:, 2] += quad_coef * quad_vs_pd
    return A


if __name__ == '__main__':
    main()

# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.
"""
Rpp modeling and inversion with linearized approximation, Aki 1980.
"""

from math import pi
import numpy as np
from optimize import optimize_l1
from utils import elapar_hs2delta


def main():
    aki1980_inv()


def aki1980_inv():
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
    t1s = np.arange(0, 60, 6)
    t1s_r = t1s / 180. * pi
    t2s_r = np.arcsin(vp2 / vp1 * np.sin(t1s_r))
    t2s = t2s_r / pi * 180.
    ave_angles = 0.5 * (t1s + t2s)

    # Calculate the coefficient matrix A in Ax=b
    A = aki1980_coe(vs_vp_ratio, ave_angles)

    # Calculate the reflection amplitude or b in Ax=b
    rpp = aki1980(vs_vp_ratio, ro_rd, vp_rd, vs_rd, ave_angles)

    model, x = optimize_l1(A, rpp)
    print('Obj: %g' % model.objVal)
    print('x =', x)

    rm = aki1980(vs_vp_ratio, x[0], x[1], x[2], ave_angles)

    m = len(t1s)
    for i in range(m):
        print(t1s[i], A[i], rpp[i], rm[i])

    print("-------------------------------")
    print("Model ro, vp, vs reldif =", ro_rd, vp_rd, vs_rd)
    print('Gurobi opt x =', x)
    print(np.linalg.lstsq(A, rpp))


def aki1980_coe(vs_vp_ratio, average_angles):
    """
    Get the Coefficient matrix A in Ax=b.

    Parameters
    ----------
    vs_vp_ratio : float
        Vs over Vp ratio of background model
    average_angles : list
        average of incident and transmission angles.
        The unit is degree. Length is m.

    Returns
    -------
    A : array
        The coefficient matrix, shape (m, 3). The three columns are for
        density, Vp, and Vs, respectively.
    """
    angles = average_angles / 180. * pi
    cs = -4 * vs_vp_ratio ** 2 * np.sin(angles) ** 2
    cd = 0.5 * (1 + cs)
    cp = 0.5 / np.cos(angles) ** 2
    A = np.column_stack((cd, cp, cs))
    return A


def aki1980(vs_vp_ratio, ro_rd, vp_rd, vs_rd, average_angles):
    """
    Calculate PP reflection amplitude using Aki and Richards (1980)
    linear approximation to the Zoeppritz equation.

    Equations refer to:
        Equation 1, Shuey, 1985, Geophysics
        Equation 2, Downton and Ursenbach, 2006, Geophysics

    Parameters
    ----------
    vs_vp_ratio : float
        Vs over Vp ratio of background model
    ro_rd : float
        relative difference of density.
        Say background (average) density is bd = (d2 + d1) / 2.
        ro_rd = (d2 - d1) / bd
    vp_rd : float
        relative difference of Vp.
    vs_rd : float
        relative difference of Vs.
    average_angles : array
        average of incident and transmission angles.
        The unit is degree. Length is m.

    Returns
    -------
    rpp : array
        P-wave reflection amplitude. Shape is (m, 1).
    """
    A = aki1980_coe(vs_vp_ratio, average_angles)
    x = np.array([ro_rd, vp_rd, vs_rd])
    R = np.matmul(A, x)
    return np.abs(R)


if __name__ == '__main__':
    main()

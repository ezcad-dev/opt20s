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
    angles = np.arange(0, 60, 6)
    ave_angles = inc2ave_angle(angles, vp_rd)

    # Calculate the coefficient matrix A in Ax=b
    A = aki1980_coe(vs_vp_ratio, ave_angles)

    # Calculate the reflection amplitude or b in Ax=b
    rpp = aki1980(vs_vp_ratio, ro_rd, vp_rd, vs_rd, ave_angles)

    model, x = optimize_l1(A, rpp)
    print('Obj: %g' % model.objVal)
    print('x =', x)

    rm = aki1980(vs_vp_ratio, x[0], x[1], x[2], ave_angles)

    m = len(angles)
    for i in range(m):
        print(angles[i], A[i], rpp[i], rm[i])

    print("-------------------------------")
    print("Model ro, vp, vs reldif =", ro_rd, vp_rd, vs_rd)
    print('Gurobi opt x =', x)
    print(np.linalg.lstsq(A, rpp))


def inc2ave_angle(inc_angles, vp_rd):
    """
    Calculate average angle from incident angle.

    Parameters
    ----------
    inc_angles : array
        incident angles in degrees.
    vp_rd : float
        relative difference of Vp.
        Say background (average) Vp is bv = (vp2 + vp1) / 2.
        vp_rd = (vp2 - vp1) / bd

    Returns
    -------
    ave_angles : array
        average angles in degrees.
    """
    r1 = (2 + vp_rd) / (2 - vp_rd)  # r1 = vp2 / vp1
    t1 = inc_angles / 180. * pi
    t2 = np.arcsin(r1 * np.sin(t1))
    ref_angles = t2 / pi * 180.
    ave_angles = 0.5 * (inc_angles + ref_angles)
    return ave_angles


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


def aki1980(vs_vp_ratio, ro_rd, vp_rd, vs_rd, average_angles,
            amp_type='real'):
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
    amp_type : str
        amplitude type, 'abs' for absolute value, 'real' for the real value.

    Returns
    -------
    rpp : array
        P-wave reflection amplitude. Shape is (m, 1).
    """
    A = aki1980_coe(vs_vp_ratio, average_angles)
    x = np.array([ro_rd, vp_rd, vs_rd])
    R = np.matmul(A, x)
    if amp_type is 'real':
        return R
    elif amp_type is 'abs':
        return np.abs(R)
    else:
        raise ValueError("Unknown amplitude type")


if __name__ == '__main__':
    main()

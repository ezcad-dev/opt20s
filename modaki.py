# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.
"""
Rpp modeling and inversion with linearized approximation, Aki 1980.
"""

from math import pi
import numpy as np


def main():
    pass


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

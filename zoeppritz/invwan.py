# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.
"""
Rpp modeling and inversion with quadratic approximation, Wang 1999.
"""

from math import pi
import numpy as np
from zoeppritz.modaki import aki1980_coe, inc2ave_angle
from zoeppritz.modwan import wang1999


def main():
    pass


def wan1itr(angles, rpp, x_ini, vs_vp_ratio=0.5):
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
    vs_vp_ratio
        Vs/Vp ratio, assumed known a priori.

    Returns
    -------
    x_new : tuple
        Updated model
    """
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

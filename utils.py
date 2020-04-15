# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.
"""
Utility functions
"""

from math import sqrt


def elapar_hs2delta(vp1, vs1, ro1, vp2, vs2, ro2):
    """
    Elastic parameterization, convert half-space to delta model.

    The half-space model has two layers, upper layers has P-wave velocity,
    S-wave velocity, density denoted by vp1, vs1, ro1, respectively.
    The lower layer has vp2, vs2, ro2.
    The delta model refers to one used in publications on linearized
    approximations to Zoeppritz equation, using delta over background.

    Parameters
    ----------
    vp1 : float
        P-wave velocity of the upper layer
    vs1 : float
        S-wave velocity of the upper layer
    ro1 : float
        Density of the upper layer
    vp2 : float
        P-wave velocity of the lower layer
    vs2 : float
        S-wave velocity of the lower layer
    ro2 : float
        Density of the lower layer

    Returns
    -------
    ro_rd : float
        Relative difference of density.
        Say background (average) density is bd = (d2 + d1) / 2.
        ro_rd = (d2 - d1) / bd
    vp_rd : float
        Relative difference of Vp.
    vs_rd : float
        Relative difference of Vs.
    vs_vp_ratio: float
        The ratio of Vs/Vp.
    """
    vp_dif = vp2 - vp1
    vs_dif = vs2 - vs1
    ro_dif = ro2 - ro1
    vp_ave = 0.5 * (vp2 + vp1)
    vs_ave = 0.5 * (vs2 + vs1)
    ro_ave = 0.5 * (ro2 + ro1)
    ro_rd = ro_dif / ro_ave
    vp_rd = vp_dif / vp_ave
    vs_rd = vs_dif / vs_ave
    vs_vp_ratio = vs_ave / vp_ave
    return ro_rd, vp_rd, vs_rd, vs_vp_ratio


def elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2):
    """
    Elastic parameterization, convert half-space to ratio model.

    The half-space model has two layers, upper layers has P-wave velocity,
    S-wave velocity, density denoted by vp1, vs1, ro1, respectively.
    The lower layer has vp2, vs2, ro2.
    The ratio model refers to one used in explicit form of Zoeppritz
    equation (Cerveny, 1977; Zhu, 2012)

    Parameters
    ----------
    vp1 : float
        P-wave velocity of the upper layer
    vs1 : float
        S-wave velocity of the upper layer
    ro1 : float
        Density of the upper layer
    vp2 : float
        P-wave velocity of the lower layer
    vs2 : float
        S-wave velocity of the lower layer
    ro2 : float
        Density of the lower layer

    Returns
    -------
    r1 : float
        vp2 / vp1
    r2 : float
        vs1 / vp1
    r3 : float
        vs2 / vp1
    r4 : float
        ro2 / ro1
    """
    r1 = vp2 / vp1
    r2 = vs1 / vp1
    r3 = vs2 / vp1
    r4 = ro2 / ro1
    return r1, r2, r3, r4


def poisson2vsvp(poisson_ratio):
    """
    Convert Poisson's ratio to Vs/Vp ratio.

    Parameters
    ----------
    poisson_ratio : float
        Poisson's ratio.

    Returns
    -------
    vsvp_ratio : float
        Vs/Vp ratio.
    """
    vpvs_ratio = poisson2vpvs(poisson_ratio)
    vsvp_ratio = 1. / vpvs_ratio
    return vsvp_ratio


def poisson2vpvs(poisson_ratio):
    """
    Convert Poisson's ratio to Vp/Vs ratio.

    Parameters
    ----------
    poisson_ratio : float
        Poisson's ratio.

    Returns
    -------
    vpvs_ratio : float
        Vp/Vs ratio.
    """
    return sqrt(2 * (poisson_ratio - 1) / (2 * poisson_ratio - 1))


def vsvp2poisson(vsvp_ratio):
    """
    Convert Vs/Vp ratio to Poisson's ratio.

    Note that if VS = 0, then Poisson's ratio equals 0.5,
    indicating either a fluid, because shear waves do not pass through fluids,
    or a material that maintains constant volume regardless of stress,
    also known as an ideal incompressible material.

    Poisson's ratio for common rocks
        sandstone   ~0.2
        carbonate   ~0.3
        shale       >0.3
        coal        ~0.4

    Parameters
    ----------
    vsvp_ratio : float
        Vs/Vp ratio.

    Returns
    -------
    poisson_ratio : float
        Poisson's ratio.
    """
    if vsvp_ratio == 0:
        return 0.5
    vpvs_ratio = 1. / vsvp_ratio
    return vpvs2poisson(vpvs_ratio)


def vpvs2poisson(vpvs_ratio):
    """
    Convert Vp/Vs ratio to Poisson's ratio.

    Parameters
    ----------
    vpvs_ratio : float
        Vp/Vs ratio.

    Returns
    -------
    poisson_ratio : float
        Poisson's ratio.
    """
    s = vpvs_ratio ** 2
    return 0.5 * (s - 2) / (s - 1)

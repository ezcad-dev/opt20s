# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.
"""
Zoeppritz partial derivatives
"""

from math import pi, sin, cos, sqrt


def main():
    test_pd()


def test_pd():
    # Two half spaces elastic model
    vp1, vp2 = 3.0, 2.0
    vs1, vs2 = 1.5, 1.0
    ro1, ro2 = 2.3, 2.0

    vp2_a = 1.2
    vp2_b = 2.5
    dx = 0.1
    nx = int((vp2_b - vp2_a) / dx)
    angle = 20.

    from utils import elapar_hs2ratio
    from cerveny1977 import cervey1977

    for i in range(nx):
        vp2 = vp2_a + (i + 1) * dx
        r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)
        fr_ana = pdr1(r1, r2, r3, r4, angle)

        delta = 0.001
        a1, _ = cervey1977(r1, r2, r3, r4, angle)
        a2, _ = cervey1977(r1+delta, r2, r3, r4, angle)
        fr_num = (a2 - a1) / delta
        res = fr_ana - fr_num

        print("pdr1 =", vp2, r1, fr_ana, fr_num, abs(res))

    # Two half spaces elastic model
    vp1, vp2 = 3.0, 2.0
    vs1, vs2 = 1.5, 1.0
    ro1, ro2 = 2.3, 2.0

    vs1_a = 1.0
    vs1_b = 2.0
    dx = 0.1
    nx = int((vs1_b - vs1_a) / dx)

    for i in range(nx):
        vs1 = vs1_a + (i + 1) * dx
        r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)
        fr_ana = pdr2(r1, r2, r3, r4, angle)

        delta = 0.001
        a1, _ = cervey1977(r1, r2, r3, r4, angle)
        a2, _ = cervey1977(r1, r2+delta, r3, r4, angle)
        fr_num = (a2 - a1) / delta
        res = fr_ana - fr_num

        print("pdr2 =", vs1, r2, fr_ana, fr_num, abs(res))

    # Two half spaces elastic model
    vp1, vp2 = 3.0, 2.0
    vs1, vs2 = 1.5, 1.0
    ro1, ro2 = 2.3, 2.0

    vs2_a = 0.8
    vs2_b = 1.4
    dx = 0.05
    nx = int((vs2_b - vs2_a) / dx)

    for i in range(nx):
        vs2 = vs2_a + (i + 1) * dx
        r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)
        fr_ana = pdr3(r1, r2, r3, r4, angle)

        delta = 0.001
        a1, _ = cervey1977(r1, r2, r3, r4, angle)
        a2, _ = cervey1977(r1, r2, r3+delta, r4, angle)
        fr_num = (a2 - a1) / delta
        res = fr_ana - fr_num

        print("pdr3 =", vs2, r3, fr_ana, fr_num, abs(res))

    # Two half spaces elastic model
    vp1, vp2 = 3.0, 2.0
    vs1, vs2 = 1.5, 1.0
    ro1, ro2 = 2.3, 2.0

    ro2_a = 1.5
    ro2_b = 2.5
    dx = 0.1
    nx = int((ro2_b - ro2_a) / dx)

    for i in range(nx):
        ro2 = ro2_a + (i + 1) * dx
        r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)
        fr_ana = pdr4(r1, r2, r3, r4, angle)

        delta = 0.001
        a1, _ = cervey1977(r1, r2, r3, r4, angle)
        a2, _ = cervey1977(r1, r2, r3, r4+delta, angle)
        fr_num = (a2 - a1) / delta
        res = fr_ana - fr_num

        print("pdr4 =", ro2, r4, fr_ana, fr_num, abs(res))


def stability_check(r1, angle):
    # require Vp1 > Vp2, i.e., r1 < 1
    if (r1 * sin(angle)) >= 1:
        raise ValueError("Cannot handle post-critical angle")
    # Physically, Vp always larger than Vs, so
    # alpha1 > beta1, alpha2 > beta2
    # if assumes alpha1 > alpha2, it follows that
    # r1 < 1, r2 < 1, r3 < 1


def getqt(r1, r2, r3, r4, angle):
    Q = 2 * sin(angle)**2 * (r4*r3**2 - r2**2)
    T0 = sin(angle) / cos(angle)
    T1 = r1 * sin(angle) / sqrt(1.0 - r1**2 * sin(angle)**2)
    T2 = r2 * sin(angle) / sqrt(1.0 - r2**2 * sin(angle)**2)
    T3 = r3 * sin(angle) / sqrt(1.0 - r3**2 * sin(angle)**2)
    return Q, T0, T1, T2, T3


def pdr1(r1, r2, r3, r4, angle):
    """
    Zoeppritz partial derivative w.r.t. r1.

    The equation is generated by Mathematica, see Zhu 2012.
    Limitation is no Vp critical angle (alpha1 > alpha2).

    Parameters
    ----------
    r1 : float
        Vp2 / Vp1
    r2 : float
        Vs1 / Vp1
    r3 : float
        Vs2 / Vp1
    r4 : float
        Ro2 / Ro1
    angle : float
        incident angle in degrees

    Returns
    -------
    grad : float
        gradient or partial derivative w.r.t. r1
    """
    if angle <= 0 or angle >= 90:
        raise ValueError("Unsupported angle value {}".format(angle))
    angle = angle * pi / 180.0
    stability_check(r1, angle)
    Q, T0, T1, T2, T3 = getqt(r1, r2, r3, r4, angle)

    # coefficient for T1 T2 T3
    tc1 = (1 + Q)**2
    tc2 = (r4 - Q)**2 
    tc3 = (r4 - Q - 1)**2

    # Notation:
    # term, page #, counting #
    # tm_p#_#
    tm_p1_1 = r4 / r1 * T1**3 * T2
    tm_p1_2 = r4 / r1 * T1 * T2
    # tm_p1_2 = r4 / r1 * T1 * T2
    # tm_p1_1 = tm2 * T1**2
    tm_p1_3 = tc2 / r1 * T1**3 * T3
    tm_p1_4 = tc2 / r1 * T1 * T3
    tm_p1_5 = tc3 / r1 * T0 * T1**3 * T2 * T3
    tm_p1_6 = tc3 / r1 * T0 * T1 * T2 * T3

    tm_b1 = tm_p1_1 + tm_p1_2 + tm_p1_3 + tm_p1_4 + tm_p1_5 + tm_p1_6

    tm_p1_7 = Q**2 + r4 * T1 * T2
    tm_p1_8 = tc2 * T1 * T3 - r4 * T0 * T3
    tm_p1_9 = - tc3 * T0 * T1 * T2 * T3
    tm_p2_1 = - tc1 * T0 * T2

    tm_b2 = tm_p1_7 + tm_p1_8 + tm_p1_9 + tm_p2_1

    tm_c1 = - tm_b1 * tm_b2
    # -- over above -- #

    tm_p2_2 = tm_p1_7
    tm_p2_3 = tc2 * T1 * T3 + r4 * T0 * T3
    tm_p2_4 = - tm_p1_9
    tm_p2_5 = - tm_p2_1

    tm_b3 = tm_p2_2 + tm_p2_3 + tm_p2_4 + tm_p2_5

    tm_c2 = tm_b3**2

    tm_c3 = tm_c1 / tm_c2
    # -- over above -- #

    tm_p2_6 = tm_p1_1
    tm_p2_7 = tm_p1_2
    tm_p2_8 = tm_p1_3
    tm_p2_9 = tm_p1_4
    tm_p2_10 = - tm_p1_5
    tm_p2_11 = - tm_p1_6

    tm_b4 = tm_p2_6 + tm_p2_7 + tm_p2_8 + tm_p2_9 + tm_p2_10 + tm_p2_11

    tm_p2_12 = tm_p1_7
    tm_p3_1 = tm_p2_3
    tm_p3_2 = - tm_p1_9
    tm_p3_3 = - tm_p2_1

    tm_b5 = tm_p2_12 + tm_p3_1 + tm_p3_2 + tm_p3_3

    tm_c4 = tm_b4 / tm_b5

    tm_c5 = tm_c3 + tm_c4

    grad = tm_c5
    return grad


def pdr2(r1, r2, r3, r4, angle):
    """
    Zoeppritz partial derivative w.r.t. r2.

    The equation is generated by Mathematica, see Zhu 2012.
    Limitation is no Vp critical angle (alpha1 > alpha2).

    Parameters
    ----------
    r1 : float
        Vp2 / Vp1
    r2 : float
        Vs1 / Vp1
    r3 : float
        Vs2 / Vp1
    r4 : float
        Ro2 / Ro1
    angle : float
        incident angle in degrees

    Returns
    -------
    grad : float
        gradient or partial derivative w.r.t. r2
    """
    if angle <= 0 or angle >= 90:
        raise ValueError("Unsupported angle value {}".format(angle))
    angle = angle * pi / 180.0
    stability_check(r1, angle)
    Q, T0, T1, T2, T3 = getqt(r1, r2, r3, r4, angle)

    # coefficient for T1 T2 T3
    tc1 = (1 + Q)
    tc2 = (r4 - Q)
    tc3 = (r4 - Q - 1)

    tc4 = 8 * r2 * sin(angle)**2

    # Notation:
    # term, page #, counting #
    # tm_p#_#
    tm_p1_1 = - tc4 * Q
    tm_p1_2 = r4 / r2 * T1 * T2**3
    tm_p1_3 = r4 / r2 * T1 * T2
    tm_p1_4 = tc2 * tc4 * T1 * T3
    tm_p1_5 = tc3 * tc4 * T0 * T1 * T2 * T3
    tm_p1_6 = tc3**2 / r2 * T0 * T1 * T2**3 * T3
    # tm_p1_7 = tm_p1_6 / T2**2
    tm_p1_7 = tc3**2 / r2 * T0 * T1 * T2 * T3
    tm_p1_8 = - tc1 * tc4 * T0 * T2
    tm_p1_9 = tc1**2 / r2 * T0 * T2**3
    tm_p1_10 = tc1**2 / r2 * T0 * T2

    tm_b1 = tm_p1_1 + tm_p1_2 + tm_p1_3 + tm_p1_4 + tm_p1_5 \
        + tm_p1_6 + tm_p1_7 + tm_p1_8 + tm_p1_9 + tm_p1_10

    tm_p1_11 = Q**2 + r4 * T1 * T2
    tm_p2_1 = tc2**2 * T1 * T3
    tm_p2_2 = - r4 * T0 * T3
    # tm_p2_3 = - tm_p1_7 * r2
    tm_p2_3 = - tc3**2 * T0 * T1 * T2 * T3
    tm_p2_4 = - tc1**2 * T0 * T2

    tm_b2 = tm_p1_11 + tm_p2_1 + tm_p2_2 + tm_p2_3 + tm_p2_4

    tm_c1 = - tm_b1 * tm_b2

    tm_p2_5 = tm_p1_11
    tm_p2_6 = tm_p2_1
    tm_p2_7 = - tm_p2_2
    tm_p2_8 = - tm_p2_3
    tm_p2_9 = - tm_p2_4

    tm_c2 = tm_p2_5 + tm_p2_6 + tm_p2_7 + tm_p2_8 + tm_p2_9

    tm_c2 = tm_c2**2

    tm_d1 = tm_c1 / tm_c2

    tm_p2_10 = tm_p1_1
    tm_p2_11 = tm_p1_2
    tm_p2_12 = tm_p1_3
    tm_p2_13 = tm_p1_4
    tm_p2_14 = - tm_p1_5
    tm_p2_15 = - tm_p1_6
    tm_p2_16 = - tm_p1_7
    tm_p3_1 = - tm_p1_8
    tm_p3_2 = - tm_p1_9
    tm_p3_3 = - tm_p1_10

    tm_c3 = tm_p2_10 + tm_p2_11 + tm_p2_12 + tm_p2_13 + tm_p2_14 \
        + tm_p2_15 + tm_p2_16 + tm_p3_1 + tm_p3_2 + tm_p3_3

    tm_p3_4 = tm_p1_11
    tm_p3_5 = tm_p2_1
    tm_p3_6 = - tm_p2_2
    tm_p3_7 = - tm_p2_3
    tm_p3_8 = - tm_p2_4

    tm_c4 = tm_p3_4 + tm_p3_5 + tm_p3_6 + tm_p3_7 + tm_p3_8

    tm_d2 = tm_c3 / tm_c4

    grad = tm_d1 + tm_d2
    return grad


def pdr3(r1, r2, r3, r4, angle):
    """
    Zoeppritz partial derivative w.r.t. r3.

    The equation is generated by Mathematica, see Zhu 2012.
    Limitation is no Vp critical angle (alpha1 > alpha2).

    Parameters
    ----------
    r1 : float
        Vp2 / Vp1
    r2 : float
        Vs1 / Vp1
    r3 : float
        Vs2 / Vp1
    r4 : float
        Ro2 / Ro1
    angle : float
        incident angle in degrees

    Returns
    -------
    grad : float
        gradient or partial derivative w.r.t. r3
    """
    if angle <= 0 or angle >= 90:
        raise ValueError("Unsupported angle value {}".format(angle))
    angle = angle * pi / 180.0
    stability_check(r1, angle)
    Q, T0, T1, T2, T3 = getqt(r1, r2, r3, r4, angle)

    # coefficient for T1 T2 T3
    tc1 = (1 + Q)
    tc2 = (r4 - Q)
    tc3 = (r4 - Q - 1)

    tc4 = 8 * r3 * r4 * sin(angle) ** 2

    # Notation:
    # term, page  # , counting #
    # tm_p  # _#
    tm_p1_1 = tc4 * Q
    tm_p1_2 = - tc2 * tc4 * T1 * T3
    tm_p1_3 = tc2 ** 2 / r3 * T1 * T3 ** 3
    tm_p1_4 = tc2 ** 2 / r3 * T1 * T3
    tm_p1_5 = r4 / r3 * T0 * T3 ** 3
    tm_p1_6 = r4 / r3 * T0 * T3
    tm_p1_7 = - tc3 * tc4 * T0 * T1 * T2 * T3
    tm_p1_8 = tc3 ** 2 / r3 * T0 * T1 * T2 * T3 ** 3
    tm_p1_9 = tc3 ** 2 / r3 * T0 * T1 * T2 * T3
    tm_p1_10 = tc1 * tc4 * T0 * T2

    tm_b1 = tm_p1_1 + tm_p1_2 + tm_p1_3 + tm_p1_4 + tm_p1_5 \
        + tm_p1_6 + tm_p1_7 + tm_p1_8 + tm_p1_9 + tm_p1_10

    tm_p2_1 = Q ** 2 + r4 * T1 * T2
    tm_p2_2 = tc2 ** 2 * T1 * T3
    tm_p2_3 = - r4 * T0 * T3
    tm_p2_4 = - tc3 ** 2 * T0 * T1 * T2 * T3
    tm_p2_5 = - tc1 ** 2 * T0 * T2

    tm_b2 = tm_p2_1 + tm_p2_2 + tm_p2_3 + tm_p2_4 + tm_p2_5

    tm_c1 = - tm_b1 * tm_b2

    tm_p2_6 = tm_p2_1
    tm_p2_7 = tm_p2_2
    tm_p2_8 = - tm_p2_3
    tm_p2_9 = - tm_p2_4
    tm_p2_10 = - tm_p2_5

    tm_c2 = tm_p2_6 + tm_p2_7 + tm_p2_8 + tm_p2_9 + tm_p2_10

    tm_c2 = tm_c2 ** 2

    tm_d1 = tm_c1 / tm_c2

    tm_p2_11 = tm_p1_1
    tm_p2_12 = tm_p1_2
    tm_p2_13 = tm_p1_3
    tm_p2_14 = tm_p1_4
    tm_p2_15 = - tm_p1_5
    tm_p2_16 = - tm_p1_6
    tm_p3_1 = - tm_p1_7
    tm_p3_2 = - tm_p1_8
    tm_p3_3 = - tm_p1_9
    tm_p3_4 = - tm_p1_10

    tm_c3 = tm_p2_11 + tm_p2_12 + tm_p2_13 + tm_p2_14 + tm_p2_15 \
        + tm_p2_16 + tm_p3_1 + tm_p3_2 + tm_p3_3 + tm_p3_4

    tm_p3_5 = tm_p2_1
    tm_p3_6 = tm_p2_2
    tm_p3_7 = - tm_p2_3
    tm_p3_8 = - tm_p2_4
    tm_p3_9 = - tm_p2_5

    tm_c4 = tm_p3_5 + tm_p3_6 + tm_p3_7 + tm_p3_8 + tm_p3_9

    tm_d2 = tm_c3 / tm_c4

    grad = tm_d1 + tm_d2
    return grad


def pdr4(r1, r2, r3, r4, angle):
    """
    Zoeppritz partial derivative w.r.t. r4.

    The equation is generated by Mathematica, see Zhu 2012.
    Limitation is no Vp critical angle (alpha1 > alpha2).

    Parameters
    ----------
    r1 : float
        Vp2 / Vp1
    r2 : float
        Vs1 / Vp1
    r3 : float
        Vs2 / Vp1
    r4 : float
        Ro2 / Ro1
    angle : float
        incident angle in degrees

    Returns
    -------
    grad : float
        gradient or partial derivative w.r.t. r4
    """
    if angle <= 0 or angle >= 90:
        raise ValueError("Unsupported angle value {}".format(angle))
    angle = angle * pi / 180.0
    stability_check(r1, angle)
    Q, T0, T1, T2, T3 = getqt(r1, r2, r3, r4, angle)

    # coefficient for T1 T2 T3
    tc1 = (1 + Q)
    tc2 = (r4 - Q)
    tc3 = (r4 - Q - 1)

    tc4 = 4 * r3 ** 2 * sin(angle) ** 2
    tc5 = 2 * (1 - 2 * r3 ** 2 * sin(angle) ** 2)

    # Notation:
    # term, page  # , counting #
    # tm_p  # _#
    tm_p1_1 = tc4 * Q
    tm_p1_2 = T1 * T2
    tm_p1_3 = tc2 * tc5 * T1 * T3
    tm_p1_4 = T0 * T3
    tm_p1_5 = tc3 * tc5 * T0 * T1 * T2 * T3
    tm_p1_6 = tc1 * tc4 * T0 * T2

    tm_b1 = tm_p1_1 + tm_p1_2 + tm_p1_3 + tm_p1_4 + tm_p1_5 + tm_p1_6

    tm_p1_7 = Q ** 2 + r4 * T1 * T2
    tm_p1_8 = tc2 ** 2 * T1 * T3
    tm_p1_9 = - r4 * T0 * T3
    tm_p1_10 = - tc3 ** 2 * T0 * T1 * T2 * T3
    tm_p1_11 = - tc1 ** 2 * T0 * T2

    tm_b2 = tm_p1_7 + tm_p1_8 + tm_p1_9 + tm_p1_10 + tm_p1_11

    tm_c1 = - tm_b1 * tm_b2

    tm_p2_1 = tm_p1_7
    tm_p2_2 = tm_p1_8
    tm_p2_3 = - tm_p1_9
    tm_p2_4 = - tm_p1_10
    tm_p2_5 = - tm_p1_11

    tm_c2 = tm_p2_1 + tm_p2_2 + tm_p2_3 + tm_p2_4 + tm_p2_5

    tm_c2 = tm_c2 ** 2

    tm_d1 = tm_c1 / tm_c2

    tm_p2_6 = tm_p1_1
    tm_p2_7 = tm_p1_2
    tm_p2_8 = tm_p1_3
    tm_p2_9 = - tm_p1_4
    tm_p2_10 = - tm_p1_5
    tm_p2_11 = - tm_p1_6

    tm_c3 = tm_p2_6 + tm_p2_7 + tm_p2_8 + tm_p2_9 + tm_p2_10 + tm_p2_11

    tm_p2_12 = tm_p1_7
    tm_p2_13 = tm_p1_8
    tm_p2_14 = - tm_p1_9
    tm_p2_15 = - tm_p1_10
    tm_p3_1 = - tm_p1_11

    tm_c4 = tm_p2_12 + tm_p2_13 + tm_p2_14 + tm_p2_15 + tm_p3_1

    tm_d2 = tm_c3 / tm_c4

    grad = tm_d1 + tm_d2
    return grad


if __name__ == '__main__':
    main()
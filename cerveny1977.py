# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.
"""
Zoeppritz equation, explicit form, Cerveny 1977.
"""

from math import pi, sin, sqrt
from cmath import phase
from utils import elapar_hs2ratio


def cervey1977(r1, r2, r3, r4, inc_angle):
    """
    Calculate Rpp using Zoeppritz equation, explicit and exact, Zhu 2012.

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
    inc_angle : float
        incident angle in degrees

    Returns
    -------
    amp : float
        Amplitude, modulus or absolute value of the complex Rpp
    pha : float
        Phase in degrees
    """
    if inc_angle < 0 or inc_angle >= 90:
        raise ValueError("Wrong angle value")

    if inc_angle == 0:
        # normal incidence, elastic Rpp reduces to acoustic
        amp = (r1 * r4 - 1) / (r1 * r4 + 1)
        pha = 0. if amp >= 0 else 180.
        return abs(amp), pha

    angle = inc_angle / 180. * pi

    r0 = 1  # dummy, Vp1 / Vp1
    CT0 = r0 * sin(angle) / complex_sqrt(r0, angle)
    CT1 = r1 * sin(angle) / complex_sqrt(r1, angle)
    CT2 = r2 * sin(angle) / complex_sqrt(r2, angle)
    CT3 = r3 * sin(angle) / complex_sqrt(r3, angle)

    Q = 2 * sin(angle) ** 2 * (r4 * r3 ** 2 - r2 ** 2)
    A = (r4 - Q) ** 2 * CT1 * CT3
    B = (r4 - Q - 1) ** 2 * CT0 * CT1 * CT2 * CT3
    C = (1 + Q) ** 2 * CT0 * CT2
    D = r4 * CT1 * CT2
    E = r4 * CT0 * CT3
    F = Q ** 2

    upp = F - E + D - C - B + A
    low = F + E + D + C + B + A
    rpp = upp / low
    amp = abs(rpp)
    pha = phase(rpp) * 180. / pi

    return amp, pha


def complex_sqrt(r, angle):
    """
    Calculate complex square root in equations 2c-f in Zhu 2012.

    Parameters
    ----------
    r : float
        A ratio, one of r0, r1, r2, r3.
    angle : float
        incident angle in radians.

    Returns
    -------
    crsr : complex
        The square root, refer to the equations.
    """
    rsin = r * sin(angle)
    rss = 1. - rsin ** 2
    rsr = sqrt(abs(rss))
    if rss >= 0:
        crsr = complex(rsr, 0.)
    else:
        crsr = complex(0., -rsr)
    return crsr


def main():
    import numpy as np

    # Two half spaces elastic model
    # vp1, vp2 = 3.0, 2.0
    # vs1, vs2 = 1.5, 1.0
    # ro1, ro2 = 2.3, 2.0
    vp1, vp2 = 2.0, 4.0
    vs1, vs2 = 0.88, 1.54
    ro1, ro2 = 2.0, 2.3

    # Change parameterization
    r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)

    # Define angles
    ans = np.arange(0, 90, 1)
    for angle in ans:
        amp, pha = cervey1977(r1, r2, r3, r4, angle)
        print("ang,amp,pha =", angle, amp, pha)


if __name__ == '__main__':
    main()

# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.
"""
Zoeppritz equation, explicit form, Cerveny 1977.
"""

from math import pi, sin, sqrt
from cmath import phase


def main():
    pass


def rps_cer1977(r1, r2, r3, r4, inc_angle, amp_type='real'):
    """
    Calculate Rps using Zoeppritz equation, explicit and exact, Zhu 2014.

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
    amp_type : str
        amplitude type, 'abs' for absolute value, 'real' for the real part
        of a complex number which preserves the sign.

    Returns
    -------
    amp : float
        Amplitude
    pha : float
        Phase in degrees
    """
    if inc_angle < 0 or inc_angle >= 90:
        raise ValueError("Wrong angle value")

    if inc_angle == 0:
        # normal incidence, no PS conversion
        amp, pha = 0, 0
        return amp, pha

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
    G = Q * (1 + Q)
    H = (r4 - Q) * (r4 - Q - 1) * CT1 * CT3

    upp = 2 * CT2 * (G + H) / r2
    low = F + E + D + C + B + A
    rpp = upp / low
    pha = phase(rpp) * 180. / pi
    if amp_type is 'real':
        return rpp.real, pha
    elif amp_type is 'abs':
        return abs(rpp), pha
    else:
        raise ValueError("Unknown amplitude type")


def rpp_cer1977(r1, r2, r3, r4, inc_angle, amp_type='real'):
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
    amp_type : str
        amplitude type, 'abs' for absolute value, 'real' for the real part
        of a complex number which preserves the sign.

    Returns
    -------
    amp : float
        Amplitude
    pha : float
        Phase in degrees
    """
    if inc_angle < 0 or inc_angle >= 90:
        raise ValueError("Wrong angle value")

    if inc_angle == 0:
        # normal incidence, elastic Rpp reduces to acoustic
        amp = (r1 * r4 - 1) / (r1 * r4 + 1)
        pha = 0. if amp >= 0 else 180.
        if amp_type is 'real':
            return amp, pha
        elif amp_type is 'abs':
            return abs(amp), pha
        else:
            raise ValueError("Unknown amplitude type")

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
    pha = phase(rpp) * 180. / pi
    if amp_type is 'real':
        return rpp.real, pha
    elif amp_type is 'abs':
        return abs(rpp), pha
    else:
        raise ValueError("Unknown amplitude type")


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


if __name__ == '__main__':
    main()

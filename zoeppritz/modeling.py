# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.

import numpy as np
from zoeppritz.utils import elapar_hs2delta, elapar_hs2ratio
from zoeppritz.modaki import aki1980, inc2ave_angle
from zoeppritz.modwan import wang1999
from zoeppritz.modcer import rpp_cer1977, rps_cer1977


def modeling(model, inc_angles, equation, reflection):
    """
    Unified API for GUI call.

    Parameters
    ----------
    model : tuple
        Two half-space elastic model (vp1, vs1, ro1, vp2, vs2, ro2)
    inc_angles : str
        Incident angles in degrees, either comma separated values, or
        1-60(2) means from 1 to 60 with step 2.
    equation : str
        modeling equation, 'linear', 'quadratic', 'zoeppritz'
    reflection : str
        reflection type, 'PP', 'PS'

    Returns
    -------
    rc : array
        amplitude and phase of the reflection coefficients at the angles.
        The array shape is mx3 of columns: incident angle, amplitude, phase.

    """
    # Change parameterization
    vp1, vs1, ro1, vp2, vs2, ro2 = model
    ro_rd, vp_rd, vs_rd, vs_vp_ratio = \
        elapar_hs2delta(vp1, vs1, ro1, vp2, vs2, ro2)
    r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)

    if '-' in inc_angles:
        a, b = inc_angles.split('-')
        c, d = b.split('(')
        e, f = d.split(')')
        a1 = float(a)
        a2 = float(c)
        ad = float(e)
        angles = np.arange(a1, a2, ad)
    else:
        angles = [float(a) for a in inc_angles.split(',')]
    ave_angles = inc2ave_angle(angles, vp_rd)

    m = len(angles)
    a, p = np.zeros(m), np.zeros(m)
    if reflection == 'PP':
        if equation == 'linear':
            a = aki1980(vs_vp_ratio, ro_rd, vp_rd, vs_rd, ave_angles)
            return np.vstack((angles, a, p)).T  # mx3 array
        elif equation == 'quadratic':
            a = wang1999(vs_vp_ratio, ro_rd, vp_rd, vs_rd, ave_angles)
            return np.vstack((angles, a, p)).T  # mx3 array
        elif equation == 'zoeppritz':
            for i in range(m):
                angle = angles[i]
                amp, pha = rpp_cer1977(r1, r2, r3, r4, angle)
                a[i], p[i] = amp, pha
            return np.vstack((angles, a, p)).T  # mx3 array
        else:
            raise NotImplementedError
    elif reflection == 'PS':
        if equation == 'linear':
            raise NotImplementedError
        elif equation == 'quadratic':
            raise NotImplementedError
        elif equation == 'zoeppritz':
            for i in range(m):
                angle = angles[i]
                amp, pha = rps_cer1977(r1, r2, r3, r4, angle)
                a[i], p[i] = amp, pha
            return np.vstack((angles, a, p)).T  # mx3 array
        else:
            raise NotImplementedError
    else:
        raise NotImplementedError

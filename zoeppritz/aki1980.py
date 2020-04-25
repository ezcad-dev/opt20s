# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.
"""
Rpp modeling and inversion with linearized approximation, Aki 1980.
"""

import numpy as np
from zoeppritz.utils import elapar_hs2delta
from zoeppritz.modaki import inc2ave_angle, aki1980_coe, aki1980


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

    from .optimize import optimize_l1
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


if __name__ == '__main__':
    main()

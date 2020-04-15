# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.
"""
Use Gurobi solve Ax=b by minimizing L1 norm of |Ax-b|
"""

import numpy as np
import gurobipy as gp
from gurobipy import GRB


def optimize_l1(A, b):
    # Create a new model
    model = gp.Model("L1 Norm")

    # Create variables
    m, n = A.shape
    x = model.addMVar(shape=n, vtype=GRB.CONTINUOUS, name="x")
    t = model.addMVar(shape=m, vtype=GRB.CONTINUOUS, name="t")

    # Set objective
    ones = np.ones(m)
    model.setObjective(ones @ t, GRB.MINIMIZE)

    # Add constraints
    model.addConstr(A @ x - b <= t, name="axb1")
    model.addConstr(A @ x - b >= -t, name="axb2")
    model.addConstr(t >= 0, name="pos")

    # Optimize model
    model.optimize()

    return model, x.X

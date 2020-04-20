#!/usr/bin/env python
# coding: utf-8

# In[1]:


import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

import sys
sys.path.append("..")
from utils import elapar_hs2delta, elapar_hs2ratio
from modaki import aki1980, aki1980_coe, inc2ave_angle
from modwan import wang1999
from invwan import wan1itr


# In[2]:


# Two half spaces elastic model
vp1, vp2 = 4.0, 2.0
vs1, vs2 = 2.0, 1.0
ro1, ro2 = 2.4, 2.0
model = vp1, vs1, ro1, vp2, vs2, ro2
angles = np.arange(1, 71, 7)  # incident angles
print("Angles", angles)
title = "Data linear, Model linear"

# Change parameterization
vp1, vs1, ro1, vp2, vs2, ro2 = model
ro_rd, vp_rd, vs_rd, vs_vp_ratio =     elapar_hs2delta(vp1, vs1, ro1, vp2, vs2, ro2)
r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)

ave_angles = inc2ave_angle(angles, vp_rd)

# Calculate the reflection amplitude or b in Ax=b
b = aki1980(vs_vp_ratio, ro_rd, vp_rd, vs_rd, ave_angles)

# Calculate the coefficient matrix A in Ax=b
A = aki1980_coe(vs_vp_ratio, ave_angles)

# Solve it by least squares
x = np.linalg.lstsq(A, b, rcond=None)[0]
print("Model parameters: ro_rd, vp_rd, vs_rd")
print("Target model: {:8.6f} {:8.6f} {:8.6f}".format(ro_rd, vp_rd, vs_rd))
print("Invert model: {:8.6f} {:8.6f} {:8.6f}".format(x[0], x[1], x[2]))

ro, vp, vs = x
bsyn = aki1980(vs_vp_ratio, ro, vp, vs, ave_angles)

d = {
    'inc_angles': angles,
    'ave_angles': ave_angles,
    'binv': b,
    'bsyn': bsyn,
}
df = pd.DataFrame(d)

# gca stands for 'get current axis'
ax = plt.gca()
df.plot(kind='scatter', x='inc_angles', y='binv', color='red',
        label='Data', ax=ax)
df.plot(kind='line', x='inc_angles', y='bsyn', color='blue',
        label='Model', ax=ax)

ax.legend(loc='upper left', frameon=False)
plt.title(title)
plt.show()


# In[3]:


# Two half spaces elastic model
vp1, vp2 = 4.0, 2.0
vs1, vs2 = 2.0, 1.0
ro1, ro2 = 2.4, 2.0
model = vp1, vs1, ro1, vp2, vs2, ro2
angles = np.arange(1, 71, 7)  # incident angles
title = "Data quadratic, Model quadratic"

# vp1, vp2 = 3.0, 3.3
# vs1, vs2 = 1.5, 1.7
# ro1, ro2 = 2.3, 2.4
# model = vp1, vs1, ro1, vp2, vs2, ro2
# angles = np.arange(0, 60, 6)  # incident angles

# Change parameterization
vp1, vs1, ro1, vp2, vs2, ro2 = model
ro_rd, vp_rd, vs_rd, vs_vp_ratio =     elapar_hs2delta(vp1, vs1, ro1, vp2, vs2, ro2)
r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)

ave_angles = inc2ave_angle(angles, vp_rd)

# Calculate the reflection amplitude or b in Ax=b
b = wang1999(vs_vp_ratio, ro_rd, vp_rd, vs_rd, ave_angles)

print("Model parameters: ro_rd, vp_rd, vs_rd")
print("Target model: {:8.6f} {:8.6f} {:8.6f}".format(ro_rd, vp_rd, vs_rd))
ro_rd_ini = ro_rd * 0.8
vp_rd_ini = vp_rd * 1.0
vs_rd_ini = vs_rd * 1.0
x_ini = (ro_rd_ini, vp_rd_ini, vs_rd_ini)
print("Initial model: {:8.6f} {:8.6f} {:8.6f}".format(ro_rd_ini, vp_rd_ini, vs_rd_ini))

for i in range(5):
    x_new = wan1itr(angles, b, x_ini, vs_vp_ratio)
    print("Itr {} model {:8.6f} {:8.6f} {:8.6f}".format(i+1,
        x_new[0], x_new[1], x_new[2]))
    x_ini = x_new
    
ro, vp, vs = x_new
bsyn = wang1999(vs_vp_ratio, ro, vp, vs, ave_angles)

d = {
    'inc_angles': angles,
    'ave_angles': ave_angles,
    'binv': b,
    'bsyn': bsyn,
}
df = pd.DataFrame(d)

# gca stands for 'get current axis'
ax = plt.gca()
df.plot(kind='scatter', x='inc_angles', y='binv', color='red',
        label='Data', ax=ax)
df.plot(kind='line', x='inc_angles', y='bsyn', color='blue',
        label='Model', ax=ax)

ax.legend(loc='upper left', frameon=False)
plt.title(title)
plt.show()


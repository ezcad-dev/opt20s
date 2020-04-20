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
from modcer import rpp_cer1977
from invcer import cer1itr


# In[2]:


# Two half spaces elastic model
vp1, vp2 = 4.0, 2.0
vs1, vs2 = 2.0, 1.0
ro1, ro2 = 2.4, 2.0
model = vp1, vs1, ro1, vp2, vs2, ro2
angles = np.arange(1, 71, 7)  # incident angles
title = "Data exact, Model exact"

# Change parameterization
vp1, vs1, ro1, vp2, vs2, ro2 = model
r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)

# Calculate the reflection amplitude or b in Ax=b
m = len(angles)
b = np.zeros(m)
for i in range(m):
    angle = angles[i]
    amp, pha = rpp_cer1977(r1, r2, r3, r4, angle)
    b[i] = amp

print("Model parameters: r1, r2, r3, r4")
print("Target model: {:6.4f} {:6.4f} {:6.4f} {:6.4f}".format(r1, r2, r3, r4))
r1_ini = 2.4 / 4.0
r2_ini = 2.2 / 4.0
r3_ini = 1.3 / 4.0
r4_ini = 1.6 / 2.4
x_ini = (r1_ini, r2_ini, r3_ini, r4_ini)
print("Initial model: {:6.4f} {:6.4f} {:6.4f} {:6.4f}".format(r1_ini, r2_ini, r3_ini, r4_ini))

b_ini = np.zeros(m)
for i in range(m):
    angle = angles[i]
    amp, pha = rpp_cer1977(r1_ini, r2_ini, r3_ini, r4_ini, angle)
    b_ini[i] = amp

for i in range(5):
    x_new = cer1itr(angles, b, x_ini)
    print("Itr {} model {:6.4f} {:6.4f} {:6.4f} {:6.4f}".format(i+1,
        x_new[0], x_new[1], x_new[2], x_new[3]))
    x_ini = x_new

r1n, r2n, r3n, r4n = x_new
bsyn = np.zeros(m)
for i in range(m):
    angle = angles[i]
    amp, pha = rpp_cer1977(r1, r2, r3, r4, angle)
    bsyn[i] = amp    

d = {
    'inc_angles': angles,
    'binv': b,
    'bsyn': bsyn,
    'bini': b_ini,
}
df = pd.DataFrame(d)

# gca stands for 'get current axis'
ax = plt.gca()
df.plot(kind='scatter', x='inc_angles', y='binv', color='red',
        label='Data', ax=ax)
df.plot(kind='line', x='inc_angles', y='bsyn', color='blue',
        label='Model', ax=ax)
df.plot(kind='line', x='inc_angles', y='bini', color='green',
        label='Model initial', ax=ax)

ax.legend(loc='upper left', frameon=False)
plt.title(title)
plt.show()


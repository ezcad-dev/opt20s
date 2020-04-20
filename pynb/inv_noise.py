#!/usr/bin/env python
# coding: utf-8

# In[21]:


import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

import sys
sys.path.append("..")
from utils import elapar_hs2delta, elapar_hs2ratio
from modcer import rpp_cer1977, rps_cer1977
from invcer import cer1itr


# In[30]:


# Two half spaces elastic model
vp1, vp2 = 4.0, 2.0
vs1, vs2 = 2.0, 1.0
ro1, ro2 = 2.4, 2.0
model = vp1, vs1, ro1, vp2, vs2, ro2
angles = np.arange(1, 60, 1)  # incident angles

# Change parameterization
vp1, vs1, ro1, vp2, vs2, ro2 = model
r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)

# Calculate the reflection amplitude or b in Ax=b
m = len(angles)
rpp = np.zeros(m)
rps = np.zeros(m)
for i in range(m):
    angle = angles[i]
    amp, pha = rpp_cer1977(r1, r2, r3, r4, angle)
    rpp[i] = amp
    amp, pha = rps_cer1977(r1, r2, r3, r4, angle, amp_type='abs')
    rps[i] = amp

# Add noise to data
ramp = np.max(rpp) - np.min(rpp)
mu, sigma = 0, 0.05 * ramp  # mean and standard deviation
noise = np.random.normal(mu, sigma, m)
rpp_noisy = rpp + noise
ramp = np.max(rps) - np.min(rps)
mu, sigma = 0, 0.05 * ramp  # mean and standard deviation
noise = np.random.normal(mu, sigma, m)
rps_noisy = rps + noise

print("Model parameters: r1, r2, r3, r4")
print("Target model: {:6.4f} {:6.4f} {:6.4f} {:6.4f}".format(r1, r2, r3, r4))
# r1_ini = 2.0 / 4.0
# r2_ini = 2.0 / 4.0
# r3_ini = 1.0 / 4.0
# r4_ini = 1.6 / 2.4
r1_ini = 2.4 / 4.0
r2_ini = 2.2 / 4.0
r3_ini = 1.3 / 4.0
r4_ini = 1.6 / 2.4
x_ini = (r1_ini, r2_ini, r3_ini, r4_ini)
print("Initial model: {:6.4f} {:6.4f} {:6.4f} {:6.4f}".format(r1_ini, r2_ini, r3_ini, r4_ini))

rpp_ini = np.zeros(m)
rps_ini = np.zeros(m)
for i in range(m):
    angle = angles[i]
    amp, pha = rpp_cer1977(r1_ini, r2_ini, r3_ini, r4_ini, angle)
    rpp_ini[i] = amp
    amp, pha = rps_cer1977(r1_ini, r2_ini, r3_ini, r4_ini, angle, amp_type='abs')
    rps_ini[i] = amp

for i in range(5):
    x_new = cer1itr(angles, rpp_noisy, x_ini, rps=rps_noisy, scale=1)
    print("Itr {} model {:6.4f} {:6.4f} {:6.4f} {:6.4f}".format(i+1,
        x_new[0], x_new[1], x_new[2], x_new[3]))
    x_ini = x_new

r1n, r2n, r3n, r4n = x_new
rpp_syn = np.zeros(m)
rps_syn = np.zeros(m)
for i in range(m):
    angle = angles[i]
    amp, pha = rpp_cer1977(r1n, r2n, r3n, r4n, angle)
    rpp_syn[i] = amp
    amp, pha = rps_cer1977(r1n, r2n, r3n, r4n, angle, amp_type='abs')
    rps_syn[i] = amp

d = {
    'inc_angles': angles,
    'rpp_clean': rpp,
    'rps_clean': rps,
    'rpp_noisy': rpp_noisy,
    'rps_noisy': rps_noisy,
    'rpp_ini': rpp_ini,
    'rps_ini': rps_ini,
    'rpp_syn': rpp_syn,
    'rps_syn': rps_syn,
}
df = pd.DataFrame(d)

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))
# ax = plt.gca()  # gca stands for 'get current axis'
df.plot(kind='scatter', x='inc_angles', y='rpp_noisy', color='red',
        label='Data', ax=axes[0])
df.plot(kind='line', x='inc_angles', y='rpp_syn', color='blue',
        label='Model final', ax=axes[0])
df.plot(kind='line', x='inc_angles', y='rpp_ini', color='green',
        label='Model initial', ax=axes[0])

df.plot(kind='scatter', x='inc_angles', y='rps_noisy', color='red',
        label='Data', ax=axes[1])
df.plot(kind='line', x='inc_angles', y='rps_syn', color='blue',
        label='Model final', ax=axes[1])
df.plot(kind='line', x='inc_angles', y='rps_ini', color='green',
        label='Model initial', ax=axes[1])

axes[0].legend(loc='upper left', frameon=False)
axes[1].legend(loc='upper left', frameon=False)
# plt.title(title)
plt.show()


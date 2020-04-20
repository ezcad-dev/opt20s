#!/usr/bin/env python
# coding: utf-8

# In[12]:


import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

import sys
sys.path.append("..")
from utils import elapar_hs2delta, elapar_hs2ratio
from modaki import aki1980, inc2ave_angle
from modwan import wang1999
from modcer import rpp_cer1977


# In[13]:


def plot(model, angles, title):
    # Change parameterization
    vp1, vs1, ro1, vp2, vs2, ro2 = model
    ro_rd, vp_rd, vs_rd, vs_vp_ratio =         elapar_hs2delta(vp1, vs1, ro1, vp2, vs2, ro2)
    r1, r2, r3, r4 = elapar_hs2ratio(vp1, vs1, ro1, vp2, vs2, ro2)

    # Define angles
    ave_angles = inc2ave_angle(angles, vp_rd)

    # Calculate the reflection amplitude or b in Ax=b
    rpp_aki1980 = aki1980(vs_vp_ratio, ro_rd, vp_rd, vs_rd, ave_angles)
    rpp_wang1999 = wang1999(vs_vp_ratio, ro_rd, vp_rd, vs_rd, ave_angles)

    m = len(angles)
    rpp_cerveny1977 = np.zeros(m)
    for i in range(m):
        angle = angles[i]
        amp, pha = rpp_cer1977(r1, r2, r3, r4, angle)
        rpp_cerveny1977[i] = amp

    d = {
        'inc_angles': angles,
        'ave_angles': ave_angles,
        'rpp_aki1980': rpp_aki1980,
        'rpp_wang1999': rpp_wang1999,
        'rpp_cerveny1977': rpp_cerveny1977,
    }
    df = pd.DataFrame(d)
    
    # gca stands for 'get current axis'
    ax = plt.gca()
    df.plot(kind='line', x='inc_angles', y='rpp_aki1980', color='blue',
            label='Linear', ax=ax)
    df.plot(kind='line', x='inc_angles', y='rpp_wang1999', color='red',
            label='Quadratic', ax=ax)
    df.plot(kind='line', x='inc_angles', y='rpp_cerveny1977', color='black',
            label='Exact', ax=ax)

    ax.legend(loc='lower left', frameon=False)
    plt.title(title)
    plt.show()


# In[14]:


# Two half spaces elastic model
vp1, vs1, ro1 = 2.2, 1.1, 2.1
vp2, vs2, ro2 = 2.0, 1.0, 2.0
model = vp1, vs1, ro1, vp2, vs2, ro2
angles = np.arange(0, 89, 1)  # incident angles
title = "Weak Contrast"
plot(model, angles, title)


# In[15]:


# Two half spaces elastic model
vp1, vs1, ro1 = 4.0, 2.0, 2.4
vp2, vs2, ro2 = 2.0, 1.0, 2.0
model = vp1, vs1, ro1, vp2, vs2, ro2
angles = np.arange(0, 89, 1)  # incident angles
title = "Strong Contrast"
plot(model, angles, title)


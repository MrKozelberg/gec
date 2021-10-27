#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 16:07:11 2021

@author: mrk
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as md
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

import numpy as np

phi0 = np.genfromtxt("../ip/90_90_0.txt")
phi1 = np.genfromtxt("../ip/90_90_1.txt")

deltaz = 0.1

ip0 = phi0[-1,1]
ip1 = phi1[-1,1]

E0 = (phi0[2:,1] - phi0[:-2,1]) / deltaz / 2 
E1 = (phi1[2:,1] - phi1[:-2,1]) / deltaz / 2

plt.plot(E0, phi0[1:-1,0])
plt.plot(E1, phi0[1:-1,0])

print((ip1 - ip0)/ ip0)
print((E1[0] - E0[0])/ E0[0])
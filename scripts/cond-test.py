#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 11:25:57 2021

@author: mrk
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as md
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

import numpy as np

data0 = np.genfromtxt('../conductivity/COND-TEST-0.txt')
data05 = np.genfromtxt('../conductivity/COND-TEST-0.5.txt')
data1 = np.genfromtxt('../conductivity/COND-TEST-1.txt')

fig = plt.figure(figsize = (6,5))

ax = fig.add_subplot(111)
ax.set_ylabel('Altitude, km')
ax.set_xlabel(r'Conductivity, $s^{-1}$')
ax.set_ylim([0,70])
ax.set_xlim([1e-4,10])
ax.set_xscale('log')
ax.plot(data0[:,1], data0[:,0], label = r'1')
ax.plot(data05[:,1], data05[:,0], label = r'2')
ax.plot(data1[:,1], data1[:,0], label = r'3')

ax.grid()
ax.legend()

fig.tight_layout()
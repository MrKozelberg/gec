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
phi025 = np.genfromtxt("../ip/90_90_025.txt")
phi05 = np.genfromtxt("../ip/90_90_05.txt")
phi075 = np.genfromtxt("../ip/90_90_075.txt")
phi1 = np.genfromtxt("../ip/90_90_1.txt")

deltaz = 0.1

E0 = (phi0[2:,1] - phi0[:-2,1]) / deltaz / 2 
E025 = (phi025[2:,1] - phi025[:-2,1]) / deltaz / 2
E05= (phi05[2:,1] - phi05[:-2,1]) / deltaz / 2 
E075 = (phi075[2:,1] - phi075[:-2,1]) / deltaz / 2
E1 = (phi1[2:,1] - phi1[:-2,1]) / deltaz / 2

print("xi {:.5} {:.5} {:.5} {:.5} {:.5}".format(0.0, 0.25, 0.5, 0.75, 1.0))
print("IP {:.5} {:.5} {:.5} {:.5} {:.5}".format(phi0[-1,1], phi025[-1,1], phi05[-1,1], phi075[-1,1], phi1[-1,1]))
print("E {:.5} {:.5} {:.5} {:.5} {:.5}".format(E0[0], E025[0], E05[0], E075[0], E1[0]))


#!/usr/bin/env python3
import numpy as np
from sys import argv

data = np.load(argv[1])
pw = data['pw']
rainc = data['rainc']
cape = data['cape']
cbot = data['cbot']
ctop = data['ctop']
rainc_new = np.zeros((49, 180, 360))
rainc_new[1:-1,:,:] = rainc[2:,:,:] - rainc[:-2,:,:]

alpha = np.zeros((49,180,360))
alpha[1:-1,:,:] = rainc_new[1:-1,:,:] / pw[1:-1,:,:]

np.savez("/home/mrk/GEC/data/"+argv[1][:-4]+"-NEW", cape = cape, ctop = ctop, cbot = cbot, alpha = alpha)

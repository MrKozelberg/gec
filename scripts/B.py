#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 08:57:17 2021

@author: mrkozelberg
"""

import numpy as np
from datetime import datetime, timedelta

a = '2015-12-31'

dt = datetime.strptime(a, '%Y-%m-%d')
data = np.load("/home/mrk/wrfmain/testdata/DATA-"+a+"-00.npz")

pw = data['pw']
rainc = data['rainc']
cape = data['cape']
cbot = data['cbot']
ctop = data['ctop']
rainc_new = np.zeros((49, 180, 360))
rainc_new[1:-1,:,:] = rainc[2:,:,:] - rainc[:-2,:,:]

alpha = np.zeros((49,180,360))
alpha[1:-1,:,:] = rainc_new[1:-1,:,:] / pw[1:-1,:,:]

np.savez("/home/mrk/GEC/data/DATA-"+a+"-00-NEW", cape = cape, ctop = ctop, cbot = cbot, alpha = alpha)

data = np.load("/home/mrk/wrfmain/testdata/DATA-"+a+"-12.npz")

pw = data['pw']
rainc = data['rainc']
cape = data['cape']
cbot = data['cbot']
ctop = data['ctop']
rainc_new = np.zeros((49, 180, 360))
rainc_new[1:-1,:,:] = rainc[2:,:,:] - rainc[:-2,:,:]

alpha = np.zeros((49,180,360))
alpha[1:-1,:,:] = rainc_new[1:-1,:,:] / pw[1:-1,:,:]

np.savez("/home/mrk/GEC/data/DATA-"+a+"-12-NEW", cape = cape, ctop = ctop, cbot = cbot, alpha = alpha)

while a != '2016-12-31':
    a = (dt + timedelta(days=1)).strftime('%Y-%m-%d')
    dt = datetime.strptime(a, '%Y-%m-%d')
    for aa in ["00","12"]:
        data = np.load("/home/mrk/wrfmain/testdata/DATA-"+a+"-"+aa+".npz")
        pw = data['pw']
        rainc = data['rainc']
        cape = data['cape']
        cbot = data['cbot']
        ctop = data['ctop']
        rainc_new = np.zeros((49, 180, 360))
        rainc_new[1:-1,:,:] = rainc[2:,:,:] - rainc[:-2,:,:]
        
        alpha = np.zeros((49,180,360))
        alpha[1:-1,:,:] = rainc_new[1:-1,:,:] / pw[1:-1,:,:]
        
        np.savez("/home/mrk/GEC/data/DATA-"+a+"-"+aa+"-NEW", cape = cape, ctop = ctop, cbot = cbot, alpha = alpha)



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 19:27:12 2021

@author: mrkozelberg 
"""

import numpy as np
from datetime import datetime, timedelta

#subtraction of linear growth while keeping the mean
def correcter(array):
    delta = array[42,:,:] - array[18,:,:]
    result = np.zeros(np.shape(array)) 
    for i in range(len(array)):
        result[i,:,:] = array[i,:,:] - ((i - (42+18)/2)/24)*delta
    return result

def correcter1(array1,array2):
    result = np.zeros((25,180,360))
    result[:18,:,:] += array1[24:42,:,:]
    result[18:,:,:] += array1[18:25,:,:]
    result[6:,:,:] += array2[18:37,:,:]
    result[:6,:,:] += array2[36:42,:,:]
    return result/2

cape = np.zeros((25, 180, 360))
cbot = np.zeros((25, 180, 360))
ctop = np.zeros((25, 180, 360))
pw = np.zeros((25, 180, 360))
rainc_new = np.zeros((25, 180, 360))
t = 0

a = '2015-12-31'
dt = datetime.strptime(a, '%Y-%m-%d')
data00 = np.load("/home/mrk/wrfmain/testdata/DATA-"+a+"-00.npz")
data12 = np.load("/home/mrk/wrfmain/testdata/DATA-"+a+"-12.npz")

cape += correcter1(correcter(data00['cape']),correcter(data12['cape']))
cbot += correcter1(correcter(data00['cbot']),correcter(data12['cbot']))
ctop += correcter1(correcter(data00['ctop']),correcter(data12['ctop']))
pw += correcter1(correcter(data00['pw']),correcter(data12['pw']))

rainc00 = data00['rainc']
rainc_new00 = np.zeros((49, 180, 360))
rainc_new00[1:-1,:,:] = rainc00[2:,:,:] - rainc00[:-2,:,:]
rainc12 = data12['rainc']
rainc_new12 = np.zeros((49, 180, 360))
rainc_new12[1:-1,:,:] = rainc12[2:,:,:] - rainc12[:-2,:,:]
rainc_new += correcter1(correcter(rainc_new00),correcter(rainc_new12))

t += 1

while a != '2016-12-31':
    a = (dt + timedelta(days=1)).strftime('%Y-%m-%d')
    dt = datetime.strptime(a, '%Y-%m-%d')
    data00 = np.load("/home/mrk/wrfmain/testdata/DATA-"+a+"-00.npz")
    data12 = np.load("/home/mrk/wrfmain/testdata/DATA-"+a+"-12.npz")
    
    cape += correcter1(correcter(data00['cape']),correcter(data12['cape']))
    cbot += correcter1(correcter(data00['cbot']),correcter(data12['cbot']))
    ctop += correcter1(correcter(data00['ctop']),correcter(data12['ctop']))
    pw += correcter1(correcter(data00['pw']),correcter(data12['pw']))
        
    rainc00 = data00['rainc']
    rainc_new00 = np.zeros((49, 180, 360))
    rainc_new00[1:-1,:,:] = rainc00[2:,:,:] - rainc00[:-2,:,:]
    
    rainc12 = data12['rainc']
    rainc_new12 = np.zeros((49, 180, 360))
    rainc_new12[1:-1,:,:] = rainc12[2:,:,:] - rainc12[:-2,:,:]
    rainc_new += correcter1(correcter(rainc_new00),correcter(rainc_new12))
    
    t += 1

cape = cape/t
cbot = cbot/t
ctop = ctop/t
pw = pw/t
rainc_new = rainc_new/t
alpha = rainc_new/pw

np.savez("/home/mrk/GEC/data/DATA-2016-FULL-NEW",cape=cape,cbot=cbot,ctop=ctop,alpha=alpha)
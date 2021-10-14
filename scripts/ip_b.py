#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 11:52:17 2021

@author: mrkozelberg
"""

import numpy as np
from datetime import datetime, timedelta

def correcter(array):
    delta = array[42] - array[18]
    result = np.zeros(len(array)) 
    for i in range(len(array)):
        result[i] = array[i] - ((i - (42+18)/2)/24)*delta
    return result 

def correcter1(array1,array2):
    result = np.zeros(25)
    result[:18] += array1[24:42]
    result[18:] += array1[18:25]
    result[6:] += array2[18:37]
    result[:6] += array2[36:42]
    return result/2

ip_sum = np.zeros(25)
t = 0

a = '2015-12-31'
ip00 = np.genfromtxt("ipne/IP-"+a+"-00.txt")
ip12 = np.genfromtxt("ipne/IP-"+a+"-12.txt")
ip_sum += correcter1(correcter(ip00[:,1]),correcter(ip12[:,1]))
t += 1

dt = datetime.strptime(a, '%Y-%m-%d')

while a != '2016-12-31':
    a = (dt + timedelta(days=1)).strftime('%Y-%m-%d')
    dt = datetime.strptime(a, '%Y-%m-%d')
    ip00 = np.genfromtxt("ipne/IP-"+a+"-00.txt")
    ip12 = np.genfromtxt("ipne/IP-"+a+"-12.txt")
    ip_sum += correcter1(correcter(ip00[:,1]),correcter(ip12[:,1]))
    t += 1

np.save("/home/mrk/GEC/ip/IP-2016-FULL-BNE.npy",ip_sum/t)
    

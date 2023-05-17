#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 12:03:40 2022

Preprocessing delle correnti: da coordinate polari a cartesiane ua, va

@author: ggronchi
"""

# this script does:
#    read currents from North Sea experiment paper (Yapa 2002) in polar coordinates 
#    transform currents in cartesian coordinates
#    south-east bound: the x-direction in my model is actually the S-E direction in the reality

import numpy as np
import pandas as pd

depth = np.loadtxt("../NORTHSEA/AMBIENT/vel-rye1996-orig.txt")[:,0]
time1 = np.loadtxt("../NORTHSEA/AMBIENT/vel-rye1996-orig.txt")[:,1]
time2 = np.loadtxt("../NORTHSEA/AMBIENT/vel-rye1996-orig.txt")[:,2]

intensity1 = np.loadtxt("../NORTHSEA/AMBIENT/vel-rye1996-orig.txt")[:,3]
direction1 = np.loadtxt("../NORTHSEA/AMBIENT/vel-rye1996-orig.txt")[:,4]

intensity2 = np.loadtxt("../NORTHSEA/AMBIENT/vel-rye1996-orig.txt")[:,5]
direction2 = np.loadtxt("../NORTHSEA/AMBIENT/vel-rye1996-orig.txt")[:,6]



ua1 = intensity1 * np.cos((direction1)/180 * np.pi)
va1 = intensity1 * np.sin((direction1)/180 * np.pi)

ua2 = intensity2 * np.cos((direction2)/180 * np.pi)
va2 = intensity2 * np.sin((direction2)/180 * np.pi)

print(depth)
print(time1)
print(ua1)
#print(va)

time = np.zeros(7)
outdata = np.stack((depth, time1, time2, ua1, va1, ua2, va2), axis=1)
print(outdata)
datafile=pd.DataFrame(data=outdata, columns=['Depth', 'Time1', 'Time2', 'ua1', 'va1', 'ua2', 'va2'])
#print(datafile)
datafile.to_csv("../NORTHSEA/AMBIENT/vel-rye1996-proc.txt", index=False, header=True, float_format=str, sep='\t', mode='w')
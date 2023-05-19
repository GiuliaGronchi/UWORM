#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 15:50:44 2022

This script plots trajectory, dynamics, mass, u and v velocity components, density difference

@author: ggronchi
"""
import numpy as np
import matplotlib.pyplot as plt


outname='test'
fname='test'

# load model output file and observations file

modelf = np.loadtxt('../OUTPUT/exp_'+fname+'/file1.txt', comments='#', delimiter='\t',skiprows=1)
paramf = np.loadtxt('../OUTPUT/exp_'+fname+'/aparams.txt', comments='#', delimiter='\t',skiprows=1)

time = modelf[:,0]  
timesec = time *60  
mass = modelf[:,1]   
u = modelf[:,2]
w = modelf[:,3]                
x = modelf[:,9]
y = modelf[:,10]
z = modelf[:,11]
r = modelf[:,8]
h = modelf[:,7]
rho = modelf[:,5]
rhoa = modelf[:,6]


time1=paramf[:,0]
alpha=paramf[:,1]
V0=paramf[:,2]
vaproj=paramf[:,3]
rhoa1 = paramf[:,4]
Qs=paramf[:,5]
Qf=paramf[:,6]

#plt.clf() #avoid overplot


# trajectory : z vs x

plt.grid(visible=None, which='major', axis='both')
#plt.plot(x,z, 'r')  #centerline trajectory
plt.plot(x-r,z, 'r', label='WORM')
plt.plot(x+r,z, 'r')
plt.xlabel('x[m]')
plt.ylabel('z[m]')
plt.title('Plume trajectory')
plt.legend()
plt.savefig("../PLOTS/exp_"+outname+"/trajectory.png", dpi=200)
plt.show()


# dynamics : z vs t

plt.grid(visible=None, which='major', axis='both')
plt.plot(time, z, 'r', label='WORM')
plt.xlabel('time[min]')
plt.ylabel('z[m]')
plt.title('Plume dynamics')
plt.legend()
plt.savefig("../PLOTS/exp_"+outname+"/dynamics.png", dpi=200)
plt.show()


# Cylinder mass : m vs t

plt.grid(visible=None, which='major', axis='both')
plt.plot(mass, time, 'r')
plt.ylabel('Element mass[kg]')
plt.xlabel('time[m]')
plt.title('Mass')
plt.savefig("../PLOTS/exp_"+outname+"/mass.png", dpi=200)
plt.show()



# vertical velocity : w vs t

plt.grid(visible=None, which='major', axis='both')
plt.plot(w, z+107, 'g', label='my_sim')
plt.xlabel('w[m/s]')
plt.ylabel('z[m]')
plt.title('Vertical velocity')
plt.savefig("../PLOTS/exp_"+outname+"/verticalvelocity.png", dpi=200)

plt.show()



# density difference between plume and ambient : rho-rhoa vs z

plt.grid(visible=None, which='major', axis='both')
plt.plot(rhoa, z, label=r'$\rho_a $')
plt.plot(rho, z, label=r'$\rho $')
plt.ylabel('z[m]')
plt.xlabel(r'$[kg/m^3]$')
plt.legend()
plt.title('Plume and Ambient Density')
plt.savefig("../PLOTS/exp_"+outname+"/densdiff.png", dpi=200)
plt.show()








#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Thu Jul 22 17:12:41 2021
@author: giulia

This code simulates a plume spill from a given depth under the ocean surface.
The plume is simulated as cylinders whose mass grows in time according to the \
     entrainment of sea water caused by turbulence.
Each cylinder is described by an array of variables, where x=(mass, u, v, w, x,y,z position, temperature and salinity)
which are updated through governing equations (RungeKutta 4 scheme)
xb is before step
xn is next step
parameters are the cylinder radius, thickness, velocity intensity, density, oil density
parameters are also ocean parameters:
ocean velocities, temperature,salinity,density from netcdf cmems files are preprocessed to csv files

'''
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import utils 

# # # Calculation of reduced gravity
# def reduced_g(params):
#     g1 = 9.8 * (params['rhoa']-params['rho'])/params['rho']
#     return g1

# # # Calculation of projected velocity differences
# def proj_vel(params):
#     ua=params['ua']
#     va=params['va']
#     u=params['u']
#     v=params['v']
#     V0=params['V0']
#     vel = np.array([u,v])
#     a_vel = np.array([ua,va])   
#     proj_vel=np.dot(vel,a_vel)/V0
#     return proj_vel

# def vdif(params):
#     vdif=np.abs(params['V0']  - utils.proj_vel(params))  
#     return vdif

# # # Computation of the shear flux
# def shear_entrain_yapa(params):
#     b=params['b']   
#     h=params['h']  
#     alpha=params['alpha'] 
#     # eq. entrainment following Yapa et al., (1997)                        
#     Qs = 2 * np.pi * b * h * alpha * utils.vdif(params)
#     return Qs

# # # Computation of the forced flux
# def forced_entrain_yapa(params):
#     ua=params['ua']
#     va=params['va']
#     b=params['b']    
#     phi=params['phi']  
#     theta=params['theta']   
#     bb=params['bb']
#     phib=params['phib']  
#     thetab=params['thetab'] 
#     ds=params['ds']
    
#     Qfx = ua * ( np.pi * b * (b-bb)* np.absolute(np.cos(theta)*np.cos(phi)) + \
#                 2 * b* ds * np.sqrt(1 - np.cos(theta)**2 * np.cos(phi)**2 ) + \
#                 np.pi * b**2 /2 * np.absolute(np.cos(theta)*np.cos(phi)-np.cos(thetab)*np.cos(phib)) )
#     Qfy = va * ( np.pi * b * (b-bb)* np.absolute(np.sin(theta)*np.cos(phi)) + \
#                 2 * b* ds * np.sqrt(1 - np.sin(theta)**2 * np.cos(phi)**2 ) + \
#                 np.pi * b**2 /2 * np.absolute(np.sin(theta)*np.cos(phi)-np.sin(thetab)*np.cos(phib)) )  
#     Qf= np.sqrt(Qfx**2 + Qfy**2) 
#     return Qf

# # # Computation of the entrainment coefficient
# def entrain_coeff_yapa(params):
#     b=params['b']   
#     phi=params['phi']    
#     b=params['b']     
#     phi=params['phi']  
#     c3=0.057
#     c4=0.554
#     c5=5
#     E=2   
#     invF1= np.sqrt(utils.reduced_g(params) * b)/(E * utils.vdif(params))
#     alpha= np.sqrt(2) *(c3 + c4 * np.sin(phi) * invF1**2 ) / (1 + c5 * utils.proj_vel(params) / utils.vdif(params))
#     return alpha

# # # GOVERNING EQUATIONS
# # # Calculation of dx, where x = ( m, um, vm, wm, cm, x, y, z, Tm, Sm ) 
# def model(x,params):
       
#     rhoa=params['rhoa']
#     ua=params['ua']
#     va=params['va']
#     ca=params['ca']
#     Ta=params['Ta']
#     Sa=params['Sa']
#     g1 = params['g1'] 

#     Qs=shear_entrain_yapa(params)
#     Qf=forced_entrain_yapa(params)
    
#     Qe = Qs + Qf  # total flux is the sum  
     
#     xdot=np.array([rhoa * Qe , rhoa * Qe * ua ,  rhoa * Qe * va , x[0] * g1, rhoa * Qe * ca, \
#                    x[1]/x[0], x[2]/x[0], x[3]/x[0], rhoa* Qe * Ta, rhoa * Qe * Sa]) 

#     return xdot


# # # RUNGE-KUTTA IV 
# def RK4(model,params,xb,dt):  
      
#     k1=dt*model(xb,params)
#     k2=dt*model(xb + k1/2, params)
#     k3=dt*model(xb + k2/2, params)
#     k4=dt*model(xb + k3, params)
#     dx=1/6 * (k1 + 2*k2 + 2*k3 + k4)
    
#     return dx

# # Read ocean interpolated variables from cmems 

df=pd.read_csv('/Users/mhoxhaj/Scrivhox/A_DEV/GG/WORM/INPUT/ocean_profiles_input.csv')
depth=df['depth']
uo=df['uo']
vo=df['vo']
thetao=df['thetao']
so=df['so']
rhoa=df['rhoa']
f_uo = interp1d(depth, uo)
f_vo = interp1d(depth, vo)
f_thetao = interp1d(depth, thetao)
f_so = interp1d(depth, so)
f_rhoa = interp1d(depth, rhoa)
plt.ylim(-200,0)
plt.grid(visible=None, which='major', axis='both')
#plt.plot(f_vo(depth),depth)
#plt.plot(f_uo(depth),depth)
plt.plot(f_rhoa(depth),depth)

# output file name
fname='test'

# total number of steps

tmax=20000

# A cylinder is generated (just one for instantaneous release)
   
for cyl in range (0,1): #range(0,tmax)

    # #  INITIAL PARAMETERS : depth, radius, vel , oil concentration, thickness
    
    z0=-107
    
    # time-step of integration in seconds
    dt=0.02
    #dt=b/V0 
    
    p = {'ca':0.001,'b':0.0508, 'bb':0.0508, 'V0':2.1, 'c0':1, \
    'phi':np.arcsin(1), 'theta':np.arctan2(0,0),\
    'phib':np.arcsin(1), 'thetab':np.arctan2(0,0), 'u':0, 'v':0}
    p['Ta'], p['Sa'], p['rhoa'] = float(f_thetao(z0)), float(f_so(z0)), float(f_rhoa(z0))
    p['ua'], p['va'] =float(f_uo(z0)), float(f_vo(z0))
    
    # Crude oil Troll API 35.8 = 843 kg/m3 at 15.5 C
    
    rho_oil_0 = 900
    p['rho_oil'] = rho_oil_0 - 0.71 * p['Ta']
    p['rho']=p['rho_oil']
    p['h']=p['V0']*dt
    p['ds']=p['V0']*dt
    
    m = p['rho_oil']* np.pi * p['b']**2 * p['h'] 

    xb = np.array([m, 0, 0, p['V0']*m, p['c0']*m, 0, 0, z0, p['Ta']*m, p['Sa']*m]) 

    Flag=True
    p['g1']= utils.reduced_g(p) if Flag==True else 0  

    p['alpha'] = utils.entrain_coeff_yapa(p)  
   
    if cyl==0 :
        # create output parameters file
        paramdata=np.array([0., p['alpha'], utils.proj_vel(p), p['V0'], p['rhoa'], utils.shear_entrain_yapa(p), utils.forced_entrain_yapa(p)])
        s = f'''
        # Initial conditions
        # Mass: {m}
        # Velocity: {p['V0']}
        # Depth: {z0}
        # dt[s]: {dt}
        # Total time [min]: {tmax*dt/60}
        # ua= {p['ua']}
        # va= {p['va']}
        # rhoa= {p['rhoa']}
        # rho_oil= {p['rho_oil']}
        # reduced g = {p['g1']}
        '''   
        print(s)
    
    # Create output variables file 
    
    outdata=np.array([0., xb[0], xb[1]/m, xb[3]/m, xb[4]/m, p['rho'],p['rhoa'], p['h'], p['b'], xb[5], xb[6], xb[7]])
    
    ###################### <<<<<<>>>>>>>>> #########################

    # TIME STARTS
      
    for t in range (0,tmax):
         
        if t < cyl:    
        # the cylinder is not yet released
      
            outdata1=np.array([(t+1)*dt/60, xb[0], xb[1]/m, xb[3]/m, xb[4]/m, p['rho'],p['rhoa'], p['h'], p['b'], xb[5], xb[6], xb[7]])       
        
        else:
        # the cylinder is released
            
            # Update variables
            
            xn = xb + utils.RK4(utils.model,p,xb,dt)
            m = xn[0]
            u,v,w = xn[1]/xn[0], xn[2]/xn[0], xn[3]/xn[0]
            c = xn[4]/xn[0]
            
            
            # Update ambient ocean data at the cyl depth up to surface
            
            if xn[7] < 0 :
                p['Ta'], p['Sa'] = float(f_thetao(xn[7])), float(f_so(xn[7]))
                p['rhoa']=  float(f_rhoa(xn[7]))
                p['ua']=float(f_uo([xn[7]]))
                p['va']=float(f_vo([xn[7]]))
            else :
                break
            
            # Update parameters
            
            p['rho_oil'] = rho_oil_0 - 0.71 * xn[8]/xn[0]
            p['rho'] =  p['rho_oil']*p['rhoa']/(p['rho_oil']*(1-c)+p['rhoa']*c)
            p['Vold']=p['V0']
            p['V0'] = np.sqrt(u**2 + v**2 + w**2)
            p['ds'] = p['V0']*dt
            p['h']= p['V0']/p['Vold'] * p['h']
            p['bb']=p['b']
            p['b']=np.sqrt(m/(p['rho']*np.pi*p['h']))
            p['thetab']=p['theta']
            p['theta']=np.arctan2(v,u)
            p['phib']=p['phi']
            p['phi']=np.arcsin(w/p['V0'])
            p['g1']= utils.reduced_g(p) if Flag==True else 0
            p['alpha']=utils.entrain_coeff_yapa(p)
            p['u']=u
            p['v']=v
            
             
            # When density is equal to ocean, print neutral buoyancy depth and set gravity force to zero with flag
            if np.abs(p['rhoa']-p['rho']) < 0.2 and Flag :
                Flag=False
                print('Neutral buoyancy reached at depth {} m and after {} mins.'.format(xn[7],(t+1)*dt/60 ))
                print('Reached vertical velocity is {:.3} m/s.'.format(w))
                print('Reached radius is {:.3} m.'.format(p['b']))
                print('Reached density is {:.6} kg/m3.'.format(p['rho']))
            
            
            if cyl == 0: 
                paramdata1=np.array([(t+1)*dt/60, p['alpha'], utils.proj_vel(p), p['V0'], p['rhoa'], utils.shear_entrain_yapa(p), utils.forced_entrain_yapa(p)])
                paramdata=np.vstack([paramdata,paramdata1])
                
                
            outdata1=np.array([((t+1))*dt/60, xb[0], xb[1]/m, xb[3]/m, xb[4]/m, p['rho'],p['rhoa'], p['h'], p['b'], xb[5], xb[6], xb[7]])
            outdata = np.vstack([outdata, outdata1]) 
            
            xb=xn
            
    # # PRINT OUTPUT

    datafile=pd.DataFrame(data=outdata, columns=['Time', 'Mass', 'U', 'W', 'C', 'Density', 'A_Density', 'Tkness', 'Radius', 'x', 'y', 'z' ])
    datafile.to_csv(r'/Users/mhoxhaj/Scrivhox/A_DEV/GG/WORM/OUTPUT/exp_'+fname+'/file'+str(cyl+1)+'.txt', index=False, header=True, float_format='%.8f', sep='\t', mode='w')


    paramfile=pd.DataFrame(data=paramdata, columns=['Time[min]', 'alpha', 'va proj', 'V0', 'rhoa', 'Qs', 'Qf'])
    paramfile.to_csv(r'/Users/mhoxhaj/Scrivhox/A_DEV/GG/WORM/OUTPUT/exp_'+fname+'/aparams.txt', index=False, header=True, float_format='%.8f', sep='\t', mode='w')


     
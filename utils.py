#!/usr/bin/env python3

import numpy as np

# # Calculation of reduced gravity
def reduced_g(params):
    g1 = 9.8 * (params['rhoa']-params['rho'])/params['rho']
    return g1


# # Calculation of projected velocity differences
def proj_vel(params):
    ua=params['ua']
    va=params['va']
    u=params['u']
    v=params['v']
    V0=params['V0']
    vel = np.array([u,v])
    a_vel = np.array([ua,va])   
    proj_vel=np.dot(vel,a_vel)/V0
    return proj_vel

def vdif(params):
    vdif=np.abs(params['V0']  - proj_vel(params))  
    return vdif


# # Computation of the shear flux
def shear_entrain_yapa(params):
    b=params['b']   
    h=params['h']  
    alpha=params['alpha'] 
    # eq. entrainment following Yapa et al., (1997)                        
    Qs = 2 * np.pi * b * h * alpha * vdif(params)
    return Qs

# # Computation of the forced flux
def forced_entrain_yapa(params):
    ua=params['ua']
    va=params['va']
    b=params['b']    
    phi=params['phi']  
    theta=params['theta']   
    bb=params['bb']
    phib=params['phib']  
    thetab=params['thetab'] 
    ds=params['ds']
    
    Qfx = ua * ( np.pi * b * (b-bb)* np.absolute(np.cos(theta)*np.cos(phi)) + \
                2 * b* ds * np.sqrt(1 - np.cos(theta)**2 * np.cos(phi)**2 ) + \
                np.pi * b**2 /2 * np.absolute(np.cos(theta)*np.cos(phi)-np.cos(thetab)*np.cos(phib)) )
    Qfy = va * ( np.pi * b * (b-bb)* np.absolute(np.sin(theta)*np.cos(phi)) + \
                2 * b* ds * np.sqrt(1 - np.sin(theta)**2 * np.cos(phi)**2 ) + \
                np.pi * b**2 /2 * np.absolute(np.sin(theta)*np.cos(phi)-np.sin(thetab)*np.cos(phib)) )  
    Qf= np.sqrt(Qfx**2 + Qfy**2) 
    return Qf

# # Computation of the entrainment coefficient
def entrain_coeff_yapa(params):
    b=params['b']   
    phi=params['phi']    
    b=params['b']     
    phi=params['phi']  
    c3=0.057
    c4=0.554
    c5=5
    E=2   
    invF1= np.sqrt(reduced_g(params) * b)/(E * vdif(params))
    alpha= np.sqrt(2) *(c3 + c4 * np.sin(phi) * invF1**2 ) / (1 + c5 * proj_vel(params) / vdif(params))
    return alpha

# # GOVERNING EQUATIONS
# # Calculation of dx, where x = ( m, um, vm, wm, cm, x, y, z, Tm, Sm ) 
def model(x,params):
       
    rhoa=params['rhoa']
    ua=params['ua']
    va=params['va']
    ca=params['ca']
    Ta=params['Ta']
    Sa=params['Sa']
    g1 = params['g1'] 

    Qs=shear_entrain_yapa(params)
    Qf=forced_entrain_yapa(params)
    
    Qe = Qs + Qf  # total flux is the sum  
     
    xdot=np.array([rhoa * Qe , rhoa * Qe * ua ,  rhoa * Qe * va , x[0] * g1, rhoa * Qe * ca, \
                   x[1]/x[0], x[2]/x[0], x[3]/x[0], rhoa* Qe * Ta, rhoa * Qe * Sa]) 

    return xdot


# # RUNGE-KUTTA IV 
def RK4(model,params,xb,dt):  
      
    k1=dt*model(xb,params)
    k2=dt*model(xb + k1/2, params)
    k3=dt*model(xb + k2/2, params)
    k4=dt*model(xb + k3, params)
    dx=1/6 * (k1 + 2*k2 + 2*k3 + k4)
    
    return dx

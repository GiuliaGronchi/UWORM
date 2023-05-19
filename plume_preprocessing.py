#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 14:35:14 2023

@author: ggronchi
"""
import glob
import xarray as xr
import numpy as np
import pandas as pd
import seawater as sw
from scipy.interpolate import RegularGridInterpolator
ds = xr.open_mfdataset(glob.glob('/Users/ggronchi/Dropbox (CMCC)/Giulia/CMEMS-MOTU/NORTHSEA/RFVL/RFVL_NWS*.nc'))
ds1 = xr.open_mfdataset(glob.glob('/Users/ggronchi/Dropbox (CMCC)/Giulia/CMEMS-MOTU/NORTHSEA/TEMPSAL/*.nc'))
lat_fix=60.17
lon_fix=2.55

time_fix=ds['time'][0]
lat=ds['latitude']
lon=ds['longitude']
depth=ds['depth']
uo=ds['uo'].sel(time=time_fix).fillna(0)
vo=ds['vo'].sel(time=time_fix).fillna(0)

thetao=ds1['thetao'].sel(time=time_fix).fillna(0)
so=ds1['so'].sel(time=time_fix).fillna(0)

uo_depth, vo_depth = np.zeros(depth.size), np.zeros(depth.size)
thetao_depth=np.zeros(depth.size)
so_depth=np.zeros(depth.size)
rhoa_depth=np.zeros(depth.size)

for k,depths in enumerate(depth.values):
    print(k,depths)
    xy_itp_uo = RegularGridInterpolator( (lat.values, lon.values), uo.isel(depth=k).values, bounds_error=False, fill_value=0.001) 
    uxy=float(xy_itp_uo([lat_fix,lon_fix]))
    uo_depth[k]= np.where(uxy!=0,uxy,uo_depth[k-1])
    xy_itp_vo = RegularGridInterpolator( (lat.values, lon.values), vo.isel(depth=k).values, bounds_error=False, fill_value=0.001) 
    vxy=float(xy_itp_vo([lat_fix,lon_fix]))
    vo_depth[k]= np.where(vxy!=0,vxy,vo_depth[k-1])
    
    xy_itp_thetao = RegularGridInterpolator( (lat.values, lon.values), thetao.isel(depth=k).values, bounds_error=False) 
    xy_itp_so = RegularGridInterpolator( (lat.values, lon.values), so.isel(depth=k).values, bounds_error=False) 
    
    txy=float(xy_itp_thetao([lat_fix,lon_fix]))
    sxy=float(xy_itp_so([lat_fix,lon_fix]))
    thetao_depth[k]= np.where(txy!=0,txy,thetao_depth[k-1])
    so_depth[k]= np.where(sxy!=0, sxy, so_depth[k-1])
    rhoa_depth[k] = sw.eos80.dens0(so_depth[k],thetao_depth[k]) 
    
#print(uo_depth)
#print(vo_depth)

data = pd.DataFrame({'depth': -depth.values, 'uo': uo_depth, 'vo': vo_depth, 'thetao': thetao_depth, 'so': so_depth, 'rhoa' : rhoa_depth})
    
oceanfile=pd.DataFrame(data)
oceanfile.to_csv(r'ocean_profiles_input.csv', index=False, header=True, float_format='%.8f', mode='w')





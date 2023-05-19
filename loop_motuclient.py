"""
This script runs a loop of motuclient command lines through a list of coordinates,
time and depth ranges via an input file (csv/xlsx).
"""
import os
import getpass
import pandas as pd
from datetime import datetime, timedelta

# User credentials
#USERNAME = input("Enter your username: ")
#PASSWORD = getpass.getpass("Enter your password: ")

USERNAME='ggronchi'
PASSWORD='Sqlupi82.'

# /!\ ONLY BLOCK TO MODIFY
# Define parameters MEDSEA
#out_dir = "MEDSEA"
#serviceID = "MEDSEA_MULTIYEAR_PHY_006_004"
#productID = "med-cmcc-cur-rean-d"

# NORTHSEA
out_dir = "NORTHSEA"
serviceID = "NWSHELF_MULTIYEAR_PHY_004_009"
#productID = "cmems_mod_nws_phy-uv_my_7km-3D_P1D-m"
#productID = "cmems_mod_nws_phy-t_my_7km-3D_P1D-m"
productID = "cmems_mod_nws_phy-s_my_7km-3D_P1D-m"


# Create output directory if does not exist
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# Boundary dates
start_date = datetime(1995, 9, 18, 12)
#print(f'start date is {start_date}')
end_date = datetime(1995, 9, 30, 12)





#coordinates - medsea
#lon = (10, 30)
#lat = (30 , 40)

#coordinates - northsea
lon = (-20, 13)
lat = (40, 65)

# variable 
#var = "--variable uo --variable vo "
var = "--variable so "

#depth
depth=(1.,6000.)



# Download loop 
while start_date <= end_date:

    # Output filename
    out_name = f'../NORTHSEA/TEMPSAL/SAL_NWS_{start_date.day:02d}{start_date.month:02d}95.nc'
    
    
    # Motuclient command line
    query = f'python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu \
    --service-id {serviceID}-TDS --product-id {productID} \
    --longitude-min {lon[0]} --longitude-max {lon[1]} --latitude-min {lat[0]} --latitude-max {lat[1]}\
    --date-min "{start_date}" --date-max "{start_date}" \
    {var} \
    --depth-min {depth[0]} --depth-max {depth[1]}\
    --out-dir {out_dir} --out-name {out_name} --user {USERNAME} --pwd {PASSWORD}'        

    print(f"============== Running request on {start_date} ==============")
    print(query[:-30])
    
    # Run the command
    os.system(query)
    
    start_date = start_date + timedelta(days=1)
    #print(f'start date plus one day is {start_date}')
    
    
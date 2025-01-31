# Calculate mean chl inside and outside of RCLVs

# Lexi Jones
# Date created: 04/18/23
# Last edited: 04/19/23

import csv,sys,warnings
import numpy as np
import xarray as xr
from matplotlib.path import Path
from matplotlib.patches import Polygon
from config import *

pol = str(sys.argv[1])
start_year = int(sys.argv[2])

chl_dir = project_output_dir + 'OCCCI_8day_v6.0_bounding_box/' 

######################## Functions ########################
def dist_btn_two_pts(lon1,lon2,lat1,lat2):
    return np.sqrt((lon2 - lon1)**2 + (lat2 - lat1)**2)

def circle_eq(radius,center_lon,center_lat):
    theta = np.linspace(0, 2*np.pi, 150)
    R = radius
    a = center_lon + R * np.cos(theta)
    b = center_lat + R * np.sin(theta)
    return a,b

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def lin_eq(x1,x2,y1,y2):
    """ Calculate slope (m) and y-intercept (b) for a linear equation """
    m = (y2-y1)/(x2-x1)
    b = y1 - (m*x1)
    return m,b

######################## Set up chl grid ########################
chl_data = xr.open_dataset(chl_dir + 'ESACCI-OC-L3S-CHLOR_A-MERGED-8D_DAILY_4km_GEO_PML_OCx-20100101-fv6.0_bounding_box.nc')
chl_lat_array = np.array(chl_data.lat)[::-1] #chl data gets read in opposite order so you have to fix it
chl_lon_array = np.array(chl_data.lon) + 360
x,y = np.meshgrid(chl_lon_array,chl_lat_array)
x, y = x.flatten(), y.flatten()
grid_points = np.vstack((x,y)).T 

######################## Read RCLV data ########################
RCLV_data = []
with open(project_output_dir + 'parcels_8day_overlap_32day_RCLVatlas/%s_RCLV_%s_%s_atlas_includes_genesis.csv'%(pol,start_year,start_year+9)) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        if i == 0:
            pass
        else:
            RCLV_data.append(row)
        i += 1
RCLV_data = np.array(RCLV_data,dtype=object)
unique_dates = np.unique([r[0] for r in RCLV_data])

local_anom_data = [] # Date,ID,Radius,Inside Mean Chl,Outside 2R Mean Chl

i = 0
for d in unique_dates:
    
    # Get the chl data
    chl_data = xr.open_dataset(chl_dir + 'ESACCI-OC-L3S-CHLOR_A-MERGED-8D_DAILY_4km_GEO_PML_OCx-%s-fv6.0_bounding_box.nc'%(d))
    chl_data_cropped = np.array(chl_data.chlor_a[0][::-1])
    
    for row in RCLV_data[np.where([r[0]==d for r in RCLV_data])]:
        if i % 1000 == 0:
            print(i)

        # RCLV coordinates
        RCLV_bounds = row[9:]
        RCLV_x_bnds = [float(coord) for coord in RCLV_bounds[0::2]]
        RCLV_y_bnds = [float(coord) for coord in RCLV_bounds[1::2]]
        center_lon,center_lat = float(row[5]),float(row[6])

        # RCLV radius
        dist1 = dist_btn_two_pts(max(RCLV_x_bnds),center_lon,np.median(RCLV_y_bnds),center_lat)
        dist2 = dist_btn_two_pts(min(RCLV_x_bnds),center_lon,np.median(RCLV_y_bnds),center_lat)
        dist3 = dist_btn_two_pts(np.median(RCLV_x_bnds),center_lon,max(RCLV_y_bnds),center_lat)
        dist4 = dist_btn_two_pts(np.median(RCLV_x_bnds),center_lon,min(RCLV_y_bnds),center_lat)
        radius = np.mean((dist1,dist2,dist3,dist4))

        # Get chl data in the area of 2*radius
        circ_x,circ_y = circle_eq(radius*2,center_lon,center_lat)
        poly_2R = Path([(circ_x[p],circ_y[p]) for p in np.arange(0,len(circ_x))])
        grid_2R = poly_2R.contains_points(grid_points)
        mask_2R = grid_2R.reshape(len(chl_lat_array),len(chl_lon_array)) # now you have a mask with points inside a polygon
        chl_2R_masked = np.where(mask_2R==False,np.nan,chl_data_cropped)

        # Get RCLV eddy data
        poly_RCLV = Path([(RCLV_x_bnds[p],RCLV_y_bnds[p]) for p in np.arange(0,len(RCLV_x_bnds))])
        grid_RCLV = poly_RCLV.contains_points(grid_points)
        mask_RCLV = grid_RCLV.reshape(len(chl_lat_array),len(chl_lon_array)) # now you have a mask with points inside a polygon
        chl_RCLV_masked = np.where(mask_RCLV==False,np.nan,chl_2R_masked) # Just the RCLV
        chl_2R_no_RCLV_masked = np.where(mask_RCLV==True,np.nan,chl_2R_masked) # Get 2R data with RCLV removed

        # Add data to array 
        with warnings.catch_warnings():
            warnings.filterwarnings(action='ignore', message='Mean of empty slice')
            local_anom_data.append([d,row[1],radius,np.nanmean(chl_RCLV_masked),np.nanmean(chl_2R_no_RCLV_masked)])

        i += 1
        
np.save(project_output_dir + '%s_RCLV_%s_%s_inside_outside_chl.npy'%(pol,start_year,start_year+9),local_anom_data)        

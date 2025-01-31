# Calculate climatological chl anomaly, seperate chl data by feature type using eddy masks, and seperate by eddy province 
# Lexi Jones
# Date created: 02/01/23
# Last edited: 02/01/23

import sys,os
import numpy as np
import xarray as xr
from config import *

################## Universal variables ##################

mask_dir = project_output_dir + 'eddy_masks_and_chl_2000_to_2020/ID_mask/'
full_domain_save_dir = project_output_dir + 'eddy_masks_and_chl_2000_to_2020/chl_feature_bins_full_domain/individual_dates/'
eddy_provinces_save_dir = project_output_dir + 'eddy_masks_and_chl_2000_to_2020/chl_feature_bins_eddy_provinces/individual_dates/'
chl_dir = project_output_dir + 'OCCCI_8day_v6.0_bounding_box/' 
clim_dir = project_output_dir + 'OCCCI_8day_v6.0_climatologies/'
LAVD_dir = project_output_dir + 'parcels_8day_overlap_32day_LAVD/'

dates = []
for filename in os.listdir(LAVD_dir):
    dates.append(filename[0:8])    
date_list = np.sort(np.unique(dates)).tolist()

months = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

lon1,lon2 = 180,230
lat1,lat2 = 15,30

sample_chl_file = 'ESACCI-OC-L3S-CHLOR_A-MERGED-8D_DAILY_4km_GEO_PML_OCx-20100101-fv6.0_bounding_box.nc'
chl_data = xr.open_dataset(chl_dir + sample_chl_file)
chl_lat1,chl_lat2 = find_nearest(chl_data.lat,lat1),find_nearest(chl_data.lat,lat2)
chl_lon1,chl_lon2 = find_nearest(chl_data.lon,lon1 - 360),find_nearest(chl_data.lon,lon2 - 360) #NOTE: this code works if all data is to the east of -180, need to fix if looking west of -180/180
chl_lat_array = np.array(chl_data.lat)[chl_lat2:chl_lat1][::-1] #chl data gets read in opposite order so you have to fix it
chl_lon_array = np.array(chl_data.lon)[chl_lon1:chl_lon2] + 360

################## Create Eddy Province Grid ##################
### KEY ###
# 0: SE
# 1: Northern lats
# 2: Upper cyclonic lee eddies
# 3: Lower anticyclonic lee eddies

def lin_eq(x1,x2,y1,y2):
    """ Calculate slope (m) and y-intercept (b) for a linear equation """
    m = (y2-y1)/(x2-x1)
    b = y1 - (m*x1)
    return m,b

# Horizontal divide
province_grid = np.zeros((len(chl_lat_array),len(chl_lon_array)))
province_grid[np.where(chl_lat_array >= 23),:] = 1

# Slope for lower lee eddies domain
anti_x1,anti_x2 = 198.5,204.5
anti_y1,anti_y2 = min(chl_lat_array),chl_lat_array[find_nearest(chl_lat_array,19)]
anti_slope,anti_intercept = lin_eq(anti_x1,anti_x2,anti_y1,anti_y2)

for y in np.arange(anti_y1,anti_y2,chl_lat_array[1]-chl_lat_array[0]):
    x_max = (y-anti_intercept)/anti_slope
    province_grid[np.where(np.isclose(chl_lat_array,y)),np.where(chl_lon_array < x_max)] = 3
    
# Vertical part for lower lee eddies at big island
for y in np.arange(anti_y2,chl_lat_array[find_nearest(chl_lat_array,19.5)],chl_lat_array[1]-chl_lat_array[0]):
    x_max = 204.5
    province_grid[np.where(np.isclose(chl_lat_array,y)),np.where(chl_lon_array < x_max)] = 3

# Slope for upper lee eddies domain
cyc_x1,cyc_x2 = 204.5,199
cyc_y1,cyc_y2 = chl_lat_array[find_nearest(chl_lat_array,20)],23
cyc_slope,cyc_intercept = lin_eq(cyc_x1,cyc_x2,cyc_y1,cyc_y2)

for y in np.arange(cyc_y1,cyc_y2,chl_lat_array[1]-chl_lat_array[0]):
    x_max = (y-cyc_intercept)/cyc_slope
    province_grid[np.where(np.isclose(chl_lat_array,y)),np.where(chl_lon_array < x_max)] = 2
    
# Vertical part for upper lee eddies at big island
for y in np.arange(chl_lat_array[find_nearest(chl_lat_array,19.5)],cyc_y1,chl_lat_array[1]-chl_lat_array[0]):
    x_max = 204.5
    province_grid[np.where(np.isclose(chl_lat_array,y)),np.where(chl_lon_array < x_max)] = 2

########################################################################

def chl_by_feature_type(d,chl_anom,eddy_mask,province,save_dir):
    
    # Categorize chl data by cell type (background, coherent cyc, coherent anti, SSH cyc, SSH coherent cyc, SSH anti, SSH coherent anti)
    if province == 'ALL':
        # Background 
        bg_chl = chl_anom[np.where((np.char.find(eddy_mask,'BG') != -1))]

        # Anticyclones
        SSH_anti_chl = chl_anom[np.where((np.char.find(eddy_mask,'AS') != -1) & (np.char.find(eddy_mask,'AR') == -1))]
        coh_anti_chl = chl_anom[np.where((np.char.find(eddy_mask,'AS') == -1) & (np.char.find(eddy_mask,'AR') != -1))]
        SSH_coh_anti_chl = chl_anom[np.where((np.char.find(eddy_mask,'AS') != -1) & (np.char.find(eddy_mask,'AR') != -1))]

        # Cyclones
        SSH_cyc_chl = chl_anom[np.where((np.char.find(eddy_mask,'CS') != -1) & (np.char.find(eddy_mask,'CR') == -1))]
        coh_cyc_chl = chl_anom[np.where((np.char.find(eddy_mask,'CS') == -1) & (np.char.find(eddy_mask,'CR') != -1))]
        SSH_coh_cyc_chl = chl_anom[np.where((np.char.find(eddy_mask,'CS') != -1) & (np.char.find(eddy_mask,'CR') != -1))]

    else: # subset mask by province grid 
        bg_chl = chl_anom[np.where((np.char.find(eddy_mask,'BG') != -1) & (province_grid==province))]
        SSH_anti_chl = chl_anom[np.where((np.char.find(eddy_mask,'AS') != -1) & (np.char.find(eddy_mask,'AR') == -1) & (province_grid==province))]
        coh_anti_chl = chl_anom[np.where((np.char.find(eddy_mask,'AS') == -1) & (np.char.find(eddy_mask,'AR') != -1) & (province_grid==province))]
        SSH_coh_anti_chl = chl_anom[np.where((np.char.find(eddy_mask,'AS') != -1) & (np.char.find(eddy_mask,'AR') != -1) & (province_grid==province))]
        SSH_cyc_chl = chl_anom[np.where((np.char.find(eddy_mask,'CS') != -1) & (np.char.find(eddy_mask,'CR') == -1) & (province_grid==province))]
        coh_cyc_chl = chl_anom[np.where((np.char.find(eddy_mask,'CS') == -1) & (np.char.find(eddy_mask,'CR') != -1) & (province_grid==province))]
        SSH_coh_cyc_chl = chl_anom[np.where((np.char.find(eddy_mask,'CS') != -1) & (np.char.find(eddy_mask,'CR') != -1) & (province_grid==province))]
       
    # Remove nan values from data
    bg_chl = bg_chl[~np.isnan(bg_chl)]
    SSH_anti_chl = SSH_anti_chl[~np.isnan(SSH_anti_chl)]
    coh_anti_chl = coh_anti_chl[~np.isnan(coh_anti_chl)]
    SSH_coh_anti_chl = SSH_coh_anti_chl[~np.isnan(SSH_coh_anti_chl)]
    SSH_cyc_chl = SSH_cyc_chl[~np.isnan(SSH_cyc_chl)]
    coh_cyc_chl = coh_cyc_chl[~np.isnan(coh_cyc_chl)]
    SSH_coh_cyc_chl = SSH_coh_cyc_chl[~np.isnan(SSH_coh_cyc_chl)]
    
    # Save as numpy arrays
    np.save(save_dir + '%s_bg_chl_clim_anom_p%s.npy'%(d,province), bg_chl)
    np.save(save_dir + '%s_SSH_anti_chl_clim_anom_p%s.npy'%(d,province), SSH_anti_chl)
    np.save(save_dir + '%s_coh_anti_chl_clim_anom_p%s.npy'%(d,province), coh_anti_chl)
    np.save(save_dir + '%s_SSH_coh_anti_chl_clim_anom_p%s.npy'%(d,province), SSH_coh_anti_chl)
    np.save(save_dir + '%s_SSH_cyc_chl_clim_anom_p%s.npy'%(d,province), SSH_cyc_chl)
    np.save(save_dir + '%s_coh_cyc_chl_clim_anom_p%s.npy'%(d,province), coh_cyc_chl)
    np.save(save_dir + '%s_SSH_coh_cyc_chl_clim_anom_p%s.npy'%(d,province), SSH_coh_cyc_chl)

for d in date_list: 
    print(d)
    
    # Load 8day chl and monthly climatology
    chl_file = 'ESACCI-OC-L3S-CHLOR_A-MERGED-8D_DAILY_4km_GEO_PML_OCx-%s-fv6.0_bounding_box.nc'%(d)
    month_ind = int(d[4:6]) - 1
    month_clim = np.load(clim_dir + 'ESACCI-OC-L3S-CHLOR_A-MERGED-8D_DAILY_4km_GEO_PML_OCx-fv6.0_bounding_box_%s_AVG.npy'%(months[month_ind])) 
    chl_data = xr.open_dataset(chl_dir + chl_file)
    
    # Chl_data and month clim have the same bounds, crop to the desired bounds
    chl_anom = np.array(chl_data.chlor_a[0,chl_lat2:chl_lat1,chl_lon1:chl_lon2][::-1]) - month_clim[chl_lat2:chl_lat1,chl_lon1:chl_lon2][::-1] 

    # Load the eddy mask
    eddy_mask = np.load('%s%s_eddy_ID_mask_includes_genesis.npy'%(mask_dir,d))
    
    chl_by_feature_type(d,chl_anom,eddy_mask,'ALL',full_domain_save_dir)
    
    # Iterate through the 4 provinces
    for p in np.arange(0,4): 
        chl_by_feature_type(d,chl_anom,eddy_mask,p,eddy_provinces_save_dir)


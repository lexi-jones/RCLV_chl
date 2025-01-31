# project_SSH_and_RCLV_atlas_on_chl_grid_v3.py

# Script creates a mask on the OC-CCI grid with SSH eddies generated from the MATLAB OceanEddies code and coherent eddies pulled from the RCLV atlas.
# Mask output will have 'BG' as background cells, AS# or CS# for SSH eddy cells where # is the ID and A/C indicates anti/cyclonic, and AR# or CR# for RCLVs.

# Updates from project_SSH_and_RCLV_atlas_on_chl_grid_v2.py include new filenames for most up-to-date dataset (OCCCI v6.0, and 2000-2020 atlases)
# Update on project_SSH_and_RCLV_atlas_on_chl_grid.py to do analysis on a full year of data, only 
# outputs the eddy mask, chl analysis moved to 'SLA_and_RCLV_eddy_chl.py'
# Update on project_SSH_and_coherent_eddies_on_grid_v3.py to load in RCLV altas rather than find them within this script.

#
# Lexi Jones-Kellett
# Date created: 01/03/23
# Last edited: 01/31/23

################################################# USER INPUT + IMPORTS  #####################################################################

import numpy as np
import xarray as xr
from matplotlib.path import Path
import csv, sys, os
from copy import copy
from config import *

################## Input arguments (Note: sys.argv[0] is the name of the script) ##################

#year = str(sys.argv[1]) # format: YYYY
lon1,lon2 = 180,230 #float(sys.argv[2]),float(sys.argv[3]) # 0 to 360
lat1,lat2 = 15,30 #float(sys.argv[4]),float(sys.argv[5])

# Output will be saved to these dirs
ID_mask_dir = project_output_dir + 'eddy_masks_and_chl_2000_to_2020/ID_mask/'
age_mask_dir = project_output_dir + 'eddy_masks_and_chl_2000_to_2020/age_mask/'

# Input dirs
RCLV_dir = project_output_dir + 'parcels_8day_overlap_32day_RCLVatlas/'
SLA_dir = project_output_dir + 'SLA_eddies/'
chl_dir = project_output_dir + 'OCCCI_8day_v6.0_bounding_box/'

anti_filename = 'anticyc_eddy_data_minage31_minsize12_tolno3_20000101_to_20191231_8day_subset_v2.csv'
cyc_filename = 'cyclonic_eddy_data_minage31_minsize12_tolno3_20000101_to_20191231_8day_subset_v2.csv'
#RCLV_filename = 'RCLV_20000101_20191227_atlas_no_genesis.csv'
RCLV_filename = 'RCLV_20000101_20191227_atlas_includes_genesis.csv'
sample_chl_file = 'ESACCI-OC-L3S-IOP-MERGED-8D_DAILY_4km_GEO_PML_OCx_QAA-20191227-fv6.0_bounding_box.nc'

########################## Get a list of files in the year of interest ############################

dates = []
for filename in os.listdir(chl_dir):
    dates.append(filename[-30:-22])    
date_list = np.sort(np.unique(dates)).tolist()

###################################### Chl grid setup ######################################

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

chl_data = xr.open_dataset(chl_dir + sample_chl_file)
chl_lat1,chl_lat2 = find_nearest(chl_data.lat,lat1),find_nearest(chl_data.lat,lat2)
chl_lon1,chl_lon2 = find_nearest(chl_data.lon,lon1 - 360),find_nearest(chl_data.lon,lon2 - 360) #NOTE: this code works if all data is to the east of -180, need to fix if looking west of -180/180
chl_lat_array = np.array(chl_data.lat)[chl_lat2:chl_lat1][::-1] #chl data gets read in opposite order so you have to fix it
chl_lon_array = np.array(chl_data.lon)[chl_lon1:chl_lon2] + 360
    
#################################### Load SLA datasets ####################################

def read_SSH_eddy_data(filename,orientation):
    """
    Reads in csv file
    
    orientation: 'A' or 'C'
    """
    eddy_data = []
    with open(SLA_dir + filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        i = 0
        for row in csv_reader:
            if i == 0:
                pass # header
            else:
                #if int(row[2]) >= 32: # uncomment this if not including genesis
                eddy_data.append([row[0]] + [orientation+'S'+row[1]] + [orientation+'S'+str(row[2])] + row[9:]) #date,ID,age,bounds
            i += 1
    return eddy_data

anti_eddy_data = read_SSH_eddy_data(anti_filename,'A')
cyc_eddy_data = read_SSH_eddy_data(cyc_filename,'C')

###################################### Load RCLV dataset ######################################

RCLV_data = []
with open(RCLV_dir + RCLV_filename) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        if i == 0:
            pass # header
       # elif row[0][0:4] == str(year): #get all the data from the desired year
        else:
            if row[2] == 'cyc':
                orientation = 'C'
            else:
                orientation = 'A' 
            RCLV_data.append([row[0]] + [orientation+'R'+row[1]] + [orientation+'R'+str(row[3])] + row[9:]) #date,ID,age,bounds
        i += 1
        
# Make all of the rows the same length
max_row_len = max([len(row) for row in RCLV_data])
for r in np.arange(0,len(RCLV_data)):
    len_diff = max_row_len - len(RCLV_data[r])
    RCLV_data[r] = np.concatenate((RCLV_data[r],['']*len_diff))

################################### CONCATENATE EDDY DATASETS ###################################

RCLV_row_len = len(RCLV_data[0])
anti_len_diff = RCLV_row_len - len(anti_eddy_data[0])
cyc_len_diff = RCLV_row_len - len(cyc_eddy_data[0])

for r in np.arange(0,len(anti_eddy_data)):
    anti_eddy_data[r] = np.concatenate((anti_eddy_data[r],['']*anti_len_diff)) 
for r in np.arange(0,len(cyc_eddy_data)):
    cyc_eddy_data[r] = np.concatenate((cyc_eddy_data[r],['']*np.abs(cyc_len_diff)))
    
full_eddy_dataset = np.concatenate((anti_eddy_data,cyc_eddy_data,RCLV_data)) 
# 'full_eddy_dataset' contains all eddies with format: [date,ID,bounds]
# ID format: orientation + type + #
# Possible combinations: 'AS###','CS###','AR###','CR###'
# A = anticyclone, C = cyclone, S = SLA eddy, R = RSLV

###################################### Project Eddies ######################################

def project_eddy_bnds_on_mask(eddy_mask,eddy_data,date):
    """
    eddy_mask: 2D array initialized with the desired grid for the projection
    eddy_data: 2D array containing eddy data
    """
    
    eddy_ID_mask = copy(eddy_mask)
    eddy_age_mask = copy(eddy_mask)
    
    for pt in eddy_data: 
        if pt[0] == date:
                        
            ID,age,bnds = pt[1],pt[2],pt[3:]

            x_bnds = [float(coord) for coord in bnds[0::2] if str(coord) != '']
            y_bnds = [float(coord) for coord in bnds[1::2] if str(coord) != '']

            # Reformat boundary points to read into matplotlib.Path
            poly_pts = []
            for pt in np.arange(0,len(x_bnds)): #iterate through each point of the polygon
                poly_pts.append((x_bnds[pt],y_bnds[pt])) #reformat with parentheses

            # Create mask from chl grid with points in the eddy are TRUE
            x,y = np.meshgrid(chl_lon_array,chl_lat_array)
            x, y = x.flatten(), y.flatten()
            grid_points = np.vstack((x,y)).T

            poly = Path(poly_pts) # make a polygon
            grid = poly.contains_points(grid_points)
            mask = grid.reshape(len(chl_lat_array),len(chl_lon_array))

            # Update mask with all eddies with this eddy's ID
            x_mask = np.where(mask == True)[0]
            y_mask = np.where(mask == True)[1]

            for j in np.arange(0,len(x_mask)):
                if eddy_ID_mask[x_mask[j],y_mask[j]] == 'BG':
                    eddy_ID_mask[x_mask[j],y_mask[j]] = ID
                else:
                    eddy_ID_mask[x_mask[j],y_mask[j]] = eddy_ID_mask[x_mask[j],y_mask[j]] + ID
                    
                if eddy_age_mask[x_mask[j],y_mask[j]] == 'BG':
                    eddy_age_mask[x_mask[j],y_mask[j]] = age
                else:
                    eddy_age_mask[x_mask[j],y_mask[j]] = eddy_age_mask[x_mask[j],y_mask[j]] + age

    # Save the eddy mask as a numpy file
    #ID_filename = '%s%s_eddy_ID_mask_no_genesis.npy'%(ID_mask_dir,d)
    #age_filename = '%s%s_eddy_age_mask_no_genesis.npy'%(age_mask_dir,d)
    ID_filename = '%s%s_eddy_ID_mask_includes_genesis.npy'%(ID_mask_dir,d)
    age_filename = '%s%s_eddy_age_mask_includes_genesis.npy'%(age_mask_dir,d)
    np.save(ID_filename, eddy_ID_mask)
    np.save(age_filename, eddy_age_mask)

######################################### RUN THE SCRIPT ###########################################

for d in date_list: 
    print(d)
    
    # Get the eddy mask
    eddy_mask = np.full((len(chl_lat_array),len(chl_lon_array)),'BG',dtype='<U15') 
    project_eddy_bnds_on_mask(eddy_mask,full_eddy_dataset,d)
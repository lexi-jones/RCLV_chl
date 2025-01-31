# Calculate Nonlin Param & Update OceanEddies Dataset

# New values: Trans speed, rotational speed, nonlinearity param, area

# LJK
# Last edit: 09/11/24

import csv,os,time,sys
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime
sys.path.append('../RCLVatlas/')
sys.path.append('../RCLVatlas/RCLVatlas/')
from mainfunctions_for_RCLV_atlas import *
from config import *

############################################

# Open SLA eddy data
eddy_dir = project_output_dir + 'SLA_eddies/'

# Old datasets that were derived from '_v221104_QC'
old_anti_df = pd.read_csv(eddy_dir + 'anticyc_eddy_data_minage31_minsize12_tolno3_20000101_to_20191231_8day_subset_v2.csv')
old_cyc_df = pd.read_csv(eddy_dir + 'cyclonic_eddy_data_minage31_minsize12_tolno3_20000101_to_20191231_8day_subset_v2.csv')

# Get the ID's of the eddies used in the project
anti_IDs = np.unique(old_anti_df['Eddy ID'])
cyc_IDs = np.unique(old_cyc_df['Eddy ID'])
    
# Full dataset to redo
full_anti_df = pd.read_csv(eddy_dir + 'anticyc_eddy_data_minage31_minsize12_tolno3_20000101_to_20191231_v221104_QC.csv')
full_cyc_df = pd.read_csv(eddy_dir + 'cyclonic_eddy_data_minage31_minsize12_tolno3_20000101_to_20191231_v221104_QC.csv')

# Filter down to only the ID's used in the project
full_anti_df = full_anti_df[full_anti_df['Eddy ID'].isin(anti_IDs)] 
full_cyc_df = full_cyc_df[full_cyc_df['Eddy ID'].isin(cyc_IDs)] 

def reformat_eddy_data(df):
    ###### 1. Drop Untrustworthy Data Columns
    df = df.drop(columns=['Ls (km)','Area (km^2)','Amplitude (cm)','Translation Speed (m/s)']) 
    
    ###### 2. Calculate Translation Speed; I manually recalculated this because the values output from OceanEddies were incorrect
    print('STEP 2')
    df.insert(3,"Trans_Speed (cm/s-1)",[np.nan]*len(df)) # insert empty column of nans
    
    def speed_from_other_pt(eddy_df,other_ind,current_row):
        """
        Function for translation speed calc 
        """

        # Get dates of eddy instances
        Tcurrent = str(int(current_row['Date']))   
        Tother = str(int(eddy_df.loc[other_ind]['Date']))

        # Convert to datetime and find difference (days)
        delT = np.abs((datetime.datetime(int(Tcurrent[0:4]),int(Tcurrent[4:6]),int(Tcurrent[6:8])) - 
                           datetime.datetime(int(Tother[0:4]),int(Tother[4:6]),int(Tother[6:8]))).days)

        # Get distance between center lat/lons (km)
        dist = distance_from_lat_lon(eddy_df.loc[other_ind]['Center Lat'],
                                        eddy_df.loc[other_ind]['Center Lon'],
                                        current_row['Center Lat'],current_row['Center Lon'])

        # Calc translation speed (cm/s)
        speed = (dist*(10**5))/(delT*(24*60*60))
        return speed
    
    # One ID at a time, get the eddy's data and calculate its translation speed at each time step
    for i in np.unique(df['Eddy ID']):
        this_eddy_df = df.where(df['Eddy ID'] == i).dropna(how='all').dropna(axis=1)
        this_eddy_df_sorted = this_eddy_df.sort_values(by=['Date']) # sort by date

        count = 0
        for index, row in this_eddy_df_sorted.iterrows():
            if (count != 0) and (count != len(this_eddy_df_sorted)-1): # can't calculate the trans speed of the 1st or last contour
                # Calc average translation speed (cm/s)
                speedminus1 = speed_from_other_pt(this_eddy_df_sorted,index-1,row)
                speedplus1 = speed_from_other_pt(this_eddy_df_sorted,index+1,row)
                avg_speed = (speedminus1+speedplus1)/2
                df.at[index,"Trans_Speed (cm/s-1)"] = avg_speed    
            count += 1

    ###### 3. Re-crop back down to 8-day subset of the data
    # Get dates to iterate through from LAVD files
    LAVD_dir = project_output_dir + 'parcels_8day_overlap_32day_LAVD/'
    dates = []
    for filename in os.listdir(LAVD_dir):
        dates.append(filename[0:8])
    date_list = np.sort(np.unique(dates)).tolist() 
    df = df[df['Date'].isin(date_list)] # subset

    ###### 4. Crop out eddies that leave the domain
    print('STEP 4')
    drop_inds = [] # store the indeces of eddies that are not within the bounds
    for index, row in df.iterrows(): # need to do all indeces in full run
        x_bounds = [c for c in row[np.where(df.columns == 'Boundary Coords')[0][0]::2].values if not(np.isnan(c))]
        y_bounds = [c for c in row[(np.where(df.columns == 'Boundary Coords')[0][0]+1)::2].values if not(np.isnan(c))]
        if ((min(x_bounds) < 180) or (max(x_bounds) > 230) or (min(y_bounds) < 15) or (max(y_bounds) > 30)):
            drop_inds.append(index)
    df = df.drop(index=drop_inds)

    ###### 5. Calculate max rotational speed and nonlinearity param along contour
    print('STEP 5')
    df.insert(4,"Max_Rot_Speed (cm/s-1)",[np.nan]*len(df)) # insert empty column of nans for Trans speed
    df.insert(5,"Nonlin_param",[np.nan]*len(df)) # insert empty column of nans for Trans speed
    
    def get_max_speed_along_contour(eddy_df,row,CMEMS_data):
        x_bounds = [c for c in row[np.where(eddy_df.columns == 'Boundary Coords')[0][0]::2].values if not(np.isnan(c))]
        y_bounds = [c for c in row[(np.where(eddy_df.columns == 'Boundary Coords')[0][0]+1)::2].values if not(np.isnan(c))]

        # Looking for the maximum speed along the contour; manually doing this because I don't trust the output from OceanEddies
        speeds = []
        for i in np.arange(0,len(x_bounds)):
            lat_ind = np.where(CMEMS_data.latitude==y_bounds[i])[0][0]
            lon_ind = np.where(CMEMS_data.longitude==x_bounds[i])[0][0]
            speed = np.sqrt(float(CMEMS_data.ugos[0,lat_ind,lon_ind])**2 + float(CMEMS_data.vgos[0,lat_ind,lon_ind])**2)*(10**2) # convert to cm/s
            speeds.append(speed)
        return np.max(speeds)
    
    for d in date_list:
        CMEMS_data = xr.open_dataset(project_output_dir + 'CMEMS_data/dt_global_allsat_phy_l4_%s.nc'%(int(d)))

        # Iterate through the eddies on this day
        df_today = df.where(df['Date'] == int(d)).dropna(how='all')
        for index, row in df_today.iterrows():
            max_rot_speed = get_max_speed_along_contour(df_today,row,CMEMS_data)
            df.at[index,"Max_Rot_Speed (cm/s-1)"] = max_rot_speed
            if (np.isnan(row["Trans_Speed (cm/s-1)"])) or (row["Trans_Speed (cm/s-1)"]==0):
                pass
            else:
                df.at[index,"Nonlin_param"] = max_rot_speed/row["Trans_Speed (cm/s-1)"]

    ###### 6. Recalculate Area to Match RCLV Data
    print('STEP 6')
    areas = []
    for index, row in df.iterrows():
        lon_bounds = [c for c in row[np.where(df.columns == 'Boundary Coords')[0][0]::2].values if not(np.isnan(c))]
        lat_bounds = [c for c in row[(np.where(df.columns == 'Boundary Coords')[0][0]+1)::2].values if not(np.isnan(c))]
        areas.append(calc_area_of_stitched_bounds(lon_bounds,lat_bounds,CMEMS_data.longitude.values,CMEMS_data.latitude.values))
    df.insert(6,"Area (km2)",areas)
    return df

print('Reformatting ANTI df...')
start_time = time.time()
new_anti_df = reformat_eddy_data(full_anti_df)
new_anti_df.to_csv(eddy_dir + 'anticyc_eddy_data_minage31_minsize12_tolno3_20000101_to_20191231_8day_subset_v3.csv')
print("--- %s seconds ---" % (time.time() - start_time))

print('Reformatting CYC df...')
start_time = time.time()
new_cyc_df = reformat_eddy_data(full_cyc_df)
new_cyc_df.to_csv(eddy_dir + 'cyc_eddy_data_minage31_minsize12_tolno3_20000101_to_20191231_8day_subset_v3.csv')
print("--- %s seconds ---" % (time.time() - start_time))



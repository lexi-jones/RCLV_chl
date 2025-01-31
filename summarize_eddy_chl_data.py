# Use data generated from 'SLA_and_RCLV_eddy_chl_v6.0_with_eddy_provinces.py' to create summary files (e.g. seasonal data & annual data)

# LJK

#### Province key
# 0: SE
# 1: Northern lats
# 2: Upper cyclonic lee eddies
# 3: Lower anticyclonic lee eddies

import os
import numpy as np
from config import *


full_data_dir = project_output_dir + 'eddy_masks_and_chl_2000_to_2020/chl_feature_bins_full_domain/individual_dates/'
full_data_save_dir = project_output_dir + 'eddy_masks_and_chl_2000_to_2020/chl_feature_bins_full_domain/'

province_data_dir = project_output_dir + 'eddy_masks_and_chl_2000_to_2020/chl_feature_bins_eddy_provinces/individual_dates/'
province_data_save_dir = project_output_dir + 'eddy_masks_and_chl_2000_to_2020/chl_feature_bins_eddy_provinces/'

LAVD_dir = project_output_dir + 'parcels_8day_overlap_32day_LAVD/'
dates = []
for filename in os.listdir(LAVD_dir):
    dates.append(filename[0:8])    
full_date_list = np.sort(np.unique(dates)).tolist()
spring_date_list = [i for i in full_date_list if (i[4:6]=='03') or (i[4:6]=='04') or (i[4:6]=='05')]
summer_date_list = [i for i in full_date_list if (i[4:6]=='06') or (i[4:6]=='07') or (i[4:6]=='08')]
fall_date_list = [i for i in full_date_list if (i[4:6]=='09') or (i[4:6]=='10') or (i[4:6]=='11')]
winter_date_list = [i for i in full_date_list if (i[4:6]=='12') or (i[4:6]=='01') or (i[4:6]=='02')]

def add_elements_to_array(sub_array,main_array):
    for e in sub_array:
        main_array.append(e)
    return main_array

def generate_summary_dataset(season,province):
    
    all_bg_chl = []
    all_SSH_anti_chl = []
    all_coh_anti_chl = []
    all_SSH_coh_anti_chl = []
    all_SSH_cyc_chl = []
    all_coh_cyc_chl = []
    all_SSH_coh_cyc_chl = []
    
    if season == 'ALL':
        date_list = full_date_list
    elif season == 'SPRING':
        date_list = spring_date_list
    elif season == 'SUMMER':
        date_list = summer_date_list
    elif season == 'FALL':
        date_list = fall_date_list
    elif season == 'WINTER':
        date_list = winter_date_list
    
    if province == 'ALL':
        data_dir = full_data_dir
        save_dir = full_data_save_dir
    else:
        data_dir = province_data_dir
        save_dir = province_data_save_dir
    
    for d in date_list:
        bg_chl = np.load(data_dir + '%s_bg_chl_clim_anom_p%s.npy'%(d,province))
        all_bg_chl = add_elements_to_array(bg_chl,all_bg_chl)
            
        SSH_anti_chl = np.load(data_dir + '%s_SSH_anti_chl_clim_anom_p%s.npy'%(d,province))
        all_SSH_anti_chl = add_elements_to_array(SSH_anti_chl,all_SSH_anti_chl)
        
        coh_anti_chl = np.load(data_dir + '%s_coh_anti_chl_clim_anom_p%s.npy'%(d,province))
        all_coh_anti_chl = add_elements_to_array(coh_anti_chl,all_coh_anti_chl)
        
        SSH_coh_anti_chl = np.load(data_dir + '%s_SSH_coh_anti_chl_clim_anom_p%s.npy'%(d,province))
        all_SSH_coh_anti_chl = add_elements_to_array(SSH_coh_anti_chl,all_SSH_coh_anti_chl)
        
        SSH_cyc_chl = np.load(data_dir + '%s_SSH_cyc_chl_clim_anom_p%s.npy'%(d,province))
        all_SSH_cyc_chl = add_elements_to_array(SSH_cyc_chl,all_SSH_cyc_chl)
        
        coh_cyc_chl = np.load(data_dir + '%s_coh_cyc_chl_clim_anom_p%s.npy'%(d,province))
        all_coh_cyc_chl = add_elements_to_array(coh_cyc_chl,all_coh_cyc_chl)
        
        SSH_coh_cyc_chl = np.load(data_dir + '%s_SSH_coh_cyc_chl_clim_anom_p%s.npy'%(d,province))
        all_SSH_coh_cyc_chl = add_elements_to_array(SSH_coh_cyc_chl,all_SSH_coh_cyc_chl)

    np.save(save_dir + '%s_bg_chl_clim_anom_p%s.npy'%(season,province), all_bg_chl)
    np.save(save_dir + '%s_SSH_anti_chl_clim_anom_p%s.npy'%(season,province), all_SSH_anti_chl)
    np.save(save_dir + '%s_coh_anti_chl_clim_anom_p%s.npy'%(season,province), all_coh_anti_chl)
    np.save(save_dir + '%s_SSH_coh_anti_chl_clim_anom_p%s.npy'%(season,province), all_SSH_coh_anti_chl)
    np.save(save_dir + '%s_SSH_cyc_chl_clim_anom_p%s.npy'%(season,province), all_SSH_cyc_chl)
    np.save(save_dir + '%s_coh_cyc_chl_clim_anom_p%s.npy'%(season,province), all_coh_cyc_chl)
    np.save(save_dir + '%s_SSH_coh_cyc_chl_clim_anom_p%s.npy'%(season,province), all_SSH_coh_cyc_chl)
    
# Full dataset

print('Summarizing full dataset...')
generate_summary_dataset('ALL','ALL')
generate_summary_dataset('SPRING','ALL')
generate_summary_dataset('SUMMER','ALL')
generate_summary_dataset('FALL','ALL')
generate_summary_dataset('WINTER','ALL')

# Each of the 4 provinces
for p in np.arange(0,4):
    print('Summarizing %s province dataset...'%(p))
    generate_summary_dataset('ALL',p)
    generate_summary_dataset('SPRING',p)
    generate_summary_dataset('SUMMER',p)
    generate_summary_dataset('FALL',p)
    generate_summary_dataset('WINTER',p)
# Bootstrapping Error Analysis

# LJK
# Date created: 01/06/25
# Last edited: 01/15/25

# 1. Get the background data 
# 2. Get the eddy type data
# 3. Re-sample the eddy type with the same # of samples as the original sample size
# 4. Compute f(dcclim) from 1 to 99% background
# 5. Repeat 1000x 
# 6. Get the 95% confidence interval at each dcclim value (200 values for full dataset, 100 values for seasonal/regional subsets)

from config import *
import csv
import numpy as np

######################## FULL DOMAIN DATA ########################

def bootstrap_full_domain():
    data_dir = project_output_dir + 'eddy_masks_and_chl_2000_to_2020/chl_feature_bins_full_domain/'
    bg_chl = np.load(data_dir + 'ALL_bg_chl_clim_anom_pALL.npy')
    
    # Background data: Get 1 to 99% quantiles of background and the density plot
    num_bins = 200
    bg_quant_min,bg_quant_max = np.quantile(bg_chl,0.01),np.quantile(bg_chl,0.99)
    bg_counts,bg_bins = np.histogram(bg_chl,bins=num_bins,density=True,range=(bg_quant_min,bg_quant_max))
    
    SLA_only_anti_chl = np.load(data_dir + 'ALL_SSH_anti_chl_clim_anom_pALL.npy')
    RCLV_only_anti_chl = np.load(data_dir + 'ALL_coh_anti_chl_clim_anom_pALL.npy')
    overlap_anti_chl = np.load(data_dir + 'ALL_SSH_coh_anti_chl_clim_anom_pALL.npy')
    SLA_only_cyc_chl = np.load(data_dir + 'ALL_SSH_cyc_chl_clim_anom_pALL.npy')
    RCLV_only_cyc_chl = np.load(data_dir + 'ALL_coh_cyc_chl_clim_anom_pALL.npy')
    overlap_cyc_chl = np.load(data_dir + 'ALL_SSH_coh_cyc_chl_clim_anom_pALL.npy')

    # Combine the type only and overlaps together
    SLA_ANTI_CHL = np.concatenate((SLA_only_anti_chl,overlap_anti_chl))
    RCLV_ANTI_CHL = np.concatenate((RCLV_only_anti_chl,overlap_anti_chl))
    SLA_CYC_CHL = np.concatenate((SLA_only_cyc_chl,overlap_cyc_chl))
    RCLV_CYC_CHL = np.concatenate((RCLV_only_cyc_chl,overlap_cyc_chl))

    eddy_datasets = [SLA_only_anti_chl,SLA_ANTI_CHL,RCLV_ANTI_CHL,SLA_only_cyc_chl,SLA_CYC_CHL,RCLV_CYC_CHL]
    eddy_labels = ['Anti_SLA_Excluding_RCLVs','Anti_SLA','Anti_RCLV','Cyc_SLA_Excluding_RCLVs','Cyc_SLA','Cyc_RCLV']

    output_dir = project_output_dir + 'eddy_masks_and_chl_2000_to_2020/chl_feature_bins_full_domain/bootstrap_error/' 
    
    # Compute confidence intervals with bootstrapping
    for d in np.arange(0,len(eddy_datasets)): 
        print(eddy_labels[d])
        
        bootstrap_yields = []
        resample_num = 1000 # resample 100 times
        for i in np.arange(0,resample_num): 
            # Resample the dataset 
            resample = np.random.choice(eddy_datasets[d], size=len(eddy_datasets[d]))
            
            # Get f(dcclim)
            counts,bins = np.histogram(resample,bins=num_bins,density=True,range=(bg_quant_min,bg_quant_max))
            y_vals = (counts-bg_counts)/bg_counts
            bootstrap_yields.append(y_vals)
        
        # Save min & max of bootstrap yields  
        np.save(output_dir+'fdclim_bootstrap_%s_resample_2.5p_%s_ALL.npy'%(resample_num,eddy_labels[d]),np.quantile(bootstrap_yields,0.025,axis=0))
        np.save(output_dir+'fdclim_bootstrap_%s_resample_97.5p_%s_ALL.npy'%(resample_num,eddy_labels[d]),np.quantile(bootstrap_yields,0.975,axis=0))

bootstrap_full_domain()

######################## PROVINCE DATA ########################

def bootstrap_provinces():
    """
    #### Province key #### 0: SE, 1: Northern lats, 2: Upper cyclonic lee eddies, 3: Lower anticyclonic lee eddies
    """
    data_dir = project_output_dir + 'eddy_masks_and_chl_2000_to_2020/chl_feature_bins_eddy_provinces/'
    output_dir = project_output_dir + 'eddy_masks_and_chl_2000_to_2020/chl_feature_bins_eddy_provinces/bootstrap_error/' 
    
    for season in ['WINTER','SPRING','SUMMER','FALL']:  
        print(season)
        for province in ['SE','N','Lee']: #### Province key #### 0: SE, 1: Northern lats, 2: Upper cyclonic lee eddies, 3: Lower anticyclonic lee eddies
            print(province)
            
            if province == 'Lee': # combine upper and lower lee eddy regions into one
                bg_chl = np.concatenate((np.load(data_dir + '%s_bg_chl_clim_anom_p%s.npy'%(season,2)),np.load(data_dir + '%s_bg_chl_clim_anom_p%s.npy'%(season,3))))
                SLA_only_anti_chl = np.concatenate((np.load(data_dir + '%s_SSH_anti_chl_clim_anom_p%s.npy'%(season,2)),np.load(data_dir + '%s_SSH_anti_chl_clim_anom_p%s.npy'%(season,3))))
                RCLV_only_anti_chl = np.concatenate((np.load(data_dir + '%s_coh_anti_chl_clim_anom_p%s.npy'%(season,2)),np.load(data_dir + '%s_coh_anti_chl_clim_anom_p%s.npy'%(season,3))))
                overlap_anti_chl = np.concatenate((np.load(data_dir + '%s_SSH_coh_anti_chl_clim_anom_p%s.npy'%(season,2)),np.load(data_dir + '%s_SSH_coh_anti_chl_clim_anom_p%s.npy'%(season,3))))
                SLA_only_cyc_chl = np.concatenate((np.load(data_dir + '%s_SSH_cyc_chl_clim_anom_p%s.npy'%(season,2)),np.load(data_dir + '%s_SSH_cyc_chl_clim_anom_p%s.npy'%(season,3))))
                RCLV_only_cyc_chl = np.concatenate((np.load(data_dir + '%s_coh_cyc_chl_clim_anom_p%s.npy'%(season,2)),np.load(data_dir + '%s_coh_cyc_chl_clim_anom_p%s.npy'%(season,3))))
                overlap_cyc_chl = np.concatenate((np.load(data_dir + '%s_SSH_coh_cyc_chl_clim_anom_p%s.npy'%(season,2)),np.load(data_dir + '%s_SSH_coh_cyc_chl_clim_anom_p%s.npy'%(season,3))))

            else:
                if province == 'SE':
                    p_key = 0
                elif province == 'N':
                    p_key = 1

                bg_chl = np.load(data_dir + '%s_bg_chl_clim_anom_p%s.npy'%(season,p_key))
                SLA_only_anti_chl = np.load(data_dir + '%s_SSH_anti_chl_clim_anom_p%s.npy'%(season,p_key))
                RCLV_only_anti_chl = np.load(data_dir + '%s_coh_anti_chl_clim_anom_p%s.npy'%(season,p_key))
                overlap_anti_chl = np.load(data_dir + '%s_SSH_coh_anti_chl_clim_anom_p%s.npy'%(season,p_key))
                SLA_only_cyc_chl = np.load(data_dir + '%s_SSH_cyc_chl_clim_anom_p%s.npy'%(season,p_key))
                RCLV_only_cyc_chl = np.load(data_dir + '%s_coh_cyc_chl_clim_anom_p%s.npy'%(season,p_key))
                overlap_cyc_chl = np.load(data_dir + '%s_SSH_coh_cyc_chl_clim_anom_p%s.npy'%(season,p_key))

            # Combine the "type only" and overlaps together
            SLA_ANTI_CHL = np.concatenate((SLA_only_anti_chl,overlap_anti_chl))
            RCLV_ANTI_CHL = np.concatenate((RCLV_only_anti_chl,overlap_anti_chl))
            SLA_CYC_CHL = np.concatenate((SLA_only_cyc_chl,overlap_cyc_chl))
            RCLV_CYC_CHL = np.concatenate((RCLV_only_cyc_chl,overlap_cyc_chl))

            # Background data: Get 1 to 99% quantiles of background and the density plot
            num_bins = 100
            bg_quant_min,bg_quant_max = np.quantile(bg_chl,0.01),np.quantile(bg_chl,0.99)
            bg_counts,bg_bins = np.histogram(bg_chl,bins=num_bins,density=True,range=(bg_quant_min,bg_quant_max))
            
            eddy_datasets = [SLA_only_anti_chl,SLA_ANTI_CHL,RCLV_ANTI_CHL,SLA_only_cyc_chl,SLA_CYC_CHL,RCLV_CYC_CHL]
            eddy_labels = ['Anti_SLA_Excluding RCLVs','Anti_SLA','Anti_RCLV','Cyc_SLA_Excluding_RCLVs','Cyc_SLA','Cyc_RCLV']

            # Compute confidence intervals with bootstrapping
            for d in np.arange(0,len(eddy_datasets)): 
                bootstrap_yields = []
                resample_num = 1000 # resample 100 times
                for i in np.arange(0,resample_num): 
                    # Resample the dataset 
                    resample = np.random.choice(eddy_datasets[d], size=len(eddy_datasets[d]))

                    # Get f(dcclim)
                    counts,bins = np.histogram(resample,bins=num_bins,density=True,range=(bg_quant_min,bg_quant_max))
                    y_vals = (counts-bg_counts)/bg_counts
                    bootstrap_yields.append(y_vals)

                # Save min & max of bootstrap yields  
                np.save(output_dir+'fdclim_bootstrap_%s_resample_2.5p_%s_%s_%s.npy'%(resample_num,eddy_labels[d],season,province),np.quantile(bootstrap_yields,0.025,axis=0))
                np.save(output_dir+'fdclim_bootstrap_%s_resample_97.5p_%s_%s_%s.npy'%(resample_num,eddy_labels[d],season,province),np.quantile(bootstrap_yields,0.975,axis=0))
                

bootstrap_provinces()

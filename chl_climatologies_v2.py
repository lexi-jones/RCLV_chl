# LJK 
# Calculate the monthly avg chl from OC-CCI data

import numpy as np
import xarray as xr
import os
from config import *


month_list = ['01','02','03','04','05','06','07','08','09','10','11','12']
months = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

chl_dir = project_output_dir + 'OCCCI_8day_v6.0_bounding_box/' #where data lives
clim_dir = project_output_dir + 'OCCCI_8day_v6.0_climatologies/' #where to save climatologies

chl_example = xr.open_dataset(chl_dir + 'ESACCI-OC-L3S-CHLOR_A-MERGED-8D_DAILY_4km_GEO_PML_OCx-20140805-fv6.0_bounding_box.nc')

for i in np.arange(0,len(months)):
    month_filenames = [f for f in os.listdir(chl_dir) if (f[-26:-24] == month_list[i])]

    month_chl = []
    for filename in month_filenames:
        temp_chl = xr.open_dataset(chl_dir + filename)
        month_chl.append(temp_chl.chlor_a[0])
    month_chl_avg = np.nanmean(month_chl,axis=0)

    np.save(clim_dir + 'ESACCI-OC-L3S-CHLOR_A-MERGED-8D_DAILY_4km_GEO_PML_OCx-fv6.0_bounding_box_%s_AVG.npy'%(months[i]), month_chl_avg)

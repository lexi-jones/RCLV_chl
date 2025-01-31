# RCLV Chlorophyll Analysis

Scripts accompany study by Jones-Kellett and Follows (https://doi.org/10.5194/egusphere-2024-3211).
RCLV tracking scripts are at github.com/lexi-jones/RCLVatlas:
    Jones-Kellett, A. (2023). RCLVatlas (Version 1.0.0) [Computer software]. https://doi.org/10.5281/zenodo.7702978

## Pipeline

1. Obtain input datasets:
	- OC-CCI chl-a (https://www.oceancolour.org/); 8-day avgs used for the paper, wrapped to center on the Pacific and cropped to domain
	- RCLVatlas v2: available for download at https://simonscmap.com/catalog/datasets/RCLV_atlas_version2 or 
	  https://zenodo.org/records/10849221. See https://github.com/lexi-jones/RCLVatlas for the source code.
	- SLA eddy atlas:  OceanEddies MATLAB software to detect and track Eulerian SLA eddy contours was obtained from 
          https://github.com/ifrenger/OceanEddies (last access: 13 October 2021). Eddies were tracked in  CMEMS Level 4, 1/4◦ SLA and geostrophic
	  velocity gridded global ocean dataset, Version 008_047 (https://data.marine.copernicus.eu/product/SEALEVEL_GLO_PHY_L4_MY_008_047/)

2. Calculate climatological chl anomalies: `chl_climatologies_v2.py`; corresponding plots: `annual_chl.ipynb`

3. Project eddy bounds onto chl-a grid to create 8-day eddy masks: `project_SSH_and_RCLV_atlas_on_chl_grid_v3.py`

4. Seperate chl data by eddy type using grid masks: `SLA_and_RCLV_eddy_chl_v6.0_with_eddy_provinces.py`

5. Seperate chl data into seasons: `summarize_eddy_chl_data.py`

6. Bootstrapping analysis to compute confidence intervals: `bootstrapping_error_analysis.py`

7. Probability density distribution plots by eddy type: `RCLV_SLA_eddy_chl_PDFs.ipynb`

8. Compute local chl anomaly: `RCLV_inside_outside_chl.py`; long-lived eddy local anom plots in `long_lived_RCLV_seasonal_province_age_OS.ipynb`

9. Nonlinearity parameter: `calc_nonlin_param.py` and figures in `nonlinearity_param.ipynb`

10. Other schematics in paper: `RCLV_SLA_eddy_chl_schematics.ipynb`

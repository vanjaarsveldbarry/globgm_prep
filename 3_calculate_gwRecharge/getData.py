from pathlib import Path
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree
import numpy as np
import subprocess
import rioxarray
import os
import shutil

dataFolder = Path('/scratch/depfg/7006713/globgm_input/_tempGLOBGM/gwRecharge/data')
dataFolder.mkdir(parents=True, exist_ok=True)
outputFolder = Path('/scratch/depfg/7006713/globgm_input/_tempGLOBGM/gwRecharge/out')
outputFolder.mkdir(parents=True, exist_ok=True)
tempFolder = Path('/scratch/depfg/7006713/globgm_input/_tempGLOBGM/gwRecharge/_temp')
tempFolder.mkdir(parents=True, exist_ok=True)



#From https://hess.copernicus.org/articles/22/2689/2018/#section2&gid=1&pid=1 Table 1 we try and reproduce and updat ethe date
#Plus some from this paper: https://www.researchgate.net/profile/Mario-Schirmer/publication/338978931_A_global-scale_dataset_of_direct_natural_groundwater_recharge_rates_A_review_of_variables_processes_and_relationships/links/64197d9c92cfd54f8418b40b/A-global-scale-dataset-of-direct-natural-groundwater-recharge-rates-A-review-of-variables-processes-and-relationships.pdf
#Precipitation                : CHELSA                           : 
#Temperature                  : CHELSA                           : 
#Potential evapotranspiration : Hargraves  pyET - from CHELSA    : 
#Number of rainy day          : NOT AVAILABLE                    : 
#Slope                        : pcr-globwb                       : DONE 
#ksat                         : pcrglobwb                        : DONE 
# soil wtare storage capacity : pcrglobwb                        : 
#Excess water                 : P - ET from chelsa                  : 
#Aridity index                : CHELSA                           : 
#Clay fraction                :                                  : DONE
#Sand fraction                :                                  : DONE 
#Silt fraction                :                                  : DONE 
#Bulk density                 :                                  : DONE
#Land use cover               :                                  : 
#Elevation                    : pcrglobwb                        : DONE
#Reccession Coeef             : pcrglobwb                        : DONE

def download_file_wget(url, outFile):
    try:
        subprocess.run(['wget', '-O', str(outFile), url], check=True)
        print('File downloaded successfully.')
    except subprocess.CalledProcessError as e:
        print(f'Failed to download file. Error: {e}')

latMin, latMax = -35.0, -30.0
lonMin, lonMax = 16.0, 33.0

gridFile='/scratch/depfg/sutan101/data/pcrglobwb_input_arise/develop/global_30sec/landSurface/topography/merit_dem_processed/version_2021-02-XX/maps_covered_with_zero/dem_average_topography_parameters_30sec_february_2021_global_covered_with_zero.nc'
_grid_lat = xr.open_dataset(gridFile).sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute().lat.values
_grid_lon = xr.open_dataset(gridFile).sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute().lon.values

elevationFile='/scratch/depfg/sutan101/data/pcrglobwb_input_arise/develop/global_30sec/landSurface/topography/merit_dem_processed/version_2021-02-XX/maps_covered_with_zero/dem_average_topography_parameters_30sec_february_2021_global_covered_with_zero.nc'
elev_30sec = xr.open_dataset(elevationFile)[['dem_average']].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
elev_30sec = elev_30sec.rename({'dem_average': 'elevation'})
elev_30sec.attrs = {}
elev_30sec.to_zarr(dataFolder / 'input.zarr', mode='w')

aiFile='/scratch/depfg/7006713/data/aridity_index/aridity_index.nc'
ai_30sec = xr.open_dataset(aiFile)[['ai']].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
ai_30sec.attrs = {}
ai_30sec.to_zarr(dataFolder / 'input.zarr', mode='a')


slopeFile='/scratch/depfg/sutan101/data/pcrglobwb_input_arise/develop/global_30sec/landSurface/topography/merit_dem_processed/version_2021-02-XX/maps_covered_with_zero/tanslope_topography_parameters_30sec_february_2021_global_covered_with_zero.nc'
slope_30sec = xr.open_dataset(slopeFile)[['tanslope']].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
slope_30sec = slope_30sec.rename({'tanslope': 'slope'})
slope_30sec.attrs = {}
slope_30sec.to_zarr(dataFolder / 'input.zarr', mode='a')

url='https://geo.public.data.uu.nl/vault-globgm/research-globgm%5B1669042611%5D/original/input/version_1.0/k_conductivity_aquifer_filled_30sec.nc'
ksatFile=tempFolder / 'temp_ksat.nc'
download_file_wget(url, ksatFile)
ksat_30sec = xr.open_dataset(ksatFile)[['k_conductivity_aquifer_filled_map']].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
ksat_30sec = ksat_30sec.rename({'k_conductivity_aquifer_filled_map': 'k_sat'})
ksat_30sec.attrs = {}
ksat_30sec.to_zarr(dataFolder / 'input.zarr', mode='a')

url='https://geo.public.data.uu.nl/vault-globgm/research-globgm%5B1669042611%5D/original/input/version_1.0/recession_coefficient_30sec.nc'
recFile=tempFolder / 'temp_rec_coeff.nc'
download_file_wget(url, recFile)
rec_coef_30sec = xr.open_dataset(recFile)[['recession_coefficient_30sec_map']].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
rec_coef_30sec = rec_coef_30sec.rename({'recession_coefficient_30sec_map': 'rec_coef'})
rec_coef_30sec.attrs = {}
rec_coef_30sec.to_zarr(dataFolder / 'input.zarr', mode='a')

test = xr.open_zarr(dataFolder / 'input.zarr')
print(test)







########################################################################################################
# # Get the clay percentage for fromSoil250m 
# def getClayData(url):
#     for layer in ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']:
#         sub_url= f'{url}/clay_{layer}cm_mean_1000.tif'
#         clay_perFile=tempFolder / f'temp_clay_perc_{layer}cm.tif'
#         download_file_wget(sub_url, clay_perFile)
#         ds = rioxarray.open_rasterio(clay_perFile, parse_coordinates=True)
#         ds = ds.where(ds != -32768, np.nan)
#         ds = ds / 10 
#         if ds.rio.crs != 'EPSG:4326':
#             ds = ds.rio.reproject('EPSG:4326')
#         ds.to_zarr(tempFolder / f'clay_perc_{layer}cm.zarr', mode='w')
#         os.remove(clay_perFile)
#     for layer in ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']:
#         ds = xr.open_zarr(tempFolder / f'clay_perc_{layer}cm.zarr').rename({'__xarray_dataarray_variable__': 'clay_perc', 
#                                                                             'band': 'layer', 
#                                                                             'x': 'lon',
#                                                                             'y': 'lat'})
#         ds = ds.drop_vars('spatial_ref')
#         if layer == '0-5':
#             clay_ds = ds
#         else:
#             clay_ds = xr.concat([clay_ds, ds], dim='layer')
#     clay_ds = clay_ds.mean(dim='layer').compute()
#     clay_ds.to_zarr(tempFolder / 'temp_clay_perc.zarr', mode='w')
    
#     for layer in ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']:
#         shutil.rmtree(tempFolder / f'clay_perc_{layer}cm.zarr')
        
# url='https://files.isric.org/soilgrids/latest/data_aggregated/1000m/clay'
# getClayData(url)

# clay_perc = xr.open_zarr(tempFolder / 'temp_clay_perc.zarr')['clay_perc']
# clay_perc = clay_perc.sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
# clay_perc = clay_perc.reindex(lat=_grid_lat, lon=_grid_lon, method='nearest')
# clay_perc.to_netcdf('/scratch/depfg/7006713/globgm_input/_tempGLOBGM/gwRecharge/_temp/test.nc')
# clay_perc.to_zarr(dataFolder / 'input.zarr', mode='a')

# # Get the sand percentage for fromSoil250m 
# def getSandData(url):
#     for layer in ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']:
#         sub_url= f'{url}/sand_{layer}cm_mean_1000.tif'
#         sand_perFile=tempFolder / f'temp_sand_perc_{layer}cm.tif'
#         download_file_wget(sub_url, sand_perFile)
#         ds = rioxarray.open_rasterio(sand_perFile, parse_coordinates=True)
#         ds = ds.where(ds != -32768, np.nan)
#         ds = ds / 10 
#         if ds.rio.crs != 'EPSG:4326':
#             ds = ds.rio.reproject('EPSG:4326')
#         ds.to_zarr(tempFolder / f'sand_perc_{layer}cm.zarr', mode='w')
#         os.remove(sand_perFile)
#     for layer in ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']:
#         ds = xr.open_zarr(tempFolder / f'sand_perc_{layer}cm.zarr').rename({'__xarray_dataarray_variable__': 'sand_perc', 
#                                                                             'band': 'layer', 
#                                                                             'x': 'lon',
#                                                                             'y': 'lat'})
#         ds = ds.drop_vars('spatial_ref')
#         if layer == '0-5':
#             sand_ds = ds
#         else:
#             sand_ds = xr.concat([sand_ds, ds], dim='layer')
#     sand_ds = sand_ds.mean(dim='layer').compute()
#     sand_ds.to_zarr(tempFolder / 'temp_sand_perc.zarr', mode='w')
    
#     for layer in ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']:
#         shutil.rmtree(tempFolder / f'sand_perc_{layer}cm.zarr')
        
# url='https://files.isric.org/soilgrids/latest/data_aggregated/1000m/sand'
# getSandData(url)
# sand_perc = xr.open_zarr(tempFolder / 'temp_sand_perc.zarr')['sand_perc']
# sand_perc = sand_perc.sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
# sand_perc = sand_perc.reindex(lat=_grid_lat, lon=_grid_lon, method='nearest')
# sand_perc.to_zarr(dataFolder / 'input.zarr', mode='a')

# #Get the silt percentage for fromSoil250m 
# def getSiltData(url):
#     for layer in ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']:
#         sub_url= f'{url}/silt_{layer}cm_mean_1000.tif'
#         silt_perFile=tempFolder / f'temp_silt_perc_{layer}cm.tif'
#         download_file_wget(sub_url, silt_perFile)
#         ds = rioxarray.open_rasterio(silt_perFile, parse_coordinates=True)
#         ds = ds.where(ds != -32768, np.nan)
#         ds = ds / 10 
#         if ds.rio.crs != 'EPSG:4326':
#             ds = ds.rio.reproject('EPSG:4326')
#         ds.to_zarr(tempFolder / f'silt_perc_{layer}cm.zarr', mode='w')
#         os.remove(silt_perFile)
#     for layer in ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']:
#         ds = xr.open_zarr(tempFolder / f'silt_perc_{layer}cm.zarr').rename({'__xarray_dataarray_variable__': 'silt_perc', 
#                                                                             'band': 'layer', 
#                                                                             'x': 'lon',
#                                                                             'y': 'lat'})
#         ds = ds.drop_vars('spatial_ref')
#         if layer == '0-5':
#             silt_ds = ds
#         else:
#             silt_ds = xr.concat([silt_ds, ds], dim='layer')
#     silt_ds = silt_ds.mean(dim='layer').compute()
#     silt_ds.to_zarr(tempFolder / 'temp_silt_perc.zarr', mode='w')
    
#     for layer in ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']:
#         shutil.rmtree(tempFolder / f'silt_perc_{layer}cm.zarr')
        
# url='https://files.isric.org/soilgrids/latest/data_aggregated/1000m/silt'
# getSiltData(url)
# silt_perc = xr.open_zarr(tempFolder / 'temp_silt_perc.zarr')['silt_perc']
# silt_perc = silt_perc.sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
# silt_perc = silt_perc.reindex(lat=_grid_lat, lon=_grid_lon, method='nearest')
# silt_perc.to_zarr(dataFolder / 'input.zarr', mode='a')










############################################################################
# simRechargeFile='/projects/0/einf4705/_data/cmip6_input/gswp3-w5e5/historical_natural/pcrglobwb_cmip6-isimip3-gswp3-w5e5_image-aqueduct_historical-natural_gwRecharge_global_monthly-total_1960_2019_basetier1.nc'
# obs_data = pd.read_csv(dataFolder / 'recharge.csv', sep=';')
# obs_data = obs_data.replace(',', '.', regex=True)
# obs_data = obs_data.apply(pd.to_numeric, errors='ignore')
# obs_data.rename(columns={'Latitude': 'lat','Longitude': 'lon','Groundwater recharge': 'obs_recharge'}, inplace=True)
# obs_data = obs_data[(obs_data['lat'] >= latMin) & (obs_data['lat'] <= latMax) & (obs_data['lon'] >= lonMin) & (obs_data['lon'] <= lonMax)]
# obs_data['obs_recharge'] = obs_data['obs_recharge'] / 1000

# latIndicies = obs_data['lat'].values
# lonIndicies = obs_data['lon'].values

# gwRech_30sec = xr.open_dataset(simRechargeFile)[['groundwater_recharge']].sel(latitude=slice(latMax, latMin), longitude=slice(lonMin, lonMax)).rename({'latitude': 'lat', 'longitude': 'lon'})
# gwRech_30sec = gwRech_30sec.resample(time='YE').sum(skipna=False)
# gwRech_30sec = gwRech_30sec.mean(dim='time', skipna=False).compute()
# data_ds = gwRech_30sec.sel(lat=xr.DataArray(latIndicies, dims="z"), lon=xr.DataArray(lonIndicies, dims="z"), method='nearest').rename({'groundwater_recharge' : 'sim_recharge'})

# ai_30sec = xr.open_dataset(dataFolder / 'ai_v3_year.nc')['Band1'].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
# ai_30sec = ai_30sec.sel(lat=xr.DataArray(latIndicies, dims="z"), lon=xr.DataArray(lonIndicies, dims="z"), method='nearest')
# ai_30sec = ai_30sec * 0.0001
# data_ds = data_ds.assign(AI=ai_30sec)

# elev_30sec = xr.open_dataset(dataFolder / 'dem_average_topography_parameters_30sec_february_2021_global_covered_with_zero.nc')['dem_average'].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
# elev_30sec = elev_30sec.sel(lat=xr.DataArray(latIndicies, dims="z"), lon=xr.DataArray(lonIndicies, dims="z"), method='nearest')
# data_ds = data_ds.assign(elev=elev_30sec)

# reces_coef_30sec = xr.open_dataset(dataFolder/ 'recession_coefficient_30sec.nc')['recession_coefficient_30sec_map'].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
# reces_coef_30sec = reces_coef_30sec.sel(lat=xr.DataArray(latIndicies, dims="z"), lon=xr.DataArray(lonIndicies, dims="z"), method='nearest')
# data_ds = data_ds.assign(J=reces_coef_30sec)
# data_ds = data_ds.to_dataframe()
# obs_coords = obs_data[['lat', 'lon']].values
# sim_coords = data_ds[['lat', 'lon']].values

# tree = cKDTree(sim_coords)
# distances, indices = tree.query(obs_coords)
# final_df = data_ds.iloc[indices]
# final_df['obs_recharge'] = obs_data['obs_recharge'].reset_index(drop=True)
# final_df.to_parquet(outputFolder / 'data_train.parquet')
# final_df.to_parquet(outputFolder / 'data_train.parquet')
# print('trainData created')
# # # # # ##Create prediction dataset
# gird= xr.open_dataset('/projects/0/einf4705/gwRecharge_correction/data/recession_coefficient_30sec.nc')
# data_ds = xr.open_dataset(simRechargeFile)[['groundwater_recharge']].sel(latitude=slice(latMax, latMin), longitude=slice(lonMin, lonMax)).rename({'latitude': 'lat', 'longitude': 'lon'})
# data_ds = data_ds.resample(time='YE').sum(skipna=False)
# data_ds = data_ds.mean(dim='time', skipna=False).compute()
# data_ds = data_ds.rename({'groundwater_recharge' : 'sim_recharge'})
# data_ds = data_ds.reindex_like(gird, method='nearest')

# elev_30sec = xr.open_dataset(dataFolder / 'dem_average_topography_parameters_30sec_february_2021_global_covered_with_zero.nc')['dem_average'].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
# elev_30sec = elev_30sec.reindex_like(data_ds, method='nearest')
# data_ds = data_ds.assign(elev=elev_30sec)

# ai_30sec = xr.open_dataset(dataFolder / 'ai_v3_year.nc')['Band1'].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
# ai_30sec = ai_30sec.reindex_like(data_ds, method='nearest')	
# ai_30sec = ai_30sec * 0.0001
# data_ds = data_ds.assign(AI=ai_30sec)

# reces_coef_30sec = xr.open_dataset(dataFolder  / 'recession_coefficient_30sec.nc')['recession_coefficient_30sec_map'].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
# reces_coef_30sec = reces_coef_30sec.reindex_like(data_ds, method='nearest')	
# data_ds = data_ds.assign(J=reces_coef_30sec)
# data_ds = data_ds.to_dataframe().reset_index()
# data_ds = data_ds.dropna()
# data_ds.to_parquet(outputFolder / 'data_predict.parquet')
# print('Predict created')

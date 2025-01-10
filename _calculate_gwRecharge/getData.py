from pathlib import Path
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree
import numpy as np
import subprocess
import rioxarray
import os
import shutil

rootFolder = Path('/scratch/depfg/7006713/gwRecharge_corrected')
tempFolder = rootFolder / '_temp'
tempFolder.mkdir(parents=True, exist_ok=True)
dataFolder = rootFolder / 'data'
dataFolder.mkdir(parents=True, exist_ok=True)


#From https://hess.copernicus.org/articles/22/2689/2018/#section2&gid=1&pid=1 Table 1 we try and reproduce and updat ethe date
#Plus some from this paper: https://www.researchgate.net/profile/Mario-Schirmer/publication/338978931_A_global-scale_dataset_of_direct_natural_groundwater_recharge_rates_A_review_of_variables_processes_and_relationships/links/64197d9c92cfd54f8418b40b/A-global-scale-dataset-of-direct-natural-groundwater-recharge-rates-A-review-of-variables-processes-and-relationships.pdf
#Precipitation                : CHELSA                           : Done
#Potential evapotranspiration : CHELSA                           : Done
#Slope                        : pcr-globwb                       : DONE 
#ksat                         : pcrglobwb                        : DONE 
#Aridity index                : NOT CHELSE BUT OTEHR HIGH RES    : DONE
#Clay fraction                : soil250m                         : DONE
#Sand fraction                : soil250m                         : DONE 
#Silt fraction                : soil250m                         : DONE 
#Bulk density                 : soil250m                         : DONE
#Land use cover               :                                  : 
#Elevation                    : pcrglobwb                        : DONE
#Reccession Coeef             : pcrglobwb                        : DONE
#30sec Recharge from chapter 1: pcrglobwb                        : DONE

def download_file_wget(url, outFile):
    try:
        subprocess.run(['wget', '-O', str(outFile), url], check=True)
        print('File downloaded successfully.')
    except subprocess.CalledProcessError as e:
        print(f'Failed to download file. Error: {e}')

latMin, latMax = -90.0, 90.0
lonMin, lonMax = -180.0, 180.0
# latMin, latMax = -44.0, -10.0
# lonMin, lonMax = 112.0, 154.0

gridFile='/scratch/depfg/sutan101/data/pcrglobwb_input_arise/develop/global_30sec/routing/surface_water_bodies/version_2020-05-XX/lddsound_30sec_version_202005XX_correct_lat.nc'
_grid_lat = xr.open_dataset(gridFile).sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute().lat.values
_grid_lon = xr.open_dataset(gridFile).sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute().lon.values

mask = xr.open_dataset(gridFile)['Band1']
mask = mask.sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
mask = mask.rename('mask')
mask = xr.where(mask != 255, 1, 0)

elevationFile='/scratch/depfg/sutan101/data/pcrglobwb_input_arise/develop/global_30sec/landSurface/topography/merit_dem_processed/version_2021-02-XX/maps_covered_with_zero/dem_average_topography_parameters_30sec_february_2021_global_covered_with_zero.nc'
elev_30sec = xr.open_dataset(elevationFile)[['dem_average']].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
elev_30sec = elev_30sec.rename({'dem_average': 'elevation'})
elev_30sec.attrs = {}
elev_30sec = elev_30sec.reindex(lat=_grid_lat, lon=_grid_lon, method='nearest')
elev_30sec = elev_30sec.where(mask==1, np.nan)
elev_30sec.to_zarr(dataFolder / 'input.zarr', mode='w')
print('Elevation data loaded.')


url='https://geo.public.data.uu.nl/vault-pcrglobwb-cmip6/research-pcrglobwb-cmip6%5B1690540205%5D/original/hypflowsci6_v1.0/output/gswp3-w5e5/historical-reference/pcrglobwb_cmip6-isimip3-gswp3-w5e5_image-aqueduct_historical-reference_gwRecharge_global_monthly-total_1960_2019_basetier1.nc'
rechFile=tempFolder / 'temp_recharge.nc'
download_file_wget(url, rechFile)
recharge = xr.open_dataset(rechFile)['groundwater_recharge']
recharge = recharge.where(recharge >= 0, 0)
recharge = recharge.resample(time='YE').sum()
recharge = recharge.mean(dim='time')
recharge = recharge.sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
recharge = recharge.reindex(lat=_grid_lat, lon=_grid_lon, method='nearest')
recharge = recharge.sortby('lat')
recharge = recharge.interpolate_na(dim="lat", method="nearest", fill_value="extrapolate")
recharge = recharge.interpolate_na(dim="lon", method="nearest", fill_value="extrapolate")
recharge = recharge.sortby('lat', ascending=False)
recharge = recharge.where(mask==1, np.nan)
recharge.to_zarr(dataFolder / 'input.zarr', mode='a')
print('Recharge data loaded.')


aiFile='/scratch/depfg/7006713/data/aridity_index/aridity_index.nc'
ai_30sec = xr.open_dataset(aiFile)[['ai']].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
ai_30sec = ai_30sec.reindex(lat=_grid_lat, lon=_grid_lon, method='nearest')
ai_30sec.attrs = {}
ai_30sec = ai_30sec.where(mask==1, np.nan)
ai_30sec.to_zarr(dataFolder / 'input.zarr', mode='a')
print('Aridity index data loaded.')

slopeFile='/scratch/depfg/sutan101/data/pcrglobwb_input_arise/develop/global_30sec/landSurface/topography/merit_dem_processed/version_2021-02-XX/maps_covered_with_zero/tanslope_topography_parameters_30sec_february_2021_global_covered_with_zero.nc'
slope_30sec = xr.open_dataset(slopeFile)[['tanslope']].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
slope_30sec = slope_30sec.reindex(lat=_grid_lat, lon=_grid_lon, method='nearest')
slope_30sec = slope_30sec.rename({'tanslope': 'slope'})
slope_30sec.attrs = {}
slope_30sec = slope_30sec.where(mask==1, np.nan)
slope_30sec.to_zarr(dataFolder / 'input.zarr', mode='a')
print('Slope data loaded.')

url='https://geo.public.data.uu.nl/vault-globgm/research-globgm%5B1669042611%5D/original/input/version_1.0/k_conductivity_aquifer_filled_30sec.nc'
ksatFile=tempFolder / 'temp_ksat.nc'
download_file_wget(url, ksatFile)
ksat_30sec = xr.open_dataset(ksatFile)[['k_conductivity_aquifer_filled_map']].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
ksat_30sec = ksat_30sec.rename({'k_conductivity_aquifer_filled_map': 'k_sat'})
ksat_30sec = ksat_30sec.reindex(lat=_grid_lat, lon=_grid_lon, method='nearest')
ksat_30sec.attrs = {}
ksat_30sec = ksat_30sec.where(mask==1, np.nan)
ksat_30sec.to_zarr(dataFolder / 'input.zarr', mode='a')
print('Ksat data loaded.')

url='https://geo.public.data.uu.nl/vault-globgm/research-globgm%5B1669042611%5D/original/input/version_1.0/recession_coefficient_30sec.nc'
recFile=tempFolder / 'temp_rec_coeff.nc'
download_file_wget(url, recFile)
rec_coef_30sec = xr.open_dataset(recFile)[['recession_coefficient_30sec_map']].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
rec_coef_30sec = rec_coef_30sec.rename({'recession_coefficient_30sec_map': 'rec_coef'})
rec_coef_30sec = rec_coef_30sec.reindex(lat=_grid_lat, lon=_grid_lon, method='nearest')
rec_coef_30sec.attrs = {}
rec_coef_30sec = rec_coef_30sec.where(mask==1, np.nan)
rec_coef_30sec.to_zarr(dataFolder / 'input.zarr', mode='a')
print('Recession coefficient data loaded.')

url='https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_pet_penman_mean_1981-2010_V.2.1.tif'
petFile=tempFolder / 'temp_pet.tif'
petncFile= tempFolder / 'temp_petTEMP.nc'
invertedPetncFile= tempFolder / 'temp_pet.nc'
download_file_wget(url, petFile)
subprocess.run(['gdal_translate', '-of', 'netCDF', str(petFile), str(petncFile)], check=True)
subprocess.run(['cdo', 'invertlat', str(petncFile), str(invertedPetncFile)], check=True)
pet_ds = xr.open_dataset(invertedPetncFile)[['Band1']].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
pet_ds = pet_ds.reindex(lat=_grid_lat, lon=_grid_lon, method='nearest')
pet_ds = (pet_ds * 0.001) * 12
pet_ds = pet_ds.rename({'Band1': 'pet'})
pet_ds = pet_ds.sortby('lat')
pet_ds = pet_ds.interpolate_na(dim="lat", method="nearest", fill_value="extrapolate")
pet_ds = pet_ds.interpolate_na(dim="lon", method="nearest", fill_value="extrapolate")
pet_ds = pet_ds.sortby('lat', ascending=False)
pet_ds.attrs = {}
pet_ds = pet_ds.where(mask==1, np.nan)
pet_ds.to_zarr(dataFolder / 'input.zarr', mode='a')
print('PET data loaded.')

url='https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio12_1981-2010_V.2.1.tif'
tpFile=tempFolder / 'temp_tp.tif'
tpncFile= tempFolder / 'temp_tpTEMP.nc'
invertedTpncFile= tempFolder / 'temp_tp.nc'
download_file_wget(url, tpFile)
subprocess.run(['gdal_translate', '-of', 'netCDF', str(tpFile), str(tpncFile)], check=True)
subprocess.run(['cdo', 'invertlat', str(tpncFile), str(invertedTpncFile)], check=True)
tp_ds = xr.open_dataset(invertedTpncFile)[['Band1']].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
tp_ds = tp_ds.reindex(lat=_grid_lat, lon=_grid_lon, method='nearest')
tp_ds = (tp_ds *0.001)
tp_ds = tp_ds.rename({'Band1': 'tp'})
tp_ds.attrs = {}
tp_ds = tp_ds.where(mask==1, np.nan)
tp_ds.to_zarr(dataFolder / 'input.zarr', mode='a')
print('TP data loaded.')

# Get the clay percentage for fromSoil250m 
def getClayData(url):
    for layer in ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']:
        sub_url= f'{url}/clay_{layer}cm_mean_1000.tif'
        clay_perFile=tempFolder / f'temp_clay_perc_{layer}cm.tif'
        download_file_wget(sub_url, clay_perFile)
        ds = rioxarray.open_rasterio(clay_perFile, parse_coordinates=True)
        ds = ds.where(ds != -32768, np.nan)
        ds = ds / 10 
        if ds.rio.crs != 'EPSG:4326':
            ds = ds.rio.reproject('EPSG:4326')
        ds.to_zarr(tempFolder / f'clay_perc_{layer}cm.zarr', mode='w')
        os.remove(clay_perFile)
    for layer in ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']:
        ds = xr.open_zarr(tempFolder / f'clay_perc_{layer}cm.zarr').rename({'__xarray_dataarray_variable__': 'clay_perc', 
                                                                            'band': 'layer', 
                                                                            'x': 'lon',
                                                                            'y': 'lat'})
        ds = ds.drop_vars('spatial_ref')
        if layer == '0-5':
            clay_ds = ds
        else:
            clay_ds = xr.concat([clay_ds, ds], dim='layer')
    clay_ds = clay_ds.mean(dim='layer').compute()
    clay_ds.to_zarr(tempFolder / 'temp_clay_perc.zarr', mode='w')
    
    for layer in ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']:
        shutil.rmtree(tempFolder / f'clay_perc_{layer}cm.zarr')
        
url='https://files.isric.org/soilgrids/latest/data_aggregated/1000m/clay'
getClayData(url)
clay_perc = xr.open_zarr(tempFolder / 'temp_clay_perc.zarr')['clay_perc']
clay_perc = clay_perc.sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
clay_perc = clay_perc.reindex(lat=_grid_lat, lon=_grid_lon, method='nearest')
clay_perc = clay_perc.sortby('lat')
clay_perc = clay_perc.interpolate_na(dim="lat", method="nearest", fill_value="extrapolate")
clay_perc = clay_perc.interpolate_na(dim="lon", method="nearest", fill_value="extrapolate")
clay_perc = clay_perc.sortby('lat', ascending=False)
clay_perc = clay_perc.where(mask==1, np.nan)
clay_perc.to_zarr(dataFolder / 'input.zarr', mode='a')
print('Clay data loaded.')

# Get the sand percentage for fromSoil250m 
def getSandData(url):
    for layer in ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']:
        sub_url= f'{url}/sand_{layer}cm_mean_1000.tif'
        sand_perFile=tempFolder / f'temp_sand_perc_{layer}cm.tif'
        download_file_wget(sub_url, sand_perFile)
        ds = rioxarray.open_rasterio(sand_perFile, parse_coordinates=True)
        ds = ds.where(ds != -32768, np.nan)
        ds = ds / 10 
        if ds.rio.crs != 'EPSG:4326':
            ds = ds.rio.reproject('EPSG:4326')
        ds.to_zarr(tempFolder / f'sand_perc_{layer}cm.zarr', mode='w')
        os.remove(sand_perFile)
    for layer in ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']:
        ds = xr.open_zarr(tempFolder / f'sand_perc_{layer}cm.zarr').rename({'__xarray_dataarray_variable__': 'sand_perc', 
                                                                            'band': 'layer', 
                                                                            'x': 'lon',
                                                                            'y': 'lat'})
        ds = ds.drop_vars('spatial_ref')
        if layer == '0-5':
            sand_ds = ds
        else:
            sand_ds = xr.concat([sand_ds, ds], dim='layer')
    sand_ds = sand_ds.mean(dim='layer').compute()
    sand_ds.to_zarr(tempFolder / 'temp_sand_perc.zarr', mode='w')
    
    for layer in ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']:
        shutil.rmtree(tempFolder / f'sand_perc_{layer}cm.zarr')       
url='https://files.isric.org/soilgrids/latest/data_aggregated/1000m/sand'
getSandData(url)
sand_perc = xr.open_zarr(tempFolder / 'temp_sand_perc.zarr')['sand_perc']
sand_perc = sand_perc.sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
sand_perc = sand_perc.reindex(lat=_grid_lat, lon=_grid_lon, method='nearest')
sand_perc = sand_perc.sortby('lat')
sand_perc = sand_perc.interpolate_na(dim="lat", method="nearest", fill_value="extrapolate")
sand_perc = sand_perc.interpolate_na(dim="lon", method="nearest", fill_value="extrapolate")
sand_perc = sand_perc.sortby('lat', ascending=False)
sand_perc = sand_perc.where(mask==1, np.nan)
sand_perc.to_zarr(dataFolder / 'input.zarr', mode='a')
print('Sand data loaded.')

# Get the silt percentage for fromSoil250m 
def getSiltData(url):
    for layer in ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']:
        sub_url= f'{url}/silt_{layer}cm_mean_1000.tif'
        silt_perFile=tempFolder / f'temp_silt_perc_{layer}cm.tif'
        download_file_wget(sub_url, silt_perFile)
        ds = rioxarray.open_rasterio(silt_perFile, parse_coordinates=True)
        ds = ds.where(ds != -32768, np.nan)
        ds = ds / 10 
        if ds.rio.crs != 'EPSG:4326':
            ds = ds.rio.reproject('EPSG:4326')
        ds.to_zarr(tempFolder / f'silt_perc_{layer}cm.zarr', mode='w')
        os.remove(silt_perFile)
    for layer in ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']:
        ds = xr.open_zarr(tempFolder / f'silt_perc_{layer}cm.zarr').rename({'__xarray_dataarray_variable__': 'silt_perc', 
                                                                            'band': 'layer', 
                                                                            'x': 'lon',
                                                                            'y': 'lat'})
        ds = ds.drop_vars('spatial_ref')
        if layer == '0-5':
            silt_ds = ds
        else:
            silt_ds = xr.concat([silt_ds, ds], dim='layer')
    silt_ds = silt_ds.mean(dim='layer').compute()
    silt_ds.to_zarr(tempFolder / 'temp_silt_perc.zarr', mode='w')
    
    for layer in ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']:
        shutil.rmtree(tempFolder / f'silt_perc_{layer}cm.zarr')
        
url='https://files.isric.org/soilgrids/latest/data_aggregated/1000m/silt'
getSiltData(url)
silt_perc = xr.open_zarr(tempFolder / 'temp_silt_perc.zarr')['silt_perc']
silt_perc = silt_perc.sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
silt_perc = silt_perc.reindex(lat=_grid_lat, lon=_grid_lon, method='nearest')
silt_perc = silt_perc.sortby('lat')
silt_perc = silt_perc.interpolate_na(dim="lat", method="nearest", fill_value="extrapolate")
silt_perc = silt_perc.interpolate_na(dim="lon", method="nearest", fill_value="extrapolate")
silt_perc = silt_perc.sortby('lat', ascending=False)
silt_perc = silt_perc.where(mask==1, np.nan)
silt_perc.to_zarr(dataFolder / 'input.zarr', mode='a')
print('Silt data loaded.')

ds = xr.open_zarr(dataFolder / 'input.zarr')
ds = ds.drop_vars('spatial_ref', errors='ignore')
ds.to_zarr(dataFolder / 'input_final.zarr', mode='w')
print('All data loaded.')
shutil.rmtree(dataFolder / 'input.zarr')
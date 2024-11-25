from pathlib import Path
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree
import numpy as np


dataFolder = Path('/projects/0/einf4705/gwRecharge_correction/data')
outputFolder = Path('/projects/0/einf4705/gwRecharge_correction/data_processed')
outputFolder.mkdir(parents=True, exist_ok=True)
simRechargeFile='/projects/0/einf4705/_data/cmip6_input/gswp3-w5e5/historical_natural/pcrglobwb_cmip6-isimip3-gswp3-w5e5_image-aqueduct_historical-natural_gwRecharge_global_monthly-total_1960_2019_basetier1.nc'

latMin, latMax = -90, 90
lonMin, lonMax = -180,180

obs_data = pd.read_csv(dataFolder / 'recharge.csv', sep=';')
obs_data = obs_data.replace(',', '.', regex=True)
obs_data = obs_data.apply(pd.to_numeric, errors='ignore')
obs_data.rename(columns={'Latitude': 'lat','Longitude': 'lon','Groundwater recharge': 'obs_recharge'}, inplace=True)
obs_data = obs_data[(obs_data['lat'] >= latMin) & (obs_data['lat'] <= latMax) & (obs_data['lon'] >= lonMin) & (obs_data['lon'] <= lonMax)]
obs_data['obs_recharge'] = obs_data['obs_recharge'] / 1000

latIndicies = obs_data['lat'].values
lonIndicies = obs_data['lon'].values

gwRech_30sec = xr.open_dataset(simRechargeFile)[['groundwater_recharge']].sel(latitude=slice(latMax, latMin), longitude=slice(lonMin, lonMax)).rename({'latitude': 'lat', 'longitude': 'lon'})
gwRech_30sec = gwRech_30sec.resample(time='YE').sum(skipna=False)
gwRech_30sec = gwRech_30sec.mean(dim='time', skipna=False).compute()
data_ds = gwRech_30sec.sel(lat=xr.DataArray(latIndicies, dims="z"), lon=xr.DataArray(lonIndicies, dims="z"), method='nearest').rename({'groundwater_recharge' : 'sim_recharge'})

ai_30sec = xr.open_dataset(dataFolder / 'ai_v3_year.nc')['Band1'].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
ai_30sec = ai_30sec.sel(lat=xr.DataArray(latIndicies, dims="z"), lon=xr.DataArray(lonIndicies, dims="z"), method='nearest')
ai_30sec = ai_30sec * 0.0001
data_ds = data_ds.assign(AI=ai_30sec)

elev_30sec = xr.open_dataset(dataFolder / 'dem_average_topography_parameters_30sec_february_2021_global_covered_with_zero.nc')['dem_average'].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
elev_30sec = elev_30sec.sel(lat=xr.DataArray(latIndicies, dims="z"), lon=xr.DataArray(lonIndicies, dims="z"), method='nearest')
data_ds = data_ds.assign(elev=elev_30sec)

reces_coef_30sec = xr.open_dataset(dataFolder/ 'recession_coefficient_30sec.nc')['recession_coefficient_30sec_map'].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
reces_coef_30sec = reces_coef_30sec.sel(lat=xr.DataArray(latIndicies, dims="z"), lon=xr.DataArray(lonIndicies, dims="z"), method='nearest')
data_ds = data_ds.assign(J=reces_coef_30sec)
data_ds = data_ds.to_dataframe()
obs_coords = obs_data[['lat', 'lon']].values
sim_coords = data_ds[['lat', 'lon']].values

tree = cKDTree(sim_coords)
distances, indices = tree.query(obs_coords)
final_df = data_ds.iloc[indices]
final_df['obs_recharge'] = obs_data['obs_recharge'].reset_index(drop=True)
final_df.to_parquet(outputFolder / 'data_train.parquet')
final_df.to_parquet(outputFolder / 'data_train.parquet')
print('trainData created')
# # # # ##Create prediction dataset
gird= xr.open_dataset('/projects/0/einf4705/gwRecharge_correction/data/recession_coefficient_30sec.nc')
data_ds = xr.open_dataset(simRechargeFile)[['groundwater_recharge']].sel(latitude=slice(latMax, latMin), longitude=slice(lonMin, lonMax)).rename({'latitude': 'lat', 'longitude': 'lon'})
data_ds = data_ds.resample(time='YE').sum(skipna=False)
data_ds = data_ds.mean(dim='time', skipna=False).compute()
data_ds = data_ds.rename({'groundwater_recharge' : 'sim_recharge'})
data_ds = data_ds.reindex_like(gird, method='nearest')

elev_30sec = xr.open_dataset(dataFolder / 'dem_average_topography_parameters_30sec_february_2021_global_covered_with_zero.nc')['dem_average'].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
elev_30sec = elev_30sec.reindex_like(data_ds, method='nearest')
data_ds = data_ds.assign(elev=elev_30sec)

ai_30sec = xr.open_dataset(dataFolder / 'ai_v3_year.nc')['Band1'].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
ai_30sec = ai_30sec.reindex_like(data_ds, method='nearest')	
ai_30sec = ai_30sec * 0.0001
data_ds = data_ds.assign(AI=ai_30sec)

reces_coef_30sec = xr.open_dataset(dataFolder  / 'recession_coefficient_30sec.nc')['recession_coefficient_30sec_map'].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).compute()
reces_coef_30sec = reces_coef_30sec.reindex_like(data_ds, method='nearest')	
data_ds = data_ds.assign(J=reces_coef_30sec)
data_ds = data_ds.to_dataframe().reset_index()
data_ds = data_ds.dropna()
data_ds.to_parquet(outputFolder / 'data_predict.parquet')
print('Predict created')

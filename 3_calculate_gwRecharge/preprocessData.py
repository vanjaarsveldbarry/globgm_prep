from pathlib import Path
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree
import numpy as np


dataFolder = Path('/scratch/depfg/7006713/globgm_input/_tempGLOBGM/gwRecharge/data')
dataFolder.mkdir(parents=True, exist_ok=True)
outputFolder = Path('/scratch/depfg/7006713/globgm_input/_tempGLOBGM/gwRecharge/out')
outputFolder.mkdir(parents=True, exist_ok=True)

simRechargeFile='/scratch/depfg/7006713/globgm_input/cmip6_input/gswp3-w5e5/historical/pcrglobwb_cmip6-isimip3-gswp3-w5e5_image-aqueduct_historical-reference_gwRecharge_global_monthly-total_1960_2019_basetier1.nc'

latMin, latMax = -35.0, -30.0
lonMin, lonMax = 16.0, 33.0

obs_data = pd.read_csv('/scratch/depfg/7006713/data/recharge/global_groundwater_recharge_moeck-et-al.csv', sep=';')
obs_data = obs_data.replace(',', '.', regex=True)
obs_data = obs_data.apply(pd.to_numeric)
obs_data.rename(columns={'Latitude': 'lat','Longitude': 'lon','Groundwater recharge': 'obs_recharge'}, inplace=True)
obs_data = obs_data[(obs_data['lat'] >= latMin) & (obs_data['lat'] <= latMax) & (obs_data['lon'] >= lonMin) & (obs_data['lon'] <= lonMax)]
obs_data['obs_recharge'] = obs_data['obs_recharge'] / 1000

latIndicies = obs_data['lat'].values
lonIndicies = obs_data['lon'].values

sim_rech = xr.open_dataset(simRechargeFile)[['groundwater_recharge']].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax))
sim_rech = sim_rech.isel(time=slice(0, 2))
sim_rech = sim_rech.resample(time='YE').sum(skipna=False)
sim_rech = sim_rech.mean(dim='time', skipna=False).compute()
data_ds = sim_rech.sel(lat=xr.DataArray(latIndicies, dims="z"), lon=xr.DataArray(lonIndicies, dims="z"), method='nearest').rename({'groundwater_recharge' : 'sim_recharge'})

ds = xr.open_zarr(dataFolder / 'input.zarr').sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax))
ds = ds.sel(lat=xr.DataArray(latIndicies, dims="z"), lon=xr.DataArray(lonIndicies, dims="z"), method='nearest')
print(ds)
data_ds = data_ds.assign(AI=ds['ai'])
data_ds = data_ds.assign(elev=ds['elevation'])
data_ds = data_ds.assign(slope=ds['slope'])
data_ds = data_ds.assign(J=ds['rec_coef'])

data_ds = data_ds.to_dataframe().reset_index()
obs_coords = obs_data[['lat', 'lon']].values
sim_coords = data_ds[['lat', 'lon']].values

tree = cKDTree(sim_coords)
distances, indices = tree.query(obs_coords)
final_df = data_ds.iloc[indices]
obs_vals = obs_data['obs_recharge'].values
final_df['obs_recharge'] = obs_data['obs_recharge'].values
final_df = final_df.dropna()
final_df.to_parquet(dataFolder / 'data_train.parquet')
print('trainData created')
#Create prediction dataset
grid_lat = xr.open_zarr(dataFolder / 'input.zarr').sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).lat.values
grid_lon = xr.open_zarr(dataFolder / 'input.zarr').sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax)).lon.values
data_ds = xr.open_dataset(simRechargeFile)[['groundwater_recharge']].sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax))
data_ds = data_ds.isel(time=slice(0, 2))
data_ds = data_ds.resample(time='YE').sum(skipna=False)
data_ds = data_ds.mean(dim='time', skipna=False).compute()
data_ds = data_ds.rename({'groundwater_recharge' : 'sim_recharge'})
data_ds = data_ds.reindex(lat=grid_lat, lon=grid_lon, method='nearest')
ds = xr.open_zarr(dataFolder / 'input.zarr')
ds = ds.sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax))
data_ds = data_ds.assign(AI=ds['ai'])
data_ds = data_ds.assign(elev=ds['elevation'])
data_ds = data_ds.assign(slope=ds['elevation'])
data_ds = data_ds.assign(J=ds['elevation'])
data_ds = data_ds.compute()
data_ds = data_ds.to_dataframe().reset_index()
data_ds = data_ds.dropna()
data_ds.to_parquet(dataFolder / 'data_predict.parquet')
print('Predict created')

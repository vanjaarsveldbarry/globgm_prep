from pathlib import Path
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree
import numpy as np


rootFolder = Path('/scratch/depfg/7006713/gwRecharge_corrected')
tempFolder = rootFolder / '_temp'
tempFolder.mkdir(parents=True, exist_ok=True)
dataFolder = rootFolder / 'data'
dataFolder.mkdir(parents=True, exist_ok=True)

latMin, latMax = -90.0, 90.0
lonMin, lonMax = -180.0, 180.0
# latMin, latMax = -44.0, -10.0
# lonMin, lonMax = 112.0, 154.0

obs_data = pd.read_csv('/scratch/depfg/7006713/data/recharge/global_groundwater_recharge_moeck-et-al.csv', sep=';')
obs_data = obs_data.replace(',', '.', regex=True)
obs_data = obs_data.apply(pd.to_numeric)
obs_data.rename(columns={'Latitude': 'lat','Longitude': 'lon','Groundwater recharge': 'obs_recharge'}, inplace=True)
obs_data = obs_data[(obs_data['lat'] >= latMin) & (obs_data['lat'] <= latMax) & (obs_data['lon'] >= lonMin) & (obs_data['lon'] <= lonMax)]
obs_data['obs_recharge'] = obs_data['obs_recharge'] / 1000

latIndicies = obs_data['lat'].values
lonIndicies = obs_data['lon'].values

data_ds = xr.open_zarr(dataFolder / 'input_final.zarr').sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax))
data_ds = data_ds.sel(lat=xr.DataArray(latIndicies, dims="z"), lon=xr.DataArray(lonIndicies, dims="z"), method='nearest')
data_ds = data_ds.to_dataframe().reset_index()
data_ds.drop(columns=['z'], inplace=True)
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
data_ds = xr.open_zarr(dataFolder / 'input_final.zarr')
data_ds = data_ds.sel(lat=slice(latMax, latMin), lon=slice(lonMin, lonMax))
data_ds = data_ds.astype('float32')
data_ds = data_ds.compute()
data_ds = data_ds.to_dataframe().reset_index()
data_ds = data_ds.dropna()
data_ds.to_parquet(dataFolder / 'data_predict.parquet')
print('Predict created')

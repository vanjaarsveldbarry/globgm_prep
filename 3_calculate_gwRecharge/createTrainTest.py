from pathlib import Path
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree
import numpy as np
import zarr

def pre_process_train(df):
    columns_to_clip = [col for col in df.columns if col not in ['lat', 'lon']]
    df[columns_to_clip] = df[columns_to_clip].clip(lower=0)
    df[['AI', 'sim_recharge' ,'J', 'elev', 'obs_recharge']] = df[['AI', 'sim_recharge' ,'J', 'elev', 'obs_recharge']].apply(np.sqrt)
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.dropna(inplace=True)
    y = df['obs_recharge'].to_numpy()
    X = df[['AI', 'sim_recharge', 'J', 'elev']].to_numpy()
    return y, X

def pre_process_predict(df):
    columns_to_clip = [col for col in df.columns if col not in ['lat', 'lon']]
    df[columns_to_clip] = df[columns_to_clip].clip(lower=0)
    df[['AI', 'sim_recharge' ,'J', 'elev',]] = df[['AI', 'sim_recharge' ,'J', 'elev']].apply(np.sqrt)
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.dropna(inplace=True)
    lat_lon = df[['lat', 'lon']]
    X = df[['AI', 'sim_recharge', 'J', 'elev']].to_numpy()
    return X, lat_lon

dataFolder = Path('/scratch/depfg/7006713/globgm_input/_tempGLOBGM/gwRecharge/data')

df = pd.read_parquet(dataFolder / 'data_train.parquet')
y, X = pre_process_train(df)
zarr.save(dataFolder / 'y_train.zarr', y)
zarr.save(dataFolder / 'X_train.zarr', X)

df = pd.read_parquet(dataFolder / 'data_predict.parquet')
X, lat_lon = pre_process_predict(df)
zarr.save(dataFolder / 'X_predict.zarr', X)
lat_lon.to_parquet(dataFolder / 'lat_lon_predict.parquet')
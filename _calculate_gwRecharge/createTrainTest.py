from pathlib import Path
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree
import numpy as np
import zarr
import seaborn as sns

def pre_process_train(df):
    columns_to_drop = ['lat', 'lon', 'groundwater_recharge']  # Replace with actual column names
    df.drop(columns=columns_to_drop, inplace=True)
    df['elevation'] = df['elevation'].clip(lower=0)
    df= df.apply(np.sqrt)
    y = df['obs_recharge'].to_numpy()
    X = df.drop(columns=['obs_recharge'])
    X = X.to_numpy()
    return y, X

def pre_process_predict(df):
    columns_to_drop = ['groundwater_recharge']  # Replace with actual column names
    df.drop(columns=columns_to_drop, inplace=True)
    df['elevation'] = df['elevation'].clip(lower=0)
    df.loc[:, df.columns.difference(['lat', 'lon'])] = df.loc[:, df.columns.difference(['lat', 'lon'])].apply(np.sqrt)
    lat_lon = df[['lat', 'lon']]
    columns_to_drop = ['lat', 'lon']  # Replace with actual column names
    df.drop(columns=columns_to_drop, inplace=True)
    X = df.to_numpy()
    return X, lat_lon

rootFolder = Path('/scratch/depfg/7006713/gwRecharge_corrected')
tempFolder = rootFolder / '_temp'
tempFolder.mkdir(parents=True, exist_ok=True)
dataFolder = rootFolder / 'data'
dataFolder.mkdir(parents=True, exist_ok=True)

df = pd.read_parquet(dataFolder / 'data_train.parquet')
y, X = pre_process_train(df)
zarr.save(dataFolder / 'y_train.zarr', y)
zarr.save(dataFolder / 'X_train.zarr', X)
# import matplotlib.pyplot as plt
# y = pd.DataFrame(y, columns=['groundwater_recharge'])
# X = pd.DataFrame(X)
# # Create a DataFrame for pairplot
# df_pairplot = pd.concat([y, X], axis=1)
# sns.pairplot(df_pairplot, y_vars=y.columns, x_vars=X.columns)
# plt.savefig('/scratch/depfg/7006713/gwRecharge_corrected/pairplot_y_vs_X.png')

df = pd.read_parquet(dataFolder / 'data_predict.parquet')
X, lat_lon = pre_process_predict(df)
zarr.save(dataFolder / 'X_predict.zarr', X)
lat_lon.to_parquet(dataFolder / 'lat_lon_predict.parquet')
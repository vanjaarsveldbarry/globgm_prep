from pygam import LinearGAM, s, l, te
from pathlib import Path
import pandas as pd
import numpy as np
import xarray as xr
import zarr
import time
import dask
from concurrent.futures import ThreadPoolExecutor, as_completed
import shutil

rootFolder = Path('/scratch/depfg/7006713/gwRecharge_corrected')
tempFolder = rootFolder / '_tempPredict'
if tempFolder.exists():
    shutil.rmtree(tempFolder)
tempFolder.mkdir(parents=True, exist_ok=True)
dataFolder = rootFolder / 'data'
dataFolder.mkdir(parents=True, exist_ok=True)
outputFolder = rootFolder / 'output'
outputFolder.mkdir(parents=True, exist_ok=True)


y = zarr.load(dataFolder / 'y_train.zarr')
X = zarr.load(dataFolder / 'X_train.zarr')
gam = LinearGAM(s(0) + s(1) + s(2) + s(3)  + s(4)  + s(5) + s(6) + s(7) + s(8) + s(9) , fit_intercept=True).gridsearch(X, y)
gam.summary()

predict_df = pd.read_parquet(dataFolder / 'lat_lon_predict.parquet')
X_pred = zarr.load(dataFolder / 'X_predict.zarr')
total_start_time = time.time()
num_rows = len(predict_df)
batch_size = 5000000
total_iterations = (num_rows + batch_size - 1) // batch_size 
print(f"Total iterations: {total_iterations} for {num_rows} rows")
def process_batch(start, end):
    predict_batch = predict_df.iloc[start:end]
    predict_batch = predict_batch.copy()
    predict_batch['predicted_recharge'] = np.square(gam.predict(X_pred[start:end]))
    predict_batch = xr.DataArray(predict_batch.groupby(['lat', 'lon'])['predicted_recharge'].first().unstack()).rename('predicted_recharge')
    predict_batch.to_zarr(tempFolder / f'predicted_recharge_{end}.zarr', mode='w')

# # Use ThreadPoolExecutor to process batches in parallel
with ThreadPoolExecutor(max_workers=22) as executor:
    futures = []
    iteration_times = []
    total_start_time = time.time()
    
    for start in range(0, num_rows, batch_size):
        end = min(start + batch_size, num_rows)
        futures.append(executor.submit(process_batch, start, end))
    
    count = 0
    for iteration, future in enumerate(as_completed(futures), 1):
        future.result()
        count += 1
        print(f"{count} / {total_iterations}")
print('Done with the parallel processing')
with dask.config.set(**{'array.slicing.split_large_chunks': False}):
    count = 0
    predicted_recharge_files = sorted(tempFolder.glob('predicted_recharge_*.zarr'))
    total = len(predicted_recharge_files)
    for file in predicted_recharge_files:
        count = count + 1
        ds = xr.open_zarr(file).compute()
        if count == 1: 
            ds_final=ds
        else:
            ds_final = ds_final.combine_first(ds).compute()
        print(f'{count} / {total}')
grid_lat = xr.open_zarr(dataFolder / 'input_final.zarr').lat.values
grid_lon = xr.open_zarr(dataFolder / 'input_final.zarr').lon.values
ds_final = ds_final.reindex(lat=grid_lat, lon=grid_lon, method='nearest')
print(ds_final)
ds_final.to_netcdf(outputFolder / 'predicted_recharge.nc', mode='w')
print('Done with merging')
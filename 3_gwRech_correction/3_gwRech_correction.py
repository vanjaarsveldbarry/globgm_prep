import xarray as xr 
from pathlib import Path
import pyinterp
import pyinterp.backends.xarray
import numpy as np 
import sys
import pyinterp.fill

saveFolder = sys.argv[1]
sim = sys.argv[2]

input_folder=Path(f'{saveFolder}')
predRechargeFile = list(input_folder.glob('*predicted_recharge*.nc'))[0]
gw_corrected_field = xr.open_dataset(f'{predRechargeFile}')['predicted_recharge']
gw_corrected_field = gw_corrected_field / 365.24
gw_corrected_field = gw_corrected_field.sortby('lat')
lat_grid = gw_corrected_field['lat'].values
lon_grid = gw_corrected_field['lon'].values

gwRecharge_file = list(input_folder.glob('*gwRecharge_global_monthly*.nc'))[0]
ds = xr.open_dataset(gwRecharge_file, chunks='auto', engine='h5netcdf')['groundwater_recharge']
ds = xr.where(ds < 0, 0.0, ds)
if 'lat' in ds.dims:
    ds = ds.rename({'lat': 'latitude'})
if 'lon' in ds.dims:
    ds = ds.rename({'lon': 'longitude'})
ds = ds.resample(time='YE').sum(skipna=False)
ds = ds.mean(dim='time', skipna=False).compute().sortby('latitude')
ds = ds / 365.25

ds_fill = pyinterp.backends.xarray.Grid2D(ds, geodetic=False)
filled = pyinterp.fill.loess(ds_fill, nx=3, ny=3)
ds.values = filled.T
ds = pyinterp.backends.xarray.Grid2D(ds, geodetic=False)
mx, my = np.meshgrid(lon_grid, lat_grid, indexing="ij")
ds = ds.bivariate(coords=dict(longitude=mx.ravel(), latitude=my.ravel()), num_threads=0)
ds = ds.reshape(mx.shape).T
ds = xr.DataArray(ds, coords=[lat_grid, lon_grid], dims=['lat', 'lon']).sortby('lat')
gw_corrected_field_min = gw_corrected_field.min().values
cf = (gw_corrected_field + gw_corrected_field_min) / (ds + gw_corrected_field_min)
cf = cf.rename('correction_factor').to_dataset()
cf = cf.sortby('lat', ascending=False)
cf['gw_corrected_field_min'] = xr.DataArray(gw_corrected_field_min)
cf.to_zarr(input_folder / 'gwRecharge_correction_factor.zarr', mode='w')
import xarray as xr

ds = xr.open_dataset('/scratch/depfg/7006713/gwRecharge_corrected/output/predicted_recharge.nc')
ds = ds.sel(lat=slice(-10, -44), lon=slice(112, 154))
ds.to_netcdf('/scratch/depfg/7006713/gwRecharge_corrected/test_out.nc')

ds = xr.open_zarr('/scratch/depfg/7006713/gwRecharge_corrected/data/input_final.zarr')
ds = ds.sel(lat=slice(-10, -44), lon=slice(112, 154))
ds.to_netcdf('/scratch/depfg/7006713/gwRecharge_corrected/test_input.nc')
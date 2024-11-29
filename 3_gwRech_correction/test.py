import xarray as xr

ds = xr.open_zarr('/scratch/depfg/7006713/globgm_input/cmip6_input/gswp3-w5e5/historical_natural/gwRecharge_correction_factor.zarr')
ds = ds.sel(lat=slice(-10, -44), lon=slice(112, 154))
ds.to_netcdf('/scratch/depfg/7006713/globgm_input/cmip6_input/gswp3-w5e5/historical_natural/test_out.nc')
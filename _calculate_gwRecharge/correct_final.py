from pathlib import Path
import xarray as xr

rootFolder = Path('/scratch/depfg/7006713/gwRecharge_corrected')
dataFolder = rootFolder / 'data'
dataFolder.mkdir(parents=True, exist_ok=True)
outputFolder = rootFolder / 'output'
outputFolder.mkdir(parents=True, exist_ok=True)


pred_gwRecharge = xr.open_dataset(outputFolder / 'predicted_recharge_pre_correction.nc', chunks={}).compute()
ds = xr.open_zarr(dataFolder / 'input_final.zarr')[['tp', 'groundwater_recharge']].compute()
pred_gwRecharge['tp'] = ds['tp']
pred_gwRecharge['groundwater_recharge'] = ds['groundwater_recharge']
pred_gwRecharge['predicted_recharge'] = xr.where(pred_gwRecharge['predicted_recharge'] > pred_gwRecharge['tp'], pred_gwRecharge['tp'], pred_gwRecharge['predicted_recharge'])
pred_gwRecharge['predicted_recharge'] = xr.where(pred_gwRecharge['groundwater_recharge'] == 0.0, 0.0, pred_gwRecharge['predicted_recharge'])
pred_gwRecharge.to_netcdf(outputFolder / 'predicted_recharge.nc')
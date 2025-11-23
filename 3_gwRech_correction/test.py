import xarray as xr
import seaborn as sns
import pandas as pd
import seaborn as sns

# ds1 = xr.open_zarr('/scratch/depfg/7006713/globgm_input/cmip6_input/gswp3-w5e5/historical_natural/NATURAL_gwRecharge_correction_factor.zarr')
# ds1 = ds1.sel({'lat': slice(-22.0, -35.0), 'lon': slice(16.0, 33.0)})
# ds1.to_netcdf('/eejit/home/7006713/projects/globgm_prep/3_gwRech_correction/test1.nc')
# print(ds1)
# ds2 = xr.open_zarr('/scratch/depfg/7006713/globgm_input/cmip6_input/gswp3-w5e5/historical/gwRecharge_correction_factor.zarr')
# ds2 = ds2.sel({'lat': slice(-22.0, -35.0), 'lon': slice(16.0, 33.0)})
# ds2.to_netcdf('/eejit/home/7006713/projects/globgm_prep/3_gwRech_correction/test2.nc')
# print(ds2)


# ds_nat = xr.open_dataarray('/eejit/home/7006713/projects/globgm_prep/3_gwRech_correction/test1.nc').values
# ds_non_nat = xr.open_dataarray('/eejit/home/7006713/projects/globgm_prep/3_gwRech_correction/test2.nc').values
# print(ds_nat)
# print(ds_non_nat)

# import matplotlib.pyplot as plt

# # Flatten the arrays
# ds_nat_flat = ds_nat.flatten()
# ds_non_nat_flat = ds_non_nat.flatten()

# # Create a DataFrame for seaborn
# data = pd.DataFrame({
#     'Natural': ds_nat_flat,
#     'Non-Natural': ds_non_nat_flat
# })

# # Plot the CDF
# sns.ecdfplot(data=data, palette=["blue", "red"])
# plt.xlabel('Value')
# plt.ylabel('ECDF')
# plt.title('CDF of Natural and Non-Natural Data')
# # plt.xlim(0, 0.2)
# plt.xscale('symlog')
# plt.savefig('/eejit/home/7006713/projects/globgm_prep/3_gwRech_correction/test.png')

import xarray as xr
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from geodatasets import get_path

# Load the data
df = pd.read_csv('/scratch/depfg/7006713/data/recharge/global_groundwater_recharge_moeck-et-al.csv', delimiter=';')
df = df.rename(columns={'Recharge [mm/yr]': 'recharge', 'Latitude': 'lat', 'Longitude': 'lon', 'Groundwater recharge': 'recharge'})	
df = df.applymap(lambda x: float(str(x).replace(',', '.')))
df = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat))
print(df)

# Plot the points on a global map
path = get_path("naturalearth.land")
world = gpd.read_file(path)
fig, ax = plt.subplots(figsize=(15, 10))
world.plot(ax=ax, color='lightgrey')
df.plot(ax=ax, markersize=10, color='red', alpha=0.4)
plt.title('Global Map with Locations of Each Point')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.savefig('/eejit/home/7006713/projects/globgm_prep/3_gwRech_correction/glyobal_map.png')
# Convert the DataFrame to a GeoDataFrame

from getpass import getpass
import json
import ssl
import sys
from pathlib import Path
from tqdm import tqdm
from irods.models import Collection
from irods.models import DataObject
from irods.session import iRODSSession
from irods.exception import DataObjectDoesNotExist
import xarray as xr
import shutil 

def setup_iRodsSession(env_config, password, ca_file):
    
    def get_irods_environment(irods_environment_file="irods_environment.json"):
        """Reads the irods_environment.json file, which contains the environment
        configuration."""
        with open(irods_environment_file, 'r') as f:
            return json.load(f)

    def setup_session(irods_environment_config,  password, require_ssl = True, ca_file = ca_file):
        """Use irods environment files to configure a iRODSSession"""

        if require_ssl:
            ssl_context = ssl.create_default_context(purpose=ssl.Purpose.SERVER_AUTH, cafile=ca_file, capath=None, cadata=None)
            ssl_settings = {'client_server_negotiation': 'request_server_negotiation',
                            'client_server_policy': 'CS_NEG_REQUIRE',
                            'encryption_algorithm': 'AES-256-CBC',
                            'encryption_key_size': 32,
                            'encryption_num_hash_rounds': 16,
                            'encryption_salt_size': 8,
                            'ssl_context': ssl_context}
            session = iRODSSession(
                irods_password=password,
                **irods_environment_config,
                **ssl_settings
            )
        else:
            session = iRODSSession(
                password=password,
                **irods_environment_config,
            )

        return session

    session = setup_session(get_irods_environment(env_config), password, ca_file)
    if session is None:
        print("Error: unable to create session.")
        sys.exit(1)
    else:
        return session
    
# Create a session
env_config = '/home/bvjaarsv/.irods/irods_environmentDAG.json'
password = 'GXqmWHNtFHmljrvGMi5ihuU7OT9ZI-4p'
ca_file = '/home/bvjaarsv/.irods/cacert.pem'

savePath = '/projects/0/einf4705/gwRecharge_correction/data/_temp'
projectDirectory='/nluu14p/home/deposit-pilot/geowat_global_output'
variables = [
    'gwRecharge_annuaTot_output.nc',
    ]

# for var in variables:
#     session = setup_iRodsSession(env_config, password, ca_file)
#     with session as s:
#         query = s.query(Collection.parent_name, Collection.name, DataObject.id, DataObject.name).filter(DataObject.name == var)
#         allFiles = [f'{result[Collection.name]}/{var}' for result in query.get_results()]
#         allFiles = [file for file in allFiles if 'initial_conditions' not in file and 'trash' not in file]
#         allFiles = [file for file in allFiles if '30min' not in file and '5min' not in file and 'europe_30sec' not in file]

#     for server_path in tqdm(allFiles, desc=f'Downloading {var}', disable=False):
#         print(server_path)
#         session = setup_iRodsSession(env_config, password ,ca_file)
#         with session as s:
#             local_path = Path(str(server_path).replace(projectDirectory, savePath))
#             local_path.parent.mkdir(parents=True, exist_ok=True)
#             s.data_objects.get(server_path, local_path)

#         def annualMean(file):
#             ds = xr.open_dataset(file, chunks=None, engine='h5netcdf')
#             if 'annua' in file.as_posix():
#                 ds = ds.mean('time', skipna=False).compute()
#             file_zarr = file.with_suffix('.zarr')
#             ds = ds.to_zarr(file_zarr, mode='w')
#             Path(file).unlink()

#         ds = annualMean(local_path)

savePath = '/projects/0/einf4705/gwRecharge_correction/data/_temp/30sec/newDownscale'
# def merge_and_average_zarr_files(input_dir, output_dir):
#     for subdir in Path(input_dir).iterdir():
#         if subdir.is_dir():
#             zarr_files = list(subdir.glob('**/*.zarr'))
#             print(zarr_files)
            
#             if len(zarr_files) == 3:
#                 ds = xr.open_mfdataset(zarr_files, engine='zarr', concat_dim='time', combine='nested')
#                 ds = ds.mean(dim='time', skipna=False)
#                 subSavePath = output_dir / subdir.name / 'geRecharge_annuaTot_output.zarr'
#                 ds.to_zarr(subSavePath, mode='w')
#                 for file in zarr_files:
#                     shutil.rmtree(file)

output_dir = Path(savePath).parent.parent / 'merged_averaged_results'
# output_dir.mkdir(parents=True, exist_ok=True)
# merge_and_average_zarr_files(savePath, output_dir)
# print('merged done')


zarr_files = list(output_dir.glob('**/*.zarr'))
count = 0
for file in tqdm(zarr_files):
    count = count + 1
    ds = xr.open_zarr(file).compute()
    if count == 1: 
        ds_final=ds.compute()
    else:
        ds_final = ds_final.combine_first(ds).compute()

grid = xr.open_dataset('/projects/0/einf4705/gwRecharge_correction/data/recession_coefficient_30sec.nc')
latGrid = grid.lat
lonGrid = grid.lon
ds_final = ds_final.reindex(lat=latGrid, lon=lonGrid, method='nearest')
ds_final.to_zarr('/projects/0/einf4705/gwRecharge_correction/data/geRecharge_annual_average_1979_2019.zarr', mode='w')
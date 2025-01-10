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
from irods.column import Like

def setup_iRodsSession(env_config, password):
    
    def get_irods_environment(irods_environment_file="irods_environment.json"):
        """Reads the irods_environment.json file, which contains the environment
        configuration."""
        with open(irods_environment_file, 'r') as f:
            return json.load(f)

    def setup_session(irods_environment_config,  password, require_ssl = True, ca_file = "/eejit/home/7006713/.irods/cacert.pem"):
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

    session = setup_session(get_irods_environment(env_config), password)
    if session is None:
        print("Error: unable to create session.")
        sys.exit(1)
    else:
        return session

saveDirectory = Path(f"{sys.argv[1]}/{sys.argv[2]}/historical/")
password = sys.argv[3]

env_config = '/eejit/home/7006713/.irods/irods_environmentDAG.json'

projectDirectory='/nluu14p/home/deposit-pilot/globgm/preprocess/corrected_recharge'

var='predicted_recharge.nc'
session = setup_iRodsSession(env_config, password)
with session as s:
    query = s.query(Collection.parent_name, Collection.name, DataObject.id, DataObject.name).filter(Like(DataObject.name, f'%{var}%'))
    allFiles = [f'{result[Collection.name]}/{result[DataObject.name]}' for result in query.get_results()]
    allFiles = [file for file in allFiles if 'trash' not in file]
    for server_path in allFiles:
        local_path = Path(str(server_path).replace(projectDirectory, str(saveDirectory)))
        s.data_objects.get(server_path, local_path, forceFlag=True)
# Downloads molecular reactions dataset and model config file

import boto3
from botocore import UNSIGNED
from botocore.client import Config
from pathlib import Path

# Set up S3 client to access public folders 
s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))

# config dictionary
data_config =   {
                    'bucket_name' : 'deepsynthesis',
                    'data_directory' : Path('public_data/MIT_mixed_augm_char/'),
                    'data_destination' : Path('OpenNMT-py/data/molecular_data/'),
                    'data_files' :  [
                                        'source_train.txt',
                                        'source_val.txt',
                                        'source_test.txt',
                                        'target_train.txt',
                                        'target_val.txt',
                                        'target_test.txt'
                                    ],
                    'config_directory' : Path('public_data/config_files/'),
                    'model_config' : 'molecular_transformer.yml',
                    'config_destination' : Path('OpenNMT-py/config/')
                }

# function to download a file if the destination does not already exist
def download_file(bucket, s3_directory, filename, destination_directory):
    file_destination = destination_directory/filename
    file_source = s3_directory/filename 
    if not file_destination.exists():
        print(f'Downloading {filename}')
        s3.download_file(bucket,
                         str(file_source),
                         str(file_destination))
    else:
        print(f'File {filename} already downloaded')


bucket = data_config['bucket_name']
data_destination = data_config['data_destination']

# create data directory
if not data_destination.exists():
    data_destination.mkdir(exist_ok=True)

print('Downloading data files')
for filename in data_config['data_files']:
    download_file(bucket, data_config['data_directory'], filename, data_destination)

print('Downloading model config')
download_file(bucket, data_config['config_directory'], 
              data_config['model_config'], data_config['config_destination'])
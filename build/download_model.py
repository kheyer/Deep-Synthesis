import boto3
from botocore import UNSIGNED
from botocore.client import Config
from pathlib import Path

s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))

bucket_name = 'deepsynthesis'
model_file = 'public_models/Molecular_transformer_step_400000.pt'
destination_file = 'models/Molecular_transformer_step_400000.pt'

path = Path('.')
file_path = path/destination_file 

if not file_path.exists():
    print('Downloading Model File')
    (path/'models').mkdir(exist_ok=True)
    s3.download_file(bucket_name, model_file, destination_file)
else:
    print('Model already downloaded')
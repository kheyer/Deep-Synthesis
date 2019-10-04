# Downloads model file to translation_function directory
# Model file must be in directory before creating the Lambda function

import boto3
from botocore import UNSIGNED
from botocore.client import Config
from pathlib import Path

s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))

bucket_name = 'deepsynthesis'
model_file = 'public_models/deep_synthesis.pt'
destination_file = 'translation_function/deep_synthesis.pt'

path = Path('.')
file_path = path/destination_file 

if not file_path.exists():
    print('Downloading Model File')
    (path/'models').mkdir(exist_ok=True)
    s3.download_file(bucket_name, model_file, destination_file)
else:
    print('Model already downloaded')
from boto3.session import Session
import boto3
import os
from pathlib import Path
import json

credentials = json.load(open('configs/credentials.json'))
access_key = credentials['access_key']
secret_key = credentials['secret_key']

bucket_name = 'deepsynthesis'
model_file = 'Molecular_transformer_step_400000.pt'

path = Path('.')
file_path = path/'models'/model_file 

if not file_path.exists():
    print('Downloading Model File')
    (path/'models').mkdir(exist_ok=True)
    session = Session(aws_access_key_id=access_key,
                    aws_secret_access_key=secret_key)

    s3 = session.resource('s3')
    s3.Bucket(bucket_name).download_file(model_file, 'models/Molecular_transformer_step_400000.pt')

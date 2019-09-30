'''
app.py

This script serves as the handler function for AWS Lambda

The code here is based off a repo written by Matt McClean 
https://github.com/mattmcclean/sam-pytorch-example

The lambda_handler function takes in an event as input,
runs predictions on the content of the event, and returns
the predictions as a base64 encoded ascii string

The input event should be structured as follows:

{
    "beam" : int,
    "n_best" : int
    "return_attention" : bool
    "data" : Tokenized SMILES data dumped to a JSON string,
    "warmup" : bool
}

This function is designed to import dependencies from a Lambda layer.
Loading in the layer is done with the unzip_requirements function
'''
  
# this import statement is needed if you want to use an AWS Lambda Layer to load dependencies
# it unzips all of the pytorch & dependency packages when the script is loaded to 
# avoid the 250 MB unpacked limit in AWS Lambda

try:
    import unzip_requirements
except ImportError:
    pass

import logging
import time
import boto3 
import requests
import json
import numpy as np
import pickle
import pdb
import gzip
import os
from io import BytesIO
import base64

from translate import *


logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Load model config, load model from config
model_description = json.load(open('model_config.json'))
logger.info('Loading model')
Translator = TranslationModel(model_description)

# Handles incoming events to lambda function
def lambda_handler(event, context):

    print("Starting event")
    logger.info(event)

    # allow ping for warmup
    if event['warmup']:
        response = {'warmup' : 'confirmed'}
    else:
        response = run_translation(event)

    print("Returning response")
    return {
        'isBase64Encoded': True,
        "statusCode": 200,
        'headers': {
            'Content-Type': 'application/json',
            'Content-Encoding': 'gzip',
            'Access-Control-Allow-Origin': '*'
        },
        "body": gzip_b64encode(response)
    }


def run_translation(event):
    # processes event info, runs translation
    beam = event['beam']
    n_best = event['n_best']
    return_attention = event['return_attention']
    data = json.loads(event['data'])

    print("Starting Prediction")
    predictions = predict(data, beam, n_best, return_attention)
    return predictions

def predict(data, beam, n_best, return_attention):
    # function runs predictions and processes outputs
    logger.info("Starting Prediction")
    logger.info(f"Prediction Scope: {len(data)} examples + {beam} beam width")
    start_time = time.time()
    scores, preds, attns = Translator.run_translation(data, 
                                        beam=beam, n_best=n_best, return_attention=return_attention)
    logger.info("--- Inference time: %s seconds ---" % (time.time() - start_time))
    response = {}
    response['predictions'] = preds 
    response['scores'] = process_scores(scores)
    response['attention'] = process_attention(attns)

    return response

def gzip_b64encode(data):
    compressed = BytesIO()
    with gzip.GzipFile(fileobj=compressed, mode='w') as f:
        json_response = json.dumps(data)
        f.write(json_response.encode('utf-8'))
    return base64.b64encode(compressed.getvalue()).decode('ascii')

def process_scores(scores):
    # scores are torch tensors, must be cast to float
    new_scores = []

    for score_set in scores:
        float_scores = [i.item() for i in score_set]
        new_scores += [float_scores]
    
    return new_scores

def process_attention(attns):
    # attention maps are numpy arrays.
    # arrays are converted to lists prior to compression
    new_attentions = []
    
    for attention_set in attns:
        list_attentions = [np.around(i, decimals=4).tolist() for i in attention_set]
        new_attentions += [list_attentions]
        
    return new_attentions

import streamlit as st 
import math
import boto3
import json
import pickle
import gzip
import numpy as np

from lambda_async import *
from postprocess import *


class LambdaInterface():
    '''
    Class to prepare data and interface with AWS Lambda function

    Data is batched in two steps.
    First the data is batched into chunks of size "payload",
    where payload is the number of samples sent to a lambda
    function during a single invocation.

    Payload sizes are hard coded based on experimenting with
    the performance of an AWS Lambda 2048 MB function

    Second, payload chunks are batched into invocation chunks
    of size "fan_size", where "fan_size" is the number of 
    concurrent lambda invocations. 
    
    Currently, "fan_size" is set to be no more than 4

    Once data is batched, running the "predict_async" function
    will use the functions in the lambda_async.py file to 
    invoke all invocation_chunks simultaneously
    '''
    def __init__(self, data, beam, n_best, attention, function, bucket):
        self.data = data
        self.beam = beam
        self.n_best = n_best
        self.attention = attention
        self.function = function
        self.bucket = bucket
        self.s3 = boto3.client('s3')
                
        self.chunksize_dict = {
                                1 : 60,
                                2 : 50,
                                3 : 40,
                                5 : 30
                              }

        self.payload_size = self.chunksize_dict[beam]

        # concurrent invocations will be up to 4
        self.fan_size = min(4, math.ceil(len(data)/self.payload_size))
        
        self.region = 'us-west-2'
        
    def chunk_data(self, data, chunksize):
        # partitions a list into chunks of size chunksize
        return [data[i:i+chunksize] for i in range(0, len(data), chunksize)]
    
    def data_to_payload(self, data):
        # writes a list of data payloads into a json format compatible
        # with the functions in lambda_async.py
        requests = []

        for data_item in data:

            data_dict = {
                            "function_name" : self.function,
                            "payload" : {"data" : json.dumps(data_item),
                                         "beam" : self.beam,
                                         "n_best" : self.n_best,
                                         "return_attention" : self.attention}
                        }

            requests.append(data_dict)

        return requests
    
    def process_response(self, response):
        # Due to the size of the response (mostly the attention maps)
        # responses are not returned directly. They are stored in a s3 bucket
        # as gzipped pickle files, and the item key is returned in the response
        s3_key = response['body']['s3_key']
        outputs = self.s3.get_object(Bucket=self.bucket, Key=s3_key)
        output_dict = pickle.loads(gzip.decompress(outputs['Body'].read()))
        
        predictions = output_dict['predictions']
        scores = output_dict['scores']
        attention = output_dict['attention']

        # delete item from bucket after it is retrieved 
        self.s3.delete_object(Bucket=self.bucket, Key=s3_key)
        
        return predictions, scores, attention
    
    def reconstruct_output(self, output):
        # Raw invocation outputs are a list of json responses 
        # containing s3 item keys pointing to the actual prediction outputs

        # This function retrieves those outputs and combines them into 
        # a single set of lists
        output = [self.process_response(i) for i in output]
    
        preds = []
        scores = []
        attns = []
        for out in output:
            preds_iter, scores_iter, attns_iter = out
            preds += preds_iter
            scores += scores_iter
            attns += attns_iter

        return preds, scores, attns

    def reconstruct_attention(self, attns):
        # Attention maps are loaded as lists of lists.
        # This function reconstructs them into numpy arrays
        reconstructed_attentions = []

        for attention_set in attns:
            array_attentions = [np.array(i) for i in attention_set]
            reconstructed_attentions += [array_attentions]

        return reconstructed_attentions
    
    def predict_async(self):
        # runs async prediction

        # break data into chunks of size payload_size
        chunked_data = self.chunk_data(self.data, self.payload_size)

        # process payload chunks into json format
        payload_chunks = self.data_to_payload(chunked_data)

        # break payloads into invocation chunks
        # each chunk contains fan_size payload chunks
        invocation_chunks = self.chunk_data(payload_chunks, self.fan_size)
        
        results = []
        # async prediction on each invocation chunk
        for invocation_chunk in invocation_chunks:
            results += invoke_all(requests=invocation_chunk, region=self.region)
            
        # reconstruct outputs
        predictions, scores, attentions = self.reconstruct_output(results)
        attentions = self.reconstruct_attention(attentions)
        return predictions, scores, attentions


@st.cache(ignore_hash=True)
def translate_data(smile_data, beam, n_best, attention, model_description):
    function = model_description['function']
    bucket = model_description['bucket']
    lambda_interface = LambdaInterface(smile_data.smiles_tokens, beam, n_best, attention, 
                                        function, bucket)

    preds, scores, attns = lambda_interface.predict_async()
    prediction = Predictions(smile_data, preds, scores, attns)
    return prediction
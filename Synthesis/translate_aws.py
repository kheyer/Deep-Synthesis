import streamlit as st 
import math
import boto3
import json
import pickle
import gzip
import numpy as np
import base64

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

    Once data is batched, running the "run_translation" function
    will use the functions in the lambda_async.py file to 
    invoke all invocation_chunks simultaneously
    '''
    def __init__(self, config):
        self.function = config['function']
        self.bucket = config['bucket']
        self.s3 = boto3.client('s3')
        self.fan_size = config['fan_size']
        self.chunksize_dict = {
                                1 : 100,
                                2 : 80,
                                3 : 60,
                                5 : 50
                              }
        self.region = 'us-west-2'
        
    def chunk_data(self, data, chunksize):
        # partitions a list into chunks of size chunksize
        return [data[i:i+chunksize] for i in range(0, len(data), chunksize)]
    
    def data_to_payload(self, data, beam, n_best, attention, warmup=False):
        # writes a list of data payloads into a json format compatible
        # with the functions in lambda_async.py
        requests = []

        for data_item in data:

            data_dict = {
                            "function_name" : self.function,
                            "payload" : {"data" : json.dumps(data_item),
                                         "beam" : beam,
                                         "n_best" : n_best,
                                         "return_attention" : attention,
                                         "warmup" : warmup}
                        }

            requests.append(data_dict)

        return requests
    
    def process_response(self, response):
        # Lambda response is converted to a json string, gzipped and base64 encoded.
        # This function reconstructs the response by reversing that process
        output_dict = json.loads(gzip.decompress(base64.b64decode(response['body'])))

        predictions = output_dict['predictions']
        scores = output_dict['scores']
        attention = output_dict['attention']

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
        if attns[0] is not None:
            reconstructed_attentions = []

            for attention_set in attns:
                array_attentions = [np.array(i) for i in attention_set]
                reconstructed_attentions += [array_attentions]
        else:
            reconstructed_attentions = attns

        return reconstructed_attentions
    
    def run_translation(self, data, beam, n_best, return_attention):
        # runs async prediction

        fan_size = self.fan_size
        payload_size = min(self.chunksize_dict[beam], math.ceil(len(data)/fan_size))

        # break data into chunks of size payload_size
        chunked_data = self.chunk_data(data, payload_size) #self.chunk_data(self.data, self.payload_size)

        # process payload chunks into json format
        payload_chunks = self.data_to_payload(chunked_data, beam, n_best, return_attention)

        # break payloads into invocation chunks
        # each chunk contains fan_size payload chunks
        invocation_chunks = self.chunk_data(payload_chunks, fan_size) #self.chunk_data(payload_chunks, self.fan_size)
        
        results = []
        # async prediction on each invocation chunk
        for invocation_chunk in invocation_chunks:
            results += invoke_all(requests=invocation_chunk, region=self.region)
            
        # reconstruct outputs
        predictions, scores, attentions = self.reconstruct_output(results)
        attentions = self.reconstruct_attention(attentions)
        return scores, predictions, attentions
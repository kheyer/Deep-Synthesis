'''
Non-blocking concurrent requests to AWS Lambda API endpoints
with Python asyncio and aiohttp

Code adapted from gist by Mathew Marcus
https://gist.github.com/byrro/439d233636ac6894f87c1df6ff1fa298

    Credits to Mathew Marcus (2019)
    https://www.mathewmarcus.com/blog/
    http://archive.is/nXkCb
'''

import asyncio
import json
import os
from typing import Dict, List
import urllib
import streamlit as st 
import time
import aiohttp
import boto3
from botocore import session, awsrequest, auth
from session_id import *

AWS_CREDENTIALS = session.Session().get_credentials()

def sign_headers(*, url: str, payload: Dict):
    '''Sign AWS API request headers'''
    segments = urllib.parse.urlparse(url).netloc.split('.')
    service = segments[0]
    region = segments[1]

    request = awsrequest.AWSRequest(
        method='POST',
        url=url,
        data=json.dumps(payload),
    )

    auth.SigV4Auth(AWS_CREDENTIALS, service, region).add_auth(request)

    return dict(request.headers.items())


async def invoke(*, url: str, payload: Dict, session):
    '''Invoke a Lambda function'''
    signed_headers = sign_headers(url=url, payload=payload)

    async with session.post(url, json=payload, headers=signed_headers) \
            as response:
        return await response.json()


def run_invocations(
        *,
        requests: List,
        base_url: str,
        session,
        is_async: bool,
        ):
    for request in requests:
        method = 'invoke-async' if is_async is True else 'invocations'

        url = os.path.join(base_url, request['function_name'], method)

        yield invoke(url=url, payload=request['payload'], session=session)


def invoke_all(
        *,
        requests: List,
        region: str = 'us-east-1',
        is_async: bool = False,
        ):
    
    base_url = f'https://lambda.{region}.amazonaws.com/2015-03-31/functions'

    async def wrapper():
        async with aiohttp.ClientSession(raise_for_status=True) as session:
            invocations = run_invocations(
                requests=requests,
                base_url=base_url,
                session=session,
                is_async=is_async,
            )

            return await asyncio.gather(*invocations)

    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)

    results = loop.run_until_complete(wrapper())

    return results


def background(f):
    from functools import wraps
    @wraps(f)
    def wrapped(*args, **kwargs):
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        if callable(f):
            return loop.run_in_executor(None, f, *args, **kwargs)
        else:
            raise TypeError('Task must be a callable')    
    return wrapped


@background
def warmup(payload, function):
    client = boto3.client('lambda', 'us-west-2')
    print('sending request')
    client.invoke(FunctionName=function,
                             InvocationType='RequestResponse',
                             Payload=json.dumps(payload))
    print("warmup completed")

# Within a given session, ping will fire every 120 seconds
@fancy_cache(unique_to_session=True, ttl=120)
def warmup_lambda(fan_size, function, seed=None):
    payload = { 'data': '""',
                'beam': '',
                'n_best': '',
                'return_attention': '',
                'warmup': True}
    for i in range(fan_size):
        print('starting warmup')
        warmup(payload, function)
        time.sleep(0.1)
    return 'warmup complete'
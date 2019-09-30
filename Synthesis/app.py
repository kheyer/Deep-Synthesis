import streamlit as st 
import argparse
import json
from utils import *
from preprocess import *
from translate_aws import *
from lambda_async import *

# Importing translate also imports OpenNMT.
# Putting this in a try/except block allows for deploying
# the app configured to run AWS Lambda predictions
# without bringing along OpenNMT
try:
    from translate import *
except ModuleNotFoundError:
    pass

# define runtime specifications on launch
parser = argparse.ArgumentParser()
parser.add_argument("runtime", help="Argument to select local or AWS runtimes ['local' or 'AWS']")
args = parser.parse_args()

# Dictionary of configuration files
config_json_opts = {
                      'AWS' : 'configs/aws_description.json',
                      'local' : 'configs/model_description.json'
                   }

runtime = args.runtime 
model_description = json.load(open(config_json_opts[runtime]))

# Define class to handle translations based on runtime
if runtime == 'local':
    translator_class = TranslationModel
elif runtime == 'AWS':
    translator_class = LambdaInterface

    # Asyncronously ping fan_size instances concurrently
    # to prevent cold start
    warmup_lambda(model_description['fan_size'], model_description['function'])
    
else:
    raise ValueError('''Please provide a valid runtime argument. Use 'local' to run predictions locally, or 
                'AWS' to run predictions on AWS (AWS inference requires permissions)''')


st.title('Deep Synthesis')

input_options = app_setup()

prediction_option = st.selectbox('Select an Input Format', input_options)

single_predict, source_param, target_param = get_data_params(input_options, prediction_option)

data_box = st.checkbox('Load Data')

if data_box:
    smile_data = load_data(single_predict, source_param, target_param)
    display_idx = display_slider(smile_data)

    display_data(smile_data, display_idx)

    st.write('Input Prediction Parameters')
    beam = int(st.selectbox('Select Beam Width', [1,2,3,5]))
    n_best = int(st.selectbox('Select Top K Predictions', [1,2,3,5]))

    prediction_box = st.checkbox('Run Prediction')

    if prediction_box:
        placeholder = st.empty()
        placeholder.text('Translation in Progress')
        prediction = translate_data(smile_data, beam, n_best, True, translator_class, model_description)
        placeholder.text('Translation Complete')

        display_prediction(prediction, display_idx)
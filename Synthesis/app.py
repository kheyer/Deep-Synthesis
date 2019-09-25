import streamlit as st 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("runtime", help="Argument to select local or AWS runtimes ['local' or 'AWS']")
args = parser.parse_args()

import json
from utils import *
from preprocess import *

if args.runtime == 'local':
    from translate import *
    model_description = json.load(open('configs/model_description.json'))
elif args.runtime == 'AWS':
    from translate_aws import *
    model_description = json.load(open('configs/aws_description.json'))
else:
    raise ValueError('''Please provide a valid runtime argument. Use 'local' to run predictions locally, or 
                'AWS' to run predictions on AWS (AWS inference requires permissions)''')

st.title('Deep Synthesis')

input_options = app_setup()

prediction_option = st.selectbox('Select an Input Format', input_options)

single_predict, source_param, target_param = get_data_params(input_options, prediction_option)

if st.checkbox('Load Data'):
    smile_data = load_data(single_predict, source_param, target_param)

    display_data(smile_data)

    st.write('Input Prediction Parameters')
    beam = int(st.selectbox('Select Beam Width', [1,2,3,5]))
    n_best = int(st.selectbox('Select Top K Predictions', [1,2,3,5]))

    if st.checkbox('Run Prediction'):
        placeholder = st.empty()
        placeholder.text('Translation in Progress')
        prediction = translate_data(smile_data, beam, n_best, True, model_description)
        placeholder.text('Translation Complete')

        display_prediction(prediction)
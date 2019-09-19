import streamlit as st 
import json

from utils import *
from preprocess import *

st.title('Deep Synthesis')

input_options = app_setup()
model_description = json.load(open('configs/model_description.json'))

prediction_option = st.selectbox('Select an Input Format', input_options)

single_predict, source_param, target_param = get_data_params(input_options, prediction_option)

if st.checkbox('Load Data'):
    smile_data = load_data(single_predict, source_param, target_param)

    display_data(smile_data)

    beam = int(st.text_input('Beam Width', '10'))
    n_best = int(st.text_input('Top K Predictions', '5'))

    if st.checkbox('Run Prediction'):
        placeholder = st.empty()
        placeholder.text('Translation in Progress')
        prediction = translate_data(smile_data, beam, n_best, model_description)
        placeholder.text('Translation Complete')
import streamlit as st 

from utils import *
from preprocess import *

st.title('Deep Synthesis')

input_options = app_setup()

prediction_option = st.selectbox('Select an Input Format', input_options)

single_predict, source_param, target_param = get_data_params(input_options, prediction_option)

if st.checkbox('Load Data'):
    smile_data = load_data(single_predict, source_param, target_param)

    display_data(smile_data)
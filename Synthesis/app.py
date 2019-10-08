import streamlit as st 
import argparse
import json
from utils import *
from preprocess import *

# define runtime specifications on launch
parser = argparse.ArgumentParser()
parser.add_argument("runtime", help="Argument to select local or AWS runtimes ['local' or 'AWS']")
args = parser.parse_args()

st.title('Deep Synthesis')

# Set up application depending on runtime environment and user specifications 
model_description, translator_class, single_predict, source_param, target_param = app_setup(args)

smile_data = load_data(single_predict, source_param, target_param)

# create slider to display data
display_idx = display_slider(smile_data)
display_data(smile_data, display_idx)

beam, n_best = prediction_params(single_predict)

st.text('Everything look right? Press the "Predict Products" button to generate predictions')
prediction = translate_data(smile_data, beam, n_best, True, translator_class, 
                            model_description)

display_prediction(prediction, display_idx)

download_data(single_predict, prediction)
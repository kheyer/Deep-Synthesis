import streamlit as st 
from utils import *

st.title('Deep Synthesis')

base_smile, input_options = app_setup()

prediction_option = st.selectbox('Select an Input Format', input_options)
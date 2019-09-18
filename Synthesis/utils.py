import streamlit as st 
import os

def app_setup():
    # starter values for prediction type 
    base_smile = 'O=Cc1cncc(Cl)c1COC1CCCCO1.OCc1c(Cl)cncc1Cl'
    input_options = ['Predict from String', 'Predict from File']

    return (base_smile, input_options) 
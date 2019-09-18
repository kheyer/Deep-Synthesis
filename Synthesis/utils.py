import streamlit as st 
import os

def app_setup():
    # starter values for prediction type 
    input_options = ['Predict from String', 'Predict from File']

    return input_options

def get_data_params(input_options, prediction_options):
    # function determines if prediction will be run on a string input by the user
    # or from a file of SMILES strings
    if prediction_options == input_options[0]:
        single_predict = True
    else:
        single_predict = False 

    if single_predict:
        # If single predict, create a text box for user input
        base_smile = 'O=Cc1cncc(Cl)c1COC1CCCCO1.OCc1c(Cl)cncc1Cl'
        smile = st.text_input('Input Smile', base_smile)
        target_smile = st.text_input('Known Target?' , '')
        # returns single_predict (bool), smile (string), target_smile (string)
        return single_predict, smile, target_smile
    else:
        # If file predict, create text boxes that let users navigate to the prediction file
        source_filename, target_filename = get_filenames()
        # returns single_predict (bool), source_filename (string to .txt filename)
        # target_filename (string to .txt filename if desired, else None)
        return single_predict, source_filename, target_filename
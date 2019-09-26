import streamlit as st 
import os
from rdkit import Chem
from rdkit.Chem import Draw
import logging
import time

from preprocess import *
from postprocess import *

logger = logging.getLogger()
logger.setLevel(logging.INFO)

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
        smile = st.text_input('Input Source SMILES (Required)', base_smile)
        target_smile = st.text_input('Input Target SMILES (Optional)' , '')
        # returns single_predict (bool), smile (string), target_smile (string)
        return single_predict, smile, target_smile
    else:
        # If file predict, create text boxes that let users navigate to the prediction file
        source_filename, target_filename = get_filenames()
        # returns single_predict (bool), source_filename (string to .txt filename)
        # target_filename (string to .txt filename if desired, else None)
        return single_predict, source_filename, target_filename

def get_filenames(path='data'):
    # Creates interface for user to select a file to load
    # file must be stored locally in a directory accessable from the streamlit app
    folder = st.text_input('Data Folder Path', path)
    source_filename = file_selector(folder_path=folder, txt='Select source file')

    if st.checkbox('Target File?'):
        # optional check box for a file of targets (not required)
        target_folder = st.text_input('Data Folder Path', folder)
        target_filename = file_selector(folder_path=target_folder, txt='Select target file')
    else:
        target_filename = None

    return source_filename, target_filename


def file_selector(folder_path='.', txt='Select a file'):
    # uses Streamlit string inputs to navigate to a file
    filenames = os.listdir(folder_path)
    selected_filename = st.selectbox(txt, filenames)
    return os.path.join(folder_path, selected_filename)

@st.cache
def load_data(single_predict, source_param, target_param):
    # triggers actually loading the data
    # loaded data is cached
    if single_predict:
        # If single_predict, load with the method for a single entry
        data = SmilesData.single_entry(source_param, target_param)
    else:
        # Else, load from file
        data = SmilesData.file_entry(source_param, target_param)

    # Returns a SilesData object
    return data 

def display_data(smile_data):
    # displays data held in a SmilesData object
    if len(smile_data) > 1:
        # If more than one SMILES is present, a streamlit slider bar is created
        display_idx = st.slider('Display Index', 0, len(smile_data)-1, 0)
        st.image(smile_data.display(idx=display_idx, img_size=(300,300)))
    else:
        st.image(smile_data.display(img_size=(300,300)))

@st.cache(ignore_hash=True)
def translate_data(smile_data, beam, n_best, attention, translator_class, model_description):
    # Important note: translator class must be instantiated within this function for 
    # Streamlit caching to work properly
    start = time.time()
    translator = translator_class(model_description)
    scores, preds, attns = translator.run_translation(smile_data.smiles_tokens, 
                                                beam=beam, n_best=n_best, return_attention=attention)
    logger.info(f'Inference Time: {time.time() - start}')
    prediction = Predictions(smile_data, preds, scores, attns)
    return prediction

@st.cache
def plot_topk(prediction_tokens, legend, img_size=(400,400)):
    mols = [Chem.MolFromSmiles(process_prediction(i)) for i in prediction_tokens]
    return Draw.MolsToGridImage(mols, legends=legend, subImgSize=img_size)

def display_prediction(prediction):
    #if prediction:
    if len(prediction) > 1:
        prediction_idx = st.slider('Prediction Index', 0, len(prediction)-1, 0)
    else:
        prediction_idx = 0

    prediction_data = display_parameters(prediction, idx=prediction_idx)

    st.write(f'Top {prediction.top_k} Predictions')
    st.image(plot_topk([i.prediction_tokens for i in prediction_data],
                        [i.legend for i in prediction_data], img_size=(300,300)))

    if len(prediction_data) > 1:
        view_idx = st.slider('View Prediction (In Order of Model Confidence)', 0, len(prediction_data)-1, 0)
    else:
        view_idx = 0

    current_prediction = prediction_data[view_idx]
    im, attn_plot = plot_prediction(current_prediction.source_tokens,
                                    current_prediction.prediction_tokens,
                                    current_prediction.attention,
                                    current_prediction.legend,
                                    img_size=(300,300))

    st.write(f'Predicted Smile: {process_prediction(current_prediction.prediction_tokens)}')
    if im:
        st.image(im)
    st.pyplot(plt.show(attn_plot), bbox_inches = 'tight', pad_inches = 0)

    st.write('\nPrediction Dataframe')
    st.dataframe(prediction.sample_df(prediction_idx))
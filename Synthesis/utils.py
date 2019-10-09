import streamlit as st 
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import logging
import time
import json
import numpy as np 
from preprocess import *
from postprocess import *
from confirm_button import *
from translate_aws import *
from lambda_async import *
from session_id import *
from landing_page import *

# Importing translate also imports OpenNMT.
# Putting this in a try/except block allows for deploying
# the app configured to run AWS Lambda predictions
# without bringing along OpenNMT
try:
    from translate import *
except ModuleNotFoundError:
    pass

logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Sample smiles for users to predict on
example_smiles = [
    'CCN.Oc1cccc2ccccc12.[Na+].[OH-].c1ccc(OP(Oc2ccccc2)Oc2ccccc2)cc1',
    'CN(C)C=O.Cn1c(=O)n(CCCI)c2ccccc21.O=C([O-])[O-].O=c1[nH]c2cc(Cl)ccc2n1C1CCNCC1.[Na+].[Na+]',
    'C1CNC1.CS(=O)(=O)OCCC#Cc1ccc2c(-c3ccc(Br)cc3)nsc2c1.O=C([O-])[O-].[Na+].[Na+]',
    'CCN(CC)CC.CO.COC(=O)c1sc(Cl)c(Cl)c1NC(C)=O.[H][H]',
    'C1CCOC1.CCC(C)(C)Cc1cn(S(=O)(=O)N(C)C)c(C(=O)Cc2ccc(Br)cc2)n1.C[Mg+].[Br-].[Ce+3].[Cl-].[Cl-].[Cl-]',
    'CC1CO1.CCCBr.CCCCCC.CCN(CC)CC.CN(C)P(=O)(N(C)C)N(C)C.COc1cc([N+](=O)[O-])c2ncccc2c1O.O',
    'C1CCOC1.CCOC(=O)CCCCCOc1c(C=CC(=O)CCc2ccccc2)occc1=O.CC[BH-](CC)CC.[Li+]',
    'C#Cc1ccccc1.CCN(CC)CC.CCN(CC)CCOc1cc(Cl)c(I)cc1Cl.[Cu]I.[Pd]',
    'Cl.I.O.O=C(c1ccccc1)c1ccc(O)cc1.[I-].[K+].[NH4+].[OH-]',
    'O=Cc1cncc(Cl)c1COC1CCCCO1.OCc1c(Cl)cncc1Cl'
]

def app_setup(args):
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

        input_options = ['Welcome to Deep Synthesis', 'How it Works', 'Deep Synthesis Tutorial',
                         'Predict from String', 'Predict from File']
        option_output = st.sidebar.selectbox('Select Page', input_options)

    elif runtime == 'AWS':
        translator_class = LambdaInterface

        # Asyncronously ping fan_size instances concurrently
        # to prevent cold start
        warmup_lambda(model_description['fan_size'], model_description['function'])
        input_options = ['Welcome to Deep Synthesis', 'How it Works',
                         'Deep Synthesis Tutorial', 'Predict from String']
        option_output = st.sidebar.selectbox('Select Page', input_options)
        
    else:
        raise ValueError('''Please provide a valid runtime argument. Use 'local' to run predictions locally, or 
                    'AWS' to run predictions on AWS (AWS inference requires permissions)''')

    # get data params during setup
    single_predict, smile, target = get_data_params(option_output, runtime)

    return model_description, translator_class, single_predict, smile, target

def get_data_params(prediction_options, runtime):
    # function determines if prediction will be run on a string input by the user
    # or from a file of SMILES strings
    if prediction_options == 'Welcome to Deep Synthesis':
        landing_page()

        return None, None, None

    elif prediction_options == 'How it Works':
        explanation_page()

        return None, None, None

    elif prediction_options == 'Deep Synthesis Tutorial':
        tutorial_page(runtime)

        return None, None, None

    elif prediction_options == 'Predict from String':
        single_predict = True
        # If single predict, create a text box for user input
        # seed UI with sample reaction
        sample_rxn = st.sidebar.selectbox('Choose sample reaction', example_smiles,
                                     format_func=format_smiles)

        message = 'Enter a SMILES string of reactants, or choose an example reaction \
                    from the sidebar on the left.\nThen press the "Load Data" button.'
        smile = st.text_input(message, sample_rxn)
        target_smile = None
        # returns single_predict (bool), smile (string), target_smile (string)
        return single_predict, smile, target_smile

    else:
        single_predict = False 
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

@cache_on_button_press('Load Data')
def load_data(single_predict, source_param, target_param):
    # triggers actually loading the data
    # loaded data is cached
    if single_predict:
        # If single_predict, load with the method for a single entry
        data = SmilesData.single_entry(source_param, target_param)
    else:
        # Else, load from file
        data = SmilesData.file_entry(source_param, target_param)

    # Returns a SimlesData object
    return data 

def display_data(smile_data, display_idx):
    # displays data held in a SmilesData object
    return st.image(smile_data.display(idx=display_idx, img_size=(500,500)), use_column_width=True)

def display_slider(data):
    if len(data) > 1:
        display_idx = st.sidebar.slider('Display Index', 0, len(data)-1, 0)
    else:
        display_idx = 0
    
    return display_idx

@cache_on_button_press('Predict Products')
def translate_data(smile_data, beam, n_best, attention, translator_class, 
                    model_description):
    # Important note: translator class must be instantiated within this function for 
    # Streamlit caching to work properly
    placeholder = st.empty()
    placeholder.text('Prediction in Progress')

    start = time.time()
    translator = translator_class(model_description)
    scores, preds, attns = translator.run_translation(smile_data.smiles_tokens, 
                                                beam=beam, n_best=n_best, return_attention=attention)
    prediction_time = time.time() - start
    prediction = Predictions(smile_data, preds, scores, attns)
    logger.info(f'Inference Time: {prediction_time}')

    placeholder.text('Prediction Complete')
    return prediction

@fancy_cache(unique_to_session=True, ttl=3600)
def plot_topk(prediction_tokens, legend, img_size=(400,400)):
    mols = [Chem.MolFromSmiles(process_prediction(i)) for i in prediction_tokens]
    return Draw.MolsToGridImage(mols, legends=legend, subImgSize=img_size)

def display_prediction(prediction, display_idx):
    prediction_data = display_parameters(prediction, idx=display_idx)

    st.write(f'Top {prediction.top_k} Predictions')
    st.image(plot_topk([i.prediction_tokens for i in prediction_data],
                        [i.legend for i in prediction_data], img_size=(400,400)), use_column_width=True)

    if len(prediction_data) > 1:
        view_idx = st.slider('View Prediction (In Order of Model Confidence)', 0, len(prediction_data)-1, 0)
    else:
        view_idx = 0

    current_prediction = prediction_data[view_idx]
    im, attn_plot = plot_prediction(current_prediction.source_tokens,
                                    current_prediction.prediction_tokens,
                                    current_prediction.attention,
                                    current_prediction.legend,
                                    img_size=(200,400))

    st.write(current_prediction.legend)
    if im:
        st.image(im, use_column_width=True)
    st.write(f'Predicted Smile: {process_prediction(current_prediction.prediction_tokens)}')
    full_rxn = process_prediction(current_prediction.source_tokens) + '>>' + process_prediction(current_prediction.prediction_tokens)
    st.write(f'Full Predicted Reaction: {full_rxn}')
    st.write('')
    st.pyplot(plt.show(attn_plot), bbox_inches = 'tight', pad_inches = 0)

    st.write('\nPrediction Dataframe')
    st.dataframe(prediction.sample_df(display_idx))

@cache_on_button_press('Save Prediction Data')
def save_data(predictions, save_folder):
    # button to save data to some input file path
    predictions.df.to_csv(save_folder, index=False)
    st.text(f'Prediction data saved to {save_folder}')

def download_data(single_predict, predictions):
    # text inpout box for path to save prediction data
    if not single_predict:
        save_folder = st.text_input('Prediction Save Destination', 'data/predictions.csv')
        save_data(predictions, save_folder)

def prediction_params(single_predict):
    # get beam search parameters depending on prediction type
    # beam/n_best parameters are only exposed for bulk prediction
    if single_predict:
        beam = 5
        n_best = 5
    else:
        st.write('Input Prediction Parameters')
        beam = int(st.selectbox('Select Beam Width', [1,2,3,5]))
        n_best = int(st.selectbox('Select Top K Predictions', [1,2,3,5]))

    return beam, n_best

def format_smiles(smile):
    # string formatting function for streamlit sidebar
    if len(smile) > 20:
        smile = smile[:20] + '...'
    return smile
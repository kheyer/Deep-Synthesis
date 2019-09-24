# Deep Synthesis
A deep learning framework for predicting chemical synthesis

## Local Setup

Follow these instructions to set up Deep Synthesis on your local machine

Clone the repository

    git clone https://github.com/kheyer/Deep-Synthesis

Run `buildEnv.sh`. This shell script will create a conda environment named `deep_synthesis`, install all required packages and download the trained model

    cd Deep-Synthesis
    chmod +x buildEnv.sh
    ./buildEnv.sh

Follow the package installation prompts.

Once package installation is complete, activate the new conda environment and start the Streamlit app

    source activate deep_synthesis`
    streamlit run Synthesis/app.py`

The Streamlit app is now running locally on port 8501

## Local Inference

Once the Streamlit app is running locally, inference can be run on local hardware. There are two ways to run inference. The "Predict from String" option predicts on a single SMILES string entered into the text entry box.

<img src="https://github.com/kheyer/Deep-Synthesis/blob/main_setup_branch/media/local_prediction_1.png" width="500" alt="local prediction from string">

Alternatively, the "Predict from File" option predicts in bulk from a text file containing multiple SMILES inputs. The inputs should be formatted with one set of reactants per line. See the example files in the `/data` directory. The "Predict from File" option expects the source file to be stored locally. The file path to the source file defaults to the `/data` directory, but can be pointed at any file path relative to the repo file.

<img src="https://github.com/kheyer/Deep-Synthesis/blob/main_setup_branch/media/local_prediction_2.png" width="500" alt="local prediction from string">

For both prediction formats, known target product SMILES can optionally be provided. If target SMILES are provided, predictions are automatically scored.

## Requisites

#### Dependencies

#### Installation

## Build Environment

## Configs

## Test

## Run Inference

## Build Model

## Serve Model

## Analysis

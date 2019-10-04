# Deep Synthesis
A deep learning framework for predicting chemical synthesis

# Running Locally

The following instructions detail how to set up the environment and run the application locally 

## Local Setup

Follow these instructions to set up Deep Synthesis on your local machine

Clone the repository

    git clone https://github.com/kheyer/Deep-Synthesis

Run `buildEnv.sh`. This shell script will create a conda environment named `deep_synthesis`, install all required packages and download the trained model

    cd Deep-Synthesis
    chmod +x build_local/buildEnv.sh
    ./build_local/buildEnv.sh

Follow the package installation prompts.

Once package installation is complete, activate the new conda environment and start the Streamlit app

    source activate deep_synthesis
    streamlit run Synthesis/app.py local

The Streamlit app is now running locally on port 8501

## Local Docker Setup

To build a Docker container running the application locally, run the following:

Clone the repository

    git clone https://github.com/kheyer/Deep-Synthesis
    
Build the container

    cd Deep-Synthesis
    docker build -f build_local/local.Dockerfile -t deep_synthesis .
    docker run -d -p 8501:8501 deep_synthesis

The Streamlit app is now running locally on port 8501


## Local Inference via Streamlit App

Once the Streamlit app is running locally, the application can be accessed at `localhost:8501`. Inference will run on local hardware. There are two ways to run inference. The "Predict from String" option predicts on a single SMILES string entered into the text entry box.

<img src="https://github.com/kheyer/Deep-Synthesis/blob/main_setup_branch/media/local_prediction_1.png" width="500" alt="local prediction from string">

Alternatively, the "Predict from File" option predicts in bulk from a text file containing multiple SMILES inputs. The inputs should be formatted with one set of reactants per line. See the example files in the `/data` directory. The "Predict from File" option expects the source file to be stored locally. The file path to the source file defaults to the `/data` directory, but can be pointed at any file path relative to the repo file.

<img src="https://github.com/kheyer/Deep-Synthesis/blob/main_setup_branch/media/local_prediction_2.png" width="500" alt="local prediction from file">

For both prediction formats, known target product SMILES can optionally be provided. If target SMILES are provided, predictions are automatically scored.

Once predictions are generated, they can be visualized using the slider bar. Source-to-prediction attention maps are also generated.

<img src="https://github.com/kheyer/Deep-Synthesis/blob/main_setup_branch/media/local_prediction_3.png" width="500" alt="local prediction evaluation">

## Local Inference

Local inference can also be run from the command line via the `translate_cli.py` file.

    python Synthesis/translate_cli.py --runtime local --source_file [source.txt] --destination_file [output_file.csv] --target_file [path/to/target.txt] --beam [beam] --n_best [n_best]
    
* `runtime` - determines if inference will be run locally or on the cloud. Pass `local` for local inference
* `source_file` - path to the file of source sequences. The source file should be a `.txt` file of SMILES inputs with one SMILES input per line.
* `destination_file` - path to the destination file where predictions will be saved. This should be a `.csv` file. Predictions are saved as a csv that includes source sequences and prediction scores in addition to predictions.
* `target_file` (optional) - path to a file of target sequences if available. If targets are provided, model predictions will be scored against the targets.
* `beam` (optional) - beam width for translation. Defaults to 5.
* `n_best` (optional) - Top K predictions returned. Defaults to 5. Must be less than or equal to `beam`

To test CLI inference with the sample data provided, run:

    python Synthesis/translate_cli.py --runtime local --source_file data/source_small.txt --destination_file data/prediction_outputs.csv --target_file data/target_small.txt


## Training Locally

If training sequence to sequence models is your thing, you can retrain the final model using the scripts in the `/train` directory. See the README in the [train directory](https://github.com/kheyer/Deep-Synthesis/tree/training/train) for full instructions.


# Running on the Cloud

Both the streamlit front end and model inference can be run on AWS. This setup is much more involved, but the general idea is straightforward. The front end is set up to run on AWS EKS, and model inference is set up to run on AWS Lambda.

For AWS Lambda setup, see the README in the `lambda_setup` directory.

For AWS EKS setup, see the README in the `eks_setup` directory.




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

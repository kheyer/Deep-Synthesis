# Local Setup

The following instructions detail how to set up Deep Synthesis on your local environment. Local setup is great for running bulk predictions from a file, or tinkering with the application.

Local inference will run on either GPU or CPU depending on what is available on the local environment. If Pytorch detects a CUDA enabled device, inference will default to GPU.

## Creating a Local Conda Environment

The following instructions show how to set up Deep Synthesis as a Conda environment on your local machine.

Clone the repository

    git clone https://github.com/kheyer/Deep-Synthesis

Run `build_local/buildEnv.sh`. This shell script will create a conda environment named `deep_synthesis`, install all required packages and download the trained model

    cd Deep-Synthesis
    chmod +x build_local/buildEnv.sh
    ./build_local/buildEnv.sh

Follow the package installation prompts.

Once package installation is complete, activate the new conda environment and start the Streamlit app

    source activate deep_synthesis
    streamlit run Synthesis/app.py local

The Streamlit app is now running locally on port 8501


## Creating a Local Docker Container

Deep Synthesis can also be built as a Docker container

Clone the repository

    git clone https://github.com/kheyer/Deep-Synthesis
    
Build the container

    cd Deep-Synthesis
    docker build -f build_local/local.Dockerfile -t deep_synthesis .
    docker run -d -p 8501:8501 deep_synthesis

The Streamlit app is now running locally on port 8501

## Running Inference Locally

Local inference can be run in two ways - from the application GUI or the command line.

### GUI Inference

The local deployment GUI supports predicting on single strings or bulk prediction from a file.

For single string prediction, choose the "Predict from String" option on the left sidebar. Paste your SMILES string into the text box and follow the prediction prompts.

To make bulk predictions from a file, create a file that has one set of SMILES reactants per line. ie

    CN(C)C(=O)c1cc(Cl)ccn1.CNC.Nc1ccc(O)cc1.O=C(Cl)c1cc(Cl)ccn1
    CCO.C[SH]=c1[nH]ccc(=O)[nH]1.NCCCN(Cc1cccc(Cl)c1)c1ccccn1.O
    C1CCOC1.CCN(CC)CC.Cc1nc2c(Cl)c(N)ccc2s1.O=C(Cl)c1cc(Cl)cc(Cl)c1
    ...

Sample data files are provided in the `/data` directory.

Choose the "Predict from File" option on the left sidebar. This opens a text input box for the folder path containing the data file, defaulting to the `/data` directory in the repo. The drop down box below the directory input lists all files in that directory. Navigate to and select the source file.

Optionally, a file of target SMILES can be provided if known products to the reactions are available. If targets are provided, the model's predictions will be automatically scored.

Load data using the button. For bulk predictions, there will now be a slider on the sidebar that lets you scroll through your data.

<img src="https://github.com/kheyer/Deep-Synthesis/blob/training/media/prediction3.png" width="500" alt="local prediction from string">

Bulk prediction mode allows you to specify the beam width and top-k results from the prediction. Once your prediction parameters are set, run predictions.

Bulk prediction mode also supports saving prediction data to a csv file based on the path inputted to the dialogue box at the bottom of the app page.

<img src="https://github.com/kheyer/Deep-Synthesis/blob/training/media/prediction4.png" width="500" alt="local prediction from string">


### Command Line Inference

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

Note: if you built Deep Synthesis inside a docker container, you should exec into the container to run CLI inference.


## Training Locally

If training sequence to sequence models is your thing, you can retrain the final model using the scripts in the `/train` directory. See the README in the [train directory](https://github.com/kheyer/Deep-Synthesis/tree/training/train) for full instructions.
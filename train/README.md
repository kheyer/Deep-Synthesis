# Training

These shell scripts allow for recreating the training process for the final model. These scripts assume you have already set up the repo environment locally as described in the main README.

To download and preprocess the full training and validation datasets, run

    chmod +x train/preprocessData.sh
    ./train/preprocessData.sh
    
This will download the dataset to the `OpenNMT-py/data/` directory and run the OpenNMT preprocessing pipeline on the data. This will create processed data files and vocabulary files needed for training, as well as the model config file that specifies the model and training protocol. The model config can be found at `OpenNMT-py/config/molecular_transformer.yml`.

To initiate training, run

    chmod +x train/trainModel.sh
    ./train/trainModel.sh
    
This will start the training protocol. The model config assumes you are running training on a GPU enabled device (you don't want to try this on a CPU).

Training the final takes around ~3 days on a single 2080 Ti GPU, but decent results can be achieved after only 16 hours.

To run local inference on a model checkpoint, the `translateData.sh` file shows how to interface with the default OpenNMT CLI for translation. Before running `translateData.sh`, it may be necessary to change the model step checkpoint (if you did not run training to completion). The default batch size may also need to be changed depending on the GPU memory available.

    chmod +x train/translateData.sh
    ./train/translateData.sh
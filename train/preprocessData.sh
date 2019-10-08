# Shell script to download and preprocess data

# Download data
python train/download_data.py 

# Enter OpenNMT
cd OpenNMT-py

# Set directory name.
# Note that this directory name must match what was used in download_data.py
# and what is listed on the model config yaml file
dataset=data/molecular_data

# Run data preprocessing
python preprocess.py -train_src ${dataset}/source_train.txt \
    -train_tgt ${dataset}/target_train.txt \
    -valid_src ${dataset}/source_val.txt \
    -valid_tgt ${dataset}/target_val.txt \
    -save_data ${dataset}/molecular_transformer_data \
    -src_seq_length 1000 \
    -tgt_seq_length 1000 \
    -src_vocab_size 1000 \
    -tgt_vocab_size 1000 \
    -share_vocab

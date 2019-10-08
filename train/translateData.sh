# Shell script for running inference on a trained model
# Currently this assumes the training process has been run for a full 500,000 iterations
# and that translation will be done with a beam width of 10 and a top k value of 10.
# The script also assumes translation will run on GPU 0

cd OpenNMT-py

dataset=data/molecular_data
model=Experiments/Checkpoints/molecular_data/Molecular_transformer_step_500000.pt

python translate.py -model ${model} \
    -src ${dataset}/source_test.txt \
    -output ${dataset}/predictions.txt \
    -batch_size 128 \
    -replace_unk \
    -max_length 200 \
    -verbose \
    -beam_size 10 \
    -n_best 10 \
    -min_length 5 \
    -gpu 0
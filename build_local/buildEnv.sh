echo "creating environment"
conda create -n deep_synthesis
echo "activating environment"
source activate deep_synthesis
echo "installing RDKit"
conda install rdkit -c rdkit
echo "installing requirements"
pip --no-cache-dir install -r requirements.txt
echo "Cloning OpenNMT"
git clone https://github.com/kheyer/OpenNMT-py
echo "Downloading Model"
python build_local/download_model.py
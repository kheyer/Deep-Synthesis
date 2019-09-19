from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd

def process_prediction(smile, canonicalize=True):
    smile = ''.join(smile.split(' '))
    mol = Chem.MolFromSmiles(smile)

    if mol is not None:
        if canonicalize:
            smile = Chem.MolToSmiles(mol, isomericSmiles=True)
    else:
        smile = ''

    return smile
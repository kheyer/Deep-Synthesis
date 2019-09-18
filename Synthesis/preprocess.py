from rdkit import Chem
from rdkit.Chem import Draw

def canonicalize_smiles(smiles, remove_stereo=True):
    # Converts input SMILES to their canonical form
    # Has the option to remove stereochemistry from SMILES if desired
    # Currently predicting with stereochemistry is not supported, but may be in the future
    mol = Chem.MolFromSmiles(smiles)

    if remove_stereo and '@' in smiles:
            Chem.rdmolops.RemoveStereochemistry(mol)
    
    assert mol is not None
    
    return Chem.MolToSmiles(mol, isomericSmiles=True)
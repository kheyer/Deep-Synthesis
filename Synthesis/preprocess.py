from rdkit import Chem
from rdkit.Chem import Draw

class SmilesData():
    def __init__(self, smiles, tokens, target=None, target_tokens=None):
        self.smiles = smiles
        self.smiles_tokens = tokens
        self.target = target
        self.target_tokens = None

    @classmethod
    def single_entry(cls, smiles, target=None):
        smiles = preprocess(smiles)
        smiles_tokens = tokenize(smiles)
        if target:
            target = [preprocess(target)]
            target_tokens = [tokenize(target[0])]
        else:
            target_tokens = None
        
        return cls([smiles], [smiles_tokens], target, target_tokens)

    @classmethod
    def file_entry(cls, source_filename, target_filename=None):
        smiles = list(open(source_filename))
        
        if target_filename:
            target = list(open(target_filename))
        else:
            target = None
        
        return cls.list_entry(smiles, target)

    @classmethod
    def list_entry(cls, source_smiles, target_smiles=None):
        smiles = [preprocess(i) for i in source_smiles]
        smiles_tokens = [tokenize(i) for i in smiles]
        
        if target_smiles:
            target = [preprocess(i) for i in target_smiles]
            target_tokens = [tokenize(i) for i in target]
        else:
            target = None
            target_tokens = None

def preprocess(smiles):
    # Function to preprocess a single SMILES text string

    # Loading from a file may add a '\n' or ' ' to the end of a SMILES string
    smiles = smiles.strip('\n')
    smiles = smiles.strip(' ')

    # if spaces are present, they are joined
    # this makes the processing pipeline robust to inputs that are already tokenized,
    # partially tokenized, tokenized with a different method or using space delimiting
    # to denote reactants from reagents
    if ' ' in smiles:
        smiles = ''.join(smiles.split(' '))

    # sometimes a '>' character is used to denote reactants from reagents
    # this convention is not supported by RDKit, and must be converted to the traditional
    # '.' delimiter
    if '>' in smiles:
        smiles = smiles.replace('>', '.')

    smiles = canonicalize_smiles(smiles, remove_stereo=True)
        
    return smiles

def tokenize(smiles):
    # tokenizes SMILES string by character
    return ' '.join([i for i in smiles])

def canonicalize_smiles(smiles, remove_stereo=True):
    # Converts input SMILES to their canonical form
    # Has the option to remove stereochemistry from SMILES if desired
    # Currently predicting with stereochemistry is not supported, but may be in the future
    mol = Chem.MolFromSmiles(smiles)

    if remove_stereo and '@' in smiles:
            Chem.rdmolops.RemoveStereochemistry(mol)
    
    assert mol is not None
    
    return Chem.MolToSmiles(mol, isomericSmiles=True)
from rdkit import Chem
from rdkit.Chem import Draw

class SmilesData():
    def __init__(self, smiles, tokens, target=None, target_tokens=None):
        # Class to hold SMILES data
        # Holds source SMILES, tokenized sources, and optionally target SMILES / target tokens
        self.smiles = smiles
        self.smiles_tokens = tokens
        self.target = target
        self.target_tokens = None

    @classmethod
    def single_entry(cls, smiles, target=None):
        # Takes as input a single text string for sources and optionally targets
        # Inputs are processed and tokenized
        smiles, smiles_tokens = process_and_tokenize(smiles)
        target, target_tokens = process_and_tokenize(target)
        # Items packaged into lists to be compatible with bulk methods
        return cls([smiles], [smiles_tokens], [target], [target_tokens])

    @classmethod
    def file_entry(cls, source_filename, target_filename=None):
        # Takes as input file paths for files containing source / target SMILES
        # SMILES are loaded as lists and passed to the list_entry method
        smiles = list(open(source_filename))
        
        if target_filename:
            target = list(open(target_filename))
        else:
            target = None
        
        return cls.list_entry(smiles, target)

    @classmethod
    def list_entry(cls, source_smiles, target_smiles=None):
        # Takes as input lists of SMILES
        # SMILES lists are processed and tokenized
        smiles, smiles_tokens = process_and_tokenize(source_smiles)
        target, target_tokens = process_and_tokenize(target_smiles)

        return cls(smiles, smiles_tokens, target, target_tokens)

def process_and_tokenize(smiles):
    if smiles is not None:
        if type(smiles) == str:
            smiles = preprocess(smiles)
            smiles_tokens = tokenize(smiles)
        else:
            smiles = [preprocess(i) for i in smiles]
            smiles_tokens = [tokenize(i) for i in smiles]
    else:
        smiles = None
        smiles_tokens = None 
    
    return (smiles, smiles_tokens)


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
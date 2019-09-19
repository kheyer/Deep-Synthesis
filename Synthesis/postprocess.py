from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd

class Predictions():
    '''
    Predictions is a class to process and manage predictions.
    The inputs are the SmilesData object used to hold data,
    as well as the raw prediction outputs (preds, scores, attns) 
    from the OpenNMT model.

    Valid SMILES predictions are canonicalized, while invalid
    SMILES predictions are set to an empty string

    If targets are provided in the data, the predictions will be scored
    in terms of top-k accuracy, where k is determined by the 'n_best'
    parameter of the prediction
    '''
    def __init__(self, data, preds, scores, attns):
        self.sources = data.smiles
        self.source_tokens = data.smiles_tokens
        self.targets = data.target
        self.target_tokens = data.target_tokens
        
        self.preds = preds
        self.scores = scores
        self.attns = attns
        self.top_k = len(scores[0])
        
        self.post_process()
        
    def post_process(self):
        # Run post processing functions
        self.do_target = bool(self.targets)
        
        self.create_prediction_df()
        
        if self.do_target:
            self.score_predictions()

        # TODO: screen predictions based on stoichiometry

    def create_prediction_df(self):
        # This function flattens predictions (provided as a list of lists from OpenNMT)
        # into a dataframe that can be easily examined.
        # The dataframe will contain both raw predictions and processed predictions
        preds_flat = []
        sources_flat = []
        scores_flat = []
        prediction_ids = []
            
        for i in range(len(self.scores)):
            preds_flat += self.preds[i]
            scores_flat += [score.item() for score in self.scores[i]]
            sources_flat += [self.sources[i] for j in range(len(self.scores[i]))]
            prediction_ids += [i for j in range(len(self.scores[i]))]
        
        prediction_data = list(zip(prediction_ids, sources_flat, preds_flat, scores_flat))
        columns = ['ID', 'Sources', 'Prediction_Tokens', 'Scores']
        self.df = pd.DataFrame(prediction_data, columns=columns)
            
        self.df['Predictions'] = self.df.Prediction_Tokens.map(lambda x: process_prediction(x))
        
    def score_predictions(self):
        # scores predictions if targets are provided
        targets_flat = []
        for i in range(len(self.scores)):
            targets_flat += [self.targets[i] for j in range(len(self.scores[i]))]

        self.df['Targets'] = targets_flat

        self.df['Correct'] = self.df.apply(lambda row: row['Targets'] == row['Predictions'], axis=1)
        gb = self.df.groupby('ID')
        self.score = (gb.Correct.mean() > 0).mean()
        print(self.score)

    def __len__(self):
        # length is the number of examples in the object
        return len(self.sources)


def process_prediction(smile, canonicalize=True):
    # Processes predictions from a list of tokens to a SMILES string
    smile = ''.join(smile.split(' '))
    mol = Chem.MolFromSmiles(smile)

    if mol is not None:
        if canonicalize:
            # if canonicalize, the predicted string is converted to its canonical form
            smile = Chem.MolToSmiles(mol, isomericSmiles=True)
    else:
        # The model may produce invalid SMILES predictions that 
        # do not correspond to a real structure.
        # In this case, mol will return None and we assign the prediction a value ''
        smile = ''

    return smile
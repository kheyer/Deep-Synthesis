from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import numpy as np
from collections import namedtuple
import seaborn as sns
import matplotlib.pyplot as plt

PredictionTuple = namedtuple('Prediction', 
                ['source_tokens', 'prediction_tokens', 'attention', 'legend'])

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
            scores_flat += [score for score in self.scores[i]]
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

    def sample_df(self, idx=0):
        sample_df = self.df[self.df.ID == idx]
        return sample_df

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

def display_parameters(prediction, idx=0):
    # Gathers all data needed to visualize a specific prediction based on idx
    # Returns source sequence tokens, prediction sequence tokens, attention arrays
    # and legends for plotting

    # Outputs are packaged into a named tuple 
    prediction_tokens = prediction.preds[idx]
    source_tokens = [prediction.source_tokens[idx] for _ in range(len(prediction_tokens))]
    scores = prediction.scores[idx]
    attentions = prediction.attns[idx]
    
    if prediction.do_target:
        correct = list(prediction.df[prediction.df.ID == idx].Correct.values)
        legend = [f'Prediction {i}, Probability {np.exp(score):.4} ({corr})' 
                  for i, (score, corr) in enumerate(zip(scores, correct))]
    else:
        legend = [f'Prediction {i}, Probability {np.exp(score):.4}' 
                  for i, score in enumerate(scores)]
        
    params = list(zip(source_tokens, prediction_tokens, attentions, legend))
    params = [PredictionTuple(*i) for i in params]
    return params

def plot_prediction(source_tokens, prediction_tokens, attention, legend, img_size=(400,400)):
    # Generated prediction plots 
    # An RDKit image of the source and prediction molecules is generated
    # Attention scores are plotted against source and prediction tokens
    ### IMPORTANT attention plot must be generated with raw prediction tokens ###
    # canonicalized predictions may be rearranged
    source_mol = Chem.MolFromSmiles(process_prediction(source_tokens))
    prediction_mol = Chem.MolFromSmiles(process_prediction(prediction_tokens))
    legends = ['Source', legend]
    
    im = Draw.MolsToGridImage([source_mol, prediction_mol], legends=legends, subImgSize=img_size,
                                 molsPerRow=2)
    
    attn_plot = plot_attention(source_tokens, prediction_tokens, attention)
    
    return im, attn_plot

def plot_attention(source, target, attention):
    # Creates a heatmap attention score plot based on source and target tokens
    # Source and target must be tokenized space delimited strings
    source_toks = source.split(' ')
    target_toks = target.split(' ')

    # Attention score is padded based on batch prediction
    # We truncate to the relevant size based on source and target tokens
    attention = attention[:len(target_toks), :len(source_toks)]
    figsize = (attention.shape[1]//4, attention.shape[0]//4)
    fig, ax1 = plt.subplots(figsize=figsize)

    ax = sns.heatmap(attention, linewidths=0.1, ax=ax1, linecolor='black',
                    xticklabels=source_toks, yticklabels=target_toks, square=True,
                    cmap=sns.color_palette("YlGn", n_colors=15), 
                    cbar_kws={"shrink": 0.5, 'pad':0.04, 'label': 'Attention Score'})

    ax.set_xlabel('Source Tokens')
    ax.set_ylabel('Target Tokens')

    loc, labels = plt.yticks()
    ax.set_yticklabels(labels, rotation=360)

    ax.tick_params(top=True, bottom=False,
                labeltop=False, labelbottom=True)

    for _, spine in ax.spines.items():
        spine.set_visible(True)

    return ax
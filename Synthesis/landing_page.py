import streamlit as st 
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

rx1 = 'C1[CH](C)C[CH](Cl)C1.[OH-]>>'
rx2 = 'C1[CH](C)C[CH](Cl)C1.[OH-]>>C1[CH](C)C[CH](O)C1.[Cl-]'
rx3 = 'CC(C)(C#N)N=NC(C)(C)C#N.CCOC(C)=O.Cc1ccn(-c2ccc(I)cc2)n1.ClC(Cl)(Cl)Cl.O=C1CCC(=O)N1Br.[N-]=[N+]=[N-]>>'
rx4 = '[N-]=[N+]=NCc1ccn(-c2ccc(I)cc2)n1'
rx5 = 'CCOC(=O)c1cnc2c(F)cc(SCc3ccccc3)cc2c1O.CO.ClC(Cl)Cl.NCc1ccc(Cl)cc1>>'
rx5_preds = ['O=C(NCc1ccc(Cl)cc1)c1cnc2c(F)cc(SCc3ccccc3)cc2c1O', 
             'CCOC(=O)c1cnc2c(F)cc(S(=O)Cc3ccccc3)cc2c1O',
             'NC(=O)c1cnc2c(F)cc(SCc3ccccc3)cc2c1O']
rx5_labels = ['Prediction 1, Probability 0.9986',
                'Prediction 2, Probability 0.00026',
                'Prediction 3, Probability 0.00019']
rx6 = 'C1COCCO1.COC(=O)c1ccc(C(C)(C)C)c(C#N)c1.[Li+].[OH-]>>CC(C)(C)c1ccc(C(=O)O)cc1C#N'

mol1 = 'COc1ccc2c(Cl)nc(Nc3cc(C)[nH]n3)cc2c1'
mol2 = 'CN(C)CCCN(C)c1ccc(C(F)(F)F)cc1N'
mol3 = 'COc1cccc(Nc2nc(Cl)nc3ccccc23)c1OC'

def rxn_to_image(rxn, img_size=(200,200)):
    return Draw.ReactionToImage(AllChem.ReactionFromSmarts(rxn, useSmiles=True), subImgSize=img_size)

def landing_page():
    st.write('Welcome to Deep Synthesis\n')
    st.write('Deep Synthesis is a deep learning tool designed ', 
            'to predict the products of an organic synthesis reaction. ',
            'Deep Synthesis was designed to help chemists explore synthesis in silico ',
            'and rapidly iterate ideas.')
    st.write('How does deep learning help chemists? Lets look at an example. ', 
            'What happens to 3-methylcyclopentylchloride ',
            'when it is exposed to a base?')
    st.image(rxn_to_image(rx1), use_column_width=True)
    st.write('For this reaction, the hydroxide would replace the chlorine to form an alcohol.')
    st.image(rxn_to_image(rx2), use_column_width=True)
    st.write('This is a fairly straightforward reaction. The SN2 displacement of the chloride ',
                'is really the only reaction possible.')
    st.write('But what if your reactant mixture was much more complicated? What if it looked something like this:')
    st.image(rxn_to_image(rx3, img_size=(150,200)), use_column_width=True)
    st.write("That's a lot of compounds. There's a lot of ways things could react. What will the ",
            'major product be?')
    st.write('Did you guess this?')
    st.image(rxn_to_image(rx3+rx4, img_size=(150,200)), use_column_width=True)
    st.write('Organic chemistry quickly becomes extremely complicated. When multiple ',
             'compounds with multiple functional groups are mixed in a single reaction, ',
             'the total combinatorial space of possible products explodes. Expert chemists ',
             'tackle this problem by relying on years of studying reaction rules to help ',
             'make educated predictions. This is a lot for someone to keep in their head!')
    st.write('Deep Synthesis makes life easier by helping chemists rapidly iterate ideas in silico. ',
             'Deep Synthesis uses a state of the art deep learning model to predict likely products ',
             'from a given set of reactants.')
    st.image(rxn_to_image(rx5), use_column_width=True)
    st.image(Draw.MolsToGridImage([Chem.MolFromSmiles(i) for i in rx5_preds], subImgSize=(400,400),
                                    legends=rx5_labels), use_column_width=True)

    st.write('## Why Deep Learning? How does Deep Synthesis work?')
    st.write('See the "How it Works" tab on the left for an explanation of the methods used')
    st.write('## How do I run predictions?')
    st.write('For a quick tutorial, see the Deep Synthesis Tutorial tab. If ',
                'you are ready to run predictions, go to the "Predict from String" tab.')
    st.write('## Can I check out the source?')
    st.write('All code is available at the [Deep Synthesis Github Repo](https://github.com/kheyer/Deep-Synthesis)')

def smiles_to_image(smiles, img_size=(400,400)):
    mols = [Chem.MolFromSmiles(i) for i in smiles]
    return Draw.MolsToGridImage(mols, subImgSize=img_size, legends=smiles)

def explanation_page():
    st.write('## How it Works')
    st.write('Deep Synthesis treats chemical reactions as a machine translation problem, ',
                'taking advantage of the fact that we can represent molecules as SMILES strings. ',
                'A SMILES string is a structured text representation of an organic compound. For ',
                'example, here are several molecules and their SMILES representations.')
    st.image(smiles_to_image([mol1, mol2, mol3]), use_column_width=True)
    st.write('A chemical reaction can be modeled as a transformation from one string to another. ',
            'The reaction:')
    st.image(rxn_to_image(rx6), use_column_width=True)
    st.write('Can be represented by the following strings:')
    st.write('#### C1COCCO1.COC(=O)c1ccc(C(C)(C)C)c(C#N)c1.[Li+].[OH-] -> CC(C)(C)c1ccc(C(=O)O)cc1C#N')
    st.write('')
    st.write('The same way we can translate text from one language to another, Deep Synthesis translates ',
             'the text of the reactants to the text of the products. Deep Synthesis was trained on ',
             'a dataset of 1.4 million SMILES reaction pairs to learn generalized rules for mapping ',
             'reactions in SMILES space.')
    st.write('## The Model')
    st.write('Inspired by [Schwaller et al.](https://arxiv.org/abs/1811.02633), Deep Synthesis ',
             'uses a sequence to sequence transformer model to translate Reactant SMILES into Product SMILES.')
    st.image('https://raw.githubusercontent.com/kheyer/Deep-Synthesis/training/media/model.png',
                width=400)

def tutorial_page(runtime):

    st.write('Deep Synthesis works by taking a reaction SMILES string as input and ',
             'generating a set of predictions for likely reaction products. ')

    st.write('To use Deep Synthesis, navigate to the "Predict from String" tab and ',
             'paste your SMILES string into the text box. If you do not have SMILES strings on hand ',
             'and you just want to try out the model, select one of the provided options from the ',
              'dropdown menu on the left.')
    st.write('When your data has been entered, click the "Load Data" button to load your data and ',
                'see a visual representation of the input SMILES string')
    st.write('Once the data is loaded, click the "Predict Products" button to run inference on your SMILES.')

    st.image('https://raw.githubusercontent.com/kheyer/Deep-Synthesis/training/media/prediction1.png',
                use_column_width=True)

    st.write('Once inference has been run, you can examine individual predictions by looking at ',
             'attention maps between the reactant SMILES and the predicted product SMILES. ',
             'Attention maps give an idea of what features the model is paying attention to for ',
             'generating the product molecule.')

    st.image('https://raw.githubusercontent.com/kheyer/Deep-Synthesis/training/media/prediction2.png',
                use_column_width=True)

    if runtime == 'local':
        st.write('## Bulk Inference from File')
        st.write("I see you're running Deep Synthesis locally. Local inference mode also supports bulk ",
                    'prediction from files. Use the "Predict from File" option in the sidebar to do this. ',
                    'You will be prompted with a text box to enter the path to your source file. ',
                    'Once you have the correct path selected, use the dropdown menu to choose the ',
                    'relevant file in that directory.')
        st.write('Your source file should be a plain text file with each reactant SMILES on a new line.')
        st.write('Optionally a text file of target sequences can be provided. If targets are provided ',
                 'the model predictions are automatically evaluated for accuracy. This is helpful in ',
                 'benchmarking model performance.')
        st.write('Once the correct files are selected, click the "Load Data" button. This loads all data ',
                 'from the selected file, and creates a slider on the left to scroll through the data.')

        st.image('https://raw.githubusercontent.com/kheyer/Deep-Synthesis/training/media/prediction3.png',
                    use_column_width=True)

        st.write('Prediction runs the same as the single string input case. Once your predictions have ',
                'finished, you can export a csv file of prediction data to a directory of your choosing. ',
                'Enter the save path into the "Prediction Save Destination" dialogue box and ',
                'click the "Save Prediction Data" button.')

        st.image('https://raw.githubusercontent.com/kheyer/Deep-Synthesis/training/media/prediction4.png',
                    use_column_width=True)

    st.write('## Extra notes for those familiar with SMILES and computational chemistry')
    st.write('Deep Synthesis does not require SMILES to be in a canonical form. Predicted ',
            'product SMILES will be returned in canonical form.')
    st.write('Reactant SMILES are expected to contain all reactants and reagents present in the ',
            'reaction to give meaningful results. Simply inputting the SMILES for a single molecule ',
            'will not result in meaningful predictions, except perhaps for a change of oxidation ',
            'of certain molecules.')
    st.write('Currently Deep Synthesis does not support stereochemistry in predictions. The intent is to ',
            'add a second stereo-specific model at a later date.')

# Deep Synthesis

Deep Synthesis is a deep learning driven application for predicting the products of an organic chemical reaction. Deep Synthesis runs a deep learning machine translation model inspired by [Schwaller et al](https://arxiv.org/abs/1811.02633) that takes as input the SMILES string representation of chemcal reactants, and "translates" them to the SMILES strings of the products.

Deep Synthesis allows chemists and chemistry enthusiasts to experiment with reactions in silico.

## How to Use

There are a number of ways to use Deep Synthesis. Choose the one that best fits your needs.

### Online Web App

Deep Synthesis is running online at [deepsynthesis.xyz](deepsynthesis.xyz). Using the web app is great if you want to play with the model or run a small number of predictions.

Input your SMILES string into the text box, or choose one of the examples in the drop down menu on the left. Clicking the "Predict Products" button generates a set of predicted reaction products.

<img src="https://github.com/kheyer/Deep-Synthesis/blob/training/media/prediction1.png" width="500" alt="prediction from string">

Predictions can be further inspected by looking at attention maps between reactant and predicted product strings.

<img src="https://github.com/kheyer/Deep-Synthesis/blob/training/media/prediction2.png" width="500" alt="prediction from string">


### Running Locally

Deep Synthesis can easily be set up on your local machine, either as a Docker container or a conda environment. Local setup is best if you want to run bulk predictions or tinker with the app.

Local setup supports bulk prediction from a text file of SMILES string, either through the app GUI or a command line interface.

For instructions on how to set up locally, see the README of the `build_local` directory.

### Running on AWS

For scalable deployments, we can set up a Kubernetes cluster on AWS to host the front end of the application and a AWS Lambda function to handle inference. If you want to replicate the setup used to host the application on [deepsynthesis.xyz](deepsynthesis.xyz), see the README of the `build_aws` directory. Note that AWS setup is much more involved than local setup, and requires an AWS IAM account with permissions for EKS, EC2 and Lambda.
# Deep Synthesis

Deep Synthesis is a deep learning driven application for predicting the products of an organic chemical reaction. Deep Synthesis runs a deep learning machine translation model inspired by [Schwaller et al](https://arxiv.org/abs/1811.02633) that takes as input the SMILES string representation of chemcal reactants, and "translates" them to the SMILES strings of the products.

Deep Synthesis allows chemists and chemistry enthusiasts to experiment with reactions in silico.

Deep Synthesis is running online at [deepsynthesis.xyz](http://deepsynthesis.xyz/).

# Repository Structure

    Deep-Synthesis
    ├── Synthesis
    │   └── Main program files
    ├── build_aws
    │   └── Instructions for AWS setup
    ├── build_local
    │   └── Instructions for local setup
    ├── configs
    ├── data
    │   └── Small sample datasets
    ├── media
    ├── train
        └── Code for retraining the model

Main application files are stored in the `Synthesis` directory

Code and instructions for setup on AWS are contained in the `build_aws` directory

Code and instructions for local setup are contained in the `build_local` directory

Code for retraining the model (assuming a local install) is contained in the `train` directory

## How to Install

The Deep Synthesis repo supports local installation and setup on AWS.

### Local Setup

Deep Synthesis can easily be set up on your local machine, either as a Docker container or a conda environment. Local setup is best if you want to run bulk predictions or tinker with the app.

Local setup supports bulk prediction from a text file of SMILES string, either through the app GUI or a command line interface.

The quickest way to get up and running is to install Deep Synthesis using Docker. If you do not have Docker installed, follow the [Docker Download Instructions](https://docs.docker.com/install/). Then run the following commands:

    git clone https://github.com/kheyer/Deep-Synthesis
    cd Deep-Synthesis
    docker build -f build_local/local.Dockerfile -t deep_synthesis .
    docker run -d -p 8501:8501 deep_synthesis

Deep Synthesis is now running locally on port 8501.

For additional local installation instructions, see the README of the `build_local` directory. [link](https://github.com/kheyer/Deep-Synthesis/tree/master/build_local). The `build_local` README details how to set up Deep Synthesis as a Conda environment as an alternative to Docker, how to run bulk predictions on your local install, and how to run predictions from the command line.

### AWS Setup

For a more scalable setup, Deep Synthesis can be run on AWS. We can set up a Kubernetes cluster on AWS to host the front end of the application and a AWS Lambda function to handle inference.

<img src="https://github.com/kheyer/Deep-Synthesis/blob/readme_updates/media/aws_setup.png" width="600" alt="AWS setup">

This is the framework being used to host the application at [deepsynthesis.xyz](http://deepsynthesis.xyz/). For full details on setting up Kubernetes and AWS Lambda, see the `build_aws` directory. [link](https://github.com/kheyer/Deep-Synthesis/tree/master/build_aws)

Note that AWS setup is much more involved than local setup, and requires an AWS IAM account with permissions for EKS, EC2 and Lambda.


## Online Web App

Deep Synthesis is running online at [deepsynthesis.xyz](http://deepsynthesis.xyz/). Using the web app is great if you want to play with the model or run a small number of predictions.

Input your SMILES string into the text box, or choose one of the examples in the drop down menu on the left. Clicking the "Predict Products" button generates a set of predicted reaction products.

<img src="https://github.com/kheyer/Deep-Synthesis/blob/master/media/prediction1.png" width="500" alt="prediction from string">

Predictions can be further inspected by looking at attention maps between reactant and predicted product strings.

<img src="https://github.com/kheyer/Deep-Synthesis/blob/master/media/prediction2.png" width="500" alt="prediction from string">

## Model Details

The model used is a sequence to sequence transformer model, implemented in [OpenNMT](https://github.com/kheyer/OpenNMT-py). This model takes as input the SMILES string representation of reactants, and "translates" them to the SMILES of the product molecule. This method was originally developed by [Schwaller et al](https://arxiv.org/abs/1811.02633). Their work is available at the [Molecular Transformer Repo](https://github.com/pschwllr/MolecularTransformer).

Compared to Schwaller, the model shown here was trained from scratch in Pytorch 1.1.0 using the expanded [Patent Reaction Dataset](https://depth-first.com/articles/2019/01/28/the-nextmove-patent-reaction-dataset/). The new model also uses character level tokenization, which reduces model size and removes the need for the chemically constrained beam search procedure used by Schwaller.

### Project Slides

For more details on the project, see the [presentation slides](https://docs.google.com/presentation/d/1YdgaQKAF6Aw3qw0qi9z3Ze6R71vwK7Lpk5uczrCd2zM/edit#slide=id.g64612c95ea_0_0)
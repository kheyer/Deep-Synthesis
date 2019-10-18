# Deep Synthesis

Deep Synthesis is a deep learning driven application for predicting the products of an organic chemical reaction. Deep Synthesis runs a deep learning machine translation model inspired by [Schwaller et al](https://arxiv.org/abs/1811.02633) that takes as input the SMILES string representation of chemcal reactants, and "translates" them to the SMILES strings of the products.

Deep Synthesis allows chemists and chemistry enthusiasts to experiment with reactions in silico.

Deep Synthesis is running online at [deepsynthesis.xyz](http://deepsynthesis.xyz/).

## How to Install

The Deep Synthesis repo supports local installation and setup on AWS.

### Local Setup

Deep Synthesis can easily be set up on your local machine, either as a Docker container or a conda environment. Local setup is best if you want to run bulk predictions or tinker with the app.

Local setup supports bulk prediction from a text file of SMILES string, either through the app GUI or a command line interface.

The quickest way to get up and running is to install Deep Synthesis using Docker. If you do not have Docker, follow the [Docker Download Instructions](https://docs.docker.com/install/). Then run the following commands:

    git clone https://github.com/kheyer/Deep-Synthesis
    cd Deep-Synthesis
    docker build -f build_local/local.Dockerfile -t deep_synthesis .
    docker run -d -p 8501:8501 deep_synthesis

Deep Synthesis is now running locally on port 8501.

For additional local installation instructions, see the README of the `build_local` directory. [link](https://github.com/kheyer/Deep-Synthesis/tree/master/build_local). The `build_local` README details how to set up Deep Synthesis as a Conda environment as an alternative to Docker, how to run bulk predictions on your local install, and how to run predictions from the command line.

### AWS Setup

For a more scalable setup, Deep Synthesis can be run on AWS. We can set up a Kubernetes cluster on AWS to host the front end of the application and a AWS Lambda function to handle inference.

** image placeholder **

This is the framework being used to host the application at [deepsynthesis.xyz](http://deepsynthesis.xyz/). For full details on setting up Kubernetes and AWS Lambda, see the `build_aws` directory. [link](https://github.com/kheyer/Deep-Synthesis/tree/master/build_aws)

Note that AWS setup is much more involved than local setup, and requires an AWS IAM account with permissions for EKS, EC2 and Lambda.


## Online Web App

Deep Synthesis is running online at [deepsynthesis.xyz](http://deepsynthesis.xyz/). Using the web app is great if you want to play with the model or run a small number of predictions.

Input your SMILES string into the text box, or choose one of the examples in the drop down menu on the left. Clicking the "Predict Products" button generates a set of predicted reaction products.

<img src="https://github.com/kheyer/Deep-Synthesis/blob/master/media/prediction1.png" width="500" alt="prediction from string">

Predictions can be further inspected by looking at attention maps between reactant and predicted product strings.

<img src="https://github.com/kheyer/Deep-Synthesis/blob/master/media/prediction2.png" width="500" alt="prediction from string">

### Project Slides

For more details on the project, see the [presentation slides](https://docs.google.com/presentation/d/1YdgaQKAF6Aw3qw0qi9z3Ze6R71vwK7Lpk5uczrCd2zM/edit#slide=id.g64612c95ea_0_0)
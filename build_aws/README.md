# AWS Setup

The contents of this directory show how to set up Deep Synthesis to run on AWS. AWS Setup is designed to run the Streamlit front end on Kubernetes via AWS EKS, while inference with the trained model is run on an AWS Lambda function.

To follow the repo instructions for setting up this deployment, you should have the following AWS things configured:
 * You are running the setup from an AWS IAM account with permissions to use EC2, AWS Lambda, ECR and EKS
 * You have configured the [AWS CLI](https://aws.amazon.com/cli/)
 * You have configured [eksctl](https://eksctl.io/), the CLI for EKS
 * You have configured the [AWS IAM authenticator](https://docs.aws.amazon.com/eks/latest/userguide/install-aws-iam-authenticator.html)
 * You have set up [kubectl](https://kubernetes.io/docs/tasks/tools/install-kubectl/)

  This is the overall architecture that we are going to creeate:

 <img src="https://github.com/kheyer/Deep-Synthesis/blob/readme_updates/media/aws_setup.png" width="600" alt="AWS setup">

 The Streamlit front end for the application will be running on an EKS cluster with a load balancer to manage incoming traffic. Inference will be run on an AWS Lambda function.

The reason for choosing this deployment framework is to create a service with maximal scalability for minimal cost. AWS Lambda is only changed per time used, which makes it the cheapest option for applications running a low number of predictions. EKS creates a scalable deployment for the front end that uses minimal resources during idle periods, but can easily scale to demand. More details on price and deployment considerations are listed in the bottom section.

## Front End Deployment on EKS

Details on deploying the Streamlit front end on AWK EKS can be found in the `EKS_setup` directory.

## Inference On AWS Lambda

Details for deploying the trained model on AWS Lambda can be found in the `lambda_setup` directory.

## Other Deployment Options

Is this the most cost effective way to deploy an application like this? It really depends on the volume of predictions you plan to handle. In a nutshell, AWS Lambda is more expensive per unit time than a comparable EC2 instance, but with EC2 you pay for overall uptime, while with AWS Lambda you only pay for used time. If the cost of time used on Lambda becomes greater than the cost of keeping an EC2 instance up for the same period, running on EC2 is likely a better option.

Putting some numbers to things:
 * t3.medium - 0.0418 $/hr
 * AWS EKS - 0.2 $/hr
 * AWS Lambda - 0.176 $/hr (of used inference time)
 * Load Balancer - 0.025 $/hr

A note on the AWS Lambda price above. This is the price per runtime of the most expensive Lambda instance (3008 MB instance). The reason this is used is that compared to a cheaper instance (1024 MB), the more expensive instance ran inference about 50% faster, and as a result was approximately the same cost to run.

### Current Deployment

The current deployment runs an AWS EKS cluster with a single t3.medium node (which can scale up to 4 nodes). The node runs up to 8 pods (each pod is allocated 250mcpus on an instance with 2 vCPUs). This results in a fixed uptime cost of 0.2668 $/hr (t3.medium price + EKS). Variable cost depends on the useage of the Lambda function.

### Running Inference on EC2

An alternative setup would be to scrap the Lambda function and run inference on a pod in an EC2 instance. At a minimum, this setup would have a fixed uptime cost of 0.2418 $/hr (t3.medium price + EKS) and no variable cost. However things get tricky when you consider CPU requirements for inference and how the application scales to multiple users. Now that inference is run on EC2, each pod requires more CPU resources allocated. CPU resources would need to be quadrupled (250mcpus to 1 CPU) to match the performance of the Lambda function. At a minimum you would need 500mcpus for slower inference. So now the pods per EC2 instance has been decreased from 8 to 2-4. This means that having the same on demand scaling as the above setup would require 2-4 nodes. This would raise the fixed cost to 0.3086 - 0.3922 $/hr.

Take the scenario where your fixed cost is 0.3922 $/hr. Compared to the EKS + Lambda setup, that's an increase of 0.1254 $/hr. So is this worth it? Converting to the costs of AWS Lambda, 0.1254 $/hr is the same cost as running AWS Lambda for 42 minutes every hour. So if your total Lambda inference time averages out to 42 minutes in every hour period (17 hours out of every 24 hour period, or 5 out of 7 days in a week), then running inference on EC2 makes more sense. But remember these usage numbers must be consistent over the lifetime of the application. 

Running the same analysis on the 0.3086 $/hr scenario, the Lambda usage equivalent is 14.2 minutes in every hour period.

### Running Inference on EC2 without EKS

We could also forego EKS entirely. This setup would carry aa minimum cost of 0.0418 $/hr for a t3.medium instance. However, this setup would not scale to multiple users. Based on some ad hoc testing of my own, this setup can handle at most two concurrent users. Any more and things start to fail.
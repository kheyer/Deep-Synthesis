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

Different deployment options come with different costs and levels of scalability. Here is a quick breakdown of different cost options and their level of scalability. First some basic numbers:

 * t3.medium - 0.0418 $/hr
 * AWS EKS - 0.2 $/hr
 * AWS Lambda - 0.176 $/hr (of used inference time for a 3008 MB instance)
 * Load Balancer - 0.025 $/hr

Different deployment options will have their cost determined by which of the above components are used.

### Simple EC2 - 1x t3.medium, 0.0418 $/hr

This is the simplest deployment case. The Streamlit app can run on a single EC2 instance on a single port. The front end and inference are run on the same machine. This is the cheapest option, but it can handle at most two simultaneous users.

### EC2 + EKS - 1x t3.mediumm + EKS + Load Balancer, 0.2668 $/hr

For increased scalability, we can use EKS to manage multiple users. EKS creates multiple pods running the application on a single t3 instance. Inference is run on the same pods that are running the front end. Because of this, each pod needs to be provisioned with sufficient CPU and memory resources to run inference. This imposes a practical limit of 4 pods running concurrently before the cluster would need to be scaled to a second t3 instance.

### EC2 + EKS + Lambda - 1x t3.medium + EKS + Load Balancer + Lambda, 0.2668 $/hr fixed + Lambda variable cost

By moving model inference to AWS Lambda, the pods on our EKS cluster no longer need the resources to run inference. This means that a single t3 instance can reasonably run up to 10 pods before needing to scale. This deployment has the same fixed cost as above, with more scalability within a single t3 instance.

### Evaluating Lambda

The two EKS examples above have the same fixed cost with different scalabilities. The first option can max out four pods before needing to scale to another instance. The second can max out 10 pods before needing to scale, with the downside that actual inference time now incurrs a variable cost. So is the added scalability worth it? It depends on your expected usage.

AWS Lambda (with the configuration used in this repo) runs a cost of 0.176 $/hr. Note that the cost is per hour of actual predictions run, not per hour of uptime. Front end uptime with no predictions incurrs no Lambda cost.

If you never expect to need more than 4 pods available, then running inference on EC2 makes sense. However, if you expect to need 8 pods, the cost of running inference on EC2 becomes 0.3086 $/hr, an additonal 0.0418 $/hr from running a second t3 instance. This additional cost would be the same as running AWS Lambda inference for 14.25 minutes per hour or equivalently 5.7 hours per day.

What is actually cheapest depends on your expected scalability needs vs expected infrence runtime. If you expect a low number of users running long inference operations, then Lambda is likely more expensive. If you expect a higher number of users running short inference operations, Lambda lets you achieve greater scalability at a much cheaper price.
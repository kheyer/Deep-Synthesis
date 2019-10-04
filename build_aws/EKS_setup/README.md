# AWS EKS Setup

These instructions show how to set up the Streamlit front end on AWS EKS. EKS deploys a Kuberentes cluster that allows the front end to be easily scalable. These instructions assumes the following setup:
 * You are running this from an AWS IAM account with permissions to use EC2, ECR and EKS
 * You have configured the [AWS CLI](https://aws.amazon.com/cli/)
 * You have configured [eksctl](https://eksctl.io/), the CLI for EKS
 * You have configured the [AWS IAM authenticator](https://docs.aws.amazon.com/eks/latest/userguide/install-aws-iam-authenticator.html)
 * You have set up [kubectl](https://kubernetes.io/docs/tasks/tools/install-kubectl/)

## Building the Dockerfile

The first step is to build a Dockerfile to run the Streamlit front end and push it to AWS ECR. If you compare this Dockerfile to the local setup Dockerfile, you will notice that the AWS Dockerfile has been slimmed down to exclude Pytorch, OpenNMT and other large dependencies. This is because the AWS front end is designed to interface with a Lambda function running inference. Consequently, none of the packages needed to run inference are required.

First build the image locally

    docker build  -f aws.Dockerfile -t \
    translation-front:1.0-SNAPSHOT .

Then log into AWS ECR

    aws ecr get-login --no-include-email

This will print out a login for ECR. Copy and paste it into the terminal

    docker login -u AWS -p [password_string] https://[aws_account_id].dkr.ecr.us-east-2.amazonaws.com

Create the repository on ECR

    aws ecr create-repository --repository-name translation-front

This returns the following output

```
{
    "repository": {
        "registryId": "[aws_account_id]",
        "repositoryName": "translation-front",
        "repositoryArn": "arn:aws:ecr:us-west-2:[aws_account_id]:repository/translation-front",
        "createdAt": 1569994192.0,
        "repositoryUri": "[aws_account_id].dkr.ecr.us-west-2.amazonaws.com/translation-front"
    }
}
```
Copy the repository URI. Then tag the local Docker image and push to ECR.

    docker tag translation-front:1.0-SNAPSHOT [aws_account_id].dkr.ecr.us-west-2.amazonaws.com/translation-front:1.0-SNAPSHOT
    docker push [aws_account_id].dkr.ecr.us-west-2.amazonaws.com/translation-front:1.0-SNAPSHOT

## Creating the Cluster

To create a cluster, run 

    eksctl create cluster --name=translation-front --nodes=1 --node-type=m4.xlarge

This will create a cluster with a single node on a m4.xlarge instance. You can change the number of nodes and instance type as needed.

Update the `kubernetes.yml` file. Change `spec.templace.spec.container.image` to match the name of the image you just pushed. Then apply the config.

    kubectl apply -f kubernetes.yml

Update your local kube-config

    aws eks --region us-west-2 update-kubeconfig --name translation-front

The config includes opening a NodePort. This is not required, but can be extremely helpful for debugging. If you want to actually use the NodePort, you need to update the security group for the nodegroup to allow access to port 31000. This is done with the following commands.

    aws ec2 describe-security-groups --filters Name=group-name,Values="*eksctl-translation-front-nodegroup*"  --query "SecurityGroups[*].{Name:GroupName,ID:GroupId}"

This will output the `[security_group]` of the nodegroup. Then run

    aws ec2 authorize-security-group-ingress --protocol tcp --port 31000 --group-id [security_group] --cidr 0.0.0.0/0

The running node can now be accessed at `[instance_public_ip]:31000`

Next we need to give the EKS cluster permissions to interface with the lambda function we are using to run inference. Find the IAM roles associated with the node group and the cluster. Apply `AWSLambdaRole` and `AWSLambdaFullAccess` to the node group and the cluster.

## Configure Autoscaling

We will configure autoscaling in two ways - Horizontal auto scaling (HSA) within a node group and cluster scaling to create more nodes.


### Horizontal Autoscaler

To configure HSA, we need to install a metrics server on our cluster. Run the following:

    curl --remote-name --location https://github.com/kubernetes-incubator/metrics-server/archive/v0.3.5.tar.gz
    tar -xzf v0.3.5.tar.gz
    kubectl apply -f metrics-server-0.3.5/deploy/1.8+/
    kubectl get deployment metrics-server -n kube-system

Next, create the following permissions policy and name it `ClusterAutoScaler`. Attach this policy to the node IAM role.

    {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Action": [
                    "autoscaling:DescribeAutoScalingGroups",
                    "autoscaling:DescribeAutoScalingInstances",
                    "autoscaling:DescribeLaunchConfigurations",
                    "autoscaling:DescribeTags",
                    "autoscaling:SetDesiredCapacity",
                    "autoscaling:TerminateInstanceInAutoScalingGroup"
                ],
                "Resource": "*"
            }
        ]
    }

Now we create the Horizontal autoscaler

    kubectl autoscale deployment translation-deployment --cpu-percent=75 --min=2 --max=10

This creates an autoscaler with between 2 and 10 pod replicas based on CPU load. When CPU utilization within a pod reaches 75% of the requested CPU resources (250 millicpu, specified in the `kubernetes.yml` file), a new replica pod is created.

### Cluster Autoscaler

To configure cluster autoscaling, go to the `Auto Scaling Groups` section of the EC2 console. You will see an existing autoscaler group for the node we created. The current default is 
 * Instances - 1
 * Desired - 1
 * Min - 1
 * Max - 1

Edit these numbers to match your expected needs. Note that creating a new node creates a new EC2 instance, which incurs charges.
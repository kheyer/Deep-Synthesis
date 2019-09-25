# Lambda Layer Setup

For deploying the model on AWS, we want to run the model in an [AWS Lambda Function](https://aws.amazon.com/lambda/). The advantage of using AWS Lambda is we only pay for the compute time we use, rather than paying for an instance that is charged 24/7 regardless of use.

The disadvantage to using Lambda is the tight constraints on the code you can deploy. A Lambda deployment is limited to 250 MB uncompressed. This is fine for our model, which is ~145 MB. The problem is loading all the python dependencies. Pytorch alone is over 1 GB. To actually get the model running on AWS Lambda, we need to do some sneaky stuff.

Specifically, we need to create a leaned out version of our dependencies, package them into a zip file, then pull in the zip file to our lambda function as a layer. This lets us pull our dependences into the `/tmp` directory of the lambda function, which has a size limit of 500 MB. To get our package of dependencies below this limit, we will remove unnecessary packages (installing torch automatically installs caffe2 and other things we don't need), as well as unnecessary files like tests.

Getting this all to work was quite difficult. Most of the code in this directory and the commands listed below are adapted from an excellent tutorial repo by [Matt McClean](https://github.com/mattmcclean/sam-pytorch-example) on creating Pytorch packages for AWS Lambda.

## Overview

To create our Lambda function, we need to do a few things:

   1. Create a zip file of our dependencies (with excess files removes)
   2. Create an AWS Lambda layer from that zip file
   3. Package and deploy our inference Lambda function with the dependency lambda layer attached

If you wish to follow these instructions to create the lambda layer and lambda function, it is assumed that you already have a functioning AWS account, as well as the AWS CLI and AWS SAM CLI installed and configured.

## Creating The Zip File

The lean dependency zip file is created with the `create_layer_file.sh` shell script in the `/layer` directory. This shell script first installs all the packages listed in the `requirements.txt` file in the `/translation_function` directory, then removes excess files. The full folder of lean dependencies comes out to about 480 MB in size. This is then zipped into a file called `OpenNMT-Pytorch-1.1.0-lambda-layer.zip`. To do this, run the following commands:

    cd lambda_setup/layer
    chmod +x create_layer_file.sh
    ./create_layer_file.sh

Once the zip file is created, we copy it to an S3 bucket

    aws s3 cp OpenNMT-Pytorch-1.1.0-lambda-layer.zip s3://{BUCKET}/lambda-layers/OpenNMT-Pytorch-1.1.0-lambda-layer.zip

## Creating the Lambda Layer

Once the zip file is on S3, we can (in theory) create our Lambda layer from the command line

    aws lambda publish-layer-version \
        --layer-name “opennmt-py” \
        --description "Lambda layer of OpenNMT Pytorch 1.1.0 zipped to be extracted with unzip_requirements file" \
        --content "S3Bucket={BUCKET},S3Key=lambda-layers/OpenNMT-Pytorch-1.1.0-lambda-layer.zip" \
        --compatible-runtimes "python3.6" 

I say in theory because I often had issues getting this command to run. If the command throws an error or hangs, go to the `Layers` tab on the AWS Lambda console and click `Create Layer`. Follow the instructions to create a layer through the AWS console.

Take note of the version ARN for the created lambda layer. ie `arn:aws:lambda:{REGION}:{ID}:layer:opennmt-py:1`. If you adapt this code to create a dependency layer of your own, update the `LambdaLayerArn` default value in the `template.yaml` file.

The Lambda layer is now functional.

## Deploying the Lambda Function

The actual directory that will become the Lambda function is the `/translation_function` directory. This directory must be a flat directory containing all the files the Lambda function will need. The actual function called when the Lambda function is invoked is the `lambda_handler` function in the `app.py` file. This is specified in the `Handler` value of the `TranslationFunction` section in `template.yaml`.

The trained model must also be stored in the `/translation_function` directory. This model is not stored on github for size reasons. It can be downloaded by running `download_model.py` from the `/lambda_setup` directory.

    cd lambda_setup
    python download_model.py

Now we can create the Lambda function. First we package our files and upload them to S3

    sam package \
        --output-template-file packaged.yaml \
        --s3-bucket {BUCKET}

Then we deploy to AWS CloudFormation

    sam deploy \
        --template-file packaged.yaml \
        --stack-name {APP_NAME} \
        --capabilities CAPABILITY_IAM \
        --parameter-overrides BucketName={BUCKET}

Once creation is complete, the function is live and operational. There's just one more thing to do.

When testing inference on AWS Lambda, I had issues with the message length limit on the response from the Lambda function. Trying to predict on a reasonable number of items and return all predictions, prediction scores and attention maps went way over the limit. To get around this, the function is designed to store all prediction outputs to S3 and return the S3 item key. To do this, the Lambda function must be given PUT permission for the S3 bucket in question.

To do this, go to the AWS IAM console. Go to `Policies`, `Create Policy`, and add the following as a JSON input.

    {
        "Version": "2012-10-17",
        "Statement": [
    {
                "Effect": "Allow",
                "Action": [
                    "s3:ListAllMyBuckets",
                    "s3:GetBucketLocation"
                ],
                "Resource": "*"
            },
            {
                "Effect": "Allow",
                "Action": "s3:*",
                "Resource": [
                    "arn:aws:s3:::{BUCKET}",
                    "arn:aws:s3:::{BUCKET}/*"
                ]
            }
        ]
    }

Then go to `Roles` and find the role associated with `{APP_NAME}`. Attach the new policy to the role.

Now everything is ready to go.

## Testing Locally

Once the function is set up and has proper S3 permissions, we can test the function locally. From the `/lambda_setup` directory, run the following.

    sam local invoke TranslationFunction -n env.json -e event.json

If everything is set up properly, this will run inference on the sample data in `event.json`, print prediction data to the console, store prediction data on S3 and return a JSON response containing the S3 object key for the stored predictions.
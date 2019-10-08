import argparse
import json
from utils import *
from preprocess import *
from translate_aws import *
from lambda_async import *

# Importing translate also imports OpenNMT.
# Putting this in a try/except block allows for deploying
# the app configured to run AWS Lambda predictions
# without bringing along OpenNMT
try:
    from translate import *
except ModuleNotFoundError:
    pass

# define arg specifications on launch
parser = argparse.ArgumentParser()
parser.add_argument("--runtime", help="Argument to select local or AWS runtimes ['local' or 'AWS']")
parser.add_argument("--source_file", help="Path to source file for translation")
parser.add_argument("--destination_file", help="Path to destination file to save translations (CSV file)")
parser.add_argument("--target_file", help="Path to file of known targets (optional)", default=None)
parser.add_argument("--beam", help="Beam search width", default=5)
parser.add_argument("--n_best", help="Top k results", default=5)
args = parser.parse_args()

# Dictionary of configuration files
config_json_opts = {
                      'AWS' : 'configs/aws_description.json',
                      'local' : 'configs/model_description.json'
                   }

runtime = args.runtime 
model_description = json.load(open(config_json_opts[runtime]))

# Define class to handle translations based on runtime
if runtime == 'local':
    translator_class = TranslationModel
elif runtime == 'AWS':
    translator_class = LambdaInterface
    
else:
    raise ValueError('''Please provide a valid runtime argument. Use 'local' to run predictions locally, or 
                'AWS' to run predictions on AWS (AWS inference requires permissions)''')

# load data
data = SmilesData.file_entry(args.source_file, args.target_file)

# run interence
translator = translator_class(model_description)
# returning attention maps for CLI inference is not supported 
scores, preds, attns = translator.run_translation(data.smiles_tokens, 
                            beam=int(args.beam), n_best=int(args.n_best), return_attention=False)

prediction = Predictions(data, preds, scores, attns)

# save prediction to destination filename
prediction.df.to_csv(args.destination_file, index=False)
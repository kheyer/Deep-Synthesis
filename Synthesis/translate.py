import sys 
sys.path.insert(0, 'OpenNMT-py')

from onmt.utils.parse import ArgumentParser
import onmt.opts as opts
from onmt.translate.translator import build_translator
import torch

import types
import io

class TranslationModel:
    '''
    Class for handing machine translation
    The TranslationModel class will load a model given a config file

    The config file should contain:
    Path to the model (string)
    Beam size (int)
    K best predictions (int)
    Minimum prediction length (int)
    Maximum prediction length (int)
    Replace unk (bool)

    This class is designed to work with models created with OpenNMT
    OpenNMT relies on an extensive argparse structure to create models
    The code here has been build to plug into that framework
    '''
    def __init__(self, config):
        self.config = config

        # Configure GPU if available
        if torch.cuda.is_available():
            self.config['opt']['gpu'] = 0
            self.gpu = True 
        else:
            self.gpu = False

        # create OpenNMT opt object
        self.opt = self.configure_opt(config)

        # use opt to load model
        self.translator = self.load_translator(self.opt)
        self.base_beam = self.translator.beam_size
        self.base_best = self.translator.n_best

    def configure_opt(self, description):
        # Uses OpenNMT's ArgumentParser class to create an object
        # That holds all the parameters needed to load the model
        parser = ArgumentParser(description='translation')
        opts.config_opts(parser)
        opts.translate_opts(parser)
        opt = {a.dest: a.default for a in parser._actions}
        opt.update(description['opt'])
        opt['models'] = [description['model']]
        opt = types.SimpleNamespace(**opt)

        return opt

    def load_translator(self, opt):
        # Loads a model based on the parameters in opt

        # OpenNMT requires some sort of output file to write to
        # Since we don't intent to write predictions to a file, we 
        # can use a StringIO object
        output = io.StringIO()

        # build_translator is an OpenNMT function
        translator = build_translator(opt, out_file=output)

        return translator

    def run_translation(self, src, beam=None, n_best=None, return_attention=True):
        # Runs src through a loaded model and generates predictions
        # src should be a list of processed and tokenized SMILES

        self.translator.out_file.seek(0)
        self.translator.out_file.truncate(0)
        
        # beam and n_best are set in opt, but we can overwrite those values if desired
        if beam:
            self.translator.beam_size = beam
        if n_best:
            self.translator.n_best = n_best
            
        # beam size must be greater than or equal to the number of final predictions
        assert self.translator.beam_size >= self.translator.n_best
        
        # sets batch size based on beam size and hardware (GPU or CPU)
        bs = self.get_batch_size(self.translator.beam_size)

        # outputs from translation
        # all outputs are lists of lists
        scores, preds, attns = self.translator.translate(src, batch_size=bs, return_attention=True)
        
        if beam or n_best:
            self.reset_params()

        if return_attention:
            return (scores, preds, attns)
        else:
            return (scores, preds)

    def reset_params(self):
        # reset beam_size and n_best parameters if desired
        self.translator.beam_size = self.base_beam
        self.translator.n_best = self.base_best
    
    def get_batch_size(self, beam):
        # Calls correct batch size function based on CPU or GPU inference
        if self.gpu:
            bs =  self.get_batch_size_gpu(beam)
        else:
            bs = self.get_batch_size_cpu(beam)

        return bs

    def get_batch_size_cpu(self, beam):
        # Optinum batch sizes for CPU inference
        # based on AWS Lambda 2048 MB instance
        if beam >= 10:
            return 10
        else:
            return 15
    
    def get_batch_size_gpu(self, beam):
        # Optimum batch sizes for GPU inference 
        # Based off K80 GPU
        if beam == 1:
            bs = 256
        if beam >=2:
            bs = 192
        if beam >= 5:
            bs = 128
        if beam >= 10:
            bs = 64
            
        return bs
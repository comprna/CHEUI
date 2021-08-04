#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 10:16:37 2021

@author: labuser
"""

import argparse

parser = argparse.ArgumentParser(prog='CHEUI_predict_model1 v0.1', description=
                                 """ 
                                 This script takes an ID+signal file generated using CHEUI_preprocess_* and predict methylation status \
                                 per read and per 9-mer.
                                 
                                 """, usage='python CHEUI_predict_model1.py -i <path_to_signlas+IDs_file> '\
                                            '-m <path_to_DL_model> -o <file_out> \nversion: %(prog)s')

OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')


REQUIRED.add_argument("-i", "--signals_input",
                      help="path to the ID+signal file",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-m", "--DL_model",
                      help="path to trainned model 1 DL model",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-o", "--file_out",
                      help="Path to the output file",
                      metavar='\b',
                      required=True)

OPTIONAL.add_argument('-v', '--version', 
                        action='version', 
                        version='%(prog)s')                  

parser._action_groups.append(OPTIONAL)

ARGS = parser.parse_args()

# required arg
signals_input = ARGS.signals_input
DL_model = ARGS.DL_model
file_out = ARGS.file_out


from tensorflow.keras import Input
from tensorflow.keras.models import Model
import _pickle as cPickle
from DL_models import build_Jasper
import pandas as pd
import numpy as np
import os
import sys

# load the trainned model
inputs = Input(shape=(100, 2))
output = build_Jasper(inputs,Deep=True)
model = Model(inputs=inputs, outputs=output)


model.load_weights(DL_model)

### load the stored data 
counter = 0
IDs_signals = {}


if os.path.isfile(file_out):
    print('WARNING: read level prediction file already exists, please delete it or change the output name')
    sys.exit()


with open(signals_input, 'rb') as signal_in:
    while True:
        try:
            counter +=1
            IDs_signals.update(cPickle.load(signal_in))
            # to avoid loading everything predict every 10k singnals
            if counter%500 == 0:
                if IDs_signals:
                    predictions = model.predict(np.array(list(IDs_signals.values())))
                    
                    predictions_df = pd.DataFrame.from_dict({'KMER': IDs_signals.keys(),
                                                            'Prediction': predictions.reshape(len(predictions)).tolist()}
                                                            )
                    
                    predictions_df.to_csv(file_out,
                                          mode='a',
                                          header=False,
                                          sep='\t', 
                                          index=False)
                    IDs_signals = {}
                    predictions_df = pd.DataFrame({'KMER': [], 'Prediction': []})
                else:
                    continue
        except:
            if IDs_signals:
                predictions = model.predict(np.array(list(IDs_signals.values())))
                predictions_df = pd.DataFrame.from_dict({'KMER': IDs_signals.keys(),
                                                        'Prediction': predictions.reshape(len(predictions)).tolist()}
                                                        )
                
                predictions_df.to_csv(file_out,
                                      mode='a',
                                      header=False,
                                      sep='\t', 
                                      index=False)
                IDs_signals = {}

                predictions_df = pd.DataFrame({'KMER': [], 'Prediction': []})
            print(file_out)
            print('All signals have been processed')
            break




















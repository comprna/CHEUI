#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 10:16:37 2021

@author: labuser
"""

import argparse

parser = argparse.ArgumentParser(prog='predict_MILONGAS v0.1', description=
                                 """ 
                                 This script takes a ID and signal files and generate predictions every A 
                                 """, usage='python predict_MILONGAS.py -i <directory_input> -o <file_out> \nversion: %(prog)s')

OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')


REQUIRED.add_argument("-i", "--directory_input",
                      help="directory input",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-o", "--file_out",
                      help="output file",
                      metavar='\b',
                      required=True)

OPTIONAL.add_argument('-v', '--version', 
                        action='version', 
                        version='%(prog)s')


OPTIONAL.add_argument("-n", "--suffix_name",
                      help='name to use for output files',
                      metavar='<str>',
                      )                      

parser._action_groups.append(OPTIONAL)

ARGS = parser.parse_args()

# required arg
nanopolish_path = ARGS.input_nanopolish
file_out = ARGS.file_out
directory_out = ARGS.out_dir


# optional arg
suffix_name = ARGS.suffix_name


from tensorflow.keras import Input
from tensorflow.keras.models import Model
import _pickle as cPickle
from DL_models import build_Jasper
import pandas as pd
import numpy as np
import sys

# load the trainned model 
inputs = Input(shape=(100, 2))
output = build_Jasper(inputs,Deep=True)
model = Model(inputs=inputs, outputs=output)
    
# test the NN trainned with only 1 A
model.load_weights('/media/labuser/Data/nanopore/m6A_classifier/results/single+double/2vectors_all/103400000model.h5')

directory_data = '/media/labuser/Data/nanopore/m6A_classifier/data/yeast/nanopolish_reads/test_MILOGNAS_preprocess'
directory_data = sys.argv[1]
directory_out = '/media/labuser/Data/nanopore/m6A_classifier/data/yeast/nanopolish_reads/test_MILOGNAS_preprocess'
directory_out = sys.argv[2]

### load the stored data 
counter = 0
IDs = []
signals = []
with open(directory_data+'/'+'signals_chr1.P', 'rb') as signal_in:
    with open(directory_data+'/'+'IDs_chr1.P', 'rb') as id_in:
        while True:
            try:
                counter +=1
                IDs.append(cPickle.load(id_in))
                signals.append(cPickle.load(signal_in))
                # to avoid loading everything predict every 10k singnals
                if counter%10000 == 0:
                    print(counter, 'signals predicted')
                    predictions = model.predict(np.array(signals))
                    predictions_df = pd.DataFrame.from_dict({'KMER': IDs,
                                                            'Prediction': predictions.reshape(len(predictions)).tolist()}
                                                            )
                    
                    predictions_df.to_csv(directory_out+'/'+file_out,
                                          mode='a',
                                          header=False,
                                          sep='\t')
                    IDs = []
                    signals = []
                    predictions_df = pd.DataFrame({'KMER': [], 'Prediction': []})
            except:
                if IDs:
                    predictions = model.predict(np.array(signals))
                    predictions_df = pd.DataFrame.from_dict({'KMER': IDs,
                                                            'Prediction': predictions.reshape(len(predictions)).tolist()}
                                                            )
                    
                    predictions_df.to_csv(directory_out+'/'+file_out,
                                          mode='a',
                                          header=False,
                                          sep='\t')
                    IDs = []
                    signals = []
                    predictions_df = pd.DataFrame({'KMER': [], 'Prediction': []})
                print('All signals have been processed')
                break




















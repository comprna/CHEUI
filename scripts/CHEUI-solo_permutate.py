#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 14:43:48 2021

@author: labuser
"""

import argparse
import numpy as np
import pandas as pd
import random

from tensorflow.keras import Input
from tensorflow.keras.models import Model
from DL_models import build_Jasper
random.seed(42)


parser = argparse.ArgumentParser(prog='CHEUI-solo_calculate_pvalues v0.1', description=
                                 """ 
                                 This script calculate empirical p-values and q-values from CHEUI predictions \
                                 
                                 """, usage='python CHEUI-solo_calculate_pvalues.py -c'\
                                            '<config_file> \nversion: %(prog)s')

OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')

REQUIRED.add_argument("-r", "--input_read_level",
                      help="path to the read level probability file",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-s", "--input_site_level",
                      help="path to the site level probability file",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-m", "--DL_model",
                      help="path to pretrainned DL model 2",
                      metavar='\b',
                      required=True)

OPTIONAL.add_argument('-v', '--version', 
                        action='version', 
                        version='%(prog)s')

REQUIRED.add_argument("-o", "--output",
                      help="Path to the output file",
                      metavar='\b',
                      required=True)


ARGS = parser.parse_args()

# required arg
input_read_level = ARGS.input_read_level
input_site_level = ARGS.input_site_level
DL_model = ARGS.DL_model
file_out_path = ARGS.output


def convert_p_to_vector(probs):
    '''
    '''
    try:
        probs = sorted(probs)
    except:
        if 'Prediction' in probs:
            return []
        probs = [float(i) for i in probs]
    prob_dist = []
    for i in range(1, 100):
        count = 0
        for j in probs:
            if j>=i/100 and j<(i+1)/100:
                count += 1
        prob_dist.append(count)
    return(prob_dist)


def get_permutations(coverages_for_permutation, list_prob):
    '''
    '''
    permutations = []
    for j in enumerate(coverages_for_permutation):
        prob = random.sample(list_prob, j[1])
        prob_vec = convert_p_to_vector(prob)
        if len(prob_vec) == 0:
            continue
        permutations.append(prob_vec)
    return permutations

# Load CNN weights
inputs = Input(shape=(99, 1))
output = build_Jasper(inputs, 1)
model = Model(inputs=inputs, outputs=output)
model.load_weights(DL_model)


# first get the coverage of all real sites 
#input_site_level = '/media/labuser/Data/nanopore/mouse_brains/results/m6A/WT_E18_read_level.txt_site_level.txt'

site_level = pd.read_csv(input_site_level, sep='\t')

#input_read_level = '/media/labuser/Data/nanopore/mouse_brains/results/m6A/test.txt'

read_level_probabilities = pd.read_csv(input_read_level,
                                       sep='\t', 
                                       usecols=[1],
                                       header=None, 
                                       dtype={0: np.float64}).iloc[:,0].tolist()

prob_permuted = get_permutations(site_level['coverage'], read_level_probabilities)

# make preditions from the permuted sites to create the null
permutations_predictions = model.predict(prob_permuted)

with open(file_out_path, 'w') as file_out:
    for i in permutations_predictions:
        print(i[0], file=file_out)

















































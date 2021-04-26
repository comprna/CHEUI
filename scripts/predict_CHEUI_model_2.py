#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 10:16:37 2021

@author: labuser
"""

import argparse

parser = argparse.ArgumentParser(prog='predict_MILONGAS v0.1', description=
                                 """ 
                                 This script takes predictions from model 1 generate predictions per site
                                 
                                 """, usage='python predict_CHEUI_model2.py -i <path_to_predictions_model_1> '\
                                            '-m <path_to_DL_model> -o <file_out> \nversion: %(prog)s')

OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')

REQUIRED.add_argument("-i", "--input",
                      help="path to the input prediction file",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-m", "--DL_model",
                      help="path to DL model",
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
input_df = ARGS.input
DL_model = ARGS.DL_model
file_out = ARGS.file_out


from tensorflow.keras import Input
from tensorflow.keras.models import Model
import _pickle as cPickle
from DL_models import build_Jasper
import pandas as pd
import numpy as np
from numba import jit


def convert_p_to_vector(probs):
    '''
    '''
    probs = sorted(probs)
    prob_dist = []
    for i in range(1, 100):
        count = 0
        for j in probs:
            if j>=i/100 and j<(i+1)/100:
                count += 1
        prob_dist.append(count)
    return(prob_dist)


@jit(nopython=True)
def convert_p_to_vector_faster(probs):
    '''
    '''
    probs = sorted(probs)
    prob_dist = []
    for i in range(1, 100):
        count = 0
        for j in probs:
            if j>=i/100 and j<(i+1)/100:
                count += 1
        prob_dist.append(count)
    return(prob_dist)


# Load CNN weights
inputs = Input(shape=(99, 1))   
output = build_Jasper(inputs, 1)
model = Model(inputs=inputs, outputs=output)
secondML_path = DL_model
model.load_weights(secondML_path)

input_df = '/media/labuser/Data/nanopore/m6A_classifier/data/ELIGOS/predictions_modifications/m7G_singleRead_predictions.sorted.txt'
input_df = '/home/pablo/lib/CHEUI/m7G_singleRead_predictions.sorted.txt'


ID = ''
predictions_site = []
predictions_dic = {}
stoichiometry_dic = {}
coverage_dic = {}
counter = 0

# code the last one

with open(input_df, 'r') as input_file:
    for line in input_file:
        line = line.strip().split('\t')
        if counter == 0:
            predictions_site.append(float(line[1]))
            ID = '_'.join(line[0].split('_')[:-1])
            counter +=1
            continue
        
        if ID != '_'.join(line[0].split('_')[:-1]):
            # if the coverage if smaller than 10 

            if len(predictions_site) < 10:
                # store the info. from the current read
                predictions_site = [float(line[1])]
                ID = '_'.join(line[0].split('_')[:-1])
                counter +=1
                
            else:
                vector_prob = convert_p_to_vector(predictions_site)
                #lr_probs = model.predict(np.array(vector_prob).reshape(1,99,1))
                predictions_dic[ID] = vector_prob
                # calculate stoichiometry         
                mod = [i for i in predictions_site if i > 0.7]
                no_mod = [i for i in predictions_site if i < 0.3]
                stoichiometry_dic[ID] = len(mod)/(len(mod)+len(no_mod))
                coverage_dic[ID] = len(predictions_site)
                predictions_site = [float(line[1])]
                ID = '_'.join(line[0].split('_')[:-1])
                counter +=1
        else:
            predictions_site.append(float(line[1]))
            ID = '_'.join(line[0].split('_')[:-1])
        
        if counter % 1000 == 0:
            print(counter,'number of lines processed')
            print(ID)







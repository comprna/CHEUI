#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 10:16:37 2021

@author: labuser
"""

import argparse

parser = argparse.ArgumentParser(prog='CHEUI_predict_model2 v0.1', description=
                                 """ 
                                 This script takes predictions from model 1 and generate RNA modification predictions per site
                                 
                                 """, usage='python CHEUI_predict_model2.py -i <path_to_predictions_model_1> '\
                                            '-m <path_to_DL_model> -o <file_out> \nversion: %(prog)s')

OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')

REQUIRED.add_argument("-i", "--input",
                      help="path to the input prediction file",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-m", "--DL_model",
                      help="path to pretrainned DL model 2",
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
file_out_path = ARGS.file_out

from tensorflow.keras import Input
from tensorflow.keras.models import Model
from DL_models import build_Jasper
import numpy as np
from numba import jit
import random
random.seed(42)



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


def biggerThan100(prob_list):
    '''
    '''
    all_probs = []
    for i in range(16):
        random_reads = random.sample(prob_list, 75)
        vector_prob = convert_p_to_vector_faster(random_reads)
        all_probs.append(model.predict(np.array(vector_prob).reshape(1,99,1)))
    high_conf_times = [i for i in all_probs if i > 0.99]
    return high_conf_times


# Load CNN weights
inputs = Input(shape=(99, 1))
output = build_Jasper(inputs, 1)
model = Model(inputs=inputs, outputs=output)
secondML_path = DL_model
model.load_weights(secondML_path)


predictions_site = []
counter = 0

# code the last one
with open(file_out_path, 'w') as file_out:
    
    print('contig'+'\t'+'position'+'\t'+'site'+'\t'+'coverage'+\
          '\t'+'stoichiometry'+'\t'+'probability', file=file_out)
        
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
    
                if len(predictions_site) < 15:
                    # do not make predictions for this site
                    # store the info. from the current read
                    predictions_site = [float(line[1])]
                    ID = '_'.join(line[0].split('_')[:-1])
                    counter +=1
                    
                else:
                    # calculate stoichiometry
                    mod = [i for i in predictions_site if i > 0.7]
                    no_mod = [i for i in predictions_site if i < 0.3]
                    stoichiometry = len(mod)/(len(mod)+len(no_mod))
                    
                    if len(predictions_site) > 100 and stoichiometry > 0.1:
                        lr_probs = biggerThan100(predictions_site)
                        if len(lr_probs) > 8:
                            lr_probs = 0.99
                        
                        else:
                            predictions_site = [float(line[1])]
                            ID = '_'.join(line[0].split('_')[:-1])
                            counter +=1
                            continue
                    else:
                        vector_prob = convert_p_to_vector_faster(predictions_site)
                        lr_probs = float(model.predict(np.array(vector_prob).reshape(1,99,1)))
                        if lr_probs < 0.99:
                            predictions_site = [float(line[1])]
                            ID = '_'.join(line[0].split('_')[:-1])
                            counter +=1
                            continue

                    # if there are more than 100 reads in a site run 11 times and get the meadian
                    coverage = len(predictions_site)
                    ID_colums = ID.split('_')
                    
                    # write results to output file
                    print(ID_colums[0]+'\t'+ID_colums[1]+'\t'+ID_colums[2]+'\t'+ \
                                   str(coverage)+'\t'+str(stoichiometry)+'\t'+ \
                                   str(lr_probs), file=file_out)
                    
                    predictions_site = [float(line[1])]
                    ID = '_'.join(line[0].split('_')[:-1])
                    counter +=1
            else:
                predictions_site.append(float(line[1]))
                ID = '_'.join(line[0].split('_')[:-1])
            
            if counter % 50000 == 0:
                print(counter,'number of lines processed')


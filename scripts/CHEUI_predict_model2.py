#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 10:16:37 2021

@author: labuser
"""

import argparse

parser = argparse.ArgumentParser(prog='CHEUI_predict_model2 v0.1', description=
                                 """ 
                                 This script takes a per-read prediction file generated using CHEUI_predict_mode1.py \
                                 and generate RNA modification predictions per-site
                                 
                                 """, usage='python CHEUI_predict_model2.py -i <path_to_predictions_model_1> '\
                                            '-m <path_to_DL_model> -o <file_out> \nversion: %(prog)s')

OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')

REQUIRED.add_argument("-i", "--input",
                      help="path to read-level prediction file from CHEUI_predict_model1.py",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-m", "--DL_model",
                      help="path to pretrainned DL model 2",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-c", "--cutoff",
                      help="model 2 probability cutoff for printing sites",
                      metavar='\b',
                      default='0'
                      )

REQUIRED.add_argument("-d", "--double_cutoff",
                      help="Model 1 probability cutoffs used to calculate the stoichiometry",
                      metavar='\b',
                      default='0.3,0.7',
                      )

REQUIRED.add_argument("-n", "--min_reads",
                      help="Minimun number of reads in a site to include in the analysis, ",
                      metavar='\b',
                      default=20
                      )

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
cutoff = float(ARGS.cutoff)
double_cutoff = ARGS.double_cutoff
lower_cutoff = float(double_cutoff.split(',')[0])
upper_cutoff = float(double_cutoff.split(',')[1])
file_out_path = ARGS.file_out
min_reads = int(ARGS.min_reads)


from tensorflow.keras import Input
from tensorflow.keras.models import Model
from DL_models import build_Jasper
import numpy as np
import random
import sys
import os
random.seed(42)


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


def biggerThan100(prob_list):
    '''
    '''
    all_probs = []
    for i in range(11):
        random_reads = random.sample(prob_list, 75)
        vector_prob = convert_p_to_vector(random_reads)
        all_probs.append(model.predict(np.array(vector_prob).reshape(1,99,1)))
    return all_probs


# Load CNN weights
inputs = Input(shape=(99, 1))
output = build_Jasper(inputs, 1)
model = Model(inputs=inputs, outputs=output)
secondML_path = DL_model
model.load_weights(secondML_path)


predictions_site = []
counter = 0
counter_predictions = 0

predictions_dic = {}


if os.path.isfile(file_out_path):
    print('WARNING: site level prediction file already exists, please delete it or change the output name')
    sys.exit()



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
    
                if len(predictions_site) < min_reads:
                    # do not make predictions for this site
                    # store the info. from the current read
                    if line[1] == 'Prediction':
                        break
                    predictions_site = [float(line[1])]
                    ID = '_'.join(line[0].split('_')[:-1])
                    counter +=1
                    
                else:
                    # calculate stoichiometry
                    try:
                        mod = [i for i in predictions_site if i > upper_cutoff]
                        no_mod = [i for i in predictions_site if i < lower_cutoff]
                        stoichiometry = len(mod)/(len(mod)+len(no_mod))
                    except:
                        print('Warning: Stoichiometry cannot compute, try less stringent double cutoff values')
                        stoichiometry = 'None'
                        
                    if len(predictions_site) > 100 and stoichiometry > 0.1:
                        lr_probs = biggerThan100(predictions_site)
                        lr_probs = np.mean(lr_probs)
                        coverage = len(predictions_site)
                        ID_colums = ['_'.join(ID.split('_')[:-2])]+ID.split('_')[-2:]
                        
                        # write results to output file
                        if lr_probs > cutoff:
                            print(ID_colums[0]+'\t'+ID_colums[1]+'\t'+ID_colums[2]+'\t'+ \
                                         str(coverage)+'\t'+str(stoichiometry)+'\t'+ \
                                         str(lr_probs), file=file_out)
                        
                        predictions_site = [float(line[1])]
                        ID = '_'.join(line[0].split('_')[:-1])
                        counter +=1
                        continue
                            
                  
                    else:
                        vector_prob = convert_p_to_vector(predictions_site)
                        # if there are more than 100 reads in a site run 11 times and get the meadian
                        coverage = len(predictions_site)
                        ID_colums = ['_'.join(ID.split('_')[:-2])]+ID.split('_')[-2:]

                        predictions_dic[ID_colums[0]+'\t'+ID_colums[1]+'\t'+ID_colums[2]+'\t'+ \
                                       str(coverage)+'\t'+str(stoichiometry)] = vector_prob
                                        
                        counter_predictions += 1

                        if counter_predictions > 2000:
                            prob_vectors = np.array(list(predictions_dic.values()))
                            ID_vectors = list(predictions_dic.keys())
                            lr_probs = model.predict(prob_vectors)
                            
                            for i in enumerate(lr_probs):
                                if i[1] > cutoff:
                                    print(ID_vectors[i[0]]+'\t'+str(i[1][0]), file=file_out)
                            counter_predictions = 0
                            predictions_dic = {}
                    
                    predictions_site = [float(line[1])]
                    ID = '_'.join(line[0].split('_')[:-1])
                    counter +=1
            else:
                predictions_site.append(float(line[1]))
                ID = '_'.join(line[0].split('_')[:-1])
            
            if counter % 50000 == 0:
                print(counter,'number of lines processed')

        if len(predictions_dic) > 0:
            prob_vectors = np.array(list(predictions_dic.values()))
            ID_vectors = list(predictions_dic.keys())
            lr_probs = model.predict(prob_vectors)
            
            for i in enumerate(lr_probs):
                if i[1] > cutoff:
                    print(ID_vectors[i[0]]+'\t'+str(i[1][0]), file=file_out)


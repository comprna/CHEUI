#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 10:05:12 2021

@author: labuser
"""

import argparse
import yaml
from scipy.stats import ranksums
import itertools
import numpy as np
from scipy import stats


parser = argparse.ArgumentParser(prog='CHEUI_compare v0.1', description=
                                 """ 
                                 This script takes all predictions at single read level and site level and reported 
                                 differential methylation sites \
                                 
                                 
                                 """, usage='python CHEUI_compare.py -c'\
                                            '<config_file> \nversion: %(prog)s')

    
OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')


REQUIRED.add_argument("-c", "--config_file",
                      help="path to the config file",
                      metavar='\b',
                      required=True)

OPTIONAL.add_argument('-v', '--version', 
                        action='version', 
                        version='%(prog)s')

ARGS = parser.parse_args()

# required arg
config_file = ARGS.config_file

with open(config_file, 'r') as stream:
    try:
        config_dic = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print('error at reading yaml file')
        print(exc)

input_file = config_dic['input']
output_file = config_dic['out']

condition1 = list(config_dic['sample_labels']['condition1'].values())
condition2 = list(config_dic['sample_labels']['condition2'].values())
lower_cutoff=config_dic['lower_cutoff']
upper_cutoff=config_dic['upper_cutoff']

def run_tests(predictions_site):
    '''
    Parameters
    ----------
    predictions_site : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    try:
        condition1_values = [predictions_site[i] for i in condition1]
        condition2_values = [predictions_site[i] for i in condition2]
    except:
        return False
    
    condition1_values_f = list(itertools.chain.from_iterable(condition1_values))
    condition2_values_f = list(itertools.chain.from_iterable(condition2_values))
    

    if len(condition1_values_f) > 20 and len(condition2_values_f) > 20:
        
        try:
            mod1 = [i for i in condition1_values_f if i > upper_cutoff]
            no_mod1 = [i for i in condition1_values_f if i < lower_cutoff]
            stoi1 = len(mod1)/(len(mod1)+len(no_mod1))
        except:
            stoi1 = 0

        try:
            mod2 = [i for i in condition2_values_f if i > upper_cutoff]
            mo_mod2 = [i for i in condition2_values_f if i < lower_cutoff]
            stoi2 = len(mod2)/(len(mod2)+len(mo_mod2))
        except:
            stoi2 = 0

        U1, p = ranksums(condition1_values_f,
                         condition2_values_f)
        
        return [len(condition1_values_f), 
                len(condition2_values_f),
                stoi1, stoi2, stoi1 - stoi2,
                U1,  p]
    else:
        return False


counter = 0
predictions_site = {}
print('reading file')
with open(input_file, 'r') as f:
    with open(output_file, 'w') as file_out:
        
        print('ID'+'\t'+'coverage_1'+\
              '\t'+'coverage_2'+'\t'+'stoichiometry_1'+\
              '\t'+'stoichiometry_2'+\
              '\t'+'stoichiometry_diff'+\
              '\t'+'statistic'+\
              '\t'+'pval_U', file=file_out)
        
        for line in f:
            counter +=1
            line_p = line.rstrip().split('\t')
            if counter == 1:
                predictions_site[line_p[-1]] = [float(line_p[1])]
                ID = '_'.join(line_p[0].split('_')[:-1])
            else:
                if ID != '_'.join(line_p[0].split('_')[:-1]):
                    results = run_tests(predictions_site)
                    if results is not False:
                         print(ID+'\t'+'\t'.join(str(x) for x in results),
                               file=file_out)
                         # Empty dictionary
                         predictions_site = {}
                         
                    if line_p[-1] in predictions_site:
                        predictions_site[line_p[-1]] += [float(line_p[1])]
                    else:
                        predictions_site[line_p[-1]] = [float(line_p[1])]
    
                    ID = '_'.join(line_p[0].split('_')[:-1])
                    counter +=1
                else:
                    if line_p[-1] in predictions_site:
                        predictions_site[line_p[-1]] += [float(line_p[1])]
                    else:
                        predictions_site[line_p[-1]] = [float(line_p[1])]
                        
                    ID = '_'.join(line_p[0].split('_')[:-1])
        if predictions_site:
            results = run_tests(predictions_site)
            if results is not False:
                print(ID+'\t'+'\t'.join(str(x) for x in results),
                      file=file_out)

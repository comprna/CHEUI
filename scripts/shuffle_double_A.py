#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 09:31:04 2020

@author: labuser
"""

import pandas as pd

path = '/media/labuser/Data/nanopore/m6A_classifier/data/Epinano/doubleAA/new_kmer_model/'

train_AA_mod = 'train_AA.csv'
train_AA_no_mod = 'train_AA_no_mod.csv'
train_A_no_mod = 'train_A_no_mod.csv'
train_A_mod = 'train_A.csv'

test_AA_mod = 'test_AA.csv.csv'
test_AA_no_mod = 'test_AA_no_mod.csv'
test_A_no_mod = 'test_A_no_mod.csv'
test_A_mod = 'test_A.csv'


train_AA_mod_df =  pd.read_csv(path+'train_AA.csv',
                               sep='\t',
                               ) 

train_AA_no_mod_df =  pd.read_csv(path+'train_AA_no_mod.csv',
                                  sep='\t',
                                  ) 

train_A_no_mod_df =  pd.read_csv(path+'train_A_no_mod.csv',
                                 sep='\t',
                                 ) 

train_A_mod_df =  pd.read_csv(path+'train_A.csv',
                              sep='\t',
                              ) 

train = train_AA_mod_df.append(train_AA_no_mod_df)
train = train.append(train_A_no_mod_df)
train = train.append(train_A_mod_df)

train.sample(frac=1).round(3).to_csv(path+'train_new_model.csv',
                                     mode='w',
                                     sep='\t',
                                     index=None)


test_AA_mod_df = pd.read_csv(path+'test_AA.csv.csv',
                             sep='\t',
                             ) 

test_AA_no_mod_df =  pd.read_csv(path+'test_AA_no_mod.csv',
                                 sep='\t',
                                 )  

test_A_no_mod_df =  pd.read_csv(path+'test_A_no_mod.csv',
                                sep='\t',
                                ) 

test_A_mod_df =  pd.read_csv(path+'test_A.csv',
                             sep='\t',
                             ) 

test = test_AA_mod_df.append(test_AA_no_mod_df)
test = test.append(test_A_no_mod_df)
test = test.append(test_A_mod_df)

test.sample(frac=1).round(3).to_csv(path+'test_new_model.csv',
                                     mode='w',
                                     sep='\t',
                                     index=None)
































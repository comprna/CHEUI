#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 10:34:12 2021

@author: labuser
"""

import numpy as np
from tensorflow.keras import Input
from tensorflow.keras.models import Model
import keras_metrics

from sklearn.metrics import precision_recall_curve
from matplotlib import pyplot

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score

import tensorflow as tf

from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession

import sys
import os

from tensorflow.keras.layers import Dropout, concatenate,\
     BatchNormalization, Softmax, Add, Dense, Activation,\
     Attention, Flatten, LayerNormalization

from tensorflow.keras.activations import relu
import keras.backend as K
import tensorflow as tf
import random

import bisect

random.seed(42)

config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)

def load_predictions(path):
    '''
    load the predictions and probabilities
    '''
    df =  pd.read_csv(path, sep='\t')
    return df


def convert_p_to_vector(probs):
    '''
    '''
    probs = sorted(probs)
    prob_dist = []
    for i in range(1, 10):
        count = 0
        for j in probs:
            if j>=i/10 and j<(i+1)/10:
                count += 1
        prob_dist.append(count)
    return(prob_dist)

        
no_modify_p = '/media/labuser/Data/nanopore/m6A_classifier/results/'+\
    'plugIN_data/MILONGAS_predictions_weight.2.1_unmodified_rep1.csv' 

modify_p = '/media/labuser/Data/nanopore/m6A_classifier/results/'+\
    'plugIN_data/MILONGAS_predictions_weight.2.1_modified_rep1_sorted.csv'

modify_df = load_predictions(modify_p)
no_modify_df = load_predictions(no_modify_p)

unique_ids = no_modify_df.iloc[:,1].unique().tolist()

# prepare the dataset
X_test = []
y_test = []

rounds = 25
for i in range(rounds):
    for id_ in unique_ids:
        # decide WT or KO
        WT = random.randint(0, 1)
        n_reads = random.randint(10, 100)
        stoichiometry = random.randint(10, 100)
        if WT == 0:
            site_df = no_modify_df[no_modify_df.iloc[:,1] == id_]
            try:
                site_df_random = site_df.sample(n=n_reads)
                X_test.append(convert_p_to_vector(site_df_random.iloc[:,2].tolist()))
                y_test.append(0)
            except:
                continue
        else:
            modify_site_df = modify_df[modify_df.iloc[:,1] == id_]
            non_modify_site_df = no_modify_df[no_modify_df.iloc[:,1] == id_]
            
            number_read_modify = int((stoichiometry/100)*n_reads)
            number_read_unmodify = n_reads-number_read_modify
            try:
                modify_site_df_random = modify_site_df.sample(n=number_read_modify)
                non_modify_site_df_random = non_modify_site_df.sample(n=number_read_unmodify)
                
                mix = modify_site_df_random.iloc[:,2].tolist()+non_modify_site_df_random.iloc[:,2].tolist()
                
                X_test.append(convert_p_to_vector(mix))
                y_test.append(1)
            except:
                continue
            

     
### build the network
inputs = Input(shape=(100, 1))

# Conv layer with stride of 2 (halves the size)
x = Dense(25)(inputs)
x = Dropout(rate=0.1)(x)
x = Dense(25)(x)
x = Dropout(rate=0.1)(x)
x = Dense(25)(x)
x = Dropout(rate=0.1)(x)
x = Dense(25)(x)
x = Dropout(rate=0.1)(x)
output = Dense(1, activation='sigmoid')(x)

model = Model(inputs=inputs, outputs=output)
  
optimizer = tf.keras.optimizers.Adam(0.001)
model.compile(optimizer=optimizer,
              loss='binary_crossentropy',
              metrics=['accuracy',
                       keras_metrics.recall(),
                       keras_metrics.precision(),
                       ])
model.summary()




#Train the model on batches
X_train = np.array(X_train)
X_train = (X_train/27).tolist()

#Train the model on batches
history = model.fit(X_train, 
                    y_train,
                    epochs=100,
                    batch_size=48)


history = model.fit(X_test, 
                    y_test, 
                    batch_size=64,
                    epochs=100)
                  
        










































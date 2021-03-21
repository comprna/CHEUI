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
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import accuracy_score

import tensorflow as tf

from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession

import sys
import os

from tensorflow.keras.layers import Dropout, concatenate,\
     BatchNormalization, Softmax, Add, Dense, Activation,\
     Attention, Flatten, LayerNormalization, Conv1D, GlobalAveragePooling1D

from tensorflow.keras.activations import relu
import keras.backend as K
import tensorflow as tf
import random
import pickle
import bisect

from DL_models import build_Jasper, build_deepbinner, build_TCN_cc


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
        # hard to classify examples 
        stoichiometry = random.randint(10, 50)
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
            

pickle.dump( X_test, open( "/media/labuser/Data/nanopore/m6A_classifier/plugIN_data/plugIN_X_train.p", "wb" ))
pickle.dump( y_test, open("/media/labuser/Data/nanopore/m6A_classifier/plugIN_data/plugIN_y_train.p", "wb" ))

with open("/media/labuser/Data/nanopore/m6A_classifier/plugIN_data/plugIN/plugIN_X_train_100_stoi100.p", 'rb') as pickle_file:
    train = pickle.load(pickle_file)
with open("/media/labuser/Data/nanopore/m6A_classifier/plugIN_data/plugIN/plugIN_y_train_100_stoi100.p", 'rb') as pickle_file:
    test = pickle.load(pickle_file)


##### try other algorithms


import xgboost as xgb

data_dmatrix = xgb.DMatrix(data=np.array(X_train_b),label=y_train_b)
xg_reg = xgb.XGBRegressor(objective ='binary:logistic', colsample_bytree = 0.3, learning_rate = 0.1,
                          max_depth = 10, alpha = 10, n_estimators = 100)

xg_reg.fit(np.array(X_train_b),y_train_b)
preds = xg_reg.predict(np.array(X_test_b))
predictions_binary = [1 if n >= 0.5 else 0 for n in preds]

print('Classifier ', classifier)
print('Accuracy', accuracy_score(y_test_b, predictions_binary))
print('Recall', recall_score(y_test_b, predictions_binary ))
print('Precision', precision_score(y_test_b, predictions_binary))
print()


from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

classifiers = [
                KNeighborsClassifier(10),
                SVC(gamma='scale', C=1),
                #GaussianProcessClassifier(1.0 * RBF(1.0)),
                RandomForestClassifier(n_estimators=500),
                MLPClassifier(hidden_layer_sizes=200, alpha=1, max_iter=1000),
                AdaBoostClassifier(n_estimators=200),
                GaussianNB(),
                QuadraticDiscriminantAnalysis()]


train_normal = StandardScaler().fit_transform(train)

X_train_b, X_test_b, y_train_b, y_test_b = train_test_split(
                                            train, test, test_size=0.2, random_state=42)


for classifier in classifiers:
    clf = classifier.fit(X_train_b, y_train_b)
    predictions = clf.predict(X_test_b)
    predictions_binary = [1 if n >= 0.5 else 0 for n in predictions]
    print('Classifier ', classifier)
    print('Accuracy', accuracy_score(y_test_b, predictions_binary))
    print('Recall', recall_score(y_test_b, predictions_binary ))
    print('Precision', precision_score(y_test_b, predictions_binary))
    print()
    
#scores = cross_validate(MLP, X_test, y_test, cv=3, scoring=['accuracy', 'precision'])

   
### build the network
inputs = Input(shape=(9))   

# Conv layer with stride of 2 (halves the size)
x = Dense(25)(inputs)
x = Dropout(rate=0.1)(x)
x = BatchNormalization()(x)
x = Dense(50)(x)
x = Dropout(rate=0.1)(x)
x = BatchNormalization()(x)
x = Dense(100)(x)
x = Dropout(rate=0.1)(x)
x = BatchNormalization()(x)
x = Dense(25)(x)
x = Dropout(rate=0.1)(x)
output = Dense(1, activation='sigmoid')(x)

inputs = Input(shape=(9, 1))   
# Conv group: 3 layers of 3-kernels
x = Conv1D(filters=48, kernel_size=2, padding='same', activation='relu')(inputs)
x = Conv1D(filters=48, kernel_size=2, padding='same', activation='relu')(x)
x = Conv1D(filters=48, kernel_size=2, padding='same', activation='relu')(x)
x = GlobalAveragePooling1D()(x)
output = Dense(1, activation='sigmoid')(x)


inputs = Input(shape=(9, 1))
output = build_Jasper(inputs, 1)
model = Model(inputs=inputs, outputs=output)

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
history = model.fit(np.array(X_train_b), 
                    np.array(y_train_b), 
                    batch_size=48,
                    validation_data=(np.array(X_test_b),
                                     np.array(y_test_b)),
                    validation_batch_size=48,
                    epochs=10)























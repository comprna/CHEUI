#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 20:29:58 2020

@author: pablo
"""

import os
import pickle

import numpy as np
#from tensorflow.keras import Input
#from tensorflow.keras import Model
#from tensorflow.keras.layers import Dense
#from tensorflow.keras.preprocessing.sequence import pad_sequences
#from keras import backend as K
from sklearn import preprocessing
#from sklearn.utils import shuffle
import itertools
import umap
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.lines import Line2D

#from tcn import TCN

#import tensorflow as tf
#assert tf.test.is_gpu_available()
#assert tf.test.is_built_with_cuda()

from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession

config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)


def load_data_events(directory):
    '''
    '''
    files = os.listdir(directory)
    dic_signal = {}
    for file in files:
        if 'events' in file:
            dic_signal[file] = pickle.load(open(directory+file,"rb" ))
    
    return dic_signal


def load_data(directory):
    '''
    '''
    files = os.listdir(directory)
    dic_signal = {}
    for file in files:
        if 'events' in file:
            continue
        else:
            dic_signal[file] = pickle.load(open(directory+file,"rb" ))
    
    return dic_signal


def recall_m(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = true_positives / (possible_positives + K.epsilon())
    return recall

def precision_m(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision

def f1_m(y_true, y_pred):
    precision = precision_m(y_true, y_pred)
    recall = recall_m(y_true, y_pred)
    return 2*((precision*recall)/(precision+recall+K.epsilon()))


def Min_Max_scaler(signal):
    '''
    Function to scale the data
    '''
    min_max_scaler = preprocessing.MinMaxScaler()
    #signal = signal.reshape((len(signal),1))
    signal = min_max_scaler.fit_transform(signal)
    return signal


def plot_signals_10(KO_signal_filtered_np, WT_signal_filtered_np, kmer, save=None):
    '''
    This function wil plot the median and sd of a bunch of signals
    '''
    # first build the whole signals with medians
    sns.set_style('darkgrid')

    KO_signal = KO_signal_filtered_np.transpose().ravel()
    KO_signal_index = np.arange(1, 51)
    KO_signal_index = np.repeat(KO_signal_index, KO_signal_filtered_np.shape[0])
    
    KO_df = pd.DataFrame({'events' : KO_signal_index, 
                          'signal' : KO_signal})
    
    
    WT_signal = WT_signal_filtered_np.transpose().ravel()
    WT_signal_index = np.arange(1, 51)
    WT_signal_index = np.repeat(WT_signal_index, WT_signal_filtered_np.shape[0])
    
    WT_df = pd.DataFrame({'events' : WT_signal_index, 
                          'signal' : WT_signal})
    plt.figure(figsize=(15, 8))
    ax = sns.lineplot(x="events", y="signal", ci='sd', estimator=np.median, data=WT_df, color='red')
    ax2 = sns.lineplot(x="events", y="signal", ci='sd', estimator=np.median, data=KO_df, color='blue')
    ax2.lines[0].set_linestyle("--")
    plt.xticks([])
    plt.yticks([])
    plt.ylabel('Signal intensity', fontsize=30)
    X = [5, 15, 25, 35, 45]
    labels = list(kmer.split('_')[5])
    # You can specify a rotation for the tick labels in degrees or with keywords.
    plt.xticks(X, labels, fontsize=30)
    
    custom_lines = [Line2D([0], [0], color='blue', lw=4, markersize=40),
                    Line2D([0], [0], color='red', lw=4, markersize=40, linestyle='--')]
    
    plt.legend(custom_lines,
               ['NO MOD','MOD'],
               fancybox=True,
               framealpha=1, 
               shadow=True, 
               borderpad=1,
               fontsize=20
               )
    if save:
        plt.title('Signal RNA IVT')
        plt.savefig(save+'.pdf',
                    format='pdf',
                    dpi=1200)
    plt.close('all')
    
    return True


mod_rep_1 = '/media/labuser/Data/nanopore/m6A_classifier/data/mod_rep1_eventalign_numpy_sites/'
mod_rep_2 = '/media/labuser/Data/nanopore/m6A_classifier/data/mod_rep2_eventalign_numpy_sites/'
no_mod_rep_1 = '/media/labuser/Data/nanopore/m6A_classifier/data/no_mod_rep1_eventalign_numpy_sites/'
no_mod_rep_2 = '/media/labuser/Data/nanopore/m6A_classifier/data/no_mod_rep2_eventalign_numpy_sites/'


mod_rep_1_signal_events = load_data_events(mod_rep_1)
#mod_rep_2_signal_events = load_data_events(mod_rep_2)
no_mod_rep_1_signal_events = load_data_events(no_mod_rep_1)
#no_mod_rep_2_signal_events = load_data_events(no_mod_rep_2)


for i in mod_rep_1_signal_events.keys():
    
    # create the diretory
    out_folder = '/media/labuser/Data/nanopore/m6A_classifier/data/plots_rep1/'
    out_folder += i.split('_')[5]
    
    
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    
    mod_list = []
    for int_list in mod_rep_1_signal_events[i]:
        mod_list.append(list(itertools.chain(*int_list)))
    
    no_mod_list = []
    for int_list in no_mod_rep_1_signal_events[i]:
        no_mod_list.append(list(itertools.chain(*int_list)))
    
    
    mod_list = np.array(mod_list[:500])
    no_mod_list = np.array(no_mod_list[:500])
    
    #### visualize data
    
    plot_signals_10(no_mod_list, mod_list, i, out_folder+'/'+i.split('_')[4]+'signal.pdf')
    #plot_signals_10(no_mod_list, mod_list, out_folder+'/'+i.split('_')[4]+'signal.pdf')
    
    
    X = Min_Max_scaler(np.concatenate([mod_list, no_mod_list]))
    
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(X)
    
    plt.figure(figsize=(15, 8))
    sns.scatterplot(embedding[:len(mod_list),0],
                embedding[:len(mod_list),1],
                color='b')
    
    sns.scatterplot(embedding[len(mod_list):,0],
                embedding[len(mod_list):,1],
                color='r')
   
    plt.savefig(out_folder+'/'+i.split('_')[4]+'UMAPembedding'+'.pdf',
                format='pdf',
                dpi=1200)
    plt.close('all')


'''


train_X = []
train_y = []

mod = 0
no_mod = 0

for i in mod_rep_2_signal.keys():
    for j in mod_rep_2_signal[i]:
        if len(j) <= 1000:    
            train_X.append(j)
            train_y.append(1)
            mod +=1

for i in no_mod_rep_2_signal.keys():
    for j in no_mod_rep_2_signal[i]:
        if len(j) <= 1000:    
            train_X.append(j)
            train_y.append(0)
            no_mod +=1

test_X = []
test_y = []

for i in mod_rep_1_signal.keys():
    for j in mod_rep_1_signal[i]:
        if len(j) <= 1000:    
            test_X.append(j)
            test_y.append(1)

for i in no_mod_rep_1_signal.keys():
    for j in no_mod_rep_1_signal[i]:
        if len(j) <= 1000:    
            test_X.append(j)
            test_y.append(0)

# pad sequences
train_X_pad = pad_sequences(train_X, padding='post')
test_X_pad = pad_sequences(test_X, padding='post')

train_X_pad, train_y = shuffle(train_X_pad, train_y)

train_y, test_y = np.array(train_y), np.array(test_y) 

i = Input(batch_shape=(None, 1000, 1)) # Batch_size, time-steps, input_dim
o = TCN(nb_filters=64,  
        kernel_size=5,
        return_sequences=True,
        dropout_rate=0.2,
        padding='same')(i)  # The TCN layers are here.

o = TCN(nb_filters=64,  
        kernel_size=5,
        return_sequences=True,
        dropout_rate=0.2,
        padding='same')(o)  # The TCN layers are here.

o = TCN(nb_filters=64,  
        kernel_size=5,
        return_sequences=True,
        dropout_rate=0.2,
        padding='same')(o)  # The TCN layers are here.

o = Dense(1, activation='sigmoid')(o)

m = Model(inputs=[i], outputs=[o])
m.compile(optimizer='adam', 
          loss='binary_crossentropy',
          metrics=['acc',f1_m,precision_m, recall_m])

m.fit(train_X_pad, train_y,
      batch_size=32,
      epochs=1,
      validation_data=(test_X_pad, test_y)
      )

m.save('./m6A_classifier.h5')


'''













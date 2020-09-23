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
from keras import backend as K
from sklearn import preprocessing
#from sklearn.utils import shuffle
import itertools
import umap
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.lines import Line2D

#from tcn import TCN

import tensorflow as tf
#assert tf.test.is_gpu_available()
#assert tf.test.is_built_with_cuda()

from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession

from scipy.spatial import distance
from scipy import signal

config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)


def load_data_smooth(directory, model_kmer_dict):
    '''
    '''
    files = os.listdir(directory)
    dic_signal = {} # dictionary of signal event values
    dic_dwell = {} # dictionary of dwelling time values
    dic_distance = {} # dictionary of distance from the non-modified to the signal
    counter = 0
    for file in files:
        counter +=1
        if 'events' in file:
            file_raw = pickle.load(open(directory+file,"rb"))
            for signal_values in file_raw:
                signal_smoothed = []
                dwell_smoothed = []
                for event in signal_values:
                    if len(event) > 1000: # maximun event lenght
                        dwell_smoothed += [1000]
                    else:
                        dwell_smoothed += [len(event)] # record the original dwelling time
                    event_smoothed = smooth_event(event, 20) # smooth the event
                    signal_smoothed += event_smoothed  # add the event to the signal
                    
                if file in dic_signal:
                    expected_smoothed = make_expected(model_kmer_dict, 
                                                    file.split('_')[-2],
                                                    20)
                    dic_distance[file] += [distance_calculator(expected_smoothed,
                                                               signal_smoothed)]
                    dic_signal[file] += [signal_smoothed]
                    dic_dwell[file] += [dwell_smoothed]
                else:
                    expected_smoothed = make_expected(model_kmer_dict, 
                                                    file.split('_')[-2],
                                                    20)
                    dic_distance[file] = [distance_calculator(expected_smoothed,
                                                         signal_smoothed)]
                    dic_signal[file] = [signal_smoothed]
                    dic_dwell[file] = [dwell_smoothed]
                    
        if counter == 5:
            return dic_signal, dic_dwell, dic_distance


def make_expected(model_kmer_dict, kmer, event_lenght):
    '''
    '''
    expected_signal = []
    for i in range(5):
        expected_signal += [model_kmer_dict[kmer[i:i+5]]]*event_lenght
    return expected_signal

def distance_calculator(signal_expected, event_smoothed):
    '''
    '''
    #vector_distance = list(np.round(abs(np.array(signal_expected) - \
    #                                    np.array(event_smoothed)), 3))
    vector_distance = list(signal.correlate(signal_expected, 
                                            event_smoothed, 
                                            mode='same',
                                            method='fft'))
    return vector_distance

def smooth_event(raw_signal, lenght_events):
    '''
    smmoth the signal 
    '''
    raw_signal_events = []
    
    if len(raw_signal) < lenght_events:
        event = top_median(raw_signal, lenght_events)
        raw_signal_events = [round(i, 3) for i in event]
        
    else:
        division = floor(len(raw_signal)/lenght_events)
        new_event = []
        for i in range(0, len(raw_signal), division):
            new_event.append(np.median(raw_signal[i:i+division]))
            if len(new_event) == lenght_events:
                break
        if len(new_event) < lenght_events:
            new_event = top_median(new_event, lenght_events)
        raw_signal_events = [round(i, 3) for i in new_event]
    return raw_signal_events


def top_median(array, lenght):
    '''
    This function top an array until some specific lenght
    '''
    extra_measure = [np.median(array)]*(lenght-len(array))
    array += extra_measure
    return array


def load_data(directory):
    '''
    '''
    files = os.listdir(directory)
    dic_signal = {}
    for file in files:
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

if __name__ == '__main__':
    
    mod_rep_1 = '/media/labuser/Data/nanopore/m6A_classifier/data/Epinano/mod_rep1_eventalign_numpy_sites_AA/'
    no_mod_rep_1 = '/media/labuser/Data/nanopore/m6A_classifier/data/Epinano/no_mod_rep1_eventalign_numpy_sites_AA/'
    
    model_kmer = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/model_kmer.csv',
                             sep=',')
    
    model_kmer_dict = dict(zip(model_kmer['model_kmer'], model_kmer['model_mean']))

    mod_rep_1_signal_events, \
    mod_rep_1_signal_dwell, \
    mod_rep_1_signal_distances = load_data_smooth(mod_rep_1, model_kmer_dict)

    
    no_mod_rep_1_signal_events, \
    no_mod_rep_1_signal_dwell, \
    no_mod_rep_1_signal_distances = load_data_smooth(no_mod_rep_1, model_kmer_dict)
    
    mod_rep_1_signal_events, \
    mod_rep_1_signal_dwell, \
    mod_rep_1_signal_distances_cross = load_data_smooth(mod_rep_1, model_kmer_dict)
    
    no_mod_rep_1_signal_events, \
    no_mod_rep_1_signal_dwell, \
    no_mod_rep_1_signal_distances_cross = load_data_smooth(no_mod_rep_1, model_kmer_dict)
    
    mod_rep_1_signal_distances_flat = [item for sublist in mod_rep_1_signal_distances.values() for item in sublist]
    mod_rep_1_signal_distances_flat = [item for sublist in mod_rep_1_signal_distances_flat for item in sublist]

    no_mod_rep_1_signal_distances_flat = [item for sublist in no_mod_rep_1_signal_distances.values() for item in sublist]
    no_mod_rep_1_signal_distances_flat = [item for sublist in no_mod_rep_1_signal_distances_flat for item in sublist]

    mod_rep_1_signal_distances_cross_flat = [item for sublist in mod_rep_1_signal_distances_cross.values() for item in sublist]
    mod_rep_1_signal_distances_cross_flat = [item for sublist in mod_rep_1_signal_distances_cross_flat for item in sublist]
    
    no_mod_rep_1_signal_distances_cross_flat = [item for sublist in no_mod_rep_1_signal_distances_cross.values() for item in sublist]
    no_mod_rep_1_signal_distances_cross_flat = [item for sublist in no_mod_rep_1_signal_distances_cross_flat for item in sublist]
    
    np.median(mod_rep_1_signal_distances_flat)
    np.std(mod_rep_1_signal_distances_flat)

    plt.title('mod_rep_1_signal_distances_flat')
    g = sns.histplot(mod_rep_1_signal_distances_flat)
    g.set(xlim=(0, 50))
    
    np.median(no_mod_rep_1_signal_distances_flat)
    np.std(no_mod_rep_1_signal_distances_flat)

    plt.title('no_mod_rep_1_signal_distances_flat')
    g = sns.histplot(no_mod_rep_1_signal_distances_flat)
    g.set(xlim=(0, 50))
    
    np.median(mod_rep_1_signal_distances_cross_flat)
    plt.title('mod_rep_1_signal_distances_cross_flat')
    g = sns.histplot(mod_rep_1_signal_distances_cross_flat)
    #g.set(xlim=(0, 50))
    
    np.median(no_mod_rep_1_signal_distances_cross_flat)
    plt.title('no_mod_rep_1_signal_distances_cross_flat')
    g = sns.histplot(no_mod_rep_1_signal_distances_cross_flat)
    # g.set(xlim=(0, 50))
    
    
    for i in mod_rep_1_signal_events.keys():
        
        # create the diretory
        out_folder = '/media/labuser/Data/nanopore/m6A_classifier/data/plots_rep1/'
        out_folder += i.split('_')[5]
        
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)
        
        i = 'cc6m_2709_t7_ecorv_60_GGACC_events'
        
        mod_list = []
        for int_list in mod_rep_1_signal_events[i]:
            mod_list.append(list(itertools.chain(*int_list)))
        
        no_mod_list = []
        for int_list in no_mod_rep_1_signal_events[i]:
            no_mod_list.append(list(itertools.chain(*int_list)))
        
        mod_list = np.array(mod_list[:500])
        no_mod_list = np.array(no_mod_list[:500])
        
        #### visualize data
        
        plot_signals_10(no_mod_list, mod_list, i)
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
    
    
    
    
    








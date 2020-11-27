#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 20:29:58 2020

@author: pablo
"""

import numpy as np
from tensorflow.keras import Input
from tensorflow.keras.models import Model
import keras_metrics
from tensorflow import keras
from sklearn.metrics import precision_recall_curve
from matplotlib import pyplot
from sklearn.neighbors import NearestNeighbors

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from sklearn.model_selection import train_test_split

from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import auc

import tensorflow as tf

from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession

from DL_models import build_Jasper, build_deepbinner, build_TCN, build_TCN_cc
import os
import sys
import pickle
import itertools


config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)

tf.random.set_seed(42)

def parse_chunk(chunk):
    '''
    '''
    events = chunk['event']
    distances = chunk['distances']
    
    events = np.array([np.array(xi) for xi in events])
    distances = np.array([np.array(xi) for xi in distances])

    combined = np.concatenate((events,
                              distances), 
                              axis=1)
    
    # reshape the vectors
    combined = combined.reshape(len(combined), 2, 100)
    
    # transpose the last two dimensions to have each channel with separate information
    combined = np.transpose(combined, (0, 2, 1))
    
    return combined


def make_expected(model_kmer_dict, kmer, event_lenght):
    '''
    '''
    expected_signal = []
    for i in range(5):
        expected_signal += [model_kmer_dict[kmer[i]]]*event_lenght
    return expected_signal


def distance_calculator(signal_expected, event_smoothed):
    '''
    '''
    vector_distance = list(np.round(abs(np.array(signal_expected) - \
                                        np.array(event_smoothed)), 3))
    
    #vector_distance_eu =    list(np.round(abs(np.array(signal_expected) - \
    #                                    np.array(event_smoothed)), 3))
    
    return vector_distance


def load_local_stored(file_path):
    '''
    '''
    data = []   
    with open(file_path, 'rb') as fr:
        try:
            while True:
                data.append(pickle.load(fr))
        except:
            pass
    return data


def number_to_sequence(sequences_n, model_kmer_dict):
    '''
    convert the numbers into sequences
    '''
    numbers_indexes = [0,20,40,60,80]
    sequences = []
    for sequence in sequences_n:
        sequences_n_min = sequence[numbers_indexes]
        sequences.append([model_kmer_dict[i] for i in sequences_n_min])

    return sequences


def gaussian_noise(X_train, 
                   y_train,
                   number_augment,
                   sequences_n, 
                   model_kmer_dict,
                   tokenizer,
                   gaussian_std):
    '''
    '''
    X_train_total = X_train
    y_train_total = y_train
    
    sequences = number_to_sequence(sequences_n, model_kmer_dict)
    expected = [make_expected(tokenizer, i, 20) for i in sequences]
    
    for i in range(number_augment):
        
        noise_plut_events = []
        # create a noise vector and add it to the original one
        noise = np.random.normal(0, gaussian_std, (X_train[:,:,0].shape))
        noise_plut_events = X_train[:,:,0] + noise
        noise_plut_events_distances = distance_calculator(expected, noise_plut_events)
        
        # create a numpy of 0s
        X2_train = np.zeros((len(X_train), 100, 2))
        # fill the array
        X2_train[:,:,0] = noise_plut_events
        X2_train[:,:,1] = noise_plut_events_distances
        
        # concatenate with the real data
        X_train_total = np.concatenate((X_train_total, X2_train), axis=0)
        y_train_total = np.concatenate((y_train_total, y_train), axis=0)
    
    return X_train_total, y_train_total


def formula_extrapolation(cj, ck, lambda_p=0.5):
    '''
    c0j = (cj − ck)λ + cj
    
    λ is the desgree of extrapolation a value in the range {0, ∞}
    λ could be drawn from a random distribution for each new sample
    '''
    return (cj - ck)*lambda_p + cj


def KNN(X_train):
    '''
    make KNN among all signals
    '''
    
    nbrs = NearestNeighbors(n_neighbors=5, algorithm='ball_tree').fit(X_train[:,:,0])
    distances, indices = nbrs.kneighbors(X_train[:,:,0])
    

def extrapolation(X_train,
                  y_train,
                  sequences_n, 
                  model_kmer_dict,
                  tokenizer):
    '''
  
    '''
    sequences = number_to_sequence(sequences_n, model_kmer_dict)
    expected = [make_expected(tokenizer, i, 20) for i in sequences]

    # create all possible combinations two by two
    element_combination = list(itertools.combinations(X_train , 2))
    extrapolated_data = []
    for pair in element_combination:
        extrapolated_data.append(formula_extrapolation(pair[0], pair[1]))
        extrapolated_data.append(formula_extrapolation(pair[1], pair[0]))
        
    noise_plut_events_distances = distance_calculator(expected, extrapolated_data)
    
    # create a numpy of 0s
    X2_train = np.zeros((len(X_train), 100, 2))
    # fill the array
    X2_train[:,:,0] = extrapolated_data
    X2_train[:,:,1] = noise_plut_events_distances
    
    # concatenate with the real data
    X_train = np.concatenate((X_train, X2_train), axis=0)
    y_train = np.concatenate((y_train, y_train), axis=0)

    return X_train, y_train
    

def bootstrap():
    '''
    '''

if __name__ == '__main__':
    
    #train_path = sys.argv[1]
    #test_path =  sys.argv[2]
    #OUT_FOLDER = sys.argv[3]
        
    test_path =  '/media/labuser/Data/nanopore/m6A_classifier/data/Epinano/singleA/test_singleA.csv'

    # parse the testing files
    test = pd.read_csv(test_path, sep='\t', converters={'event': eval,
                                                        'distances': eval})
    X_test = parse_chunk(test)
    y_test = test['label']
    
    for gaussian_std in [0, 0.5, 1, 1.5, 2]:
        
        multiply_data_times = 10
        
        # process the training data
        directory = '/media/labuser/Data/nanopore/m6A_classifier/data/Epinano/singleA/'
        mod_rep_2_signal_events = load_local_stored(directory+'mod_rep2_eventalign_numpy_sitesevent_30.npy')
        mod_rep_2_signal_distances = load_local_stored(directory+'mod_rep2_eventalign_numpy_sitesdistances_euclidean_30.npy')
        mod_rep_2_signal_sequences = load_local_stored(directory+'mod_rep2_eventalign_numpy_sitessequences_30.npy')
    
        no_mod_rep_2_signal_events = load_local_stored(directory+'no_mod_rep2_eventalign_numpy_sitesevent_30.npy')
        no_mod_rep_2_signal_distances = load_local_stored(directory+'no_mod_rep2_eventalign_numpy_sitesdistances_euclidean_30.npy')
        no_mod_rep_2_signal_sequences = load_local_stored(directory+'mod_rep2_eventalign_numpy_sitessequences_30.npy')
        
        # concatenate input vectors
        no_mod_input = np.concatenate([no_mod_rep_2_signal_events,
                                       no_mod_rep_2_signal_distances],
                                       axis=1)
        
        mod_input = np.concatenate([mod_rep_2_signal_events,
                                    mod_rep_2_signal_distances
                                    ],
                                    axis=1)
        
        no_mod_input = no_mod_input.reshape(len(no_mod_input), 2, 100)
        mod_input = mod_input.reshape(len(mod_input), 2, 100)
        
        # transpose the last two dimensions to have each channel with separate information
        no_mod_input =  np.transpose(no_mod_input, (0, 2, 1))
        mod_input = np.transpose(mod_input, (0, 2, 1))
        
        
        total_data = np.concatenate([no_mod_input, mod_input],
                                    axis=0)
        
        # Create y
        target = np.concatenate([np.zeros(len(no_mod_input)),
                                 np.ones(len(mod_input))],
                                 axis=0)
        
        # Use all replicate 2 as training data
        X_train, y_train = total_data, target
        
    
        
        
        #################### Data Agumentation #####################
        
        # I need to know the sequences of nucleotides to calculate the expected sequences
        model_kmer = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/model_kmer.csv',
                                 sep=',')
        
        # Create a column with the numbers that we will use as tokens from 0 to 1
        model_kmer['number'] = [round((i+1)/1025, 5) for i in range(len(model_kmer))]
    
        model_kmer_dict = dict(zip(model_kmer['number'] , model_kmer['model_kmer']))
        tokenizer = dict(zip(model_kmer['model_kmer'], model_kmer['model_mean']))
    
        sequences_n = np.concatenate([no_mod_rep_2_signal_sequences,
                                      mod_rep_2_signal_sequences])
        
        # gaussian noisse
   
        X_train, y_train = gaussian_noise(X_train,
                                          y_train,
                                          multiply_data_times,
                                          sequences_n,
                                          model_kmer_dict,
                                          tokenizer,
                                          gaussian_std)
        
        # Define model Input
        inputs = Input(shape=(100, 2))
        output = build_Jasper(inputs,Deep=True)
        model = Model(inputs=inputs, outputs=output)
        
        # load pre-train weights
        model.load_weights('/media/labuser/Data/nanopore/m6A_classifier/results/doubleA/jasper2vectors/93600000model.h5')
        
        # freeze layers until the last 3 convolutions
        for i in enumerate(model.layers[:-9]):
            model.layers[i[0]].trainable = False
        
        model.compile(optimizer='adam',
                      loss='binary_crossentropy',
                      metrics=['accuracy',
                               keras_metrics.recall(),
                               keras_metrics.precision(),
                               ])
        
        checkpoint_filepath = '/media/labuser/Data/nanopore/m6A_classifier/results/singleA/singleA30TL_aug_noise_100_'+str(gaussian_std)+'/'
        
        if not os.path.exists(checkpoint_filepath):
            os.makedirs(checkpoint_filepath)
        
        
        model_checkpoint_callback = tf.keras.callbacks.ModelCheckpoint(
                                                                       filepath=checkpoint_filepath,
                                                                       save_weights_only=True,
                                                                       monitor='val_acc',
                                                                       mode='max',
                                                                       save_best_only=True)
        #Train the model on batches
        history = model.fit(X_train, 
                            y_train, 
                            batch_size=125,
                            validation_data=(X_test,
                                             y_test),
                            validation_batch_size=125,
                            epochs=50,
                            callbacks=[model_checkpoint_callback]
                            )
        acc = []
        val_acc = []
        loss = []
        val_loss = []
        precision = []
        val_precision = []
        recall = []
        val_recall = []
        
        acc = history.history['accuracy']
        val_acc = history.history['val_accuracy']
        loss = history.history['loss']
        val_loss = history.history['val_loss']
        precision = history.history['precision']
        val_precision = history.history['val_precision']
        recall = history.history['recall']
        val_recall = history.history['val_recall']
        
        # summarize history for accuracy
        f, ax = plt.subplots( figsize=(13,9))
        sns.lineplot(x=np.arange(len(acc)), y=np.array(acc), palette="tab10", linewidth=2.5, label='Accuracy')
        sns.lineplot(x=np.arange(len(acc)), y=np.array(val_acc), palette="tab10", linewidth=2.5, label='Val accuracy')
        plt.ylabel('Accuracies', fontsize=17)
        plt.savefig(checkpoint_filepath+'/accuracy.pdf',
                    format='pdf',
                    dpi=1200,
                    bbox_inches='tight', pad_inches=0)
        
        f, ax = plt.subplots( figsize=(13,9))
        sns.lineplot(x=np.arange(len(acc)), y=np.array(loss), palette="tab10", linewidth=2.5, label='loss')
        sns.lineplot(x=np.arange(len(acc)), y=np.array(val_loss), palette="tab10", linewidth=2.5, label='Val loss')
        plt.ylabel('Loss', fontsize=17)
        plt.savefig(checkpoint_filepath+'/Loss.pdf',
                    format='pdf',
                    dpi=1200,
                    bbox_inches='tight', pad_inches=0)
        
        f, ax = plt.subplots( figsize=(13,9))
        sns.lineplot(x=np.arange(len(acc)), y=np.array(precision), palette="tab10", linewidth=2.5, label='Precision')
        sns.lineplot(x=np.arange(len(acc)), y=np.array(val_precision), palette="tab10", linewidth=2.5, label='Val precision')
        plt.ylabel('Precision', fontsize=17)
        plt.savefig(checkpoint_filepath+'/Precision.pdf',
                    format='pdf',
                    dpi=1200,
                    bbox_inches='tight', pad_inches=0)
        
        f, ax = plt.subplots( figsize=(13,9))
        sns.lineplot(x=np.arange(len(acc)), y=np.array(recall), palette="tab10", linewidth=2.5, label='Recall')
        sns.lineplot(x=np.arange(len(acc)), y=np.array(val_recall), palette="tab10", linewidth=2.5, label='Val Recall')
        plt.ylabel('Recall', fontsize=17)
        plt.savefig(checkpoint_filepath+'/Recall.pdf',
                    format='pdf',
                    dpi=1200,
                    bbox_inches='tight', pad_inches=0)
        
        
        lr_probs = model.predict(X_test)
       
        lr_auc = roc_auc_score(y_test, lr_probs)
        # calculate roc curves
        lr_fpr, lr_tpr, _ = roc_curve(y_test, lr_probs)
        # plot the roc curve for the model
        f, ax = plt.subplots( figsize=(13,9))
        pyplot.plot(lr_fpr, lr_tpr, marker='.', label='Logistic')
        plt.ylabel('True Positive Rate', fontsize=17)
        plt.xlabel('False Positive Rate', fontsize=17)
        plt.text(0.7, 0.2, 'AUC: '+str(round(lr_auc, 3)), fontsize=27)
        plt.savefig(checkpoint_filepath+'/ROC_AUC.pdf',
                    format='pdf',
                    dpi=1200,
                    bbox_inches='tight', pad_inches=0)
    
        # predict probabilities
        # predict class values
        lr_precision, lr_recall, _ = precision_recall_curve(y_test, lr_probs)
        
        RP_auc = auc(lr_recall, lr_precision)
        
        f, ax = plt.subplots( figsize=(13,9))
        pyplot.plot(lr_recall, lr_precision, marker='.', label='Logistic')  
        plt.ylabel('Precision', fontsize=17)
        plt.xlabel('Recall', fontsize=17)
        plt.text(0.2, 0.7, 'RP_AUC: '+str(round(RP_auc, 3)), fontsize=27)
        plt.savefig(checkpoint_filepath+'/precision_recall_curve.pdf',
                    format='pdf',
                    dpi=1200,
                    bbox_inches='tight', pad_inches=0)
        
        
        
        

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
from scipy import stats
from random import shuffle

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
    vector_distance = np.round(abs(np.array(signal_expected) - \
                                        np.array(event_smoothed)), 3)
    
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


def add_permutation(signal,
                   number_augment,
                   sequences_n, 
                   model_kmer_dict,
                   tokenizer,
                  ):
    '''
    
    '''
    sequences = number_to_sequence(np.array(sequences_n), model_kmer_dict)
    expected = [make_expected(tokenizer, i, 20) for i in sequences]
    sequences_np = np.array(sequences_n)
    
    temp_same_motif = []
    temp_expected = []
    
    augmented_signals = signal
    augmented_signals_expected = expected
    
    prev = sequences_np[0]
    for i in enumerate(sequences_np):
        if (i[1] == prev).all() == True:
            temp_same_motif.append(signal[i[0]])
            temp_expected.append(expected[i[0]])
            prev = i[1]
        else:
            break
            # first clone the same array n times            
            temp_expected_aug =  np.tile(temp_expected,(number_augment,1))

            augmented_permutation = []
            for j in range(0,100,20):
                temp_augmented_permutation = []
                for k in range(number_augment):
                    ori_shape = np.array(temp_same_motif)[:,j:j+20].shape
                    temp_signals_event = np.array(temp_same_motif)[:,j:j+20]
                    temp_signals_event = temp_signals_event.reshape(len(temp_signals_event)*temp_signals_event.shape[1])
                    shuffle(temp_signals_event)
                    temp_signals_event = temp_signals_event.reshape(ori_shape)
                    
                    if len(temp_augmented_permutation) == 0:
                        temp_augmented_permutation = temp_signals_event
                    else:
                        temp_augmented_permutation = np.concatenate((temp_augmented_permutation,
                                                                     temp_signals_event), axis=0)
                if len(augmented_permutation) ==0:
                    augmented_permutation = temp_augmented_permutation
                else:
                    augmented_permutation = np.concatenate((augmented_permutation,
                                                            temp_augmented_permutation), axis=1)
                
            augmented_signals += temp_augmented_permutation.tolist()
            augmented_signals_expected += temp_expected_aug.tolist()
            
            
            temp_same_motif = []
            temp_expected = []
            temp_same_motif.append(signal[i[0]])
            temp_expected.append(expected[i[0]])

            prev = i[1]
        
    return augmented_signals, distance_calculator(augmented_signals_expected, augmented_signals)



def bootstrap():
    '''
    '''

if __name__ == '__main__':
    
    #train_path = sys.argv[1]
    #test_path =  sys.argv[2]
    #OUT_FOLDER = sys.argv[3]
    
    multiply_data_times = 10
    
    test_path =  '/media/labuser/Data/nanopore/m6A_classifier/data/Epinano/singleA/test_singleA.csv'
    
    # process the training data
    directory = '/media/labuser/Data/nanopore/m6A_classifier/data/Epinano/singleA/'
    mod_rep_2_signal_events = load_local_stored(directory+'mod_rep2_eventalign_numpy_sitesevent_30.npy')
    mod_rep_2_signal_distances = load_local_stored(directory+'mod_rep2_eventalign_numpy_sitesdistances_euclidean_30.npy')
    mod_rep_2_signal_sequences = load_local_stored(directory+'mod_rep2_eventalign_numpy_sitessequences_30.npy')

    no_mod_rep_2_signal_events = load_local_stored(directory+'no_mod_rep2_eventalign_numpy_sitesevent_30.npy')
    no_mod_rep_2_signal_distances = load_local_stored(directory+'no_mod_rep2_eventalign_numpy_sitesdistances_euclidean_30.npy')
    no_mod_rep_2_signal_sequences = load_local_stored(directory+'mod_rep2_eventalign_numpy_sitessequences_30.npy')
    
    
    # parse the testing files
    test = pd.read_csv(test_path, sep='\t', converters={'event': eval,
                                                        'distances': eval})
    X_test = parse_chunk(test)
    y_test = test['label']
    
    #################### Data Agumentation #####################
    
    # I need to know the sequences of nucleotides to calculate the expected sequences
    model_kmer = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/model_kmer.csv',
                             sep=',')
    
    # Create a column with the numbers that we will use as tokens from 0 to 1
    model_kmer['number'] = [round((i+1)/1025, 5) for i in range(len(model_kmer))]

    model_kmer_dict = dict(zip(model_kmer['number'] , model_kmer['model_kmer']))
    tokenizer = dict(zip(model_kmer['model_kmer'], model_kmer['model_mean']))
    
    # gaussian noisse
    augment_mod_signals, augment_mod_distances  = add_permutation(mod_rep_2_signal_events,
                                                                  multiply_data_times,
                                                                  mod_rep_2_signal_sequences,
                                                                  model_kmer_dict,
                                                                  tokenizer,
                                                                  )

    augment_no_mod_signals, augment_no_mod_distances  = add_permutation(no_mod_rep_2_signal_events,
                                                                        multiply_data_times,
                                                                        no_mod_rep_2_signal_sequences,
                                                                        model_kmer_dict,
                                                                        tokenizer,
                                                                        )
    # concatenate input vectors
    no_mod_input = np.concatenate([augment_no_mod_signals,
                                   augment_no_mod_distances],
                                   axis=1)
    
    mod_input = np.concatenate([augment_mod_signals,
                                augment_mod_distances
                                ],
                                axis=1)
    
    no_mod_input = no_mod_input.reshape(len(no_mod_input), 2, 100)
    mod_input = mod_input.reshape(len(mod_input), 2, 100)
    
    # transpose the last two dimensions to have each channel with separate information
    no_mod_input =  np.transpose(no_mod_input, (0, 2, 1))
    mod_input = np.transpose(mod_input, (0, 2, 1))
    
    
    X_train = np.concatenate([no_mod_input, mod_input],
                                axis=0)
    
    y_train =np.concatenate([np.zeros(len(no_mod_input)),
                             np.ones(len(mod_input))])
    
    # Define model Inout
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
    
    checkpoint_filepath =   '/media/labuser/Data/nanopore/m6A_classifier/results/singleA/singleA30TL_aug_permutation_eventWise_150epochs/'
    
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
                        epochs=150,
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
    
        
        
        

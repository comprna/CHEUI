#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 20:29:58 2020

@author: pablo
"""

import os
import pickle

import numpy as np
from tensorflow.keras import Input
from tensorflow.keras import Model
from tensorflow.keras.layers import Dense
from tensorflow.keras.preprocessing.sequence import pad_sequences
from keras import backend as K
from sklearn import preprocessing
from sklearn.utils import shuffle
from sklearn.model_selection import train_test_split
from tensorflow.keras.models import Model
from tensorflow.keras.callbacks import ModelCheckpoint
import keras_metrics

from sklearn.metrics import precision_recall_curve
from matplotlib import pyplot

from itertools import product
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

from DL_models import build_Jasper, build_deepbinner

#import ray
#ray.init()


config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)


def plot_signals_10(KO_signal_filtered_np, WT_signal_filtered_np, kmer, save=None):
    '''
    This function wil plot the median and sd of a bunch of signals
    '''
    # first build the whole signals with medians
    sns.set_style('darkgrid')

    KO_signal = KO_signal_filtered_np.transpose().ravel()
    KO_signal_index = np.arange(1, 101)
    KO_signal_index = np.repeat(KO_signal_index, KO_signal_filtered_np.shape[0])
    
    KO_df = pd.DataFrame({'events' : KO_signal_index, 
                          'signal' : KO_signal})
    
    
    WT_signal = WT_signal_filtered_np.transpose().ravel()
    WT_signal_index = np.arange(1, 101)
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
    X = [10, 30, 50, 70, 90]
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
    plt.show()
    if save:
        plt.title('Signal RNA IVT')
        plt.savefig(save+'.pdf',
                    format='pdf',
                    dpi=1200)
    plt.close('all')
    
    return True


def plot_UMAP(no_mod_list, mod_list, plot=None):
    '''
    '''
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(np.concatenate([mod_list, no_mod_list]))
    
    plt.figure(figsize=(15, 8))
    sns.scatterplot(embedding[:len(mod_list),0],
                    embedding[:len(mod_list),1],
                    color='b')
        
    sns.scatterplot(embedding[len(mod_list):,0],
                    embedding[len(mod_list):,1],
                    color='r')
    if plot:
        plt.savefig(plot+'UMAPembedding'+'.pdf',
                    format='pdf',
                    dpi=1200)
    plt.show()
    plt.close('all')
    return True
    
#@ray.remote
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


if __name__ == '__main__':
    
    # load the local stored data
    directory = '/media/labuser/Data/nanopore/m6A_classifier/data/Epinano/doubleAA/'
    mod_rep_1_signal_dwell_id = load_local_stored(directory+'mod_rep1_eventalign_numpy_sites_AAdwell_50.npy')
    mod_rep_1_signal_events_id = load_local_stored(directory+'mod_rep1_eventalign_numpy_sites_AAevent_50.npy')
    mod_rep_1_signal_distances_id = load_local_stored(directory+'mod_rep1_eventalign_numpy_sites_AAdistances_euclidean_50.npy')
    #mod_rep_1_signal_sequences_id = load_local_stored.remote(directory+'mod_rep1_eventalign_numpy_sites_AAsequences.npy')
    
    no_mod_rep_1_signal_dwell_id = load_local_stored(directory+'no_mod_rep1_eventalign_numpy_sites_AAdwell_50.npy')
    no_mod_rep_1_signal_events_id = load_local_stored(directory+'no_mod_rep1_eventalign_numpy_sites_AAevent_50.npy')
    no_mod_rep_1_signal_distances_id = load_local_stored(directory+'no_mod_rep1_eventalign_numpy_sites_AAdistances_euclidean_50.npy')
    #no_mod_rep_1_signal_sequences_id = load_local_stored.remote(directory+'no_mod_rep1_eventalign_numpy_sites_AAsequences.npy')
    
    mod_rep_2_signal_dwell_id = load_local_stored(directory+'mod_rep2_eventalign_numpy_sites_AAdwell_1000.npy')
    mod_rep_2_signal_events_id = load_local_stored(directory+'mod_rep2_eventalign_numpy_sites_AAevent_1000.npy')
    mod_rep_2_signal_distances_id = load_local_stored(directory+'mod_rep2_eventalign_numpy_sites_AAdistances_euclidean_1000.npy')
    #mod_rep_1_signal_sequences_id = load_local_stored.remote(directory+'mod_rep1_eventalign_numpy_sites_AAsequences.npy')

    no_mod_rep_2_signal_dwell_id = load_local_stored(directory+'no_mod_rep2_eventalign_numpy_sites_AAdwell_1000.npy')
    no_mod_rep_2_signal_events_id = load_local_stored(directory+'no_mod_rep2_eventalign_numpy_sites_AAevent_1000.npy')
    no_mod_rep_2_signal_distances_id = load_local_stored(directory+'no_mod_rep2_eventalign_numpy_sites_AAdistances_euclidean_1000.npy')
    #no_mod_rep_1_signal_sequences_id = load_local_stored.remote(directory+'no_mod_rep1_eventalign_numpy_sites_AAsequences.npy')
    
    '''
    mod_rep_1_signal_dwell = ray.get(mod_rep_1_signal_dwell_id)
    mod_rep_1_signal_events = ray.get(mod_rep_1_signal_events_id)
    mod_rep_1_signal_distances = ray.get(mod_rep_1_signal_distances_id)
    #mod_rep_1_signal_sequences = ray.get(mod_rep_1_signal_sequences_id)
    
    # check the previous veriable 
    no_mod_rep_1_signal_dwell  = ray.get(no_mod_rep_1_signal_dwell_id)
    no_mod_rep_1_signal_events = ray.get(no_mod_rep_1_signal_events_id)
    no_mod_rep_1_signal_distances = ray.get(no_mod_rep_1_signal_distances_id)
    #no_mod_rep_1_signal_sequences = ray.get(no_mod_rep_1_signal_sequences_id)

    mod_rep_2_signal_dwell = ray.get(mod_rep_2_signal_dwell_id)
    mod_rep_2_signal_events = ray.get(mod_rep_2_signal_events_id)
    mod_rep_2_signal_distances = ray.get(mod_rep_2_signal_distances_id)
    #mod_rep_1_signal_sequences = ray.get(mod_rep_1_signal_sequences_id)

    no_mod_rep_2_signal_dwell  = ray.get(no_mod_rep_2_signal_dwell_id)
    no_mod_rep_2_signal_events = ray.get(no_mod_rep_2_signal_events_id)
    no_mod_rep_2_signal_distances = ray.get(no_mod_rep_2_signal_distances_id)
    #no_mod_rep_1_signal_sequences = ray.get(no_mod_rep_1_signal_sequences_id)
    '''

    # concatenate vectors test
    no_mod_rep1 = np.concatenate([no_mod_rep_1_signal_events_id,
                                  no_mod_rep_1_signal_dwell_id,
                                  no_mod_rep_1_signal_distances_id,
                                 ],
                                   axis=1)
    
    mod_input_rep1 = np.concatenate([mod_rep_1_signal_events_id,
                                     mod_rep_1_signal_dwell_id,
                                     mod_rep_1_signal_distances_id,
                                    ],
                                    axis=1)
    
    # reshape the vectors for test
    no_mod_rep1 = no_mod_rep1.reshape(len(no_mod_rep1), 3, 100)
    mod_input_rep1 = mod_input_rep1.reshape(len(mod_input_rep1), 3, 100)
    
    # transpose the last two dimensions to have each channel with separate information
    no_mod_rep1 =  np.transpose(no_mod_rep1, (0, 2, 1))
    mod_input_rep1 = np.transpose(mod_input_rep1, (0, 2, 1))
    
    X_test = np.concatenate([no_mod_rep1, mod_input_rep1],
                                axis=0)
    
    # Create y_test
    y_test = np.concatenate([np.zeros(len(no_mod_rep1)),
                             np.ones(len(mod_input_rep1))],
                            axis=0)
    
    # shuffle test data
    X_test , y_test = shuffle(X_test, y_test, random_state=42)
    
    
    # concatenate training data
    no_mod_rep2 = np.concatenate([no_mod_rep_2_signal_events_id,
                                  no_mod_rep_2_signal_dwell_id,
                                  no_mod_rep_2_signal_distances_id,
                                 ],
                                   axis=1)
    
    mod_input_rep2 = np.concatenate([mod_rep_2_signal_events_id,
                                     mod_rep_2_signal_dwell_id,
                                     mod_rep_2_signal_distances_id,
                                    ],
                                    axis=1)
    
    # reshape the vectors for test
    no_mod_rep2 = no_mod_rep2.reshape(len(no_mod_rep2), 3, 100)
    mod_input_rep2 = mod_input_rep2.reshape(len(mod_input_rep2), 3, 100)
    
    # transpose the last two dimensions to have each channel with separate information
    no_mod_rep2 =  np.transpose(no_mod_rep2, (0, 2, 1))
    mod_input_rep2 = np.transpose(mod_input_rep2, (0, 2, 1))
    
    X_train = np.concatenate([no_mod_rep2, mod_input_rep2],
                                axis=0)
    
    # Create y_train
    y_train = np.concatenate([np.zeros(len(no_mod_rep2)),
                             np.ones(len(mod_input_rep2))],
                            axis=0)
    
    # shuffle test data
    X_train , y_train = shuffle(X_train, y_train, random_state=42)

    # Define model Inout
    inputs = Input(shape=(100, 3))
    
    output = build_Jasper(inputs,Deep=True)
    #output = build_deepbinner(inputs, 1)

    model = Model(inputs=inputs, outputs=output)
    
    model.compile(optimizer='adam',
                  loss='binary_crossentropy',
                  metrics=['accuracy',
                           keras_metrics.recall(),
                           keras_metrics.precision(),
                           ])
    
    #checkpoint = ModelCheckpoint('/media/labuser/Data/nanopore/m6A_classifier/scripts/models/m6A_classifier_DB.h5',
    #                             monitor='val_acc', verbose=1,
    #                             save_best_only=True, mode='max')

    OUT_FOLDER = '/media/labuser/Data/nanopore/m6A_classifier/plots/'
    
    model.summary()

    epochs =1
    histories = []
    
    for e in range(1,epochs+1):
        
        print("Epoch %d" %e)
        
        batch_size = 100000
        
        for i in range(0, len(X_train), batch_size):
            # data augmentation or noise
            if i+batch_size > len(X_train):
                X_train_batch = X_train[i:len(X_train)]
                y_train_batch = y_train[i:len(X_train)]
            
            else:
                X_train_batch = X_train[i:i+batch_size]
                y_train_batch =  y_train[i:i+batch_size]
        
            #Train the model on batches
            history = model.fit(X_train_batch, 
                                y_train_batch, 
                                batch_size=125,
                                validation_data=(X_test, y_test))
            
            histories.append(history)
            model.save(OUT_FOLDER+str(e)+str(i)+"model.h5")
    
    acc = []
    val_acc = []
    loss = []
    val_loss = []
    precision = []
    val_precision = []
    recall = []
    val_recall = []
    
    for i in histories:
        acc.append(i.history['accuracy'][0])
        val_acc.append(i.history['val_accuracy'][0])
        loss.append(i.history['loss'][0])
        val_loss.append(i.history['val_loss'][0])
        precision.append(i.history['precision'][0])
        val_precision.append(i.history['val_precision'][0])
        recall.append(i.history['recall'][0])
        val_recall.append(i.history['val_recall'][0])
    
    
    # summarize history for accuracy 
    f, ax = plt.subplots( figsize=(13,9))
    sns.lineplot(x=np.arange(len(acc)), y=np.array(acc), palette="tab10", linewidth=2.5, label='Accuracy')
    sns.lineplot(x=np.arange(len(acc)), y=np.array(val_acc), palette="tab10", linewidth=2.5, label='Val accuracy')
    plt.ylabel('Accuracies', fontsize=17)
    plt.savefig(OUT_FOLDER+'/accuracy.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    f, ax = plt.subplots( figsize=(13,9))
    sns.lineplot(x=np.arange(len(acc)), y=np.array(loss), palette="tab10", linewidth=2.5, label='loss')
    sns.lineplot(x=np.arange(len(acc)), y=np.array(val_loss), palette="tab10", linewidth=2.5, label='Val loss')
    plt.ylabel('Loss', fontsize=17)
    plt.savefig(OUT_FOLDER+'/Loss.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    f, ax = plt.subplots( figsize=(13,9))
    sns.lineplot(x=np.arange(len(acc)), y=np.array(precision), palette="tab10", linewidth=2.5, label='Precision')
    sns.lineplot(x=np.arange(len(acc)), y=np.array(val_precision), palette="tab10", linewidth=2.5, label='Val precision')
    plt.ylabel('Precision', fontsize=17)
    plt.savefig(OUT_FOLDER+'/Precision.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    f, ax = plt.subplots( figsize=(13,9))
    sns.lineplot(x=np.arange(len(acc)), y=np.array(recall), palette="tab10", linewidth=2.5, label='Recall')
    sns.lineplot(x=np.arange(len(acc)), y=np.array(val_recall), palette="tab10", linewidth=2.5, label='Val Recall')
    plt.ylabel('Recall', fontsize=17)
    plt.savefig(OUT_FOLDER+'/Recall.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    
    from sklearn.metrics import roc_curve
    from sklearn.metrics import roc_auc_score
    
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
    plt.savefig(OUT_FOLDER+'/ROC_AUC.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)

    # predict probabilities
    # predict class values
    lr_precision, lr_recall, _ = precision_recall_curve(y_test, lr_probs)

    f, ax = plt.subplots( figsize=(13,9))
    pyplot.plot(lr_recall, lr_precision, marker='.', label='Logistic')
    plt.ylabel('Precision', fontsize=17)
    plt.xlabel('Recall', fontsize=17)
    plt.savefig(OUT_FOLDER+'/precision_recall_curve.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    
    
    
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

from DL_models import build_Jasper, build_deepbinner, build_TCN_cc
import sys
import os

config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)



def parse_chunk(chunk):
    '''
    '''
    events = chunk['event']
    #dwell = chunk['dwell']
    distances = chunk['distances']
    
    events = np.array([np.array(xi) for xi in events])
    #dwell = np.array([np.array(xi) for xi in dwell])
    distances = np.array([np.array(xi) for xi in distances])

    combined = np.concatenate((events,
                               #dwell,
                              distances), 
                              axis=1)
    
    # reshape the vectors
    combined = combined.reshape(len(combined), 2, 100)
    
    # transpose the last two dimensions to have each channel with separate information
    combined = np.transpose(combined, (0, 2, 1))
    
    return combined


if __name__ == '__main__':
    
    # Define model Inout
    #inputs = Input(shape=(100, 2))
    
    #output = build_Jasper(inputs, Deep=True)
    #output = build_deepbinner(inputs, 1)
    
    inputs = Input(shape=(100, 2))
    output = build_Jasper(inputs, 1)

    model = Model(inputs=inputs, outputs=output)
    
    model.compile(optimizer='adam',
                  loss='binary_crossentropy',
                  metrics=['accuracy',
                           keras_metrics.recall(),
                           keras_metrics.precision(),
                           ])
    
    train_path = sys.argv[1]
    test_path = sys.argv[2]
    OUT_FOLDER = sys.argv[3]

    if not os.path.exists(OUT_FOLDER):
        os.makedirs(OUT_FOLDER)
    
    #model.summary()
    
    # parse the testing files
    test = pd.read_csv(test_path, sep='\t', converters={'event': eval,
                                                        #'dwell': eval,
                                                        'distances': eval})
    X_test = parse_chunk(test)
    y_test = test['label']
    
    
    epochs = 10
    histories = []
    
    for e in range(1,epochs+1):
        
        print("Epoch %d" %e)
        
        for chunk in pd.read_csv(train_path, chunksize=200000, sep='\t', converters={'event': eval,
                                                                                     #'dwell': eval,
                                                                                     'distances': eval}):
            X_train = parse_chunk(chunk)
            
            y_train = chunk['label']
                    
            #Train the model on batches
            history = model.fit(X_train, 
                                y_train, 
                                batch_size=125,
                                validation_data=(X_test,
                                                 y_test),
                                validation_batch_size=125)
            
            histories.append(history)
            model.save(OUT_FOLDER+str(e)+str(chunk.index[0])+"model.h5")
    
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
    
    


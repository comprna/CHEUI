#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 20:29:58 2020

@author: pablo
"""

import pandas as pd
import matplotlib.pyplot as plt
from six.moves import urllib
import zipfile
from scipy import stats
import numpy as np
import seaborn as sns



config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)

tf.random.set_seed(42)

def parse_chunk(chunk):
    '''
    '''
    events = chunk['event']
    distances = chunk['distances']
    sequences = chunk['sequences']
    
    events = np.array([np.array(xi) for xi in events])
    distances = np.array([np.array(xi) for xi in distances])
    sequences = np.array([np.array(xi) for xi in sequences])
    
    return events, sequences


if __name__ == '__main__':
    
    
    train_path = sys.argv[1]
    test_path =  sys.argv[2]
    OUT_FOLDER = sys.argv[3]
    
    train_path = '/media/labuser/Data/nanopore/m6A_classifier/data/Epinano/singleA/train_singleA_2000.csv'
    test_path = '/media/labuser/Data/nanopore/m6A_classifier/data/Epinano/singleA/test_singleA.csv'
    
    if not os.path.exists(OUT_FOLDER):
        os.makedirs(OUT_FOLDER)
    
    
    # parse the testing files
    test = pd.read_csv(test_path, sep='\t', converters={'event': eval,
                                                        'distances': eval,
                                                        'sequences' : eval})
    events, sequences = parse_chunk(test)
   
    unique_sequences = np.unique(sequences, axis=0)
    
    sequences_ = {}
    control_ = {}
    counter = 0
    for i in unique_sequences:
        counter +=1
        for j in enumerate(sequences):
            if (i == j[1]).all() == True:
                if counter in sequences_:
                    sequences_[counter] += [events[j[0]]]
                    control_[counter] += [sequences[j[0]]]
                else:
                    sequences_[counter] = [events[j[0]]]
                    control_[counter] = [sequences[j[0]]]
    
    
    
    # parse the testing files
    train = pd.read_csv(train_path, sep='\t', converters={'event': eval,
                                                        'distances': eval,
                                                        'sequences' : eval})
    events_t, sequences_t = parse_chunk(train)
   
    unique_sequences_t = np.unique(sequences_t, axis=0)
    
    sequences_t = {}
    control_t = {}
    countert = 0
    for i in unique_sequences_t:
        countert +=1
        for j in enumerate(sequences_t):
            if (i == j[1]).all() == True:
                if counter in sequences_:
                    sequences_t[countert] += [events_t[j[0]]]
                    control_t[countert] += [sequences_t[j[0]]]
                else:
                    sequences_t[countert] = [events_t[j[0]]]
                    control_t[countert] = [sequences_t[j[0]]]
    
    list_of_dists = ['alpha','beta','chi2','cosine','expon',\
                     'f','gamma',\
                     'logistic',\
                     'lognorm','maxwell','norm','pareto','powerlaw',\
                    ]
    
    distributions = {}
    for signals_id in sequences_.keys():
        signals_id = 100
        temp_signals = sequences_[signals_id]
        temp_signals = np.array(temp_signals)
        for i in range(0,100,20):
            temp_signals_event = temp_signals[:,i:i+20]
            temp_signals_event = temp_signals_event.reshape(len(temp_signals_event)*temp_signals_event.shape[1])
            for distribution in list_of_dists:
                dist = getattr(stats, distribution)
                parameters = dist.fit(temp_signals_event)
                if distribution in distributions:
                    distributions[distribution] += [stats.kstest(temp_signals_event, distribution, parameters)[1]]
                else:
                    distributions[distribution] = [stats.kstest(temp_signals_event, distribution, parameters)[1]]
    
    
    
    # what is the median of these distributions
    for name in distributions.keys():
        print(name)
        print(np.median(distributions[name]))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

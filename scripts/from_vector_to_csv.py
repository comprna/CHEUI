#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 09:35:47 2020

@author: labuser
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

import ray
ray.init()


config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)


    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
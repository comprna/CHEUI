#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 13:36:21 2021

@author: labuser
"""

import os
import pickle as cPickle

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



# first load m6A and unmodify using dictionaries

m6A = {}
unmodify = {}
m1A = {}

path = '/media/labuser/Data/nanopore/m6A_classifier/data/ELIGOS/'
m1A_id_f = 'IVT_m1A/all_oligos.fasq.gz_IVT_m1A_fast5_nanopolish_IDs.p'
m1A_signal_f = 'IVT_m1A/all_oligos.fasq.gz_IVT_m1A_fast5_nanopolish_signals.p'

unmodify_dir_p = '/media/labuser/Data/nanopore/m6A_classifier/data/ELIGOS/IVT_normalA_Eligos_csv/'
m6A_p = '/media/labuser/Data/nanopore/m6A_classifier/data/ELIGOS/IVT_m6A_Eligos_csv/'

out = '/media/labuser/Data/nanopore/m6A_classifier/results/signal_plots/' 

### load the stored data from m1A
IDs = []
signals = []
with open(path+m1A_signal_f, 'rb') as signal_in:
    with open(path+m1A_id_f, 'rb') as id_in:
        while True:
            try:
                IDs.append(cPickle.load(id_in))
                signals.append(cPickle.load(signal_in))
            except:
                break


### load the stored data from m1A
IDs_new_para = {}
with open(path+'all_oligos.fasq.gz_IVT_m1A_fast5_nanopolish_signals+IDS.p', 'rb') as signal_in:
    while True:
        try:
            IDs_new_para.update(cPickle.load(signal_in))
            # to avoid loading everything predict every 10k singnals
        except:
            break


m1A = {}
# parse signals 
for i in enumerate(IDs_new_para.keys()):
    if '_'.join(i[1].split('_')[:3]) in m1A:
        m1A['_'.join(i[1].split('_')[:3])] += [IDs_new_para[i[1]][:,0]]
    else:
        m1A['_'.join(i[1].split('_')[:3])] = [IDs_new_para[i[1]][:,0]]

# try make df 

#IDs_ = ['_'.join(i.split('_')[:3]) for i in IDs]
#signals_ = [i[:,0] for i in signals]

# load other modification
#m1A_df = pd.DataFrame({'ID': IDs, 
#                       'signals': signals_})

'''
### load the stored data from m1A
IDs = []
signals = []
with open('/media/labuser/Data/nanopore/m6A_classifier/data/ELIGOS/m1A/all_oligos.fasq.gz_IVT_m1A_fast5_nanopolish_signals.p', 'rb') as signal_in:
    with open('/media/labuser/Data/nanopore/m6A_classifier/data/ELIGOS/m1A/all_oligos.fasq.gz_IVT_m1A_fast5_nanopolish_IDs.p', 'rb') as id_in:
        while True:
            try:
                IDs.append(cPickle.load(id_in))
                signals.append(cPickle.load(signal_in))
                # to avoid loading everything predict every 10k singnals
            except:
                break
'''          
m1A = {}
# parse signals 
for i in enumerate(IDs):
    if '_'.join(i[1].split('_')[:3]) in m1A:
        m1A['_'.join(i[1].split('_')[:3])] += [signals[i[0]][:,0]]
    else:
        m1A['_'.join(i[1].split('_')[:3])] = [signals[i[0]][:,0]]


# load m6A signals
for file in os.listdir(m6A_p):
    df = pd.read_csv(m6A_p+file, sep='\t', converters={'event': eval,
                                                      'distances': eval})
    for row in df.iterrows():
        if '_'.join(file.split('_')[:3]) in m6A:
            m6A['_'.join(file.split('_')[:3])] += [row[1][0]]
        else:
             m6A['_'.join(file.split('_')[:3])] = [row[1][0]]


# load unmodify signals
for file in os.listdir(unmodify_dir_p):
    df = pd.read_csv(unmodify_dir_p+file, sep='\t', converters={'event': eval,
                                                      'distances': eval})
    for row in df.iterrows():
        if '_'.join(file.split('_')[:3]) in unmodify:
            unmodify['_'.join(file.split('_')[:3])] += [row[1][0]]
        else:
             unmodify['_'.join(file.split('_')[:3])] = [row[1][0]]



# check the key that have in commonn
m1A_keys = set(list(m1A.keys()))
m6A_keys = set(list(m6A.keys()))
unmodify_keys =set(list(unmodify.keys()))
 
# common IDs
common_ids = m6A_keys & unmodify_keys & m1A_keys

def repeat(arr, count):
    return np.stack([arr for _ in range(count)], axis=0)

counter = 0
for i in common_ids:
    m6A_signals = np.array(m6A[i])
    m1A_signals = np.array(m1A[i])
    unmo_signals = np.array(unmodify[i])
    index = np.arange(1, 101)
    
    m1A_index = np.repeat(index, m1A_signals.shape[0])
    m6A_index = np.repeat(index, m6A_signals.shape[0])
    unmod_index = np.repeat(index, unmo_signals.shape[0]) 
    
    m6A_df = pd.DataFrame({'signal' : m6A_signals.transpose().ravel(), 
                          'events' : m6A_index.tolist()})
                        
    m1A_df = pd.DataFrame({'signal' : m1A_signals.transpose().ravel(), 
                          'events': m1A_index.tolist()})
                             
    unmod_df = pd.DataFrame({'signal' : unmo_signals.transpose().ravel(), 
                            'events' : unmod_index.tolist()})
    
    
    sns.set_context("talk")
    sns.set(font_scale = 3.5)
    sns.set_style("ticks", {"xtick.major.size": 50, "ytick.major.size": 50})
    f, ax = plt.subplots( figsize=(13,9))
    
    ax = sns.lineplot(x="events", y="signal", ci='sd', estimator=np.median, data=m6A_df, color='red')
    ax2 = sns.lineplot(x="events", y="signal", ci='sd', estimator=np.median, data=unmod_df, color='blue')
    ax3 = sns.lineplot(x="events", y="signal", ci='sd', estimator=np.median, data=m1A_df, color='green')
    
    #ax = sns.lineplot(x="events", y="signal", data=m6A_df, color='red')
    #ax2 = sns.lineplot(x="events", y="signal", data=unmod_df, color='blue')
    #ax3 = sns.lineplot(x="events", y="signal",  data=m1A_df, color='green')

    ax2.lines[0].set_linestyle("--")
    ax3.lines[0].set_linestyle("--")
    plt.xticks([])
    plt.yticks([])
    plt.ylabel('Signal intensity', fontsize=30)
    X = [0, 25, 50, 75, 100]
    labels = list(i.split('_')[-1][2:-2])
    # You can specify a rotation for the tick labels in degrees or with keywords.
    plt.xticks(X, labels, fontsize=30)
    
    custom_lines = [Line2D([0], [0], color='blue', lw=4, markersize=40),
                    Line2D([0], [0], color='red', lw=4, markersize=40, linestyle='--'),
                    Line2D([0], [0], color='green', lw=4, markersize=40, linestyle='--')]
    
    plt.legend(custom_lines,
               ['UNMOD','M6A','M1A'],
               fancybox=True,
               framealpha=1, 
               shadow=True, 
               borderpad=1,
               fontsize=20
               )
    
    plt.title('Signal RNA IVT')
    plt.show()
    plt.savefig(out+i+'.pdf',
                format='pdf',
                dpi=1200)
    plt.close('all')
    
    counter +=1
    if counter == 5:
        break
    
    
    
    
###### print individual signals 

counter = 0
for i in common_ids:
    m6A_signals = np.array(m6A[i])
    m1A_signals = np.array(m1A[i])
    unmo_signals = np.array(unmodify[i])
    index = np.arange(1, 101)
    
    sns.set_context("talk")
    sns.set(font_scale = 3.5)
    sns.set_style("ticks", {"xtick.major.size": 50, "ytick.major.size": 50})
    f, ax = plt.subplots( figsize=(13,9))
    '''
    for j in enumerate(m7G_signals):     
        ax = sns.lineplot(x=index, y=j[1], color='red')
        if j[0] == 50:
                break
    '''     
    for j in enumerate(m6A_signals):     
        ax = sns.lineplot(x=index, y=j[1], color='red')
        if j[0] == 50:
            break
    
    for j in enumerate(m1A_signals):     
        ax = sns.lineplot(x=index, y=j[1], color='green')
        if j[0] == 50:
            break
    
    
    for j in enumerate(unmo_signals):     
        ax = sns.lineplot(x=index, y=j[1], color='blue')
        if j[0] == 50:
            break
    
    
    #ax2.lines[0].set_linestyle("--")
    plt.xticks([])
    plt.yticks([])
    plt.ylabel('Signal intensity', fontsize=30)
    X = [0, 25, 50, 75, 100]
    labels = list(i.split('_')[-1][2:-2])
    # You can specify a rotation for the tick labels in degrees or with keywords.
    plt.xticks(X, labels, fontsize=30)
    
    custom_lines = [Line2D([0], [0], color='blue', lw=4, markersize=40),
                    Line2D([0], [0], color='red', lw=4, markersize=40, linestyle='--'),
                    Line2D([0], [0], color='green', lw=4, markersize=40, linestyle='--')]
    
    plt.legend(custom_lines,
               ['UNMOD','M6A','M1A'],
               fancybox=True,
               framealpha=1, 
               shadow=True, 
               borderpad=1,
               fontsize=20
               )
    
    plt.title('Signal RNA IVT')
    plt.show()
    #plt.savefig(out+i+'.sinlge_signals.pdf',
    #            format='pdf',
    #            dpi=1200)
    plt.close('all')
    
    counter +=1
    if counter == 5:
        break
    
'''
from sklearn.manifold import TSNE
X = np.array(m7G_signals)
X_embedded = TSNE(n_components=2).fit_transform(X)
sns.scatterplot(X_embedded[:,0], X_embedded[:,1])
'''
# Compare normalA with different preprocessing

# new
path_normalA = '/media/labuser/Data/nanopore/m6A_classifier/data/ELIGOS/normalA/'
### load the stored data from m1A
IDs = []
signals = []
with open(path_normalA+'IVT_normalA_Eligos_signals.p', 'rb') as signal_in:
    with open(path_normalA+'IVT_normalA_Eligos_IDs.p', 'rb') as id_in:
        while True:
            try:
                IDs.append(cPickle.load(id_in))
                signals.append(cPickle.load(signal_in))
                # to avoid loading everything predict every 10k singnals
            except:
                break


normalA = {}
# parse signals 
for i in enumerate(IDs):
    if '_'.join(i[1].split('_')[:3]) == 'A2_951_TTTTATTTT':
        if '_'.join(i[1].split('_')[:3]) in normalA:
            normalA['_'.join(i[1].split('_')[:3])] += [signals[i[0]][:,0]]
        else:
            normalA['_'.join(i[1].split('_')[:3])] = [signals[i[0]][:,0]]



# check the key that have in commonn
normalA_keys = set(list(normalA.keys()))

 
# common IDs
common_ids = normalA_keys & unmodify_keys

    
###### print individual signals 

counter = 0
for i in common_ids:
    NormalA = np.array(normalA[i])
    unmo_signals = np.array(unmodify[i])
    m1A_signals = np.array(m1A[i])
    index = np.arange(1, 101)
    
    sns.set_context("talk")
    sns.set(font_scale = 3.5)
    sns.set_style("ticks", {"xtick.major.size": 50, "ytick.major.size": 50})
    f, ax = plt.subplots( figsize=(13,9))

    for j in enumerate(NormalA):     
        ax = sns.lineplot(x=index, y=j[1], color='red')
        if j[0] == 50:
                break
      
    for j in enumerate(unmo_signals):     
            ax = sns.lineplot(x=index, y=j[1], color='blue')
            if j[0] == 50:
                break
    
    for j in enumerate(m1A_signals):     
        ax = sns.lineplot(x=index, y=j[1], color='green')
        if j[0] == 50:
            break
    '''
    for j in enumerate(unmo_signals):     
        ax = sns.lineplot(x=index, y=j[1], color='blue')
        if j[0] == 50:
            break
    '''
    
    ax2.lines[0].set_linestyle("--")
    plt.xticks([])
    plt.yticks([])
    plt.ylabel('Signal intensity', fontsize=30)
    X = [0, 25, 50, 75, 100]
    labels = list(i.split('_')[-1][2:-2])
    # You can specify a rotation for the tick labels in degrees or with keywords.
    plt.xticks(X, labels, fontsize=30)
    
    custom_lines = [Line2D([0], [0], color='blue', lw=4, markersize=40),
                    Line2D([0], [0], color='red', lw=4, markersize=40, linestyle='--'),
                    Line2D([0], [0], color='green', lw=4, markersize=40, linestyle='--')]
    
    plt.legend(custom_lines,
               ['UNMOD','M6A','M1A'],
               fancybox=True,
               framealpha=1, 
               shadow=True, 
               borderpad=1,
               fontsize=20
               )
    
    plt.title('Signal RNA IVT')
    plt.show()
    #plt.savefig(out+i+'.sinlge_signals.pdf',
    #            format='pdf',
    #            dpi=1200)
    plt.close('all')
    
    counter +=1
    if counter == 5:
        break










### load the stored data from m1A
signals_counter = 0
signals = []
with open('/media/labuser/Data/nanopore/m6A_classifier/data/ELIGOS/m7G/all_oligos.fasq.gz_IVT_m7G_fast5_nanopolish_signals.p', 'rb') as signal_in:
    while True:
        try:
            signals_counter +=1
            signals.append(cPickle.load(signal_in))
            # to avoid loading everything predict every 10k singnals
        except:
            break

### load the stored data from m1A
IDs = []
id_counter = 0
with open('/media/labuser/Data/nanopore/m6A_classifier/data/ELIGOS/m7G/all_oligos.fasq.gz_IVT_m7G_fast5_nanopolish_IDs.p', 'rb') as id_in:
    while True:
        try:
            id_counter +=1
            IDs.append(cPickle.load(id_in))
            # to avoid loading everything predict every 10k singnals
        except:
            break            
            
m7G = {}
# parse signals 
for i in enumerate(IDs):
    if '_'.join(i[1].split('_')[:3]) in m7G:
        m7G['_'.join(IDs[i[0]].split('_')[:3])] += [signals[i[0]][:,0]]
    else:
        m7G['_'.join(IDs[i[0]].split('_')[:3])] = [signals[i[0]][:,0]]




def compare_sublists(l, lol):
    for sublist in lol:
        temp = [i for i in sublist if i in l]
        if sorted(temp) == sorted(l):
            return True
    return False


LOL_normal = normalA['A2_951_TTTTATTTT']
LOL_unmod = unmodify['A2_951_TTTTATTTT']

LOL_normal = [i.tolist() for i in LOL_normal]


for list_signal in LOL_normal:
    print(compare_sublists(list_signal, LOL_unmod))



#############################################3

### ######################################### AJ 


### load the stored data from m1A
IDs_new_para = {}
counter = 0
with open('/media/labuser/project_NSUN2/AJ/NX_signals+IDS.p', 'rb') as signal_in:
    while True:
        try:
            counter +=1
            if counter == 10000:
                break
            IDs_new_para.update(cPickle.load(signal_in))
            # to avoid loading everything predict every 10k singnals
        except:
            break

NX = {}
# parse signals 
for i in enumerate(IDs_new_para.keys()):
    if '_'.join(i[1].split('_')[:3]) in NX:
        NX['_'.join(i[1].split('_')[:3])] += [IDs_new_para[i[1]][:,0]]
    else:
        NX['_'.join(i[1].split('_')[:3])] = [IDs_new_para[i[1]][:,0]]





### load the stored data from m1A
IDs_new_para = {}
counter = 0
with open('/media/labuser/project_NSUN2/AJ/SXL_signals+IDS.p', 'rb') as signal_in:
    while True:
        try:
            counter +=1
            if counter == 10000:
                break
            IDs_new_para.update(cPickle.load(signal_in))
            # to avoid loading everything predict every 10k singnals
        except:
            break


SXL = {}
# parse signals 
for i in enumerate(IDs_new_para.keys()):
    if '_'.join(i[1].split('_')[:3]) in SXL:
        SXL['_'.join(i[1].split('_')[:3])] += [IDs_new_para[i[1]][:,0]]
    else:
        SXL['_'.join(i[1].split('_')[:3])] = [IDs_new_para[i[1]][:,0]]


# check the key that have in commonn
NX_keys = set(list(NX.keys()))
SXL_keys = set(list(SXL.keys()))
 
# common IDs
common_ids = NX_keys & SXL_keys


counter = 0
for i in common_ids:
    SXL_signals = np.array(SXL[i])
    NX_signals = np.array(NX[i])
    index = np.arange(1, 101)
    
    index_SXL = np.repeat(index, SXL_signals.shape[0])
    index_NX = np.repeat(index, NX_signals.shape[0])
    
    SXL_df = pd.DataFrame({'signal' : SXL_signals.transpose().ravel(), 
                          'events' : index_SXL.tolist()})
                        
    NX_df = pd.DataFrame({'signal' : NX_signals.transpose().ravel(), 
                          'events': index_NX.tolist()})
    
    sns.set_context("talk")
    sns.set(font_scale = 3.5)
    sns.set_style("ticks", {"xtick.major.size": 50, "ytick.major.size": 50})
    f, ax = plt.subplots( figsize=(13,9))
    
    ax = sns.lineplot(x="events", y="signal", ci='sd', estimator=np.median, data=SXL_df, color='red')
    ax2 = sns.lineplot(x="events", y="signal", ci='sd', estimator=np.median, data=NX_df, color='blue')

    ax2.lines[0].set_linestyle("--")
    plt.xticks([])
    plt.yticks([])
    plt.ylabel('Signal intensity', fontsize=30)
    X = [0, 25, 50, 75, 100]
    labels = list(i.split('_')[-1][2:-2])
    # You can specify a rotation for the tick labels in degrees or with keywords.
    plt.xticks(X, labels, fontsize=30)
    
    custom_lines = [Line2D([0], [0], color='blue', lw=4, markersize=40),
                    Line2D([0], [0], color='red', lw=4, markersize=40, linestyle='--')
                    ]
    plt.legend(custom_lines,
               ['SXL','NX'],
               fancybox=True,
               framealpha=1, 
               shadow=True, 
               borderpad=1,
               fontsize=20
               )
    
    plt.title('AJ intense guys please do not hear him')
    plt.show()
    plt.close('all')
    
    counter +=1
    if counter == 5:
        break
    
    
    
    






#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 20:29:58 2020

@author: pablo
"""

import os
import pickle

import numpy as np

import umap
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.lines import Line2D
from sklearn.utils import shuffle
from math import floor

import sys

#@ray.remote
def load_data_smooth(directory, model_kmer_dict, lenght_event, label, df):
    '''
    '''
    files = os.listdir(directory)
 
    for file in files:
        counter = 0
        file_raw = pickle.load(open(directory+file,"rb"))
        for signal_values in file_raw:
            signal_smoothed = []
            dwell_smoothed = []
            sequence_smothed = []
            for event in signal_values:
                if len(event) > 1000: # maximun event lenght
                    dwell_smoothed += [1000/1000]*lenght_event
                else:
                    dwell_smoothed += [len(event)/1000]*lenght_event # record the original dwelling time
                event_smoothed = smooth_event(event, lenght_event) # smooth the event
                signal_smoothed += event_smoothed  # add the event to the signal
            
            expected_smoothed = make_expected(model_kmer_dict, 
                                              file.split('_')[-2],
                                              lenght_event)
            
            sequence = file.split('_')[-2]
            
            sequence_smothed = make_sequences(sequence, lenght_event)
            
            try: # fix the equal kmers bug in script Epinano_site_parse_noA_all.py
                distance_vector = distance_calculator(expected_smoothed,
                                                      signal_smoothed)
                counter +=1
            except:
                break
            
            df_temp = pd.DataFrame({'event': [signal_smoothed],
                                        'distances': [distance_vector],
                                        'dwell': [dwell_smoothed],
                                        'sequences': [sequence_smothed],
                                        'label':label \
                                        })
            df = df.append(df_temp)
            
            if counter == 100: # if there is already that many signal go to the next file
                break
                #return [dic_signal, dic_dwell, dic_distance]
    return df


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
    vector_distance = list(np.round(abs(np.array(signal_expected) - \
                                        np.array(event_smoothed)), 3))
    
    #vector_distance_eu =    list(np.round(abs(np.array(signal_expected) - \
    #                                    np.array(event_smoothed)), 3))
    
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


def make_sequences(nucleotides, lenght):
    '''
    '''
    sequence_to_number = []
    for i in range(len(nucleotides)-4):
        sequence_to_number += [tokenizer[nucleotides[i:i+5]]]*lenght
    return(sequence_to_number)
    # return list(map(int, ''.join([OneHotDNA[i] for i in nucleotides])))


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
    
    directory = sys.argv[1]
    
    model_kmer = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/model_kmer.csv',
                             sep=',')
    
    # create a dictionary with each kmer and its current value
    model_kmer_dict = dict(zip(model_kmer['model_kmer'], model_kmer['model_mean']))
    
    # Create a column with the numbers that we will use as tokens from 0 to 1
    model_kmer['number'] = [round((i+1)/1025, 5) for i in range(len(model_kmer))]
    
    tokenizer = dict(zip(model_kmer['model_kmer'], model_kmer['number']))
    
    df = pd.DataFrame({'event':[],
                       'distances':[],
                       'dwell':[],
                       'sequences':[],
                       'label':[]})
    
    # start the load the raw signal files and preprocess the signals
    # also creating the vectors for distances, sequences and dwelling time
    df = load_data_smooth(directory, model_kmer_dict, 20, 1, df)
    
    #directory = '/media/labuser/Data/nanopore/m6A_classifier/data/Epinano/doubleAA/no_mod_rep2_eventalign_numpy_sites_AA/'
    directory = sys.argv[2]
    
    df = load_data_smooth(directory, model_kmer_dict, 20, 0, df)
    
    path_csv = '/'.join(directory.split('/')[:-2])
    
    df = shuffle(df)
    
    df.round(3).to_csv(path_csv+'/test_s.csv',
               mode='w', 
               sep='\t',
               index=None)




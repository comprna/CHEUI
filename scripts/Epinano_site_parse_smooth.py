
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 10:58:32 2020

@author: labuser
"""

import pandas as pd
import sys
import os
import pickle
from sklearn import preprocessing
from math import floor
import numpy as np

position_single_A = {'cc6m_2244_t7_ecorv': [103,
  122,
  191,
  212,
  610,
  638,
  726,
  802,
  877,
  973,
  1100,
  1139,
  1207,
  1294,
  1372,
  1408,
  1726,
  1880,
  1886,
  2011,
  2064,
  2097,
  2108],
 'cc6m_2459_t7_ecorv': [40,
  67,
  73,
  108,
  155,
  283,
  358,
  369,
  465,
  524,
  578,
  622,
  653,
  824,
  843,
  932,
  1008,
  1045,
  1282,
  1335,
  1343,
  1495,
  1505,
  1558,
  1730,
  1814,
  1821,
  1884,
  1962,
  1983,
  2082,
  2242,
  2248,
  2274,
  2424],
 'cc6m_2595_t7_ecorv': [21,
  99,
  243,
  368,
  389,
  407,
  445,
  508,
  634,
  651,
  782,
  841,
  870,
  881,
  920,
  944,
  997,
  1010,
  1018,
  1027,
  1091,
  1122,
  1152,
  1188,
  1194,
  1238,
  1485,
  1498,
  1600,
  1607,
  1621,
  1734,
  1894,
  1907,
  1916,
  1923,
  2003,
  2082,
  2112,
  2151,
  2195,
  2226,
  2257,
  2581],
 'cc6m_2709_t7_ecorv': [60,
  111,
  136,
  210,
  331,
  501,
  531,
  604,
  614,
  641,
  667,
  676,
  770,
  800,
  806,
  822,
  830,
  889,
  1016,
  1057,
  1137,
  1154,
  1247,
  1307,
  1659,
  1791,
  1880,
  2079,
  2124,
  2175,
  2234,
  2342,
  2442]}


def min_max_scale(signal):
    '''
    Function to scale the data
    '''
    min_max_scaler = preprocessing.MinMaxScaler()
    signal = signal.reshape((len(signal),1))
    signal = min_max_scaler.fit_transform(signal)
    return signal


def top_median(array, lenght):
    '''
    This function top an array until some specific lenght
    '''
    extra_measure = [np.median(array)]*(lenght-len(array))
    array += extra_measure
    return array


def smooth_event(raw_signal, lenght_events):
    '''
    smmoth the signal 
    '''
    raw_signal_events = []
    if len(raw_signal) < lenght_events:
        event = top_median(raw_signal, lenght_events)
        raw_signal_events = raw_signal_events + [event]
    else:
        division = floor(len(raw_signal)/lenght_events)
        new_event = []
        for i in range(0, len(raw_signal), division):
            new_event.append(np.median(raw_signal[i:i+1]))
            if len(new_event) == 10:
                break
            
        if len(new_event) < lenght_events:
            new_event = top_median(new_event, lenght_events)
        raw_signal_events = raw_signal_events + [new_event]

    return raw_signal_events



eventalign = sys.argv[1]


out_folder = eventalign[:-4]+'_numpy_sites'

dic_site = {}

# read pandas in chunks
for chunk in pd.read_csv(eventalign, chunksize=1000, sep='\t'):
    # check how many different contigs there are in the chunk
    contigs = chunk['contig'].unique()
    # go through every different chunk
    for contig in contigs:
        chunk_temp = chunk[chunk['contig'] == contig]
        for pos in position_single_A[contig]:
            pos_5_mer = [pos, pos-1, pos-2, pos-3, pos-4]
            chunk_pos = chunk_temp[chunk_temp['position'].isin(pos_5_mer)]
            for read in chunk_pos['read_name'].unique():
                chunk_read = chunk_pos[chunk_pos['read_name'] == read]
                chunk_read = chunk_read[chunk_read['model_kmer'] != 'NNNNN']
                if len(set(chunk_read['position'])) == 5: # check there is no NNNNNs in any of the kemr model
                    kmer_name = chunk_read[chunk_read['position'] == pos-2]['reference_kmer'].unique()[0]
                    temp_signal = []
                    for kmer in chunk_read['reference_kmer'].unique():
                        
                        kmer_sample_signal = list(chunk_read[chunk_read['reference_kmer'] == kmer]['samples'])
                        
                        kmer_sample_signal = [float(i) for i in ','.join(kmer_sample_signal).split(',')]
                        
                        temp_signal += smooth_event(kmer_sample_signal, 10)
                        
                    if contig+'_'+str(pos)+'_'+kmer_name in dic_site: # be careful with the reference and model kmer
                        dic_site[contig+'_'+str(pos)+'_'+kmer_name] += [temp_signal]
                    else:
                        dic_site[contig+'_'+str(pos)+'_'+kmer_name] = [temp_signal]
    
                


if not os.path.exists(out_folder):
    os.makedirs(out_folder)
'''
for site in dic_site.keys():
    site_array = []
    for indiv_signal in dic_site[site]:
        indiv_signal_entire = ','.join(indiv_signal)
        float_signal = [float(i) for i in indiv_signal_entire.split(',')]
        site_array.append(float_signal)
    with open(out_folder+'/'+site, "wb") as f:
        pickle.dump(site_array, f)
'''
for site in dic_site.keys():
    with open(out_folder+'/'+site+'_events', "wb") as f:
        pickle.dump(dic_site[site], f)






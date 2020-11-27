#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 08:35:11 2020

This script 


@author: labuser
"""

import pandas as pd
import sys
import os
import h5py


def process_sample(sample, transcript, pos, kmer, 
                   median_values, probabilities, samples_WT_KO,
                   index_reads, read_dic):
    '''
    This function goes through one of the transcripts scanning all the reads and
    getting information from a specific position, storing probabilitites and read ID
    '''

    for i in enumerate(sample[transcript].keys()):
        for pos_inf in sample[transcript+'/'+i[1]+'/events']:
            # check that position, kmer and median value are same as the model values
            if pos_inf[2] == pos and pos_inf[3].decode('utf-8') == kmer \
                and pos_inf[1].decode('utf-8') == transcript and \
                round(pos_inf[-1],2) == round(median_values[index_reads],2):
                
                # Add the position in a dictionary with transcripts and positions as keys and reads 
                # and probabilities as values
                if pos_inf[1].decode('utf-8')+'_'+str(pos_inf[2])+'_'+pos_inf[3].decode('utf-8') in read_dic:
                    lista.append(pos_inf)
                    read_dic[pos_inf[1].decode('utf-8')+'_'+str(pos_inf[2])+'_'+pos_inf[3].decode('utf-8')] += \
                            [[pos_inf[0].decode('utf-8'), 
                             probabilities[index_reads], samples_WT_KO[index_reads].decode('utf-8')]]
                else:
                    read_dic[pos_inf[1].decode('utf-8')+'_'+str(pos_inf[2])+'_'+pos_inf[3].decode('utf-8')] = \
                            [[pos_inf[0].decode('utf-8'),
                             probabilities[index_reads], samples_WT_KO[index_reads].decode('utf-8')]]
                index_reads +=1
        
    return read_dic, index_reads



def readID_prob(transcript, pos, kmer, median_values,
                probabilities, samples_WT_KO, read_dic, 
                *argv):
    '''
    This function will collect the information to assign to each read id for each site
    a probability of baing modify
    '''
    index_reads = 0
    
    for sample in argv:
        read_dic, index_reads = process_sample(sample, transcript, pos, kmer, 
                                               median_values, probabilities, samples_WT_KO,
                                               index_reads, read_dic
                                               )
    
    return read_dic, index_reads
    
if __name__ == '__main__':

    sys.argv[1]
    
    # import the orthogonal method
    Schwartz = pd.read_csv('/media/labuser/Data/nanopore/yeast_Schwartz.csv', sep='\t')
    
    Schwartz.dropna(axis=0, subset=['confGroup'], inplace=True)
    
    Schwartz = Schwartz[Schwartz['confGroupOri'] >=1 ]
    
    # import diffmos filtered files
    diffmod = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/xpore_results/all_intersection_r.csv')
    
    diffmod = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/xpore_results/all_intersection_f.csv',
                          sep='\t')
    
    
    # import the models
    model_path = '/home/labuser/lib/xpore/demo/out_no_genome/models'
    models = os.listdir(model_path)
    
    
    # Get the eventalign files from each sample
    WT_rep1 = '/home/labuser/lib/xpore/demo/data/HEK293T-WT-rep1/dataprep_without_genome/eventalign.hdf5'
    KO_rep1 = '/home/labuser/lib/xpore/demo/data/HEK293T-METTL3-KO-rep1/dataprep_without_genome/eventalign.hdf5'
    
    WT_rep1_r = h5py.File(WT_rep1, 'r')
    KO_rep1_r = h5py.File(KO_rep1, 'r')
    
    read_dic = {}
    lista = []
    # go through the significant sites in diffmode.table
    for row in diffmod.itertuples():
        transcript = row[1]
        pos = row[2]
        kmer = row[3]
        
        if kmer[2] != 'A':
            continue
        
        for model in models:
            if model[:-5] == transcript:
                model_temp = model
    
        # open the model corresponding to the transcript in models
        read_model = h5py.File(model_path+'/'+model_temp, 'r')
        
        try:
            if read_model[transcript+'/'+str(pos)].attrs['kmer'].decode('utf-8') == kmer:
                pass
            else:
                print('Kmer from xpore-diffmod table does not match the information in the models')
        except:
            continue
        
        
        # take mean values from all the samples to match with the reads
        median_values = list(read_model[transcript+'/'+str(pos)+'/nodes/y/data'])
    
        # take the information from the samples, WT or KO
        samples_WT_KO = list(read_model[transcript+'/'+str(pos)+'/nodes/y_run_names/data'])
        
        # take the information from the 
        probabilities = list(read_model[transcript+'/'+str(pos)+'/nodes/z/prob'])
        
        read_dic, index_reads = readID_prob(transcript, pos, kmer, median_values,
                                            probabilities, samples_WT_KO, read_dic,
                                            KO_rep1_r, WT_rep1_r)
        
        out_file = '/media/labuser/Data/nanopore/m6A_classifier/parse_models.txt'
        
        with open(out_file, 'w') as f_out:
            print('transcript'+'\t'+'position'+'\t'+'kmer'+'\t'+'read_id'+'\t'+'p'+'\t'+'q'+'\t'+'sample', file=f_out)
            for site in read_dic.keys():
                for read in read_dic[site]:
                    print(site.split('_')[0]+\
                          '\t'+site.split('_')[1]+\
                          '\t'+site.split('_')[2]+\
                          '\t'+read[0]+\
                          '\t'+str(read[1][0])+\
                          '\t'+str(read[1][1])+\
                          '\t'+read[2],
                          file=f_out)
    













#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 08:35:11 2020

This script 


@author: labuser
"""

import pandas as pd
import numpy as np
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
                #print('asdasd')
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


def find_model(model_dir, transcript, pos):
    '''
    '''
    
    for chromosome in os.listdir(model_dir):
        read_model = h5py.File(model_dir+chromosome, 'r')
        try:
            read_model[transcript+'/'+str(pos)].attrs['kmer'].decode('utf-8')
            #print('Model'+ model_dir +' YES contain information about '+chromosome+'/'+str(pos))
            return read_model

        except:
            #print('Model'+ model_dir +' does not contain information about '+chromosome+'/'+str(pos))
            continue
    
    return False


if __name__ == '__main__':

    sys.argv[1]
    
    # import the intersec sites from orthogonal methods and xpore
    #diffmod = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/xpore_results/all_intersection_r.csv',
    #                      sep='\t')
    
    diffmod = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/xpore_results/all_intersection_f.csv',
                          sep='\t')
    
    # import the models
    #model_path = '/home/labuser/lib/xpore/demo/out_no_genome/models'
    #models = os.listdir(model_path)
    
    model_paths = ['/media/labuser/Data/nanopore/m6A_classifier/yeast/models/diffmod_out_m6A_pos_f_1/models/',
                   '/media/labuser/Data/nanopore/m6A_classifier/yeast/models/diffmod_out_m6A_pos_f_2/models/',
                   '/media/labuser/Data/nanopore/m6A_classifier/yeast/models/diffmod_out_m6A_pos_f_3/models/']
    
    # Get the eventalign files from each sample from the fordward strand
    path_eventalign = '/media/labuser/Data/nanopore/m6A_classifier/yeast/eventaligns/'
    
    path_1 = 'dataprep_filtered_m6A_pos_f_1/eventalign.hdf5'
    path_2 = 'dataprep_filtered_m6A_pos_f_2/eventalign.hdf5'
    path_3 = 'dataprep_filtered_m6A_pos_f_3/eventalign.hdf5'
    
    KO1_RNAAA024588_f1 = h5py.File(path_eventalign+'KO1_RNAAA024588/'+path_1, 'r')
    KO1_RNAAA024588_f2 = h5py.File(path_eventalign+'KO1_RNAAA024588/'+path_2, 'r')
    KO1_RNAAA024588_f3 = h5py.File(path_eventalign+'KO1_RNAAA024588/'+path_3, 'r')
    
    KO2_RNAAB058843_f1 = h5py.File(path_eventalign+'KO2_RNAAB058843/'+path_1, 'r')
    KO2_RNAAB058843_f2 = h5py.File(path_eventalign+'KO2_RNAAB058843/'+path_2, 'r')
    KO2_RNAAB058843_f3 = h5py.File(path_eventalign+'KO2_RNAAB058843/'+path_3, 'r')
    
    KO3_RNAAB059882_f1 = h5py.File(path_eventalign+'KO3_RNAAB059882/'+path_1, 'r')
    KO3_RNAAB059882_f2 = h5py.File(path_eventalign+'KO3_RNAAB059882/'+path_2, 'r')
    KO3_RNAAB059882_f3 = h5py.File(path_eventalign+'KO3_RNAAB059882/'+path_3, 'r')

    WT1_RNAAA023484_f1 = h5py.File(path_eventalign+'WT1_RNAAA023484/'+path_1, 'r')
    WT1_RNAAA023484_f2 = h5py.File(path_eventalign+'WT1_RNAAA023484/'+path_2, 'r')
    WT1_RNAAA023484_f3 = h5py.File(path_eventalign+'WT1_RNAAA023484/'+path_3, 'r')
    
    WT2_RNAAB056712_f1 = h5py.File(path_eventalign+'WT2_RNAAB056712/'+path_1, 'r')
    WT2_RNAAB056712_f2 = h5py.File(path_eventalign+'WT2_RNAAB056712/'+path_2, 'r')
    WT2_RNAAB056712_f3 = h5py.File(path_eventalign+'WT2_RNAAB056712/'+path_3, 'r')
    
    WT3_RNAAB057791_f1 = h5py.File(path_eventalign+'WT3_RNAAB057791/'+path_1, 'r')
    WT3_RNAAB057791_f2 = h5py.File(path_eventalign+'WT3_RNAAB057791/'+path_2, 'r')
    WT3_RNAAB057791_f3 = h5py.File(path_eventalign+'WT3_RNAAB057791/'+path_3, 'r')
    
    #WT_rep1 = '/home/labuser/lib/xpore/demo/data/HEK293T-WT-rep1/dataprep_without_genome/eventalign.hdf5'
    #KO_rep1 = '/home/labuser/lib/xpore/demo/data/HEK293T-METTL3-KO-rep1/dataprep_without_genome/eventalign.hdf5'
    #WT_rep1_r = h5py.File(WT_rep1, 'r')
    #KO_rep1_r = h5py.File(KO_rep1, 'r')
    
    '''
    I have to include these 
    chr12 839845 118.0737185145632 [126.64876572106172, 117.89343541561658]
    chr15  946537  127.166509  123.010404 (models 122.976, 127.176)
    chr7 316299 94.71927  101.599916  (102.167, 94.711)
    
    # only to recuse the mean problem ones 
    diffmod = diffmod[diffmod['position'].isin([839845,946537,316299])]
    
    # import the previous version of the parsed models
    df = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/results/select_intersect_reads/parse_models_f_p.txt',
                     sep='\t')
    
    # get rid of the position where the probabilities are wrong
    df = df[df['position'] != 839845 ]
    df = df[df['position'] != 946537 ]
    df = df[df['position'] != 316299 ]
    
    df.to_csv('/media/labuser/Data/nanopore/m6A_classifier/results/select_intersect_reads/parse_models_f_p.txt',
                     sep='\t')
    '''
    
    read_dic = {}
    lista = []
    
    # go through the significant sites in diffmode.table
    for row in diffmod.itertuples():
        transcript = row[1]
        pos = row[2]
        kmer = row[3]
        mu_mod = row[20]
        
     
        # Open the model that contain the position we are parsing
        for model_dir in model_paths:
            read_model = find_model(model_dir, transcript, pos)
            if read_model == False:
                continue
            else:
                break
    
        try:
            # take mean values from all the samples to match with the reads
            median_values = list(read_model[transcript+'/'+str(pos)+'/nodes/y/data'])

            # take the information from the samples, WT or KO
            samples_WT_KO = list(read_model[transcript+'/'+str(pos)+'/nodes/y_run_names/data'])
            
            # check the total means of KO and WT to see which of the columns contain 
            # the probability of being modify
            global_means = list(read_model[transcript+'/'+str(pos)+'/nodes/y/mean'])
            
            # take the information from the 
            probabilities_array = np.array(read_model[transcript+'/'+str(pos)+'/nodes/z/prob'])
            
            # if the mu_mod is the first integer then take the first probability column.
            # if the mu_mod is the second number then take the second column of probablity numbers.
            if int(global_means[0]) == int(mu_mod):
                probabilities = list(probabilities_array[:,0])
            elif int(global_means[1]) == int(mu_mod):
                probabilities = list(probabilities_array[:,1])
            else:
                if int(round(global_means[0])) == round(mu_mod):
                    probabilities = list(probabilities_array[:,0])
                elif round(global_means[1]) == round(mu_mod):
                    probabilities = list(probabilities_array[:,1])
                else:
                    print('error, global medians do not match with models medians')
                    print(transcript, pos, mu_mod, global_means)
                    probabilities = []
            
            
            read_dic, index_reads = readID_prob(transcript, pos, kmer, median_values,
                                                probabilities, samples_WT_KO, read_dic,
                                                KO1_RNAAA024588_f1, KO1_RNAAA024588_f2, KO1_RNAAA024588_f3,
                                                KO2_RNAAB058843_f1, KO2_RNAAB058843_f2, KO2_RNAAB058843_f3,
                                                KO3_RNAAB059882_f1, KO3_RNAAB059882_f2, KO3_RNAAB059882_f3,
                                                WT1_RNAAA023484_f1, WT1_RNAAA023484_f2, WT1_RNAAA023484_f3,
                                                WT2_RNAAB056712_f1, WT2_RNAAB056712_f2, WT2_RNAAB056712_f3,
                                                WT3_RNAAB057791_f1, WT3_RNAAB057791_f2, WT3_RNAAB057791_f3
                                                )
            
            out_file = '/media/labuser/Data/nanopore/m6A_classifier/results/select_intersect_reads/parse_models_f_p.txt'
            
            with open(out_file, 'a') as f_out:
                print('transcript'+'\t'+'position'+'\t'+'kmer'+'\t'+'read_id'+'\t'+'p'+'\t'+'sample', file=f_out)
                for site in read_dic.keys():
                    for read in read_dic[site]:
                        print(site.split('_')[0]+\
                              '\t'+site.split('_')[1]+\
                              '\t'+site.split('_')[2]+\
                              '\t'+read[0]+\
                              '\t'+str(read[1])+\
                              '\t'+read[2],
                              file=f_out)
        
        except:
            print('The information from transcript '+transcript+' position '+str(pos)+' could not be found')
            













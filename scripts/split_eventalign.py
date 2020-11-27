#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 30 15:12:42 2020

@author: pablo
"""

import sys
import pandas as pd

eventalign = sys.argv[1]
summary = sys.argv[2]

index_chr = {'chr10':766161,
             'chr11':678764,
             'chr12':1071616,
             'chr13':919118,
             'chr14':799305,
             'chr15':1104414,
             'chr16':960530,
             'chr1':258084,
             'chr2':804360,
             'chr3':336578,
             'chr4':1530282,
             'chr5':591864,
             'chr6':273963,
             'chr7':1117752,
             'chr8':572889,
             'chr9':457059,
            }


eventalign = '/media/labuser/Data/nanopore/S_cerevisiae_2_3/KO2/nanopolish/KO2_nanopolish.txt'


read_indexes_f_1 = set()
read_indexes_f_2 = set()
read_indexes_f_3 = set()

read_indexes_r_1 = set()
read_indexes_r_2 = set()
read_indexes_r_3 = set()


def split_chunk(chunk_read, index_chr):
    '''
    '''
    chunk_read_1 = pd.DataFrame()
    chunk_read_2 = pd.DataFrame()
    chunk_read_3 = pd.DataFrame()
    
    for chromosome in index_chr.keys():
        chunk_read_chr = chunk_read[chunk_read['contig'] == chromosome]
        chunk_read_1 = pd.concat([chunk_read_1, \
                                  chunk_read_chr[chunk_read_chr['position'] < index_chr[chromosome]/3]],\
                                  axis=0)
        
        chunk_read_2 =  pd.concat([chunk_read_2, \
                                   chunk_read_chr[(chunk_read_chr['position'] > index_chr[chromosome]/3) &\
                                                  (chunk_read_chr['position'] < (index_chr[chromosome]/3)*2)]],\
                                  axis=0)
        
        chunk_read_3 = pd.concat([chunk_read_3, \
                                   chunk_read_chr[(chunk_read_chr['position'] > (index_chr[chromosome]/3)*2)]],\
                                  axis=0)

    return chunk_read_1, chunk_read_2, chunk_read_3
    

counter_f = 0
counter_r = 0
# read pandas in chunks
for chunk in pd.read_csv(eventalign, chunksize=200000, sep='\t'):
    chunk = chunk[chunk['reference_kmer'].str.contains('^[A-Z]{2}A[A-Z]{2}$', regex=True)]
    for read_id in chunk['read_name'].unique():
        chunk_read = chunk[chunk['read_name'] == read_id]
        for i in enumerate(chunk_read['model_kmer']):
            if i[1] != 'NNNNN':
                if list(chunk_read['model_kmer'])[i[0]] == list(chunk_read['reference_kmer'])[i[0]]:
                    
                    chunk_read_1, chunk_read_2, chunk_read_3 = split_chunk(chunk_read, index_chr)
                    
                    counter_f +=1
                    if counter_f == 1:
                        chunk_read_1.to_csv(eventalign[:-4]+'filtered_f_1.txt', sep='\t', \
                                          index=False, float_format='%g')
                        chunk_read_2.to_csv(eventalign[:-4]+'filtered_f_2.txt', sep='\t', \
                                          index=False, float_format='%g')
                        chunk_read_3.to_csv(eventalign[:-4]+'filtered_f_3.txt', sep='\t', \
                                          index=False, float_format='%g')
                    else:
                        chunk_read_1.to_csv(eventalign[:-4]+'filtered_f_1.txt', sep='\t', \
                                          index=False, header=False, float_format='%g', mode='a')
                        chunk_read_2.to_csv(eventalign[:-4]+'filtered_f_2.txt', sep='\t', \
                                          index=False, header=False, float_format='%g', mode='a')
                        chunk_read_3.to_csv(eventalign[:-4]+'filtered_f_3.txt', sep='\t', \
                                          index=False, header=False, float_format='%g', mode='a')
                    
                    read_indexes_f_1.update(set(list(chunk_read_1['read_name'])))
                    read_indexes_f_2.update(set(list(chunk_read_2['read_name'])))
                    read_indexes_f_3.update(set(list(chunk_read_3['read_name'])))

                else:
                    chunk_read_1, chunk_read_2, chunk_read_3 = split_chunk(chunk_read, index_chr)
                    counter_r +=1
                    if counter_r == 1:
                        chunk_read_1.to_csv(eventalign[:-4]+'filtered_r_1.txt', sep='\t', \
                                          index=False, float_format='%g')
                        
                        chunk_read_2.to_csv(eventalign[:-4]+'filtered_r_2.txt', sep='\t', \
                                          index=False, float_format='%g')
                        
                        chunk_read_3.to_csv(eventalign[:-4]+'filtered_r_3.txt', sep='\t', \
                                          index=False, float_format='%g')
                        
                    else:
                        chunk_read_1.to_csv(eventalign[:-4]+'filtered_r_1.txt', sep='\t', \
                                          index=False, header=False, float_format='%g', mode='a')
                        
                        chunk_read_2.to_csv(eventalign[:-4]+'filtered_r_2.txt', sep='\t', \
                                          index=False, header=False, float_format='%g', mode='a')
                        
                        chunk_read_3.to_csv(eventalign[:-4]+'filtered_r_3.txt', sep='\t', \
                                          index=False, header=False, float_format='%g', mode='a')
                    
                    read_indexes_r_1.update(set(list(chunk_read_1['read_name'])))
                    read_indexes_r_2.update(set(list(chunk_read_2['read_name'])))
                    read_indexes_r_3.update(set(list(chunk_read_3['read_name'])))
                break
            else:
                continue


df_summary_f = pd.read_csv(summary, sep='\t')
df_summary_f_1 = df_summary_f[df_summary_f['read_index'].isin(list(read_indexes_f_1))]
df_summary_f_1.to_csv(summary[:-4]+'indexes_f_1.txt', sep='\t', float_format='%.3f', index=None)

df_summary_f = pd.read_csv(summary, sep='\t')
df_summary_f_2 = df_summary_f[df_summary_f['read_index'].isin(list(read_indexes_f_2))]
df_summary_f_2.to_csv(summary[:-4]+'indexes_f_2.txt', sep='\t', float_format='%.3f', index=None)

df_summary_f = pd.read_csv(summary, sep='\t')
df_summary_f_3 = df_summary_f[df_summary_f['read_index'].isin(list(read_indexes_f_3))]
df_summary_f_3.to_csv(summary[:-4]+'indexes_f_3.txt', sep='\t', float_format='%.3f', index=None)

df_summary_r = pd.read_csv(summary, sep='\t')
df_summary_r_1 = df_summary_r[df_summary_r['read_index'].isin(list(read_indexes_r_1))]
df_summary_r_1.to_csv(summary[:-4]+'indexes_r_1.txt', sep='\t', float_format='%.3f', index=None)

df_summary_r = pd.read_csv(summary, sep='\t')
df_summary_r_2 = df_summary_r[df_summary_r['read_index'].isin(list(read_indexes_r_2))]
df_summary_r_2.to_csv(summary[:-4]+'indexes_r_2.txt', sep='\t', float_format='%.3f', index=None)

df_summary_r = pd.read_csv(summary, sep='\t')
df_summary_r_3 = df_summary_r[df_summary_r['read_index'].isin(list(read_indexes_r_3))]
df_summary_r_3.to_csv(summary[:-4]+'indexes_r_3.txt', sep='\t', float_format='%.3f', index=None)


























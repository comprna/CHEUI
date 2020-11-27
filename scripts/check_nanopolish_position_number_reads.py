#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 12:05:34 2020

@author: labuser
"""
import pandas as pd 
import sys

### check the number of reads per site in a 

# get sites from excel file
#m6A_seq = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/orthogonal_sites/mmc1.csv', sep='\t')
#m6A_seq = m6A_seq[['Peak chr', 'Peak genomic coordinate']]

m6A_seq = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/orthogonal_sites/m6A_Schwartz_ConfFiltered.csv', sep='\t')
m6A_seq =m6A_seq[m6A_seq['strand'] == '+']
#m6A_seq = m6A_seq[m6A_seq['confGroup'] >= 0]
m6A_seq = m6A_seq[['chr', 'start']]
m6A_seq['start'] = m6A_seq['start'].astype(int)

m6A_seq.rename(columns={'chr':'Peak chr', 'start':'Peak genomic coordinate'}, inplace=True)

# the position are shifted between nanopolish and the csv file, have to substract 3
m6A_seq['Peak genomic coordinate'] = m6A_seq['Peak genomic coordinate']-3

files = [sys.argv[1],
         sys.argv[2],
         sys.argv[3]]

intersection_pd = pd.DataFrame()
for file in files:
    for chunk in pd.read_csv(file, chunksize=200000, sep='\t'):
        for chromosome in chunk['contig'].unique():
            chunk_chr = chunk[chunk['contig'] == chromosome]
            m6A_seq_pos = m6A_seq[m6A_seq['Peak chr'] == chromosome]['Peak genomic coordinate']
            chunk_chr_pos = chunk_chr[chunk_chr['position'].isin(m6A_seq_pos)]
            intersection_pd = pd.concat([intersection_pd,chunk_chr_pos], axis=0)



intersection_pd.to_csv(sys.argv[1][:-5]+'intersection_eventaling_'+sys.argv[4]+'.txt',
                       sep='\t', float_format='%g')

####### count the number of reads in eventaligns
'''
folder_eventaligns = '/media/labuser/Data/nanopore/m6A_classifier/data/yeast'

WT_intersection = {}
KO_intersection = {}

for file in os.listdir(folder_eventaligns):
    if 'WT' in file:
        WT_temp = pd.read_csv(folder_eventaligns+'/'+file, sep='\t')
        for chromosome in m6A_seq['Peak chr'].unique():
            WT_temp_chr = WT_temp[WT_temp['contig'] == chromosome]
            for position in WT_temp_chr['position'].unique():
                WT_temp_chr_pos = WT_temp_chr[WT_temp_chr['position'] == position]
                if str(WT_temp_chr_pos['contig'].iloc[0])+'_'+str(WT_temp_chr_pos['position'].iloc[0]) in WT_intersection:
                    WT_intersection[str(WT_temp_chr_pos['contig'].iloc[0])+'_'+str(WT_temp_chr_pos['position'].iloc[0])] += \
                        len(WT_temp_chr_pos)
                else:
                    WT_intersection[str(WT_temp_chr_pos['contig'].iloc[0])+'_'+str(WT_temp_chr_pos['position'].iloc[0])] = \
                        len(WT_temp_chr_pos)
    else:
        break
        KO_temp = pd.read_csv(folder_eventaligns+'/'+file, sep='\t')
        for chromosome in m6A_seq['Peak chr'].unique():
            KO_temp_chr = KO_temp[KO_temp['contig'] == chromosome]
            for position in KO_temp_chr['position'].unique():
                KO_temp_ch_pos = KO_temp_chr[KO_temp_chr['position'] == position]
                if str(KO_temp_ch_pos['contig'].iloc[0])+'_'+str(KO_temp_ch_pos['position'].iloc[0]) in KO_intersection:
                    KO_intersection[str(KO_temp_ch_pos['contig'].iloc[0])+'_'+str(KO_temp_ch_pos['position'].iloc[0])] += \
                        len(KO_temp_ch_pos)
                else:
                    KO_intersection[str(KO_temp_ch_pos['contig'].iloc[0])+'_'+str(KO_temp_ch_pos['position'].iloc[0])] = \
                        len(KO_temp_ch_pos)
    
'''
        
        
        
        
        
        
        
        
        
        
        
        
        
        


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 08:52:32 2020

@author: labuser

This will check the number of reads in the positions in samfiles, eventalign 
diffmod and diffmod significant
"""
import pandas as pd
import numpy as np
import os 
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


### check the diffmod table intersection with the m6A fodrward sites

main_f = '/media/labuser/Data/nanopore/m6A_classifier/data/xpore_results/all'

diffmod_f = pd.DataFrame()
diffmod_r = pd.DataFrame()

for folder in os.listdir(main_f):
    if '_f_' in folder:
        diffmod_new = pd.read_csv(main_f+'/'+folder+'/'+'diffmod.table', sep=',')
        diffmod_f = pd.concat([diffmod_f, diffmod_new], axis=0)
    else:
        diffmod_new = pd.read_csv(main_f+'/'+folder+'/'+'diffmod.table', sep=',')
        diffmod_r = pd.concat([diffmod_r, diffmod_new], axis=0)

# significant sites 
diffmod_f = diffmod_f[diffmod_f['pval_KO_vs_WT'] < 0.05]
diffmod_r = diffmod_r[diffmod_r['pval_KO_vs_WT'] < 0.05]

print(diffmod_f.shape)
print(diffmod_r.shape)

# get rid of modifications that are in discordant with the majority of mod_assignment
diffmod_f_grouped = pd.DataFrame()

for kmer in diffmod_f['kmer'].unique():
    temp = diffmod_f[diffmod_f['kmer'] == kmer]
    if len(temp) == 1:
        diffmod_f_grouped = pd.concat([diffmod_f_grouped, temp])
    elif len(temp[temp['mod_assignment'] == 'lower']) > len(temp[temp['mod_assignment'] == 'higher']):
        diffmod_f_grouped = pd.concat([diffmod_f_grouped, temp[temp['mod_assignment'] == 'lower']])
    elif len(temp[temp['mod_assignment'] == 'lower']) < len(temp[temp['mod_assignment'] == 'higher']):
        diffmod_f_grouped = pd.concat([diffmod_f_grouped, temp[temp['mod_assignment'] == 'higher']])

diffmod_r_grouped = pd.DataFrame()
for kmer in diffmod_r['kmer'].unique():
    temp = diffmod_r[diffmod_r['kmer'] == kmer]
    if len(temp) == 1:
        diffmod_r_grouped = pd.concat([diffmod_r_grouped, temp])
    elif len(temp[temp['mod_assignment'] == 'lower']) > len(temp[temp['mod_assignment'] == 'higher']):
        diffmod_r_grouped = pd.concat([diffmod_r_grouped, temp[temp['mod_assignment'] == 'lower']])
    elif len(temp[temp['mod_assignment'] == 'lower']) < len(temp[temp['mod_assignment'] == 'higher']):
        diffmod_r_grouped = pd.concat([diffmod_r_grouped, temp[temp['mod_assignment'] == 'higher']])

print(diffmod_f_grouped.shape)
print(diffmod_r_grouped.shape)


## positive strand
m6A_seq = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/orthogonal_sites/m6A_Schwartz.csv', sep='\t')

m6A_seq = m6A_seq[m6A_seq['confGroup'] >= 0]
m6A_seq = m6A_seq[m6A_seq['strand'] == '+']

#diffmod_f_ACA = diffmod_r[diffmod_r['kmer'].str.contains('^TGT[A-Z]{2}$', regex=True)]
diffmod_f_ACA = diffmod_f_grouped[diffmod_f_grouped['kmer'].str.contains('^[A-Z]{2}ACA$', regex=True)]

m6A_seq['position'] = m6A_seq['start'].astype(int) - 1

intersection_diffmod_mazter = []
for chromosome in m6A_seq['chr'].unique():
    m6A_seq_pos =   list(m6A_seq[m6A_seq['chr'] == chromosome]['position'])
    diffmod_f_chr = list(diffmod_f_ACA[diffmod_f_ACA['id'] == chromosome]['position'])
    overlap = set(m6A_seq_pos) & set([int(i) for i in diffmod_f_chr])
    print(m6A_seq[m6A_seq['position'].isin(list(overlap))]['stoichiometryModelFit'])
    intersection_diffmod_mazter += list(overlap)

print(len(intersection_diffmod_mazter))

# high confident mazter-seq sites 

m6A_seq = m6A_seq[m6A_seq['confGroup'] >=1]

m6A_seq_intersect_f = pd.DataFrame()
diffmod_f_ACA_intersect = pd.DataFrame()

intersection_diffmod_mazter = []
for chromosome in m6A_seq['chr'].unique():
    m6A_seq_pos = list(m6A_seq[m6A_seq['chr'] == chromosome]['position'])
    diffmod_f_chr = list(diffmod_f_ACA[diffmod_f_ACA['id'] == chromosome]['position'])
    overlap = set(m6A_seq_pos) & set([int(i) for i in diffmod_f_chr])
    intersection_diffmod_mazter += list(overlap)
    print(m6A_seq[m6A_seq['position'].isin(list(overlap))]['stoichiometryModelFit'])
    m6A_seq_intersect_f = pd.concat([m6A_seq_intersect_f, m6A_seq[(m6A_seq['chr'] == chromosome) & (m6A_seq['position'].isin(list(overlap)))]])
    diffmod_f_ACA_intersect =  pd.concat([diffmod_f_ACA_intersect, diffmod_f_ACA[(diffmod_f_ACA['id'] == chromosome) & (diffmod_f_ACA['position'].isin(list(overlap)))]])
    
print(len(intersection_diffmod_mazter))


# save the intersection sites
m6A_seq_intersect_f.to_csv('/media/labuser/Data/nanopore/m6A_classifier/data/xpore_results/m6A_seq_intersect_f.csv',
                           sep='\t', index=False)
diffmod_f_ACA_intersect.to_csv('/media/labuser/Data/nanopore/m6A_classifier/data/xpore_results/diffmod_f_ACA_intersect.csv',
                 sep='\t', index=False)

######################## Negative strand ################################

m6A_seq = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/orthogonal_sites/m6A_Schwartz.csv', sep='\t')

m6A_seq = m6A_seq[m6A_seq['confGroup'] >= 0]
m6A_seq = m6A_seq[m6A_seq['strand'] == '-']

# significant sites 

diffmod_r_TGT = diffmod_r_grouped[diffmod_r_grouped['kmer'].str.contains('^TGT[A-Z]{2}$', regex=True)]

m6A_seq['position'] = m6A_seq['start'].astype(int) - 1

intersection_diffmod_mazter = []
for chromosome in m6A_seq['chr'].unique():
    m6A_seq_pos =   list(m6A_seq[m6A_seq['chr'] == chromosome]['position'])
    diffmod_r_chr = list(diffmod_r_TGT[diffmod_r_TGT['id'] == chromosome]['position'])
    overlap = set(m6A_seq_pos) & set([int(i) for i in diffmod_r_chr])
    print(m6A_seq[m6A_seq['position'].isin(list(overlap))]['stoichiometryModelFit'])
    intersection_diffmod_mazter += list(overlap)

print(len(intersection_diffmod_mazter))

# high confident mazter-seq sites

m6A_seq = m6A_seq[m6A_seq['confGroup'] >=1]

m6A_seq_intersect_r = pd.DataFrame()
diffmod_r_TGT_intersect = pd.DataFrame()

intersection_diffmod_mazter = []
for chromosome in m6A_seq['chr'].unique():
    m6A_seq_pos = list(m6A_seq[m6A_seq['chr'] == chromosome]['position'])
    diffmod_r_chr = list(diffmod_r_TGT[diffmod_r_TGT['id'] == chromosome]['position'])
    overlap = set(m6A_seq_pos) & set([int(i) for i in diffmod_r_chr])
    intersection_diffmod_mazter += list(overlap)
    m6A_seq_intersect_r = pd.concat([m6A_seq_intersect_r, m6A_seq[(m6A_seq['chr'] == chromosome) & (m6A_seq['position'].isin(list(overlap)))]])
    print(m6A_seq[m6A_seq['position'].isin(list(overlap))]['stoichiometryModelFit'])
    diffmod_r_TGT_intersect =  pd.concat([diffmod_r_TGT_intersect, diffmod_r_TGT[(diffmod_r_TGT['id'] == chromosome) & (diffmod_r_TGT['position'].isin(list(overlap)))]])

print(len(intersection_diffmod_mazter))


# Saving dataframes negative strand

m6A_seq_intersect_r.to_csv('/media/labuser/Data/nanopore/m6A_classifier/data/xpore_results/m6A_seq_intersect_r.csv',
                           sep='\t', index=False)

diffmod_r_TGT_intersect.to_csv('/media/labuser/Data/nanopore/m6A_classifier/data/xpore_results/diffmod_r_TGT_intersect.csv',
                 sep='\t', index=False)

#######################  check the other confident places. ###############################

mmc1 = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/orthogonal_sites/mmc1.csv', sep='\t')
mmc1.columns

# positive strand 
mmc1['position'] = mmc1['Peak genomic coordinate'] - 1

mmc1_intersect_f = pd.DataFrame()
diffmod_f_ACA_intersect_mmc1 = pd.DataFrame()

intersection_diffmod_mazter = []
for chromosome in mmc1['Peak chr'].unique():
    mmc1_pos = list(mmc1[mmc1['Peak chr'] == chromosome]['position'])
    diffmod_f_chr = list(diffmod_f_grouped[diffmod_f_grouped['id'] == chromosome]['position'])
    overlap = set(mmc1_pos) & set(diffmod_f_chr)
    intersection_diffmod_mazter += list(overlap)
    mmc1_intersect_f = pd.concat([mmc1_intersect_f, mmc1[(mmc1['Peak chr'] == chromosome) & (mmc1['position'].isin(list(overlap)))]])
    diffmod_f_ACA_intersect_mmc1 =  pd.concat([diffmod_f_ACA_intersect_mmc1, diffmod_f_grouped[(diffmod_f_grouped['id'] == chromosome) &\
                                                                                               (diffmod_f_grouped['position'].isin(list(overlap)))]])
    
print(len(intersection_diffmod_mazter))

mmc1_intersect_f.to_csv('/media/labuser/Data/nanopore/m6A_classifier/data/xpore_results/mmc1_intersect_f.csv',
                           sep='\t', index=False)

diffmod_f_ACA_intersect_mmc1.to_csv('/media/labuser/Data/nanopore/m6A_classifier/data/xpore_results/diffmod_f_ACA_intersect_mmc1.csv',
                 sep='\t', index=False)


# negative strand 
mmc1['position'] = mmc1['Peak genomic coordinate'] - 1
intersection_diffmod_mazter = []

mmc1_intersect_r = pd.DataFrame()
diffmod_r_TGT_intersect_mmc1 = pd.DataFrame()

intersection_diffmod_mazter = []
for chromosome in mmc1['Peak chr'].unique():
    mmc1_pos = list(mmc1[mmc1['Peak chr'] == chromosome]['position'])
    diffmod_r_chr = list(diffmod_r_grouped[diffmod_r_grouped['id'] == chromosome]['position'])
    overlap = set(mmc1_pos) & set(diffmod_r_chr)
    intersection_diffmod_mazter += list(overlap)
    mmc1_intersect_r = pd.concat([mmc1_intersect_r, mmc1[(mmc1['Peak chr'] == chromosome) & (mmc1['position'].isin(list(overlap)))]])
    diffmod_r_TGT_intersect_mmc1 =  pd.concat([diffmod_r_TGT_intersect_mmc1, diffmod_r_grouped[(diffmod_r_grouped['id'] == chromosome) &\
                                                                                               (diffmod_r_grouped['position'].isin(list(overlap)))]])
print(len(intersection_diffmod_mazter))

mmc1_intersect_r.to_csv('/media/labuser/Data/nanopore/m6A_classifier/data/xpore_results/mmc1_intersect_r.csv',
                           sep='\t', index=False)

diffmod_r_TGT_intersect_mmc1.to_csv('/media/labuser/Data/nanopore/m6A_classifier/data/xpore_results/diffmod_r_TGT_intersect_mmc1.csv',
                 sep='\t', index=False)


#################### concatenate all of the diffmod sites #####################
pd.concat([diffmod_f_ACA_intersect, diffmod_f_ACA_intersect_mmc1]).drop_duplicates().to_csv('/media/labuser/Data/nanopore/m6A_classifier/data/xpore_results/all_intersection_f.csv',
                 sep='\t', index=False)

pd.concat([diffmod_r_TGT_intersect, diffmod_r_TGT_intersect_mmc1]).drop_duplicates().to_csv('/media/labuser/Data/nanopore/m6A_classifier/data/xpore_results/all_intersection_r.csv',
                 sep='\t', index=False)











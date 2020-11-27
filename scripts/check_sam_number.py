#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 12:05:34 2020

@author: labuser
"""
import pandas as pd 
import os
import seaborn as sns

folder = '/media/labuser/Data/nanopore/m6A_classifier/data/orthogonal_sites/coverages_samfiles_mmc1_pos'


dic_counts = {}
for file in os.listdir(folder):
    file_read = folder+'/'+file
    with open(file_read, 'r') as file_in:
        for line in file_in:
            line = line.rstrip()
            if line.split('\t')[0]+'_'+line.split('\t')[1] in dic_counts:
                dic_counts[str(line.split('\t')[0])+'_'+str(line.split('\t')[1])] += int(line.split('\t')[2])
            else:
                dic_counts[str(line.split('\t')[0])+'_'+str(line.split('\t')[1])] = int(line.split('\t')[2])

distribution = []
for site in dic_counts.keys():
    distribution.append(dic_counts[site])
    
plt.xlim(0,100)
sns.distplot(distribution)
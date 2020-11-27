#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 15:35:57 2020

@author: labuser
"""

import pandas
import sys

# import the orthogonal method
Schwartz = pd.read_csv('/media/labuser/Data/nanopore/m6A_sSchwartz.csv', sep='\t')

dic_pos = {}

for i in Schwartz.itertuples():
    if i[2] in dic_pos:
        dic_pos[i[2]] += [i[4]]
    else:
        dic_pos[i[2]] = [i[4]]

mapped_reads = sys.argv[1]
mapped_reads = '/media/labuser/Data/nanopore/m6A_classifier/data/WT2_mapped.sam'

with open(mapped_reads[:-4]+'_coordinates.sam', 'w') as file_out:
    with open(mapped_reads, 'r') as file_in:
        for line in file_in: 
            line = line.rstrip().split('\t')
            try:
                if len(list(set(dic_pos[line[2]]).intersection(set(list(range(int(line[3]),int(line[3])+ len(line[9]))))))) > 0:
                    print('\t'.join(line), file=file_out)
            except:
                continue
                    
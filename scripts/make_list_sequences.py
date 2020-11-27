#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 14:03:12 2020

@author: labuser
"""

import os

folder = '/media/labuser/Data/nanopore/m6A_classifier/data/Epinano/singleA/mod_rep1_eventalign_numpy_sites'

files = os.listdir(folder)

sequences = []
for file in files:
    sequences.append(file.split('_')[-1])


with open('/media/labuser/Data/nanopore/m6A_classifier/sequences/m6A_classifier_sequences.txt', 'w') as file_out:
    counter = 0
    for sequence in set(sequences):
        print('>'+str(counter), file=file_out)
        print(sequence, file=file_out)
        counter +=1



##### check the model sequences in mazter seq

m6A_seq = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/orthogonal_sites/m6A_Schwartz.csv', sep='\t')

m6A_seq = m6A_seq[m6A_seq['confGroup'] >= 1]


for row in m6A_seq.itertuples():
    if row[8][4:7] in [i[3:-3] for i in sequences]:
        print(row)
    
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 10:58:32 2020

@author: labuser
"""


Epinano_ref = '/media/labuser/Data/nanopore/Epinanot_IVT/GSE124309_FASTA_sequences_of_Curlcakes.txt'

from Bio import SeqIO

position_single_A = {}

for record in SeqIO.parse(Epinano_ref, "fasta"):
    for i in enumerate(record.seq):
        try:
            if i[1] == 'A':
                neighbor = record.seq[i[0]-4:i[0]+4]
                if 'A' in neighbor:
                    continue
                else:
                    if record.id in position_single_A:
                        print(i[0])
                        #position_single_A[record.id] += [i[0]]
                    #else:
                        #position_single_A[record.id] = [i[0]]
                     #   print(i[0])



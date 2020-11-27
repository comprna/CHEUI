#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 12:08:09 2020

@author: labuser
"""

import pandas as pd

file = '/media/labuser/Data/nanopore/m6A_classifier/data/m6A_Schwartz.csv'

df = pd.read_csv(file, sep='\t')

df.dropna(subset =["stoichiometryModelFit"], inplace=True, axis=0)

df = df[df['confGroup'] >=0]

df.to_csv('/media/labuser/Data/nanopore/m6A_classifier/data/m6A_Schwartz_ConfFiltered.csv', sep='\t')

## make .bed file

df_bed = df['chr']+'\t'+df['start'].astype(int).astype(str)+'\t'+(df['start']+1).astype(int).astype(str)

file = '/media/labuser/Data/nanopore/m6A_classifier/data/mmc1.csv'

df = pd.read_csv(file, sep='\t')
df_bed2 = df['Peak chr']+'\t'+df['Peak genomic coordinate'].astype(int).astype(str)+'\t'+(df['Peak genomic coordinate']+1).astype(int).astype(str)

bed_final = pd.concat([df_bed, df_bed2], axis=0, ignore_index=True)

bed_final.to_csv('/media/labuser/Data/nanopore/m6A_classifier/data/6mA_total.bed', 
                 index=False, 
                 header=False)






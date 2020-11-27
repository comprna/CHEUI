#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 09:31:04 2020

@author: labuser
"""

import pandas as pd
 
df_mod_2 = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/Epinano/doubleAA/train_AA_mod_2.csv',
                       sep='\t',
                       )

df_no_mod_2 = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/Epinano/doubleAA/train_AA_no_mod_2.csv',
                          sep='\t',
                          )

df = df_mod_2.append(df_no_mod_2)

df.sample(frac=1).round(3).to_csv('/media/labuser/Data/nanopore/m6A_classifier/data/Epinano/doubleAA/train_AA_shuffle.csv',
                                   mode='w', 
                                   sep='\t',
                                   index=None)



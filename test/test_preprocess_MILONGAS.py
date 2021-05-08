#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 13:32:57 2021

@author: labuser
"""

import pytest  # pylint: disable=unused-import
import _pickle as cPickle
import subprocess
import os

import sys
sys._run_from_cmdl = True  # pylint: disable=protected-access

def test_preprocess_MILONGAS():
    
    # Delete the previous generated files
    if os.path.exists('./pytest_signals.p') is True:
        os.remove('./pytest_signals.p')
        os.remove('./pytest_IDs.p')
  
    # run preprocess_MILONGAS
    bashCmd = ['python', 'preprocess_CHEUI_parallel.py' ,'-i' ,'chr1_human_ivt_test_head.txt',\
               '-m', '../KMER_models/model_kmer.csv','-o', './', '-n', '2', '-s', 'pytest'
               ]
    
    process = subprocess.Popen(bashCmd, stdout=subprocess.PIPE)
    
    output, error = process.communicate()
    
    # load ground thruth files
    IDs = []
    signals = []
    with open('./test_old_IDs.p', 'rb') as id_in:
        with open('./test_old_signals.p', 'rb') as signal_in:
            while True:
                try:
                    IDs.append(cPickle.load(id_in))
                    signals.append(cPickle.load(signal_in))
                except:
                    break
    
    # load tested thruth files
    IDs_pytest = []
    signals_pytest = []
    with open('./pytest_signals.p', 'rb') as signal_in:
        with open('./pytest_IDs.p', 'rb') as id_in:
            while True:
                try:
                    IDs_pytest.append(cPickle.load(id_in))
                    signals_pytest.append(cPickle.load(signal_in))
                except:
                    break
    
    assert len(IDs) == len(IDs_pytest)
    assert len(signals) == len(signals_pytest)





















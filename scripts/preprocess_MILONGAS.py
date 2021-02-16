#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 16:43:33 2021

@author: labuser
"""
import sys


def _check_line(line,
                contig_idx,
                position_idx,
                reference_kmer_idx,
                model_kmer_idx,
                samples_idx):
    """
    Do a regex search against all defined regexes and
    return the key and match result of the first matching regex

    """
    line_split = line.rstrip().split('\t')
    
    # check model and reference are the same
    if line_split[reference_kmer_idx] != line_split[model_kmer_idx]:
        return None
    # check the model is not NNNNN
    if line_split[model_kmer_idx] == 'NNNNN':
        return None
    # check if there is any A in the model
    if 'A' not in line_split[model_kmer_idx]:
        return None
    
    return line_split


def _parse_line(line,
                contig_idx,
                position_idx,
                reference_kmer_idx,
                model_kmer_idx,
                samples_idx):
    """
    Do a regex search against all defined regexes and
    return the key and match result of the first matching regex

    """
    
    
    return line_split


def parse_nanopolish(filepath):
    """
    Parse nanopolish

    Parameters
    ----------
    filepath : str
        Filepath 

    Returns
    -------
    data :

    """
    data = []  # create an empty list to collect the data
    # open the file and read through it line by line
    grab_kmer = None
    line_counter = 0
    with open(filepath, 'r') as file_object:
        line = file_object.readline()
        # check the header
        header = line.rstrip().split('\t')
        if header[-1] != 'samples':
            print('nanopolish samples are not found, please run nanopolish with flag --samples')
            sys.exit()
        # get the column
        try:
            contig_idx = header.index('contig')
            position_idx = header.index('position')
            reference_kmer_idx = header.index('reference_kmer')
            model_kmer_idx = header.index('model_kmer')
            samples_idx = header.index('samples')
        except:
            print('Some nanopolish columns are not found')
    
        while line != '':  # The EOF char is an empty string
            line = file_object.readline()
            
            checked_line = _check_line(line, 
                                       contig_idx,
                                       position_idx,
                                       reference_kmer_idx,
                                       model_kmer_idx,
                                       samples_idx)
            if checked_line:
                if checked_line[model_kmer_idx][-1] == 'A':
                    
                    
                    
                    
                    
                   
                    
                
        
if __name__ == '__main__':
    
    
    nanopolish = sys.argv[1]
    filepath = '/media/labuser/Data/nanopore/m6A_classifier/data/yeast/nanopolish_reads/'\
                 'KO1_nanopolish_samples_nm6A_pos_complete_names_filtered_HighCong.txt'
                 
    filepath =  '/home/pablo/lib/MILONGAS/m6Aclassifier/nanopolish_test.txt'
                 
    data = parse_nanopolish(filepath)
    
    
    
    
    
    
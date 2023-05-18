#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 11:57:37 2023

@author: Akanksha Srivastava
"""
import os
import pickle
import argparse

parser = argparse.ArgumentParser(prog='combine_split_files', description=
                                 """ 
                                 This script takes the path of folder with individual split binary files from CHEUI C++ preprocessing script
                                 
                                 """, usage='combine_split_files.py  -i <path of directory with  individual binary files > '\
                                            ' -o < output file path > ')

OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')


REQUIRED.add_argument("-i", "--input_dir",
                      help="path of directory with  individual binary files",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-o", "--output_file",
                      help="output file path",
                      metavar='\b',
                      required=True)

ARGS = parser.parse_args()

# required arg
input_path = ARGS.input_dir

file_out = ARGS.output_file
# Define the folder containing the pickle files


# Define a function that yields the contents of each file one at a time
def load_data(folder):
    for file in os.listdir(folder):
        if file.endswith('.p'):
            with open(os.path.join(folder, file), 'rb') as f:
                while True:
                    try:
                        yield pickle.load(f)
                    except:
                        break


# Save the combined data to a new pickle file
if os.path.exists(file_out):
    print("The file already exist")
else:
    with open(file_out, 'ab') as f:
        for dataset in load_data(input_path):
            pickle.dump(dataset, f)

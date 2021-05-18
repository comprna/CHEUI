#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 09:21:05 2021

@author: labuser
"""
"""
Created on Tue Feb 16 16:43:33 2021

@author: Pablo Acera
"""

import argparse

parser = argparse.ArgumentParser(prog='preprocess_MILONGAS v0.1', description=
                                 """ 
                                 This script takes a nanopolish file, parse it and extract 
                                 information about signals and read iDs, creating two files. 
                                 
                                 """, usage='python preprocess_MILONGAS.py -i <nanopolish_file.txt> '\
                                     '-m <kmer_model> -o <out_dir> \nversion: %(prog)s')

    
OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')

#Inputs
## CHANGE -m  -t -l -f to OPTIONAL and CREATE RELATIVE PATHS FOR THESE FILES

REQUIRED.add_argument("-i", "--input_nanopolish",
                      help="Nanopolish file. Run nanopolish with the following flags: "\
                          " nanopolish eventalign --reads <in.fasta>" \
                         "--bam <in.bam> --genome <genome.fa> --print-read-names"\
                         "--scale-events --samples  > <out.txt>",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-m", "--kmer_model",
                      help="model with k-mers and expected means",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-o", "--out_dir",
                      help="output directory",
                      metavar='\b',
                      required=True)

OPTIONAL.add_argument('-v', '--version', 
                        action='version', 
                        version='%(prog)s')

OPTIONAL.add_argument("-s", "--suffix_name",
                      help='name to use for output files',
                      metavar='<str>',
                      )

OPTIONAL.add_argument("-n", "--cpu",
                      help='Number of cores to use',
                      type=int,
                      default=1
                      )                      

parser._action_groups.append(OPTIONAL)

ARGS = parser.parse_args()

# required arg
nanopolish_path = ARGS.input_nanopolish
model_kmer_path = ARGS.kmer_model
directory_out = ARGS.out_dir

# optional arg
suffix_name = ARGS.suffix_name
cpu_number = ARGS.cpu
lenght_event = 20

import sys
import os
from math import floor
import numpy as np
import _pickle as cPickle
import pandas as pd
from multiprocessing import Pool, cpu_count


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
    if 'C' not in line_split[model_kmer_idx]:
        return None
    
    return line_split

def _recycle_kmers(kmer_dict):
    '''
    '''
    keep_kmers = []
    reverse_list = list(kmer_dict.keys())[::-1]
    for i in enumerate(reverse_list):
        try:
            if i[1]-1 == reverse_list[i[0]+1]:
                keep_kmers.append(i[1])
            else:
                break
        except:
            break
    
    keep_kmers_dict = {k: v for k, v in kmer_dict.items() if k in keep_kmers}
    
    # now check if the first item has an C at the end if not, delete item
    for i in enumerate(sorted(keep_kmers_dict.keys())):
        if keep_kmers_dict[i[1]][0][-1] != 'C':
            del keep_kmers_dict[i[1]]
        
    return keep_kmers_dict



def _parse_kmers(checked_line,
                 contig_idx,
                 position_idx,
                 reference_kmer_idx,
                 model_kmer_idx,
                 samples_idx,
                 file_object,
                 parsed_kmer,
                 counter):
    """
    """
    if parsed_kmer:
        kmer_lines = parsed_kmer
        position_1 = sorted(parsed_kmer)[0]
    else:
        kmer_lines = {}
        position_1 = int(checked_line[position_idx])

    positions = [position_1, 
                 position_1+1,
                 position_1+2,
                 position_1+3,
                 position_1+4]
   
    while len(kmer_lines.keys()) < 5:
        if checked_line:
            samples = [float(i) for i in checked_line[samples_idx].split(',')]
            
            if int(checked_line[position_idx]) in kmer_lines:
                kmer_lines[int(checked_line[position_idx])] = [checked_line[model_kmer_idx], 
                                                              kmer_lines[int(checked_line[position_idx])][1] +\
                                                              samples,
                                                              checked_line[contig_idx],checked_line[3]]
            else:
                kmer_lines[int(checked_line[position_idx])] = [checked_line[model_kmer_idx], 
                                                               samples,
                                                               checked_line[contig_idx],checked_line[3]]
                
            line = file_object.readline()
            counter +=1
            if line == '':
                return {}, file_object, counter, None

            checked_line = _check_line(line, 
                                       contig_idx,
                                       position_idx,
                                       reference_kmer_idx,
                                       model_kmer_idx,
                                       samples_idx)
            try:
                if int(checked_line[1]) > max(positions):
                    break
            except:
                continue

        else:
            line = file_object.readline()
            counter +=1
            if line == '':
                return {}, file_object, counter, None

            checked_line = _check_line(line, 
                                       contig_idx,
                                       position_idx,
                                       reference_kmer_idx,
                                       model_kmer_idx,
                                       samples_idx)
            try:
                if int(checked_line[1]) > max(positions):
                    break
            except:
                continue
              
    if list(kmer_lines.keys()) == positions:
        return kmer_lines, file_object, counter, checked_line
    else:
        # check if its the end of the file
        if len(list(kmer_lines.keys())) < 5:
            return {}, file_object, counter, checked_line
        
        # grab the index of the mistmatch between the expected kmers and the 
        # ones extracted without fail
        matches =[]
        for i in enumerate(kmer_lines.keys()):
            if i[1] == positions[i[0]]:
                matches.append(i[0])
                
        # Delete before the fail, to recicly signals 
        delete_kmers = []
        for i in enumerate(matches):
            delete_kmers.append(sorted(kmer_lines.keys())[i[1]])
            
        for i in delete_kmers:
            del kmer_lines[i]
        
        # check if numbers are consecutive and take the last consecutive ones ending in A
        kmer_lines = _recycle_kmers(kmer_lines)
        
        if kmer_lines:
            return kmer_lines, file_object, counter, checked_line
        else:
            return {}, file_object, counter, checked_line



def _smooth_kmer(parsed_kmer):
    '''
    Smooth the signals to fix the lenght
    '''
    
    kmer_name = sorted(parsed_kmer.items())[0][1][0]+\
                sorted(parsed_kmer.items())[1][1][0][-1]+\
                sorted(parsed_kmer.items())[2][1][0][-1]+\
                sorted(parsed_kmer.items())[3][1][0][-1]+\
                sorted(parsed_kmer.items())[4][1][0][-1]
    
    id_kmer = list(parsed_kmer.values())[0][-2]+'_'+\
              str(sorted(parsed_kmer)[0])+'_'+kmer_name+\
                  '_'+sorted(parsed_kmer.items())[0][1][-1]
    
    
    signal_smoothed = []
    
    for pos in sorted(parsed_kmer):
        event = parsed_kmer[pos][1]
        event_smoothed = smooth_event(event) # smooth the event
        signal_smoothed += event_smoothed  # add the event to the signal
    # create an expected signal according to the kmer
    
    
    expected_smoothed = make_expected(kmer_name)
    
    # claculate distance between expected and actual signal
    distance_vector = distance_calculator(expected_smoothed,
                                          signal_smoothed)
    
    return signal_smoothed, distance_vector, id_kmer
    

def make_expected(kmer):
    
    model_kmer = pd.read_csv(model_kmer_path,sep=',')

    # create a dictionary with each kmer and its current value
    model_kmer_dict = dict(zip(model_kmer['model_kmer'], model_kmer['model_mean']))
    expected_signal = []
    for i in range(5):
        expected_signal += [model_kmer_dict[kmer[i:i+5]]]*20
    return expected_signal


def distance_calculator(signal_expected, event_smoothed):
    '''
    '''
    vector_distance = list(np.round(abs(np.array(signal_expected) - \
                                        np.array(event_smoothed)), 3))    
    return vector_distance

def smooth_event(raw_signal):
    '''
    smmoth the signal 
    '''
    raw_signal_events = []
    
    if len(raw_signal) < 20:
        event = top_median(raw_signal)
        raw_signal_events = [round(i, 3) for i in event]
        
    else:
        division = floor(len(raw_signal)/20)
        new_event = []
        for i in range(0, len(raw_signal), division):
            new_event.append(np.median(raw_signal[i:i+division]))
            if len(new_event) == 20:
                break
        if len(new_event) < 20:
            new_event = top_median(new_event)
        raw_signal_events = [round(i, 3) for i in new_event]
    return raw_signal_events


def top_median(array):
    '''
    This function top an array until some specific lenght
    '''
    extra_measure = [np.median(array)]*(20-len(array))
    array += extra_measure
    return array


def _combine_vectors(smooth_signal,
                     smooth_distance,
                     ):
    '''
    combine signals and distance vectors
    '''
    events = np.array(smooth_signal).reshape(len(smooth_signal), 1)
    distances = np.array(smooth_distance).reshape(len(smooth_distance), 1)

    combined = np.concatenate((events,
                               distances), 
                               axis=1)
    return combined

def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

def split_file(nanopolish_path, num_file):
    
    with open(nanopolish_path, 'r') as file_object:
        headerline = file_object.readline()
        headerline=headerline.rstrip()
        # check the header
        header = headerline.split('\t')
        if header[-1] != 'samples':
            print('nanopolish samples are not found, please run nanopolish with flag --samples')
            sys.exit()
        # get the column index
        
        line_count = 0
        for line in file_object:
            if line != "\n":
                line_count += 1
        line_count=line_count+1
   
    
    with open(nanopolish_path, 'r') as file_object:
        counter_reads=0
        chunk_size=line_count//num_file
        i=0
        header=None
        N_lines_print = []
        for line in file_object:
            counter_reads=counter_reads+1
            if counter_reads<chunk_size:
                
                if header:
                        N_lines_print.append(headerline)
                        N_lines_print.append('\n')
                        header=None
                        
                N_lines_print.append(line)
                
                if len(N_lines_print) >= 200000:
                    with open(directory_out+'/temp_'+str(i)+'.tmp','a') as f:
                        for row in N_lines_print:
                            f.write(row)
                        N_lines_print = []          
            else:
                if len(N_lines_print) > 0:
                    with open(directory_out+'/temp_'+str(i)+'.tmp','a') as f:
                        for row in N_lines_print:
                            f.write(row)
                        N_lines_print = []          
                            
                if 'C' in line.split('\t')[2]:
                    with open(directory_out+'/temp_'+str(i)+'.tmp','a') as f:
                        f.write(line)
                else:
                    with open(directory_out+'/temp_'+str(i)+'.tmp','a') as f:
                        f.write(line)
                    i=i+1
                    counter_reads=0
                    header=True
                    
        # in the end it could be reads in the list without being printed
        with open(directory_out+'/temp_'+str(i)+'.tmp','a') as f:
            for row in N_lines_print:
                f.write(row)
            N_lines_print = []
    return True


def parse_nanopolish(file):
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
    
    counter = 0
    parsed_kmer = {}
    stored_line = None
    signals_IDs = {}

    with open(file, 'r') as file_object:
        
        line = file_object.readline()
        
        # check the header
        header = line.rstrip().split('\t')
        if header[-1] != 'samples':
            print('nanopolish samples are not found, please run nanopolish with flag --samples')
            sys.exit()
        # get the column index
        try:
            contig_idx = header.index('contig')
            position_idx = header.index('position')
            reference_kmer_idx = header.index('reference_kmer')
            model_kmer_idx = header.index('model_kmer')
            samples_idx = header.index('samples')
        except:
            print('Some nanopolish columns are not found')
        
        while line != '':  # The EOF char is an empty string
                
            if not stored_line:
                
                line = file_object.readline()
                counter +=1
                
                if line == '':
                    break
                # check line is fordward and does not have NNNNN in the model
                checked_line = _check_line(line, 
                                           contig_idx,
                                           position_idx,
                                           reference_kmer_idx,
                                           model_kmer_idx,
                                           samples_idx)
            else:
                checked_line = stored_line
                stored_line = None
               
            if checked_line:
           
                if checked_line[model_kmer_idx][-1] == 'C' or parsed_kmer:
                   
                    parsed_kmer,\
                    file_object, \
                    counter, \
                    stored_line = _parse_kmers(checked_line,
                                               contig_idx,
                                               position_idx,
                                               reference_kmer_idx,
                                               model_kmer_idx,
                                               samples_idx,
                                               file_object,
                                               parsed_kmer,
                                               counter)

                    ### if parsed kmer fail the below code does not get executed and new_parser_kmers never get recycle
                    if parsed_kmer:
                        
                        if len(parsed_kmer) < 5:
                            continue

                        smooth_signal, smooth_distance, ID = _smooth_kmer(parsed_kmer)
                        
                        combined_signals = _combine_vectors(smooth_signal,
                                                            smooth_distance
                                                            )
                                                
                        signals_IDs[ID] = combined_signals
                        if counter%10000==0:
                        
                            with open(name_out, "ab") as sig_out:
                                cPickle.dump(signals_IDs, sig_out)
                                signals_IDs = {}
                          
                        # check if there is other As in the nine-mer to re-use lines
                        try:
                            index_ = find(ID.split('_')[-2][5:],'C')[0]-4
                            
                            use_positions = sorted(parsed_kmer.keys())[index_:]
                                
                            # Delete elements in dictionary that have a different kmer
                            # use position to select these KMERS
                            all_positions = sorted(parsed_kmer.keys()) # create all kmers to avoid error for iterating
                            for j in all_positions:
                                if j not in use_positions:
                                    del parsed_kmer[j]
                        except:
                            parsed_kmer = {}
                            continue# recover the information with the kmers
                            
            if counter%500000==0:
                print(counter, 'processed lines')
        
        if signals_IDs:
            with open(name_out, "ab") as sig_out:
                cPickle.dump(signals_IDs, sig_out)
                signals_IDs = {}

#    return True       

    
if __name__ == '__main__':

    #nanopolish_path='./chr1_human_ivt_test_head.txt'
    #model_kmer_path = './model_kmer.csv'
    #directory_out = '.'
    #suffix_name =False
    
    model_kmer = pd.read_csv(model_kmer_path,sep=',')
    
    if directory_out[0] != '.' and directory_out[0] != '/':
        directory_out = './'+directory_out
    
    if directory_out[-1] == '/':
        directory_out = directory_out[:-1]
    
    # create directory if it does not exits
    if os.path.exists(directory_out) is False:
        os.makedirs(directory_out)
    
    if suffix_name:
        name_out = directory_out+'/'+suffix_name+'_signals+IDS.p'
    else:
        name_out = directory_out+'/'+os.path.split(nanopolish_path)[1][:-4]+'_signals+IDS.p'
        
    if os.path.exists(name_out) is True:
        print()
        print('WARNING: Signal file already exist, deleting the previous generated signals+IDs files' )
        print()
        os.remove(name_out)

    if cpu_number == 1:
        parse_nanopolish(nanopolish_path)
    else:
        split_file(nanopolish_path, cpu_number)
        pathlist=[directory_out+'/'+i for i in os.listdir(directory_out) if i[-4:]=='.tmp']
        p=Pool(cpu_number)
        p.map_async(parse_nanopolish, pathlist)
        p.close()
        p.join()
        if cpu_number > 1:
            for i in pathlist:
                os.remove(i)
    



















        


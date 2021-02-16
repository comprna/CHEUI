#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 16:43:33 2021

@author: labuser
"""
import sys



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
    line_split = line.rstrip().split('\t')
    
    # check model and reference are the same
    if line_split[reference_kmer_idx] != line_split[model_kmer_idx]:
        return None 
    # check the model is not NNNNN
    if line_split[model_kmer_idx] = 
        return None 
    
    
    # if there are no matches
    return None, None


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
    line_counter = 0
    with open(filepath, 'r') as file_object:
        line_counter +=1
        line = file_object.readline()
        
        while line:
        
            if line_counter == 1: # check the header
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
            else:
                _parse_line(line, 
                            contig_idx,
                            position_idx,
                            reference_kmer_idx,
                            model_kmer_idx,
                            samples_idx)
            
                
            
       
            
            # at each line check for a match with a regex
            key, match = _parse_line(line)

            # extract school name
            if key == 'school':
                school = match.group('school')

            # extract grade
            if key == 'grade':
                grade = match.group('grade')
                grade = int(grade)

            # identify a table header 
            if key == 'name_score':
                # extract type of table, i.e., Name or Score
                value_type = match.group('name_score')
                line = file_object.readline()
                # read each line of the table until a blank line
                while line.strip():
                    # extract number and value
                    number, value = line.strip().split(',')
                    value = value.strip()
                    # create a dictionary containing this row of data
                    row = {
                        'School': school,
                        'Grade': grade,
                        'Student number': number,
                        value_type: value
                    }
                    # append the dictionary to the data list
                    data.append(row)
                    line = file_object.readline()

            line = file_object.readline()

        # create a pandas DataFrame from the list of dicts
        data = pd.DataFrame(data)
        # set the School, Grade, and Student number as the index
        data.set_index(['School', 'Grade', 'Student number'], inplace=True)
        # consolidate df to remove nans
        data = data.groupby(level=data.index.names).first()
        # upgrade Score from float to integer
        data = data.apply(pd.to_numeric, errors='ignore')
    return data

if __name__ == '__main__':
    
    
    nanopolish = sys.argv[1]
    filepath = '/media/labuser/Data/nanopore/m6A_classifier/data/yeast/nanopolish_reads/'\
                 'KO1_nanopolish_samples_nm6A_pos_complete_names_filtered_HighCong.txt'
                 
    data = parse(filepath)
    print(data)
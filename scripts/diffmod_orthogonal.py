#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 08:52:32 2020

@author: labuser

This will check the number of reads in the positions in samfiles, eventalign 
diffmod and diffmod significant
"""
import pandas as pd
import numpy as np
import os 
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


'''
# make a plot between the coverage and pvalue
diffmod_f_pval_sorted = diffmod_f_p.sort_values(by='pval_KO_vs_WT', ascending=True)

pval_f = list(diffmod_f_pval_sorted['pval_KO_vs_WT'])
#change NA with 0 for calculating the coverage
diffmod_f_pval_sorted.fillna(0, inplace=True)

WT_c_f = list(diffmod_f_pval_sorted['coverage_WT-rep1']+\
              diffmod_f_pval_sorted['coverage_WT-rep2']+\
              diffmod_f_pval_sorted['coverage_WT-rep3'])
             
KO_c_f = list(diffmod_f_pval_sorted['coverage_KO-rep1']+\
              diffmod_f_pval_sorted['coverage_KO-rep2']+\
              diffmod_f_pval_sorted['coverage_KO-rep3'])

plt.title('P_values and coverages')
plt.figure(figsize=(15,10))
sns.lineplot(np.arange(len(diffmod_f_pval_sorted)), pval_f, label='p_val')
plt.plot(np.arange(len(diffmod_f_pval_sorted)), WT_c_f,label='WT coverage')
sns.lineplot(np.arange(len(diffmod_f_pval_sorted)), KO_c_f, label='KO coverage', sort=False, ci=None)
'''

# see if there is overlap with orthogonal methods
#m6A_seq = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/orthogonal_sites/mmc1.csv', sep='\t')
#m6A_seq_pos = m6A_seq[['Peak chr', 'Peak genomic coordinate']]

#m6A_seq_pos = m6A_seq_pos.rename(columns={'Peak chr':'chr',
#                                  'Peak genomic coordinate': 'position'},
#                                )

m6A_seq = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/orthogonal_sites/m6A_Schwartz_ConfFiltered.csv', sep='\t')
m6A_seq = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/orthogonal_sites/m6A_Schwartz.csv', sep='\t')

m6A_seq = m6A_seq[m6A_seq['confGroup'] >=0]
m6A_seq =m6A_seq[m6A_seq['strand'] == '-']

#m6A_seq[m6A_seq['chr'] == 'chr7'].to_csv('/media/labuser/Data/nanopore/m6A_classifier/data/orthogonal_sites/m6A_Schwartz_ConfFiltered_chr7.bed', 
#                                     sep='\t',
#                                     index=False, 
#                                     header=False)
m6A_seq_bed = m6A_seq[['chr', 'start']]

m6A_seq_bed['start'] = m6A_seq_bed['start'].astype(int)

# the start and end position are the same, sum 1 to be able to use samtools with the bed
#m6A_seq_bed['end'] = m6A_seq_bed['start']+1
#m6A_seq_bed.to_csv('/media/labuser/Data/nanopore/m6A_classifier/data/orthogonal_sites/m6A_Schwartz_ConfFiltered.bed', 
#                   sep='\t',
#                   index=False, 
#                   header=False)

##### Check reads coverage for m6A sites in Schwartz

folder = '/media/labuser/Data/nanopore/m6A_classifier/data/orthogonal_sites/coverages_samfiles_mmc1_pos'

dic_counts_WT = {}
dic_counts_KO = {}

for file in os.listdir(folder):
    file_read = folder+'/'+file
    if 'KO' in file and 'Schwartz' in file:
        with open(file_read, 'r') as file_in:
            for line in file_in:
                line = line.rstrip()
                if line.split('\t')[0]+'_'+line.split('\t')[1] in dic_counts_KO:
                    dic_counts_KO[str(line.split('\t')[0])+'_'+str(line.split('\t')[1])] += int(line.split('\t')[2])
                else:
                    dic_counts_KO[str(line.split('\t')[0])+'_'+str(line.split('\t')[1])] = int(line.split('\t')[2])
    elif 'WT' in file and 'Schwartz' in file:
        with open(file_read, 'r') as file_in:
            for line in file_in:
                line = line.rstrip()
                if line.split('\t')[0]+'_'+line.split('\t')[1] in dic_counts_WT:
                    dic_counts_WT[str(line.split('\t')[0])+'_'+str(line.split('\t')[1])] += int(line.split('\t')[2])
                else:
                    dic_counts_WT[str(line.split('\t')[0])+'_'+str(line.split('\t')[1])] = int(line.split('\t')[2])

len(dic_counts_WT)
len(dic_counts_KO)

distribution_WT = []
for site in dic_counts_WT.keys():
    distribution_WT.append(dic_counts_WT[site])

distribution_KO = []
for site in dic_counts_KO.keys():
    distribution_KO.append(dic_counts_KO[site])


plt.figure(figsize=(15,10))
plt.xlim(0, 200)
plt.title('Bamfile coverage of m6A fodward sites')
sns.histplot(distribution_WT, color='blue', bins=3000)
sns.histplot(distribution_KO, color='red', bins=3000) 
plt.xlabel('Coverage', fontsize=15)
plt.xticks( fontsize=15)
plt.legend(title='Transcripts', loc='upper right', labels=['WT', 'KO'])
plt.savefig('/media/labuser/Data/nanopore/m6A_classifier/plots/sam_coverage_f_Schwartz.svg', format='svg', dpi=1200)



#### check the coverage in eventalign files

'''
When looking at position of the MAZTE-seq one , the position in eventalign are shifted, so 

For the forward strand pos_eventalign = MAZTER-seq -3 (to have the methylated A in the middle of the 5-mer)
'''

folder_eventaligns = '/media/labuser/Data/nanopore/m6A_classifier/data/yeast'

WT_intersection = {}
KO_intersection = {}

for file in os.listdir(folder_eventaligns):
    if 'WT' in file:
        WT_temp = pd.read_csv(folder_eventaligns+'/'+file, sep='\t')
        for chromosome in m6A_seq['Peak chr'].unique():
            WT_temp_chr = WT_temp[WT_temp['contig'] == chromosome]
            for position in WT_temp_chr['position'].unique():
                WT_temp_chr_pos = WT_temp_chr[WT_temp_chr['position'] == position]
                if str(WT_temp_chr_pos['contig'].iloc[0])+'_'+str(WT_temp_chr_pos['position'].iloc[0]) in WT_intersection:
                    WT_intersection[str(WT_temp_chr_pos['contig'].iloc[0])+'_'+str(WT_temp_chr_pos['position'].iloc[0])] += \
                        len(WT_temp_chr_pos)
                else:
                    WT_intersection[str(WT_temp_chr_pos['contig'].iloc[0])+'_'+str(WT_temp_chr_pos['position'].iloc[0])] = \
                        len(WT_temp_chr_pos)

    else:
        KO_temp = pd.read_csv(folder_eventaligns+'/'+file, sep='\t')
        for chromosome in m6A_seq['Peak chr'].unique():
            KO_temp_chr = KO_temp[KO_temp['contig'] == chromosome]
            for position in KO_temp_chr['position'].unique():
                KO_temp_ch_pos = KO_temp_chr[KO_temp_chr['position'] == position]
                if str(KO_temp_ch_pos['contig'].iloc[0])+'_'+str(KO_temp_ch_pos['position'].iloc[0]) in KO_intersection:
                    KO_intersection[str(KO_temp_ch_pos['contig'].iloc[0])+'_'+str(KO_temp_ch_pos['position'].iloc[0])] += \
                        len(KO_temp_ch_pos)
                else:
                    KO_intersection[str(KO_temp_ch_pos['contig'].iloc[0])+'_'+str(KO_temp_ch_pos['position'].iloc[0])] = \
                        len(KO_temp_ch_pos)

len(WT_intersection)
len(KO_intersection)

eventalign_WT = []
for site in WT_intersection.keys():
    eventalign_WT.append(WT_intersection[site])

eventalign_KO = []
for site in KO_intersection.keys():
    eventalign_KO.append(KO_intersection[site])


plt.figure(figsize=(15,10))
plt.xlim(0, 200)
plt.title('Eventalign coverage of m6A fodward sites')
sns.histplot(eventalign_WT, color='blue', bins=10000)
sns.histplot(eventalign_KO, color='r', bins=10000) 
plt.xlabel('Coverage', fontsize=15)
plt.xticks( fontsize=15)
plt.legend(title='Transcripts', loc='upper right', labels=['WT', 'KO'])
plt.savefig('/media/labuser/Data/nanopore/m6A_classifier/plots/eventalign_coverage_f_Schwartz.svg', format='svg', dpi=1200)



### check the diffmod table intersection with the m6A fodrward sites

main_f = '/media/labuser/Data/nanopore/m6A_classifier/data/xpore_results/all'

diffmod_f = pd.DataFrame()
diffmod_r = pd.DataFrame()

for folder in os.listdir(main_f):
    if '_f_' in folder:
        diffmod_new = pd.read_csv(main_f+'/'+folder+'/'+'diffmod.table', sep=',')
        diffmod_f = pd.concat([diffmod_f, diffmod_new], axis=0)
    else:
        diffmod_new = pd.read_csv(main_f+'/'+folder+'/'+'diffmod.table', sep=',')
        diffmod_r = pd.concat([diffmod_r, diffmod_new], axis=0)

#diffmod_f.to_csv('/media/labuser/Data/nanopore/m6A_classifier/data/xpore_results/diffmod_f.csv',
#                 sep='\t', index=False)

diffmod_f = pd.read_csv('/media/labuser/Data/nanopore/m6A_classifier/data/xpore_results/diffmod_r_pval.csv', sep=',',
                        )

diffmod_f_ACA = diffmod_f[diffmod_f['kmer'].str.contains('^TGT[A-Z]{2}$', regex=True)]
diffmod_f_ACA = diffmod_f[diffmod_f['kmer'].str.contains('^[A-Z]{2}ACA$', regex=True)]

diffmod_f_ACA['pval_KO_vs_WT']  = [float(i) for i in diffmod_f_ACA['pval_KO_vs_WT'] ]
diffmod_f_ACA_pval = diffmod_f_ACA[diffmod_f_ACA['pval_KO_vs_WT'] < 0.05]

m6A_seq['position'] = m6A_seq['start'] - 1

intersection_diffmod_mazter = []
for chromosome in m6A_seq['chr'].unique():
    m6A_seq_pos = list(m6A_seq[m6A_seq['chr'] == chromosome]['position'])
    diffmod_f_chr = list(diffmod_f_ACA[diffmod_f_ACA['id'] == chromosome]['position'])
    overlap = set(m6A_seq_pos) & set([int(i) for i in diffmod_f_chr])
    intersection_diffmod_mazter += list(overlap)

print(len(intersection_diffmod_mazter))
### check the diffmod ONLY SIGNIFICANT SITES intersection with the m6A fodrward sites

diffmod_f = diffmod_f[diffmod_f['pval_KO_vs_WT'] < 0.05]
m6A_seq = m6A_seq[m6A_seq['confGroup'] >=1]

intersection_diffmod_mazter = []
for chromosome in m6A_seq['chr'].unique():
    m6A_seq_pos = list(m6A_seq[m6A_seq['chr'] == chromosome]['position'])
    diffmod_f_chr = list(diffmod_f_ACA_pval[diffmod_f_ACA_pval['id'] == chromosome]['position'])
    overlap = set(m6A_seq_pos) & set([int(i) for i in diffmod_f_chr])
    intersection_diffmod_mazter += list(overlap)
print(len(intersection_diffmod_mazter))


'''
plot ven diagrams
venn2(subsets = (10, 5, 2), set_labels = ('Group A', 'Group B'))
plt.show()
 
MAZTER-seq    VS     diffmod(complete)
23155	      392	        14120

MAZTER-seq    VS     diffmod (significant sites)
23155         160  	        4581

MAZTER-seq(confgroup >=1)    VS     diffmod (significant sites)
673                          41         4581
'''

venn2(subsets = (23155, 893, 392), set_labels = ('MAZTER-seq', 'diffmod(all sites)'))
plt.savefig('/media/labuser/Data/nanopore/m6A_classifier/plots/venn_complete_ACA.svg', format='svg', dpi=1200)

plt.figure(figsize=(15,10))
venn2(subsets = (23155, 358, 160), set_labels = ('MAZTER-seq', 'diffmod significant sites'))
plt.savefig('/media/labuser/Data/nanopore/m6A_classifier/plots/venn_diffmod_significants_ACA.svg', format='svg', dpi=1200)

venn2(subsets = (673, 358, 41), set_labels = ('MAZTER-seq (confident sites)', 'diffmod significant sites'))
plt.tight_layout()
plt.savefig('/media/labuser/Data/nanopore/m6A_classifier/plots/venn_significant_conf_group_ACA.svg', format='svg', dpi=1200)


### select the significant ones
diffmod_f_p = diffmod_f[diffmod_f['pval_KO_vs_WT'] < 0.05]
diffmod_r_p = diffmod_r[diffmod_r['pval_KO_vs_WT'] < 0.05]

## filter sites with less than 30 reads per site in total

diffmod_f_p['coverage_WT-rep1'].fillna(0, inplace=True)
diffmod_f_p['coverage_WT-rep2'].fillna(0, inplace=True)
diffmod_f_p['coverage_WT-rep3'].fillna(0, inplace=True)

diffmod_f_p['coverage_KO-rep1'].fillna(0, inplace=True)
diffmod_f_p['coverage_KO-rep2'].fillna(0, inplace=True)
diffmod_f_p['coverage_KO-rep3'].fillna(0, inplace=True)

diffmod_r_p['coverage_WT-rep1'].fillna(0, inplace=True)
diffmod_r_p['coverage_WT-rep2'].fillna(0, inplace=True)
diffmod_r_p['coverage_WT-rep3'].fillna(0, inplace=True)

diffmod_r_p['coverage_KO-rep1'].fillna(0, inplace=True)
diffmod_r_p['coverage_KO-rep2'].fillna(0, inplace=True)
diffmod_r_p['coverage_KO-rep3'].fillna(0, inplace=True)

WT_c_f = diffmod_f_p['coverage_WT-rep1']+\
         diffmod_f_p['coverage_WT-rep2']+\
         diffmod_f_p['coverage_WT-rep3']
             
KO_c_f = diffmod_r_p['coverage_KO-rep1']+\
         diffmod_r_p['coverage_KO-rep2']+\
         diffmod_r_p['coverage_KO-rep3']

WT_c_f_30 = WT_c_f > 30
KO_c_f_30 = KO_c_f > 30

diffmod_f_p_30r = diffmod_f_p[WT_c_f_30]
diffmod_r_p_30r = diffmod_r_p[KO_c_f_30]


#m6A_seq_bed['end'] = m6A_seq_bed['Peak genomic coordinate']+1

#m6A_seq_bed.to_csv('/media/labuser/Data/nanopore/m6A_classifier/data/orthogonal_sites/mmc1.bed', sep='\t',
#                   header=None, index=None)

for chromosome in m6A_seq['chr'].unique():
    
    m6A_seq_pos = list(m6A_seq[m6A_seq['chr'] == chromosome]['start'])
    
    m6A_seq_pos_plus = [i+1 for i in m6A_seq_pos]
    m6A_seq_pos_minus = [i-1 for i in m6A_seq_pos]
    
    m6A_seq_c = m6A_seq_pos + m6A_seq_pos_plus + m6A_seq_pos_minus
    
    #diffmod_f_pos = diffmod_f[(diffmod_f['id'] == chromosome) & \
    #                          (diffmod_f['position'].isin(m6A_seq_pos))]
    
    #diffmod_f_p_30r_chr = list(diffmod_f_p_30r[diffmod_f_p_30r['id'] == chromosome]['position'])
    #diffmod_r_p_30r_chr = list(diffmod_r_p_30r[diffmod_r_p_30r['id'] == chromosome]['position'])
    
    diffmod_f_chr = list(diffmod_f_p[diffmod_f_p['id'] == chromosome]['position'])
    diffmod_r_chr = list(diffmod_r_p[diffmod_r_p['id'] == chromosome]['position'])
    
    total_diffmod = diffmod_f_chr + diffmod_r_chr
    
    overlap = set(m6A_seq_c) & set(total_diffmod)
    
    print(overlap)
    


## I have to filter out same kmers with different mod_assignment
for kmer in diffmod_f_p['kmer'].unique():
    diffmod_f_p_temp = diffmod_f_p[diffmod_f_p['kmer'] == kmer]
    print(len(diffmod_f_p_temp[diffmod_f_p_temp['mod_assignment'] == 'higher']))
    print(len(diffmod_f_p_temp[diffmod_f_p_temp['mod_assignment'] == 'lower']))
    print()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
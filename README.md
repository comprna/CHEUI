# Methylation estimation using electrical current (CHUEI) <img src="https://github.com/comprna/CHEUI/blob/master/misc/CHEUI_logo.png" width="280" height="250">


**About CHEUI**

CHEUI (Methylation (CH3) Estimation Using Ionic Current) is a two-stage deep learning method able to detect m6A and m5C transcriptome-wide at the level of individual reads and individual sites. 


------------------------------------------
# Dependencies
------------------------------------------
```
numpy==1.19.2
pandas==1.2.2
tensorflow==2.4.1
keras==2.4.3
```

------------------------------------------
# Outline of CHEUI scripts 
------------------------------------------
 <img src="https://github.com/comprna/CHEUI/blob/master/misc/pipeline_CHEUI_github.png" width="560" height="500">

## Preprocessing

CHEUI starting point is a Nanopolish output file (https://nanopolish.readthedocs.io/en/latest/).
First, a bamfile sorted and indexed with samtools is needed to run Nanopolish. 
We provide an example of how to run Nanopolish with the right flags:  

```
nanopolish index -s <sequencing_summary.txt> -d <fast5_folder> <fastq>

nanopolish eventalign -t 20 \
--reads <fastq> \
--bam <sorted.bam> \
--genome <references.fasta> \
--scale-events --signal-index  --samples --print-read-names > \
nanopolish_out.txt
```
## Download CHEUI
```
git clone https://github.com/comprna/CHEUI.git
cd CHEUI/test
```

----------------------------
# Detect m6A RNA modifications
----------------------------

## To run CHEUI to detect m6A modifications, first parse the signals centred in A nucleotides
```
python3 ../scripts/CHEUI_preprocess_m6A.py \
-i nanopolish_output_test.txt -m ../kmer_models/model_kmer.csv -o out_test_signals+IDs.p -n 15
```
## Run CHEUI model 1 for m6A to obtain m6A methylation probability per read and per 9-mer
```
python ../scripts/CHEUI_predict_model1.py -i out_test_signals+IDs.p/nanopolish_output_test_signals+IDS.p \
-m ../CHEUI_trained_models/CHEUI_m6A_model1.h5 -o ./read_level_predictions.txt
```

## We have to sort the prediction file to group all the signals from the same site
```sort -k1  --parallel=15  ./read_level_predictions.txt > ./read_level_predictions_sorted.txt```

## Now run CHEUI model 2 for m6A to get the methylation status per site
```
python3 ../scripts/CHEUI_predict_model2.py -i read_level_predictions_sorted.txt \
-m  ../CHEUI_trained_models/CHEUI_m6A_model2.h5 -o site_level_predictions.txt
```
----------------------------
# Detect m5C RNA modifications
----------------------------

To run CHEUI to detect m5C modifications, first parse the signals centred in C nucleotides
```
python3 ../scripts/CHEUI_preprocess_m5C.py \
-i nanopolish_output_test.txt -m ../kmer_models/model_kmer.csv -o out_test_signals+IDs.p -n 15
```
##  Run CHEUI model 1 for m5C to obtain m5C methylation probability per read and per 9-mer
```
python ../scripts/CHEUI_predict_model1.py -i out_test_signals+IDs.p/nanopolish_output_test_signals+IDS.p \
-m ../CHEUI_trained_models/CHEUI_m5C_model1.h5 -o ./read_level_predictions.txt
```

## We have to sort the prediction file to group all the signals from the same site
```sort -k1  --parallel=15  ./read_level_predictions.txt > ./read_level_predictions_sorted.txt```

## Now run CHEUI model 2 for m5C to get the methylation status per site
```
python3 ../scripts/CHEUI_predict_model2.py -i read_level_predictions_sorted.txt \
-m  ../CHEUI_trained_models/CHEUI_m5C_model2.h5 -o site_level_predictions.txt
```

----------------------------
Output file examples
----------------------------

An example of the read-level prediction file can be found in test/read_level_predictions.txt.
It contains 2 columns, the first column contains information about chromosome_location_k-mer_readID.
Second column contains the probability of the middle A/C of the k-mer is methylated.
```
chr10_343786_AGTACTAAG_06685ba0-c2f9-4540-9805-3e1746df432f     0.8907145261764526
chr10_343794_GAAACGGCG_06685ba0-c2f9-4540-9805-3e1746df432f     0.039006322622299194
chr10_343797_ACGGCGATG_18dad7fd-796a-4f1a-a242-27e4c5226234     0.33043956756591797
chr10_343797_ACGGCGATG_6a96d4a8-4fea-4117-b114-0e4b15537a65     0.28788983821868896
chr10_343803_ATGACAATG_06685ba0-c2f9-4540-9805-3e1746df432f     0.010302156209945679
chr10_343803_ATGACAATG_18dad7fd-796a-4f1a-a242-27e4c5226234     0.19777730107307434
chr10_343803_ATGACAATG_2d1bb785-22d1-4eeb-945d-97a9b79247d4     0.7025286555290222
```
An example of the site-level prediction file can be found in test/site_level_predictions.txt
This file is a tab separated file containing; contig, position, site, coverage, stoichiometry of the site and probability of the site being methylated.
```
contig  position        site    coverage        stoichiometry   probability
chr14   571401  TGGGCCTCC       189     0.83    0.9988387
chr12   366758  TTAACAAGA       34      1.0     0.9977418
chr12   366878  AAAACAAGA       35      1.0     0.99779636
chr12   839596  TGCACACCG       26      1.0     0.9992101
chr12   839646  ATTGCTATT       25      1.0     0.9903585
chr12   839756  AACTCGGCT       31      0.95    0.99048686
```

----------------------------
WARNINGS
----------------------------
Please follow the instructions carefully. Notice that to detect m6A or m5C a different preprocessing script is needed (CHEUI_preprocess_m6A.py/CHEUI_preprocess_m5C.py), and then the appropriate matching trained models for m6A or m5C (CHEUI_m5C_model1.h5/CHEUI_m6A_model1.h5...etc).


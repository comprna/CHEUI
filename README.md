# Methylation estimation using electrical current (CHUEI) <img src="https://github.com/comprna/CHEUI/blob/master/misc/CHEUI_logo.png" width="280" height="250">


**About CHEUI**

CHEUI (Methylation (CH3) Estimation Using Ionic Current) is a two-stage deep learning method able to detect m6A and m5C transcriptome-wide at the level of individual reads and individual sites. 


------------------------------------------
# Dependencies
------------------------------------------
```
python=<3.5
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
python3 ../scripts/CHEUI_preprocess_m6A.py -i nanopolish_output_test.txt -m ../kmer_models/model_kmer.csv -o out_test_signals+IDs.p -n 15
```
## Run CHEUI model 1 for m6A to obtain m6A methylation probability per read and per 9-mer
```
python ../scripts/CHEUI_predict_model1.py -i out_test_signals+IDs.p/nanopolish_output_test_signals+IDS.p -m ../CHEUI_trained_models/CHEUI_m6A_model1.h5 -o ./read_level_predictions.txt
```

## We have to sort the prediction file to group all the signals from the same site
```sort -k1  --parallel=15  ./read_level_predictions.txt > ./read_level_predictions_sorted.txt```

## Now run CHEUI model 2 for m6A to get the methylation status per site
```
python3 ../scripts/CHEUI_predict_model2.py -i read_level_predictions_sorted.txt -m  ../CHEUI_trained_models/CHEUI_m6A_model2.h5 -o site_level_predictions.txt -c 0.5
```
----------------------------
# Detect m5C RNA modifications
----------------------------

To run CHEUI to detect m5C modifications, first parse the signals centred in C nucleotides
```
python3 ../scripts/CHEUI_preprocess_m5C.py -i nanopolish_output_test.txt -m ../kmer_models/model_kmer.csv -o out_test_signals+IDs.p -n 15
```
##  Run CHEUI model 1 for m5C to obtain m5C methylation probability per read and per 9-mer
```
python ../scripts/CHEUI_predict_model1.py -i out_test_signals+IDs.p/nanopolish_output_test_signals+IDS.p -m ../CHEUI_trained_models/CHEUI_m5C_model1.h5 -o ./read_level_predictions.txt
```

## We have to sort the prediction file to group all the signals from the same site
```sort -k1  --parallel=15  ./read_level_predictions.txt > ./read_level_predictions_sorted.txt```

## Now run CHEUI model 2 for m5C to get the methylation status per site
```
python3 ../scripts/CHEUI_predict_model2.py -i read_level_predictions_sorted.txt -m  ../CHEUI_trained_models/CHEUI_m5C_model2.h5 -o site_level_predictions.txt
```

----------------------------
Example data files
----------------------------

An example of the read-level prediction for m6A file can be found in test/read_level_predictions.txt.
It contains 2 columns, the first column contains information about chromosome_location_k-mer_readID.
Second column contains the probability of the middle A/C of the k-mer is methylated.
```
chr10_444122_TTGTAGATG_3386eb53-8805-4c11-a721-02a23fc73cb4     0.39094510674476624
chr10_445500_TTGCAGAAA_87b56740-d8db-4d17-8cd1-aa5019d4750b     0.58213871717453
chr10_344399_TGTGAAGAA_06685ba0-c2f9-4540-9805-3e1746df432f     0.08152690529823303
chr10_343799_GGCGATGAC_18dad7fd-796a-4f1a-a242-27e4c5226234     0.5041788816452026
chr10_445260_TTTAAGAGT_a760a4ac-597e-4f57-892e-37eb0a6e1c56     0.19357600808143616
chr10_444385_CTTTACGAG_397c1862-b29c-4dd1-93a7-753da410535b     0.42834094166755676
chr10_444122_TTGTAGATG_5b66ec6a-1f4f-4638-bc69-63b54786cd6d     0.17935198545455933
```
An example of the site-level prediction file for m6A can be found in test/site_level_predictions.txt
This file is a tab separated file containing; contig, position, site, coverage, stoichiometry of the site and probability of the site being methylated.
```
contig  position        site    coverage        stoichiometry   probability
chr10   344099  TGTTAATAA       15      0.3     0.6243764
chr10   344100  GTTAATAAA       16      0.5     0.8859474
chr10   344130  AATCATAAG       15      0.6     0.6003279
chr10   344157  GAAGAGTAT       17      0.3571  0.8588939
chr10   344160  GAGTATGGG       17      0.33    0.5015969
chr10   344168  GGAAACAAC       16      0.214   0.80923474

```

----------------------------
IMPORTANT
----------------------------
Please follow the instructions carefully. Notice that to detect m6A or m5C a different preprocessing script is needed (CHEUI_preprocess_m6A.py and CHEUI_preprocess_m5C.py, respectively), and then the appropriate matching trained models for m6A or m5C must be used (CHEUI_m5C_model1.h5, CHEUI_m6A_model1.h5...etc).


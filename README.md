# CHEUI: Methylation (CH<sub>3</sub>) Estimation Using Ionic current <img src="https://github.com/comprna/CHEUI/blob/master/misc/CHEUI_logo.png" width="280" height="250">


**About CHEUI**

CHEUI (Methylation (CH<sub>3</sub>) Estimation Using Ionic current) is an RNA modification detection software. CHEUI can be used to detect m6A and m5C in individual reads at single-nucleotide resolution from any sample (e.g. single condition), or detect differential m6A or m5C between any two conditions. CHEUI uses a two-stage deep learning method to detect m6A and m5C transcriptome-wide at single-read and single-site resolution in any sequence context (i.e. without any sequence constrains).

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
# Outline of CHEUI-solo and CHEUI-diff
------------------------------------------
 <img src="https://github.com/comprna/CHEUI/blob/master/misc/pipeline_CHEUI-solo+diff_github.png" width="900" height="500">

## Before running CHEUI (IMPORTANT!)

Before running CHEUI:
1. fast5 files should be base-called, we recommend guppy version 4 or higher. 
2. Fastqs should be mapped to a reference TRANSCRIPTOME. e.g.```minimap2 -ax map-ont -k14 <reference transcript> <fastq>```
3. Run Nanopolish (https://nanopolish.readthedocs.io/en/latest/). We provide an example of how to run Nanopolish with the right flags:  
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
# IMPORTANT
----------------------------
Please follow the instructions below carefully.

1. Notice that for detecting m6A or m5C, the nanopolish output files require different preprocessing scripts: ``CHEUI_preprocess_m6A.py`` for m6A and ```CHEUI_preprocess_m5C.py``` for m5C.

2. CHEUI model 1 (read level predictions) and model 2 (site level predictions) use different predictive models for m6A and m5C that have to be specified using the --DL_model flag: 
        for m6A: ```../CHEUI_trained_models/CHEUI_m6A_model1.h5``` and ```../CHEUI_trained_models/CHEUI_m6A_model2.h5``` 
        For m5C: ```../CHEUI_trained_models/CHEUI_m5C_model1.h5``` and ```../CHEUI_trained_models/CHEUI_m5C_model2.h5```


----------------------------
# Detect m6A RNA modifications in one condition    
----------------------------

----------------------------
## CHEUI_preprocess_m6A
----------------------------

This script takes the output of nanopolish and creates a file containing signals corresponding to 9-mers centered in As and IDs.  
```
../scripts/CHEUI_preprocess_m6A.py --help

required arguments:
  -i, --input_nanopolish 
                        Nanopolish file. Run nanopolish with the following
                        flags: nanopolish eventalign --reads <in.fasta>--bam
                        <in.bam> --genome <genome.fa> --print-read-names--
                        scale-events --samples > <out.txt>
  -m, --kmer_model  file containing all the expected signal k-mer means
  -o, --out_dir     output directory

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -s <str>, --suffix_name <str>
                        name to use for output files
  -n CPU, --cpu CPU     Number of cores to use
```
```
python3 ../scripts/CHEUI_preprocess_m6A.py -i nanopolish_output_test.txt -m ../kmer_models/model_kmer.csv -o out_test_signals+IDs.p -n 15
```
----------------------------
## CHEUI_predict_model1
----------------------------
CHEUI m6A model 1 takes the previous preprocess signals and calculates m6A methylation probability per individual signal.
```
 ../scripts/CHEUI_predict_model1.py --help 
 required arguments:
  -i, --signals_input 
                        path to the ID+signal file
  -m, --DL_model    path to the m6A trainned model 1 
  -l LABEL, --label LABEL
                        label of the condition of the sample, e.g. WT_rep1
  -o, --file_out    Path to the output file
  -r, --resume          Continue running predictions

optional arguments:
  -h, --help            show help message and exit
  -v, --version         show program's version number and exit
```
```
python ../scripts/CHEUI_predict_model1.py -i out_test_signals+IDs.p/nanopolish_output_test_signals+IDS.p -m ../CHEUI_trained_models/CHEUI_m6A_model1.h5 -o ./read_level_predictions.txt -l WT_rep1
```

### Sort the predictions to group all the predictions from the same site
#### Using several replicates 
In case there are more than 1 replicate combine the two read_level_predictions files: ```cat ./read_level_predictions.txt ./read_level_predictions2.txt >> ./read_level_predictions_combined.txt``` then sort this file. 
```sort -k1  --parallel=15  ./read_level_predictions_combined.txt > ./read_level_predictions_sorted.txt```
#### Using only 1 replicate
```sort -k1  --parallel=15  ./read_level_predictions.txt > ./read_level_predictions_sorted.txt```

----------------------------
## CHEUI_predict_model2
----------------------------
CHEUI model 2 for m6A calculates the methylation probability per site
```
 ../scripts/CHEUI_predict_model2.py
 -i, --input       path to read-level prediction file from CHEUI_predict_model1.py 
 -m, --DL_model    path to pretrainned model 2
 -c, --cutoff      model 2 probability cutoff for printing sites. Default value: 0
 -d, --double_cutoff 
                       Model 1 probability cutoffs used to calculate the
                       stoichiometry. Default values: 0.3 and 0.7
 -n, --min_reads   Minimun number of reads in a site to include in the analysis. Default value: 20
 -o, --file_out    Path to the output file

optional arguments:
 -h, --help            show this help message and exit
 -v, --version         show program's version number and exit
```
```
python3 ../scripts/CHEUI_predict_model2.py -i read_level_predictions_sorted.txt -m  ../CHEUI_trained_models/CHEUI_m6A_model2.h5 -o site_level_predictions.txt -c 0.5
```
----------------------------
# Detect m5C RNA modifications in one condition
----

----------------------------
## CHEUI_preprocess_m6A
----------------------------
This script takes the output of nanopolish and creates a file containing signals corresponding to 9-mers centered in Cs and IDs.  
```
../scripts/CHEUI_preprocess_m5C.py --help

required arguments:
  -i, --input_nanopolish 
                        Nanopolish file. Run nanopolish with the following
                        flags: nanopolish eventalign --reads <in.fasta>--bam
                        <in.bam> --genome <genome.fa> --print-read-names--
                        scale-events --samples > <out.txt>
  -m, --kmer_model  file containing all the expected signal k-mer means
  -o, --out_dir     output directory

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -s <str>, --suffix_name <str>
                        name to use for output files
  -n CPU, --cpu CPU     Number of cores to use
```
```
python3 ../scripts/CHEUI_preprocess_m5C.py -i nanopolish_output_test.txt -m ../kmer_models/model_kmer.csv -o out_test_signals+IDs.p -n 15
```

----------------------------
## CHEUI_predict_model1
----------------------------
CHEUI m5C model 1 takes the previous preprocess signals and calculates m5C methylation probability per individual signal.
```
 ../scripts/CHEUI_predict_model1.py --help 
 required arguments:
  -i, --signals_input 
                        path to the ID+signal file
  -m, --DL_model    path to the m6A trainned model 1 
  -l LABEL, --label LABEL
                        label of the condition of the sample, e.g. WT_rep1
  -o, --file_out    Path to the output file
  -r, --resume          Continue running predictions

optional arguments:
  -h, --help            show help message and exit
  -v, --version         show program's version number and exit
```
```
python ../scripts/CHEUI_predict_model1.py -i out_test_signals+IDs.p/nanopolish_output_test_signals+IDS.p -m ../CHEUI_trained_models/CHEUI_m5C_model1.h5 -o ./read_level_predictions.txt -l WT_rep1
```

### Sort the predictions to group all the predictions from the same site
#### Using several replicates 
In case there are more than 1 replicate combine the two read_level_predictions files: ```cat ./read_level_predictions.txt ./read_level_predictions2.txt >> ./read_level_predictions_combined.txt``` then sort this file. 
```sort -k1  --parallel=15  ./read_level_predictions_combined.txt > ./read_level_predictions_sorted.txt```
#### Using only 1 replicate
```sort -k1  --parallel=15  ./read_level_predictions.txt > ./read_level_predictions_sorted.txt```

----------------------------
## CHEUI_predict_model2
----------------------------

CHEUI model 2 for m5C calculates the methylation probability per site
```
 ../scripts/CHEUI_predict_model2.py
 -i, --input       path to read-level prediction file from CHEUI_predict_model1.py 
 -m, --DL_model    path to pretrainned model 2
 -c, --cutoff      model 2 probability cutoff for printing sites. Default value: 0
 -d, --double_cutoff 
                       Model 1 probability cutoffs used to calculate the
                       stoichiometry. Default values: 0.3 and 0.7
 -n, --min_reads   Minimun number of reads in a site to include in the analysis. Default value: 20
 -o, --file_out    Path to the output file

optional arguments:
 -h, --help            show this help message and exit
 -v, --version         show program's version number and exit
```
```
python3 ../scripts/CHEUI_predict_model2.py -i read_level_predictions_sorted.txt -m  ../CHEUI_trained_models/CHEUI_m5C_model2.h5 -o site_level_predictions.txt -c 0.5
```
----------------------------
Example data files for CHEUI-solo
----------------------------

An example of the read-level prediction for m6A file generated using ../scripts/CHEUI_predict_model1.py 
can be found in test/read_level_predictions.txt.
It contains 2 columns, the first column contains information about contig/transcriptome_location_k-mer_readID.
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
An example of the site-level prediction file for m6A generated using ../scripts/CHEUI_predict_model2.py can be found in test/site_level_predictions.txt
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
# Identify differential m6A RNA modifications between two conditions, A and B
----------------------------

## Run CHEUI_predict_model1 for m6A, for both X and Y conditions

First use [CHEUI_preprocess_m6A](#CHEUI_preprocess_m6A) to preprocess signals centerd in A for both replicates.

```
python3 ../scripts/CHEUI_preprocess_m6A.py -i nanopolish_output_test.txt -m ../kmer_models/model_kmer.csv -o condition_X_signals+IDs.p -n 15

```
```
python3 ../scripts/CHEUI_preprocess_m6A.py -i nanopolish_output_test.txt -m ../kmer_models/model_kmer.csv -o condition_Y_signals+IDs.p -n 15
```

----------------------------
### IMPORTANT
----------------------------
Please when using ../scripts/CHEUI_predict_model1.py choose the correct --label for each condition. Later the label name will be used for the differential methylation config file.

Run [CHEUI_predict_model1](#CHEUI_predict_model1), that takes the previous preprocess signals and calculates m6A methylation probability per individual signal. For the two conditions. 
Please notice that the --label will be used later to run the differential m6A modification.

```
python ../scripts/CHEUI_predict_model1.py -i condition_X_signals+IDs.p/nanopolish_output_test_signals+IDS.p -m ../CHEUI_trained_models/CHEUI_m6A_model1.h5 -o  condition_X_read_level_predictions.txt -l X_rep1
```
```
python ../scripts/CHEUI_predict_model1.py -i condition_Y_signals+IDs.p/nanopolish_output_test_signals+IDS.p -m ../CHEUI_trained_models/CHEUI_m6A_model1.h5 -o  condition_Y_read_level_predictions.txt -l Y_rep1
```

### combine read-level probability results and sort them
```
cat condition_X_read_level_predictions.txt condition_Y_read_level_predictions.txt > CHEUI_read_level_X_Y.txt

sort -k1 --parallel=20 CHEUI_read_level_X_Y.txt > CHEUI_read_level_X_Y.sorted.txt 
```

## Run CHEUI-diff
First write a config.yml file to provide information about condition and replicates

config.yml:
```
# path to input
input: ./CHEUI_read_level_X_Y.sorted.txt 

# sample labels used to run /scripts/CHEUI_predict_model1.py
sample_labels:
    condition1:
        rep1: X_rep1
    condition2:
        rep1: Y_rep1

# cutoff used to classify methylated reads
upper_cutoff: 0.7

# cutoff used to classify unmethylated reads
lower_cutoff: 0.3
```

```
python3 ../scripts/CHEUI_diffenrentialRNAMod.py -c config.yml
```

----------------------------
Example data files for CHEUI-diff
----------------------------
```
ID                     coverage_1      coverage_2      stoichiometry_1        stoichiometry_2       stoichiometry_diff      statistic   pval_U
chr10_344212_TGGCAAATT  21               21              0.06666666666666667     0.06666666666666667     0.0                    0.0     1.0
chr10_344237_TGCTATGCC  24               24              0.0                     0.0                     0.0                    0.0     1.0
chr10_344255_CGGGACTTT  23               23              0.0                     0.0                     0.0                    0.0     1.0
chr10_344262_TTTGAAGAA  24               24              0.0                     0.0                     0.0                    0.0     1.0
```


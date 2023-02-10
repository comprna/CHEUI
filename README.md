# CHEUI: Methylation (CH<sub>3</sub>) Estimation Using Ionic current <img src="https://github.com/comprna/CHEUI/blob/master/misc/CHEUI_logo.png" width="280" height="250">


**About CHEUI**

CHEUI (Methylation (CH<sub>3</sub>) Estimation Using Ionic current) is an RNA modification detection software for Oxford Nanopore direct RNA sequencing data. CHEUI can be used to detect m6A and m5C in individual reads at single-nucleotide resolution from any sample (e.g. single condition), or detect differential m6A or m5C between any two conditions. CHEUI uses a two-stage deep learning method to detect m6A and m5C transcriptome-wide at single-read and single-site resolution in any sequence context (i.e. without any sequence constrains).

CHEUI is open source and freely available under an Academic Public License (see copy of the license in this repository).


------------------------------------------
# Dependencies

```
python=3.7
numpy==1.21.5
pandas==1.3.4
tensorflow-gpu==2.4.1
keras-preprocessing==1.1.2
```

------------------------------------------
# Outline of CHEUI-solo and CHEUI-diff

 <img src="https://github.com/comprna/CHEUI/blob/master/misc/pipeline_CHEUI-solo+diff_github.png" width="900" height="500">

## Preprocessing data before running CHEUI:

Before running CHEUI:
1. Raw signal data (fast5) should be basecalled using Guppy 4.0.11+ (4.0.11 or later) (https://community.nanoporetech.com/downloads/guppy/)(basecaller model used template_rna_r9.4.1_70bps*) 
2. Basecalled sequences (fastq) should be aligned to a reference transcriptome using minimap2 and primary, positive strand alignments should be selected, e.g.
```
minimap2 -ax map-ont -k14 <transcriptome fasta> <read fastq> | samtools view -F 2324 -b > <sorted-bam-file>
samtools index <sorted-bam-file>
```

3. Signal data should be resquiggled to aligned sequences using Nanopolish (https://nanopolish.readthedocs.io/en/latest/), ensuring that events are rescaled, e.g.
 
```
nanopolish index -s <sequencing_summary.txt> -d <fast5_folder> <read fastq>

nanopolish eventalign -t 48 \
--reads <read fastq> \
--bam <sorted-bam-file> \
--genome <transcriptome fasta> \
--scale-events --signal-index  --samples --print-read-names > nanopolish_out.txt
```
## Install CHEUI
Installation can be performed manually or using Conda (recommended).

Manual installation: 
```
git clone https://github.com/comprna/CHEUI.git
cd CHEUI/test
```

Conda installation with manual CUDA installation (recommended): 
```
conda create --name cheui python=3.7 tensorflow-gpu=2.4.1 pandas=1.3.4 -y && conda activate cheui
git clone https://github.com/comprna/CHEUI.git
cd CHEUI/test
```

Conda installation with integrated CUDA installation (not recommended): 
```
conda create --name cheui python=3.7 tensorflow-gpu=2.4.1 pandas=1.3.4 conda-forge::cudatoolkit-dev -y && conda activate cheui
git clone https://github.com/comprna/CHEUI.git
cd CHEUI/test
```

----------------------------
# IMPORTANT
----------------------------
Please follow the instructions below carefully.

1. Notice that for detecting m6A or m5C, the nanopolish output files require different preprocessing scripts: ``CHEUI_preprocess_m6A.py`` for m6A and ``CHEUI_preprocess_m5C.py`` for m5C.

2. CHEUI model 1 (read level predictions) and model 2 (site level predictions) use different predictive models for m6A and m5C that have to be specified using the --DL_model flag: 

        for m6A: 
        ```../CHEUI_trained_models/CHEUI_m6A_model1.h5``` and ```../CHEUI_trained_models/CHEUI_m6A_model2.h5``` 
        For m5C: 
        ```../CHEUI_trained_models/CHEUI_m5C_model1.h5``` and ```../CHEUI_trained_models/CHEUI_m5C_model2.h5```


----------------------------
# Detect m6A and m5C modifications in one condition    
----------------------------

----------------------------
## CHEUI preprocessing step
----------------------------

This script takes the output from nanopolish and creates a file containing signals corresponding to 9-mers centered in As and IDs.  
```
../scripts/CHEUI_preprocess_m6A.py --help

required arguments:
  -i, --input_nanopolish  Nanopolish output file. Nanopolish should be run with the following flags:
                          nanopolish eventalign --reads <in.fasta>--bam
                          <in.bam> --genome <genome.fa> --print-read-names--
                          scale-events --samples > <out.txt>
  -m, --kmer_model        file containing the expected signal k-mer means
                          (available at CHEUI/kmer_models/model_kmer.csv)
  -o, --out_dir           output directory

optional arguments:
  -h, --help              show this help message and exit
  -v, --version           show program's version number and exit
  -s <str>, --suffix_name <str>
                          name to use for output files
  -n CPU, --cpu CPU       Number of CPUs (threads) to use
```

Example command of the preprocessing step for m6A:
```
python3 ../scripts/CHEUI_preprocess_m6A.py -i nanopolish_output_test.txt -m ../kmer_models/model_kmer.csv -o out_A_signals+IDs.p -n 15
```

The processing of the Nanopolish output for m5C is very similar:
```
../scripts/CHEUI_preprocess_m5C.py --help

required arguments:
  -i, --input_nanopolish  Nanopolish output file. Nanopolish should be run with the following flags:
                          nanopolish eventalign --reads <in.fasta>--bam
                          <in.bam> --genome <genome.fa> --print-read-names--
                          scale-events --samples > <out.txt>
  -m, --kmer_model        file containing the expected signal k-mer means
                          (available at CHEUI/kmer_models/model_kmer.csv)
  -o, --out_dir           output directory

optional arguments:
  -h, --help              show this help message and exit
  -v, --version           show program's version number and exit
  -s <str>, --suffix_name <str>
                          name to use for output files
  -n CPU, --cpu CPU       Number of cores to use
```


Example command of the preprocessing step for m5C:
```
python3 ../scripts/CHEUI_preprocess_m5C.py -i nanopolish_output_test.txt -m ../kmer_models/model_kmer.csv -o out_C_signals+IDs.p -n 15
```

----------------------------
### CHEUI preprocessing step -- C++ version
----------------------------
A faster method to run the CHEUI preprocessing step. The C++ version is 2-10x times faster than the python version. 

Installation
```
cd ../scripts/preprocessing_CPP/
./build.sh
```

Parameters of the program
```
$ ./CHEUI -h
required arguments:
  -i, --input-nanopolish  Nanopolish output file. Nanopolish should be run with the following flags:
                          nanopolish eventalign --reads <in.fasta>--bam
                          <in.bam> --genome <genome.fa> --print-read-names--
                          scale-events --samples > <out.txt>
  -m, --kmer-model        file containing the expected signal k-mer means
                          (available at CHEUI/kmer_models/model_kmer.csv)
  -o, --out-dir           output directory
  --m6A/--m6C             preprocessing type

optional arguments:
  -h, --help              show this help message and exit
  -s <str>, --suffix_name <str>
                          name to use for output files
  -n CPU, --cpu CPU       Number of cores to use
```


Example command of the preprocessing step for m6A:
```
./CHEUI -i ../../test/nanopolish_output_test.txt -o ../../test/out_A_signals+IDs.p/ -m ../../kmer_models/model_kmer.csv -n 16 --m6A
```

Example command of the preprocessing step for m5C:
```
./CHEUI -i ../../test/nanopolish_output_test.txt -o ../../test/out_A_signals+IDs.p/ -m ../../kmer_models/model_kmer.csv -n 16 --m5C
```

----------------------------
## CHEUI model 1: prediction of modifications in individual reads
----------------------------

CHEUI model 1 takes the output from the previous step of preprocessing signals and calculates with the m6A_model_1
the m6A methylation probability in individual read signals at each individual A nucleotide, and with the m5C_model_1
the m5C methylation probability in individual read signals at each individual C nucleotide.
```
 ../scripts/CHEUI_predict_model1.py --help 
 
 required arguments:
  -i, --signals_input      path to the ID+signal file
  -m, --DL_model           path to the trained (m6A or m5C) model 1 
  -l LABEL, --label LABEL  label of the condition of the sample, e.g. WT_rep1
  -o, --file_out           Path to the output file
  -r, --resume             Continue running predictions

optional arguments:
  -h, --help               show help message and exit
  -v, --version            show program's version number and exit
```

Example command for the prediction of m6A in individual reads at each A nucleotide:
```
python ../scripts/CHEUI_predict_model1.py -i out_A_signals+IDs.p/nanopolish_output_test_signals+IDS.p -m ../CHEUI_trained_models/CHEUI_m6A_model1.h5 -o ./read_level_m6A_predictions.txt -l WT_rep1
```

Example command for the prediction of m5C in individual reads at each C nucleotide:
```
python ../scripts/CHEUI_predict_model1.py -i out_C_signals+IDs.p/nanopolish_output_test_signals+IDS.p -m ../CHEUI_trained_models/CHEUI_m5C_model1.h5 -o ./read_level_m5C_predictions.txt -l WT_rep1
```

### Sort the predictions to group all the predictions from the same site

The output from model 1 needs to be sorted before using model 2 to predict m6A (or m5C) at 
transcriptomic sites. We provide below two ways of doing this depending of whether you have
one or more replicates:

#### Using only 1 replicate

```
sort -k1  --parallel=15  ./read_level_m6A_predictions.txt > ./read_level_m6A_predictions_sorted.txt
sort -k1  --parallel=15  ./read_level_m5C_predictions.txt > ./read_level_m5C_predictions_sorted.txt
```

#### Using several replicates 

In case there are N replicates and you want to combine them for model 2:

```
cat ./read_level_m6A_predictions_1.txt ... ./read_level_m6A_predictions_N.txt >> ./read_level_m6A_predictions_combined.txt
cat ./read_level_m5C_predictions_1.txt ... ./read_level_m5C_predictions_N.txt >> ./read_level_m5C_predictions_combined.txt
``` 

then sort these files: 
```
sort -k1  --parallel=15  ./read_level_m6A_predictions_combined.txt > ./read_level_m6A_predictions_sorted.txt
sort -k1  --parallel=15  ./read_level_m5C_predictions_combined.txt > ./read_level_m5C_predictions_sorted.txt
```

----------------------------
## CHEUI model 2: prediction of stoichiometry and modification probability at transcriptomic sites
----------------------------
CHEUI model 2 calculates the probability and stoichiometry for m6A (or m5C) at each transcriptomic site: 
each site from the transcriptome reference used for mapping above: 
```
 ../scripts/CHEUI_predict_model2.py
 -i, --input            path to read-level prediction file from CHEUI_predict_model1.py 
 -m, --DL_model         path to pretrained model 2
 -c, --cutoff           model 2 probability cutoff for printing tested sites. Default value: 0
 -d, --double_cutoff    Model 1 probability double cutoff used to calculate the
                        stoichiometry. Default values: 0.3 and 0.7 (<0.3 is considered not-modified,
                        >0.7 is considered modified, and all other reads are ignored).
 -n, --min_reads        Minimun number of reads in a site to include in the analysis. Default value: 20
 -o, --file_out         Path to the output file

optional arguments:
 -h, --help            show this help message and exit
 -v, --version         show program's version number and exit
```

Example command for the prediction of m6A probability and stoichiometry at every A nucleotide site in the reference transcriptome:
```
python3 ../scripts/CHEUI_predict_model2.py -i read_level_m6A_predictions_sorted.txt -m  ../CHEUI_trained_models/CHEUI_m6A_model2.h5 -o site_level_m6A_predictions.txt
```

Example command for the prediction of m5C probability and stoichiometry at every C nucleotide site in the reference transcriptome:
```
python3 ../scripts/CHEUI_predict_model2.py -i read_level_m5C_predictions_sorted.txt -m  ../CHEUI_trained_models/CHEUI_m5C_model2.h5 -o site_level_m5C_predictions.txt
```


----------------------------
Example output files for CHEUI models 1 and 2 (CHEUI solo outputs)
----------------------------

Example of the read-level prediction for m6A file generated using ../scripts/CHEUI_predict_model1.py 

The file contains 3 columns, the first column contains information about contig/transcriptome_location_k-mer_readID.
Second column contains the read levele and k-mer probability of the middle A (in this example) to be methylated, 
the third column contains the label of the sample:
```
ENST00000000233.10_1003_TTGTAGATG_3386eb53-8805-4c11-a721-02a23fc73cb4    0.39094510674476624 WT_1
ENST00000000233.10_1007_TTGCAGAAA_87b56740-d8db-4d17-8cd1-aa5019d4750b    0.58213871717453    WT_1
ENST00000000233.10_133_TGTGAAGAA_06685ba0-c2f9-4540-9805-3e1746df432f     0.08152690529823303 WT_1 
ENST00000000412.8_2120_GGCGATGAC_18dad7fd-796a-4f1a-a242-27e4c5226234     0.5041788816452026  WT_1
ENST00000000412.8_2120_GGCGATGAC_a760a4ac-597e-4f57-892e-37eb0a6e1c56     0.19357600808143616 WT_1
ENST00000000412.8_2120_GGCGATGAC_397c1862-b29c-4dd1-93a7-753da410535b     0.42834094166755676 WT_1
```
An example of the site-level prediction file for m6A generated using ../scripts/CHEUI_predict_model2.py

This file is a tab separated file containing; contig, position, site, coverage, stoichiometry of the site and probability of the site being methylated:

```
contig              position    site    coverage        stoichiometry   probability
ENST00000000233.10 1003     CTTGAGTAA   648              0.10132158     0.11857438
ENST00000000233.10 1003     CTTGAGTAA   648              0.10132158     0.11857438
ENST00000000233.10 1007     AGTAATAAA   628              0.32236842     0.54757184
ENST00000000233.10 1333     AAGCAGATG   467              0.25602409     0.3631263
ENST00000000412.8  2120     AGAAACCTG   63               0.05           0.038466703
ENST00000000412.8  2126     CTGGACTGA   68               0.92424242     0.9999988
ENST00000000412.8  2130     ACTGATCTT   67               0.05           0.08036163
```

----------------------------
# Identify differential  RNA modifications between two conditions
----------------------------

## Run CHEUI_predict_model1 for m6A (or m5C) for both conditions, e.g. X and Y

First use [CHEUI_preprocess_m6A](#CHEUI_preprocess_m6A) to preprocess signals centerd in A for both replicates.

For m6A:
```
python3 ../scripts/CHEUI_preprocess_m6A.py -i nanopolish_X_output_test.txt -m ../kmer_models/model_kmer.csv -o condition_X_m6A_signals+IDs.p -n 15
python3 ../scripts/CHEUI_preprocess_m6A.py -i nanopolish_Y_output_test.txt -m ../kmer_models/model_kmer.csv -o condition_Y_m6A_signals+IDs.p -n 15
```
For m5C:
```
python3 ../scripts/CHEUI_preprocess_m5C.py -i nanopolish_X_output_test.txt -m ../kmer_models/model_kmer.csv -o condition_X_m5C_signals+IDs.p -n 15
python3 ../scripts/CHEUI_preprocess_m5C.py -i nanopolish_Y_output_test.txt -m ../kmer_models/model_kmer.csv -o condition_Y_m5C_signals+IDs.p -n 15
```

----------------------------
### IMPORTANT
----------------------------
When using ``../scripts/CHEUI_predict_model1.py`` please choose the correct **label** for each condition. Later the label name will be used for the differential methylation config file.

Run [CHEUI_predict_model1](#CHEUI_predict_model1) that takes the preprocessed signals and calculates m6A (m5C) probability per individual signal for the two conditions. 

Example for m6A for the two conditions (note the use of different labels):

```
python ../scripts/CHEUI_predict_model1.py -i condition_X_m6A_signals+IDs.p/nanopolish_output_test_signals+IDS.p -m ../CHEUI_trained_models/CHEUI_m6A_model1.h5 -o  condition_X_m6A_read_level_predictions.txt -l X_rep1

python ../scripts/CHEUI_predict_model1.py -i condition_Y_m6A_signals+IDs.p/nanopolish_output_test_signals+IDS.p -m ../CHEUI_trained_models/CHEUI_m6A_model1.h5 -o  condition_Y_m6A_read_level_predictions.txt -l Y_rep1
```

### combine the read-level probability results and sort them
```
cat condition_X_m6A_read_level_predictions.txt condition_Y_m6A_read_level_predictions.txt > CHEUI_m6A_read_level_X_Y.txt

sort -k1 --parallel=20 CHEUI_m6A_read_level_X_Y.txt > CHEUI_m6A_read_level_X_Y.sorted.txt 
```

## Run CHEUI-diff
First write a **config.yml** file to provide the information about the conditions and replicates

Example of **config.yml**:
```
# path to input
input: ./CHEUI_m6A_read_level_X_Y.sorted.txt 

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

Example command to run CHEUI-diff
```
python3 ../scripts/CHEUI_diffenrentialRNAMod.py -c config.yml
```

----------------------------
Example data files for CHEUI-diff
----------------------------
```
ID                     coverage_1      coverage_2      stoichiometry_1        stoichiometry_2       stoichiometry_diff      statistic   pval_U
txt1_212_TGGCAAATT  21               21              0.06666666666666667     0.06666666666666667     0.0                    0.0     1.0
txt1_237_TGCTATGCC  24               24              0.0                     0.0                     0.0                    0.0     1.0
txt1_255_CGGGACTTT  23               23              0.0                     0.0                     0.0                    0.0     1.0
txt1_262_TTTGAAGAA  24               24              0.0                     0.0                     0.0                    0.0     1.0
```


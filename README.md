(add default values and examples of the code when running help)
# CHEUI: Methylation (CH<sub>3</sub>) Estimation Using Ionic current <img src="https://github.com/comprna/CHEUI/blob/master/misc/CHEUI_logo.png" width="280" height="250">


**About CHEUI**

CHEUI (Methylation (CH<sub>3</sub>) Estimation Using Ionic current) is a two-stage deep learning method capable of detecting m6A and m5C transcriptome-wide at nucleotide resolution, without any sequence constraint, at the level of individual reads and individual transcript sites, and in single conditions.


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

## Use Nanopolish for re-squiggling 

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
# IMPORTANT
----------------------------
Please follow the instructions below carefully. Notice that detecting m6A or m5C requires different preprocessing scripts, CHEUI_preprocess_m6A.py for m6A and CHEUI_preprocess_m5C.py for m5C, and then the appropriate matching trained models for m6A or m5C must be used (CHEUI_m5C_model1.h5, CHEUI_m6A_model1.h5...etc).


----------------------------
# Detect m6A RNA modifications in one condition    
----------------------------

### - [Run preprocessing step for m6A](#preprocessing_m6A)
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
### Run CHEUI_predict_model1 for m6A
CHEUI m6A model 1 takes the previous preprocess signals and calculates m6A methylation probability per individual signal.
```
 ../scripts/CHEUI_preprocess_m6A.py --help 
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
```sort -k1  --parallel=15  ./read_level_predictions.txt > ./read_level_predictions_sorted.txt```

### Run CHEUI_predict_model2 for m6A
CHEUI model 2 for m6A calculates the methylation probability per site
```
 ../scripts/CHEUI_predict_model2.py
 -i, --input       path to read-level prediction file from CHEUI_predict_model1.py 
 -m, --DL_model    path to pretrainned model 2
 -c, --cutoff      model 2 probability cutoff for printing sites
 -d, --double_cutoff 
                       Model 1 probability cutoffs used to calculate the
                       stoichiometry
 -n, --min_reads   Minimun number of reads in a site to include in the
                       analysis,
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
### Run preprocessing step for m5C

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
### Run CHEUI_predict_model1 for m5C
CHEUI m5C model 1 takes the previous preprocess signals and calculates m5C methylation probability per individual signal.
```
 ../scripts/CHEUI_preprocess_m5C.py --help 
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
```sort -k1  --parallel=15  ./read_level_predictions.txt > ./read_level_predictions_sorted.txt```

### Run CHEUI_predict_model2 for m5C
CHEUI model 2 for m5C calculates the methylation probability per site
```
 ../scripts/CHEUI_predict_model2.py
 -i, --input       path to read-level prediction file from CHEUI_predict_model1.py 
 -m, --DL_model    path to pretrainned model 2
 -c, --cutoff      model 2 probability cutoff for printing sites
 -d, --double_cutoff 
                       Model 1 probability cutoffs used to calculate the
                       stoichiometry
 -n, --min_reads   Minimun number of reads in a site to include in the
                       analysis,
 -o, --file_out    Path to the output file

optional arguments:
 -h, --help            show this help message and exit
 -v, --version         show program's version number and exit
```
```
python3 ../scripts/CHEUI_predict_model2.py -i read_level_predictions_sorted.txt -m  ../CHEUI_trained_models/CHEUI_m5C_model2.h5 -o site_level_predictions.txt -c 0.5
```
----------------------------
Example data files
----------------------------

An example of the read-level prediction for m6A file generated using ../scripts/CHEUI_predict_model1.py 
can be found in test/read_level_predictions.txt.
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

### Run CHEUI_predict_model1 for m6A, for both conditions

First use [CHEUI_preprocess_m6A.py](#preprocessing_m6A) to preprocess signals centerd in A for both replicates.

Run CHEUI m6A model 1, that takes the previous preprocess signals and calculates m6A methylation probability per individual signal. For the two conditions. 
Please notice that the --label will be used later to run the differential m6A modification.

```
python ../scripts/CHEUI_predict_model1.py -i condition_A_signals+IDS.p -m ../CHEUI_trained_models/CHEUI_m6A_model1.h5 -o  condition_A_read_level_predictions.txt -l A_rep1
```
```
python ../scripts/CHEUI_predict_model1.py -i condition_B_signals+IDS.p -m ../CHEUI_trained_models/CHEUI_m6A_model1.h5 -o  condition_B_read_level_predictions.txt -l B_rep1
```

### combine read-level probability results and sort them
```
cat condition_A_read_level_predictions.txt condition_B_read_level_predictions.txt > CHEUI_read_level_A_B.txt

sort -k1 --parallel=20 CHEUI_read_level_A_B.txt > CHEUI_read_level_A_B.sorted.txt 
```

### Run CHEUI-diff
First write a config.yml file to provide information about condition and replicates

config.yml:
```
# path to input
input: ./CHEUI_read_level_A_B.sorted.txt 

# sample labels used to run /scripts/CHEUI_predict_model1.py
sample_labels:
    condition1:
        rep1: A_rep1
    condition2:
        rep1: B_rep1

# cutoff used to classify methylated reads
upper_cutoff: 0.7

# cutoff used to classify unmethylated reads
lower_cutoff: 0.3
```

```
python3 ../scripts/CHEUI_diffenrentialRNAMod.py -c config.yml
```

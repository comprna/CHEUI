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
numba==0.53.1
```

------------------------------------------
# Outline of CHEUI scripts 
------------------------------------------
 <img src="https://github.com/comprna/CHEUI/blob/master/misc/pipeline_CHEUI_github.png" width="560" height="500">

## Preprocessing

CHEUI starting point is a Nanopolish output file (https://nanopolish.readthedocs.io/en/latest/).
First, a sorted and indexed bamfile with samtools is needed to run Nanopolish. 
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

## Prepare signals to detect m6A
```
python3 ../script/CHEUI_preprocess_m6A.py \
-i nanopolish_output_test.txt -m ./kmer_models/model_kmer.csv -o out_test_signals+IDs.p -n 15
```

## Run CHEUI model 1 to obtain m6A methylation probability per read per 9-mer
```
python ../scripts/CHEUI_predict_model1.py -i out_test_signals+IDs.p/nanopolish_output_test_signals+IDS.p \
-m ../CHEUI_trainned_models/CHEUI_m6A_model1.h5 -o ./read_level_predictions.txt
```

## We have to sort the prediction file to group all the signals from the same site
```sort -k1  --parallel=15  ./read_level_predictions.txt > ./read_level_predictions_sorted.txt```

## Now run CHEUI model 2 to get methylation status per site
```
python3 ../scripts/CHEUI_predict_model2.py -i read_level_predictions_sorted.txt \
-m  ../CHEUI_trainned_models/CHEUI_m6A_model2.h5 -o site_level_predictions.txt
```








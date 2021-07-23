# Methylation estimation using electrical current (CHUEI) <img src="https://github.com/comprna/CHEUI/blob/master/msc/CHEUI_logo.png" width="280" height="250">


**About CHEUI**

CHEUI (Methylation (CH3) Estimation Using Ionic Current) is a two-stage deep learning method able to detect m6A and m5C transcriptome-wide at the level of individual reads and individual sites. 

------------------------------------------
# Outline of CHEUI scripts 
------------------------------------------
 <img src="https://github.com/comprna/CHEUI/blob/master/msc/pipeline_CHEUI_github.png" width="560" height="500">

## Preprocessing

A sorted and indexed bamfile with samtools is needed to run Nanopolish (https://nanopolish.readthedocs.io/en/latest/). 
We provide an example of how to run Nanopolish with the right flags:  

```
nanopolish index -s <sequencing_summary.txt> -d <fast5_folder> <fastq>

/g/data/xc17/pm1122/lib/nanopolish/nanopolish eventalign -t 20 \
--reads ${path}/fastq/${fasta} \
--bam ${path}/${fasta}_${fast5_folder}_sorted.bam \
--genome /g/data/xc17/pm1122/datasets/ELIGOS_data/IVT_references.txt \
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
-i <Nanopolish_file.txt> \
-m ./KMER_models/model_kmer.csv \
-o out_signals+IDs.p \
-n 20
```

## Run CHEUI model 1 to obtain m6A methylation probability per read per 9-mer
```
python CHEUI_predict_model1.py \
-i test_signals_IDs.p \
-m ./CHEUI_trainned_models/MILONGAS_model1_.2.1.h5 \
-o ./read_level_predictions.txt
```

## We have to sort the prediction file to group all the signals from the same site
```sort -k1  --parallel=15  ./test_predict_model1.txt > ./test_predict_model1_sorted.txt ```

## Now run CHEUI model 2 to get methylation status per site
```python3 predict_CHEUI_model_2.py -i test_predict_model1.txt -m  CHEUI_trainned_models/CHEUI_m6A_model2.h5 -o site_level_predictions.txt```








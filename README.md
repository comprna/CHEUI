# Methylation estimation using electrical current (CHUEI) <img src="https://github.com/comprna/CHEUI/blob/master/msc/CHEUI_logo.png" width="280" height="250">


**About CHEUI**

CHEUI (Methylation (CH3) Estimation Using Ionic Current) is a two-stage deep learning method able to detect m6A and m5C transcriptome-wide at the level of individual reads and individual sites. 

------------------------------------------
# Tutorial 
------------------------------------------
 <img src="https://github.com/comprna/CHEUI/blob/master/msc/pipeline_CHEUI_github.png" width="280" height="250">

## First run Guppy, minimap2 and nanopolish 
```
fasta=all_oligos.fasq.gz
path=/g/data/xc17/pm1122/datasets/ELIGOS_data/m5C
fast5_folder=IVT_m5C_fast5

minimap2 -ax map-ont -k5 /g/data/xc17/pm1122/datasets/ELIGOS_data/IVT_references.txt \
${path}/fastq/${fasta} > \
${path}/${fasta}_${fast5_folder}.sam

samtools view -S -b ${path}/${fasta}_${fast5_folder}.sam > \
${path}/${fasta}_${fast5_folder}.bam

samtools sort ${path}/${fasta}_${fast5_folder}.bam -o \
${path}/${fasta}_${fast5_folder}_sorted.bam

samtools index ${path}/${fasta}_${fast5_folder}_sorted.bam

/g/data/xc17/pm1122/lib/nanopolish/nanopolish index -s \
${path}/fastq/sequencing_summary.txt \
-d ${path}/${fast5_folder} \
${path}/fastq/${fasta}

/g/data/xc17/pm1122/lib/nanopolish/nanopolish eventalign -t 20 \
--reads ${path}/fastq/${fasta} \
--bam ${path}/${fasta}_${fast5_folder}_sorted.bam \
--genome /g/data/xc17/pm1122/datasets/ELIGOS_data/IVT_references.txt \
--scale-events --signal-index  --samples --print-read-names > \
${path}/${fasta}_${fast5_folder}_nanopolish.txt
```

## then run preprocess_CHEUI_parallel.py
```
python3 preprocess_CHEUI_parallel_faster.py \
-i <Nanopolish_file.txt> \
-m ./KMER_models/model_kmer.csv \
-o test_signals_IDs.p \
-n 20
```
## Now run CHEUI model 1 to get m6A methylation probability per 9-mer centered in Adenine
```
python predict_CHEUI_model_1.py \
-i test_signals_IDs.p \
-m ./CHEUI_trainned_models/MILONGAS_model1_.2.1.h5 \
-o ./read_level_predictions.txt
```

## We have to sort the prediction file to group all the signals from the same site
```sort -k1  --parallel=15  ./test_predict_model1.txt > ./test_predict_model1_sorted.txt ```

## Now run CHEUI model 2 to get methylation status per site
```python3 predict_CHEUI_model_2.py -i test_predict_model1.txt -m  CHEUI_trainned_models/CHEUI_m6A_model2.h5 -o site_level_predictions.txt```








# Methylation estimation using electrical current (CHUEI) <img src="https://github.com/comprna/CHEUI/blob/master/msc/CHEUI_logo.png" width="280" height="250">


**About CHEUI**

CHEUI 

------------------------------------------
# Tutorial 
------------------------------------------

## First run Guppy, minimap2 and nanopolish 

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



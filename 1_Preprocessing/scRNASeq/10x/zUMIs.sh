#!/bin/bash
#PBS -l nodes=1:ppn=36

threadN=36
sample=10X_LL_T
sample2=HCC1395
cell_num=2684
fastq_path=/genomics/1_Projects/FDA_QC/Fastq/10x/merge

code_path=/home/genomics/bioinformatics/app/zumi/zUMIs
input_1=$fastq_path/$sample2/$sample*R1"_001.fastq.gz"
input_2=$fastq_path/$sample2/$sample*R2"_001.fastq.gz"
ref_genome=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/star/
gtf=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf
output=/genomics/1_Projects/FDA_QC/Results/10x/merge/zUMI/$sample2

mkdir -p $output

bash $code_path/zUMIs-master.sh -f $input_1 -r $input_2 -n $sample -g $ref_genome -a $gtf -c 1-16 -m 17-26 -l 101 \
-w Counting -q 10 -Q 10 -p $threadN -d 5000,10000,25000,50000,100000,150000,200000,250000 -x "--quantMode TranscriptomeSAM" -o $output









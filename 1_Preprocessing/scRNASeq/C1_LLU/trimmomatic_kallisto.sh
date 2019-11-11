#!/bin/bash
#PBS -l nodes=1:ppn=16
input_1=/genomics/1_Projects/FDA_QC/Fastq_C1/HCC1395_H3_S73_R1_001.fastq.gz
input_2=/genomics/1_Projects/FDA_QC/Fastq_C1/HCC1395_H3_S73_R2_001.fastq.gz
output_1=/genomics/1_Projects/FDA_QC/Results/C1/Trimmomatic/Pair/HCC1395_H3_S73_R1_001_paired.fq.gz
output_2=/genomics/1_Projects/FDA_QC/Results/C1/Trimmomatic/Unpair/HCC1395_H3_S73_R1_001_unpaired.fq.gz
output_3=/genomics/1_Projects/FDA_QC/Results/C1/Trimmomatic/Pair/HCC1395_H3_S73_R2_001_paired.fq.gz
output_4=/genomics/1_Projects/FDA_QC/Results/C1/Trimmomatic/Unpair/HCC1395_H3_S73_R2_001_unpaired.fq.gz
output=/genomics/1_Projects/FDA_QC/Results/C1/Trimmomatic/

mkdir -p $output/Pair
mkdir -p $output/Unpair

java -jar /home/genomics/bioinformatics/app/trimmomatic/0.35/trimmomatic-0.35.jar PE -threads 4 -phred33 $input_1 $input_2 $output_1 $output_2 $output_3 $output_4 ILLUMINACLIP:/home/genomics/bioinformatics/app/trimmomatic/0.35/adapters/NexteraPE-PE.fa:2:36:10 LEADING:10 TRAILING:10 MAXINFO:50:0.97 MINLEN:20



threadN=12
sample=HCC1395_H3_S73
sample2=HCC1395_H3
anno_category=Ensembl_GRCh38
transcript_index=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/kallisto/transcripts.idx
trim_path=/genomics/1_Projects/FDA_QC/Results/C1/Trimmomatic/Pair/
kallisto_path=/genomics/1_Projects/FDA_QC/Results/C1/Kallisto/$anno_category/Gtf/$sample2

mkdir -p $kallisto_path

kallisto quant -i $transcript_index -o $kallisto_path -t $threadN $trim_path/$sample"_R1_001_paired.fq.gz" $trim_path/$sample"_R2_001_paired.fq.gz"

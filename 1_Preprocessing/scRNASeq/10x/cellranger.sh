#!/bin/bash
#PBS -l nodes=1:ppn=24

sample=HCC1395
sample2=10X_LL_T
threadN=24

input=/genomics/1_Projects/FDA_QC/Fastq/10x/merge
output=/genomics/1_Projects/FDA_QC/Results/10x/merge/cellranger
ref_path=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0

mkdir -p $output
cd $output

cellranger count --id=$sample --transcriptome=$ref_path --fastqs=$input/$sample --sample=$sample2 --localcores=$threadN --localmem=120










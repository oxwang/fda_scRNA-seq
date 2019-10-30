#!/bin/bash
#PBS -l nodes=1:ppn=12

sample=CFG1396
anno_category=NCBI_GRCh38
threadN=12

## 1. trimming process
output=/genomics/1_Projects/FDA_QC/Results/EPIC/merge/trim/$sample
input=/genomics/1_Projects/FDA_QC/Fastq/EPIC/merge/

mkdir -p ${output}
trim_galore -q 20 --phred33 --fastqc --gzip -o $output --paired $input/$sample"_R1".fastq.gz $input/$sample"_R2".fastq.gz

## 2. bismark alignment
output=/genomics/1_Projects/FDA_QC/Results/EPIC/merge/bismark/$anno_category/$sample/Alignment
input=/genomics/1_Projects/FDA_QC/Results/EPIC/merge/trim/$sample
ref=/genomics/Ref_genome/Human/$anno_category
mkdir -p ${output}
bismark --multicore 2 --phred33-quals -o $output --gzip $ref -1 $input/$sample"_R1_val_1.fq.gz" -2 $input/$sample"_R2_val_2.fq.gz"

## 3. methylation_extractor
output=/genomics/1_Projects/FDA_QC/Results/EPIC/merge/bismark/$anno_category/$sample/Alignment/Methy_extractor
input=/genomics/1_Projects/FDA_QC/Results/EPIC/merge/bismark/$anno_category/$sample/Alignment
mkdir -p $output
bismark_methylation_extractor -o $output --comprehensive --multicore $threadN --bedGraph --buffer_size 20G --cytosine_report --genome_folder $ref $input/$sample"_R1_val_1_bismark_bt2_pe".bam

## 4. extract methylation call files
R_code=/genomics/1_Projects/Epi_clock/Code/get_meth_call.R
log_file=/genomics/1_Projects/FDA_QC/Results/EPIC/merge/meth_call/$anno_category/log/$sample"_log.txt"
log_path=/genomics/1_Projects/FDA_QC/Results/EPIC/merge/meth_call/$anno_category/log/
mkdir -p $log_path
R CMD BATCH --no-save --no-restore '--args /genomics/1_Projects/FDA_QC/Results/EPIC/merge/bismark/NCBI_GRCh38/CFG1396/Alignment/Methy_extractor/ /genomics/1_Projects/FDA_QC/Results/EPIC/merge/meth_call/NCBI_GRCh38/CFG1396_cpg.txt' $R_code $log_file






#!/bin/bash
#PBS -l nodes=1:ppn=4

threadN=4
sample=HCC1395-COL01_S1_ROW01
anno_category=Ensembl_GRCh38
transcript_index=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/kallisto/transcripts.idx
trim_path=/genomics/1_Projects/FDA_QC/Results/C1_ceber/cutadapt/$sample
kallisto_path=/genomics/1_Projects/FDA_QC/Results/C1_ceber/Kallisto/$anno_category/Gtf/$sample

mkdir -p $kallisto_path

kallisto quant -i $transcript_index -o $kallisto_path -t $threadN --single -l 500 -s 120 $trim_path/$sample"_R2_trimmed.fq.gz"

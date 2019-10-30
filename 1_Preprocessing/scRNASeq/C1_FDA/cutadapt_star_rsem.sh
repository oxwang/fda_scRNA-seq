#!/bin/bash
#PBS -l nodes=1:ppn=8

threadN=8
sample=HCC1395-COL01_S1_ROW01
anno_category=Ensembl_GRCh38
gtf=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf
ref_genome=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/star
fasta_genome=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/fasta/genome.fa
trim_path=/genomics/1_Projects/FDA_QC/Results/C1_ceber/cutadapt/$sample
star_path=/genomics/1_Projects/FDA_QC/Results/C1_ceber/Star_featureCounts/$anno_category/Gtf/$sample
star_path2=/genomics/1_Projects/FDA_QC/Results/C1_ceber/Star_RSEM/$anno_category/Gtf/$sample

mkdir -p ${output}

STAR --runThreadN $threadN --genomeDir $ref_genome --sjdbGTFfile $gtf --readFilesIn $trim_path/$sample --readFilesCommand zcat  \
--outFileNamePrefix $star_path --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts \

ref_rsem=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/genes/rsem
mkdir -p $star_path2/rsem
rsem-calculate-expression --alignments -p $threadN --output-genome-bam --sort-bam-by-coordinate \
--single-cell-prior $star_path/Aligned.toTranscriptome.out.bam $ref_rsem $star_path2/rsem/rsem_exprs


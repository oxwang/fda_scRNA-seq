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



output=/genomics/1_Projects/FDA_QC/Results/C1/Star_RSEM/Ensembl_GRCh38/Gtf/HCC1395_H3/
gtf=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf
ref_genome=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/star
fasta_genome=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/fasta/genome.fa
input_1=/genomics/1_Projects/FDA_QC/Results/C1/Trimmomatic/Pair/HCC1395_H3_S73_R1_001_paired.fq.gz
input_2=/genomics/1_Projects/FDA_QC/Results/C1/Trimmomatic/Pair/HCC1395_H3_S73_R2_001_paired.fq.gz

mkdir -p ${output}

STAR --runThreadN 12 --genomeDir $ref_genome --sjdbGTFfile $gtf --readFilesIn $input_1 $input_2 --readFilesCommand zcat  \
--outFileNamePrefix $output --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts \

bed_file=/genomics/Ref_genome/Human/NCBI_GRCh38/RefSeq2.bed

read_distribution.py -i $output/Aligned.sortedByCoord.out.bam -r $bed_file > $output/exon_intron_summary.txt

grep -Eo '[0-9]{7,10}' $output/align_summary.txt > $output/tophat_align_summary.txt

ref_rsem=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/genes/rsem
mkdir -p $output/rsem
rsem-calculate-expression --paired-end --alignments -p 12 --output-genome-bam --sort-bam-by-coordinate \
$output/Aligned.toTranscriptome.out.bam $ref_rsem $output/rsem/rsem_exprs


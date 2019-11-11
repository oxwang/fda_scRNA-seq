#!/bin/bash
#PBS -l nodes=1:ppn=8

threadN=8
sample=HCC1395-COL01_S1_ROW01
anno_category=Ensembl_GRCh38
gtf=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf
ref_genome=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/star
fasta_genome=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/fasta/genome.fa
fastq_path=/genomics/1_Projects/FDA_QC/Fastq/C1_ceber/
trim_path=/genomics/1_Projects/FDA_QC/Results/C1_ceber/cutadapt/$sample
star_path=/genomics/1_Projects/FDA_QC/Results/C1_ceber/Star_featureCounts/$anno_category/Gtf/$sample

mkdir -p $trim_path
mkdir -p $star_path
mkdir -p $star_path/featureCounts

cutadapt -a "A{10}" -m 20 -q 10,10 -o $trim_path/$sample"_R2_trimmed.fq.gz" $fastq_path/$sample"_R2".fastq.gz

STAR --runThreadN $threadN --genomeDir $ref_genome --sjdbGTFfile $gtf --readFilesIn $trim_path/$sample"_R2_trimmed.fq.gz" \
--readFilesCommand zcat --outFileNamePrefix $star_path/ --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts \

featureCounts -T $threadN -p -a $gtf -o $star_path/featureCounts/counts.txt $star_path/Aligned.sortedByCoord.out.bam
samtools sort $star_path/featureCounts/Aligned.sortedByCoord.out.bam.featureCounts.bam -o $star_path/featureCounts/assigned_sorted.bam
samtools index $star_path/featureCounts/assigned_sorted.bam


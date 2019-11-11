#!/bin/bash
#PBS -l nodes=1:ppn=12

threadN=12
sample=HCC1395_H3_S73
sample2=HCC1395_H3
anno_category=Ensembl_GRCh38
gtf=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf
ref_genome=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/star
fasta_genome=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/fasta/genome.fa
fastq_path=/genomics/1_Projects/FDA_QC/Fastq/C1/
trim_path=/genomics/1_Projects/FDA_QC/Results/C1/Trim_galore/$sample2
star_path=/genomics/1_Projects/FDA_QC/Results/C1/Star_featureCounts/$anno_category/Gtf/$sample2

mkdir -p $trim_path
mkdir -p $star_path
mkdir -p $star_path/featureCounts

trim_galore --nextera --fastqc --gzip -o $trim_path --paired $fastq_path/$sample"_R1_001".fastq.gz $fastq_path/$sample"_R2_001".fastq.gz

STAR --runThreadN $threadN --genomeDir $ref_genome --sjdbGTFfile $gtf --readFilesIn $trim_path/$sample"_R1_001_val_1.fq.gz" $trim_path/$sample"_R2_001_val_2.fq.gz" \
--readFilesCommand zcat --outFileNamePrefix $star_path/ --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts \

featureCounts -T $threadN -p -a $gtf -o $star_path/featureCounts/counts.txt $star_path/Aligned.sortedByCoord.out.bam
samtools sort $star_path/featureCounts/Aligned.sortedByCoord.out.bam.featureCounts.bam -o $star_path/featureCounts/assigned_sorted.bam
samtools index $star_path/featureCounts/assigned_sorted.bam


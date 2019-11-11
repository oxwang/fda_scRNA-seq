#!/bin/bash
#PBS -l nodes=1:ppn=16

threadN=16
sample=10X_LL_T_S2_L001
sample2=HCC1395
cell_num=2684
fastq_path=/genomics/1_Projects/FDA_QC/Fastq/10x/merge

input_1=$fastq_path/$sample2/$sample"_R1_001".fastq.gz
input_2=$fastq_path/$sample2/$sample"_R2_001".fastq.gz
ref_genome=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/star/
gtf=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf
output=/genomics/1_Projects/FDA_QC/Results/10x/merge/umi_tools/$sample2

mkdir -p $output
mkdir -p $output/STAR
mkdir -p $output/featureCounts

umi_tools whitelist -I $input_1 -p CCCCCCCCCCCCCCCCNNNNNNNNNN --plot-prefix=$output/$sample2"_BC" \
-S $output/whitelist.txt -L $output/log_whitelist.txt --log2stderr

umi_tools extract -p CCCCCCCCCCCCCCCCNNNNNNNNNN -I $input_1 -S $output/$sample"_R1_001_"extracted.fastq.gz \
--read2-in=$input_2 --read2-out=$output/$sample"_R2_001_"extracted.fastq.gz --quality-filter-threshold=10 --quality-encoding=phred33 \
--filter-cell-barcode -L $output/log_extract.txt --whitelist=$output/whitelist.txt --log2stderr

STAR --runThreadN $threadN --genomeDir $ref_genome --sjdbGTFfile $gtf --readFilesIn $output/$sample"_R2_001_"extracted.fastq.gz \
--readFilesCommand zcat --outFileNamePrefix $output/STAR/ --outFilterMultimapNmax 1 --outReadsUnmapped Fastx \
--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts

featureCounts -T $threadN -R BAM -a $gtf -o $output/featureCounts/ $output/STAR/Aligned.sortedByCoord.out.bam
samtools sort $output/featureCounts/Aligned.sortedByCoord.out.bam.featureCounts.bam -o $output/featureCounts/assigned_sorted.bam
samtools index $output/featureCounts/assigned_sorted.bam

umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I $output/featureCounts/assigned_sorted.bam -S $output/counts.tsv.gz









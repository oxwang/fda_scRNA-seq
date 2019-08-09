#!/bin/bash
#PBS -l nodes=1:ppn=4

threadN=4
sample=wta432_1
anno_category=Ensembl_GRCh38
gtf=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf
ref_genome=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/star
fasta_genome=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/fasta/genome.fa
trim_path=/genomics/1_Projects/FDA_QC/Results/WaferGen/Second_run/cutadapt/$sample
star_path=/genomics/1_Projects/FDA_QC/Results/WaferGen/Second_run/Star_featureCounts/$anno_category/Gtf/$sample
star_path2=/genomics/1_Projects/FDA_QC/Results/WaferGen/Second_run/Star_RSEM/$anno_category/Gtf/$sample

mkdir -p ${output}

#STAR --runThreadN $threadN --genomeDir $ref_genome --sjdbGTFfile $gtf --readFilesIn $trim_path/$sample --readFilesCommand zcat  \
#--outFileNamePrefix $star_path --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts \

#bed_file=/genomics/Ref_genome/Human/NCBI_GRCh38/RefSeq2.bed

#read_distribution.py -i $output/Aligned.sortedByCoord.out.bam -r $bed_file > $output/exon_intron_summary.txt

#grep -Eo '[0-9]{7,10}' $output/align_summary.txt > $output/tophat_align_summary.txt

ref_rsem=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/genes/rsem
mkdir -p $star_path2/rsem
rsem-calculate-expression --alignments -p $threadN --output-genome-bam --sort-bam-by-coordinate \
--single-cell-prior $star_path/Aligned.toTranscriptome.out.bam $ref_rsem $star_path2/rsem/rsem_exprs


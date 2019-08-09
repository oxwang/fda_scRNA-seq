#!/bin/bash
#PBS -l nodes=1:ppn=12

sample=HCC1395_1
anno_category=NCBI_GRCh38
fastq_path=/genomics/1_Projects/FDA_QC/Fastq/bulk_RNA/merge
results_path=/genomics/1_Projects/FDA_QC/Results/RNAseq/merge

input_1=$fastq_path/$sample/$sample"_R1".fastq.gz
input_2=$fastq_path/$sample/$sample"_R2".fastq.gz
output=$results_path/Trimmomatic/
output_1=$output/$sample/$sample"_R1_paired".fq.gz
output_2=$output/$sample/$sample"_R1_unpaired".fq.gz
output_3=$output/$sample/$sample"_R2_paired".fq.gz
output_4=$output/$sample/$sample"_R2_unpaired".fq.gz

mkdir -p $output/$sample

java -jar /home/genomics/bioinformatics/app/trimmomatic/0.35/trimmomatic-0.35.jar PE -threads 4 -phred33 $input_1 $input_2 $output_1 $output_2 $output_3 $output_4 ILLUMINACLIP:/home/genomics/bioinformatics/app/trimmomatic/0.35/adapters/NexteraPE-PE.fa:2:36:10 LEADING:10 TRAILING:10 MAXINFO:50:0.97 MINLEN:20


output=$results_path/Star/$anno_category/Gtf/$sample/
gtf=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf
ref_genome=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/star
input_1=$results_path/Trimmomatic/$sample/$sample"_R1_paired".fq.gz
input_2=$results_path/Trimmomatic/$sample/$sample"_R2_paired".fq.gz

mkdir -p $output

STAR --runThreadN 12 --genomeDir $ref_genome --sjdbGTFfile $gtf --readFilesIn $input_1 $input_2 --readFilesCommand zcat  \
--outFileNamePrefix $output --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts

bed_file=/genomics/Ref_genome/Human/NCBI_GRCh38/RefSeq2.bed
read_distribution.py -i $output/Aligned.sortedByCoord.out.bam -r $bed_file > $output/exon_intron_summary.txt

ref_rsem=/genomics/Ref_genome/10x_genomics/refdata-cellranger-GRCh38-1.2.0/genes/rsem
mkdir -p $output/rsem
rsem-calculate-expression --paired-end --alignments -p 12 --output-genome-bam --sort-bam-by-coordinate \
$output/Aligned.toTranscriptome.out.bam $ref_rsem $output/rsem/rsem_exprs




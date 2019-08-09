# ATAC-seq data,charles Wang, LLu
# breast cancer and B-cell raw data
cd /ddn/gs1/home/xuxiao/work/SEQCII/ATAC-seq/
mkdir HCC1395
cd HCC1395

mkdir fastq fastqc  bowtie1 fastq_cutadaptor

# fastqc 
/ddn/gs1/home/xuxiao/tools/FastQC/fastqc
for f in `seq 1 3`
do
/ddn/gs1/home/xuxiao/tools/FastQC/fastqc  -o  ./fastqc -t 6   ./fastq/HCC1395_${f}_R1.fastq.gz   ./fastq/HCC1395_${f}_R2.fastq.gz &
done

for f in `seq 1 3`
do
/ddn/gs1/home/xuxiao/tools/FastQC/fastqc  -o  ./fastqc -t 6   ./fastq/HCC1395BL_${f}_R1.fastq.gz   ./fastq/HCC1395BL_${f}_R2.fastq.gz &
done

# cutadaptor
 # Cut the adaptor
for f in `seq 1 3`
do
~/tools/trim_galore.4/trim_galore  --paired   --stringency 1   --dont_gzip --length 20 -q 20 --trim1 -o ./fastq_cutadaptor      ./fastq/HCC1395_${f}_R1.fastq.gz    ./fastq/HCC1395_${f}_R2.fastq.gz  &
done

for f in `seq 1  3`
do
~/tools/trim_galore.4/trim_galore  --paired   --stringency 1   --dont_gzip --length 20 -q 20 --trim1 -o ./fastq_cutadaptor      ./fastq/HCC1395BL_${f}_R1.fastq.gz    ./fastq/HCC1395BL_${f}_R2.fastq.gz  &
done


# align read to GRCh38 using bowtie1

 ##best hit
for f in `seq 1 3`
    do
bowtie  -p 20 -m 1 -v 2 --best --strata -X 1500  -S ~/work/Reference_genome/GRCH38/indexBowtie_GRCh38/GRCh38.p12  -1 ./fastq_cutadaptor/HCC1395_${f}_R1_val_1.fq  -2  ./fastq_cutadaptor/HCC1395_${f}_R2_val_2.fq   ./bowtie1/HCC1395_${f}.sam &
    done


for f in `seq 1 3`
do
bowtie  -p 20 -m 1 -v 2 --best --strata -X 1500  -S ~/work/Reference_genome/GRCH38/indexBowtie_GRCh38/GRCh38.p12  -1 ./fastq_cutadaptor/HCC1395BL_${f}_R1_val_1.fq  -2  ./fastq_cutadaptor/HCC1395BL_${f}_R2_val_2.fq   ./bowtie1/HCC1395BL_${f}.sam &
done

#sort
for f in `seq 1 3`
do
java -Xmx8g -jar ~/tools/picard-tools-1.115/SortSam.jar  TMP_DIR=./TMP/C_${f}    INPUT=./bowtie1/HCC1395_${f}.sam    OUTPUT=./bowtie1/HCC1395_${f}_Sorted.bam   SORT_ORDER=coordinate &
done

for f in  `seq 1  3`
do
java -Xmx8g -jar ~/tools/picard-tools-1.115/SortSam.jar  TMP_DIR=./TMP/BL${f}    INPUT=./bowtie1/HCC1395BL_${f}.sam    OUTPUT=./bowtie1/HCC1395BL_${f}_Sorted.bam   SORT_ORDER=coordinate &
done


# merge same treatment
samtools merge  ./bowtie1/Combined_HCC1395_Sorted.bam  ./bowtie1/HCC1395_1_Sorted.bam ./bowtie1/HCC1395_2_Sorted.bam  ./bowtie1/HCC1395_3_Sorted.bam&
samtools merge  ./bowtie1/Combined_HCC1395BL_Sorted.bam  ./bowtie1/HCC1395BL_1_Sorted.bam ./bowtie1/HCC1395BL_2_Sorted.bam  ./bowtie1/HCC1395BL_3_Sorted.bam&


############

java -Xmx8g -jar ~/tools/picard-tools-1.115/MarkDuplicates.jar  TMP_DIR=./TMP/HCC1395 INPUT=./bowtie1/Combined_HCC1395_Sorted.bam 	OUTPUT=./bowtie1/Combined_HCC1395_Sorted_dedup.bam   REMOVE_DUPLICATES=TRUE  METRICS_FILE=./bowtie1/summary/DeDUPLICATE_Combined_HCC1395_summary.txt  &

java -Xmx8g -jar ~/tools/picard-tools-1.115/MarkDuplicates.jar  TMP_DIR=./TMP/HCC1395BL INPUT=./bowtie1/Combined_HCC1395BL_Sorted.bam     OUTPUT=./bowtie1/Combined_HCC1395BL_Sorted_dedup.bam   REMOVE_DUPLICATES=TRUE  METRICS_FILE=./bowtie1/summary/DeDUPLICATE_Combined_HCC1395BL_summary.txt  &


# make open chromatin
bamToBed -i Combined_HCC1395_Sorted_dedup.bam  |awk '{if($6=="+") {$3=$2+9;} else $2=$3-9;}1' OFS="\t"  >./BedFile/Combined_HCC1395_Sorted_dedup_OpenChromatin.bed   &

bamToBed -i Combined_HCC1395BL_Sorted_dedup.bam  |awk '{if($6=="+") {$3=$2+9;} else $2=$3-9;}1' OFS="\t"  >./BedFile/Combined_HCC1395BL_Sorted_dedup_OpenChromatin.bed   &

# clean the chromosome

egrep "chr"  ./BedFile/Combined_HCC1395_Sorted_dedup_OpenChromatin.bed  >./BedFile/Combined_HCC1395_Sorted_dedup_OpenChromatin_clean.bed&
egrep "chr"  ./BedFile/Combined_HCC1395BL_Sorted_dedup_OpenChromatin.bed  >./BedFile/Combined_HCC1395BL_Sorted_dedup_OpenChromatin_clean.bed &

## peak calling
## peak calling
#/ddn/gs1/biotools/python/bin/macs2.1.1


EFFECTIVE_GENOME_SIZE=3.2e9
OUTDIR="/ddn/gs1/home/xuxiao/work/SEQCII/ATAC-seq/HCC1395/Peaks"
FDR=0.0001

macs2 callpeak -t ./bowtie1/BedFile/Combined_HCC1395_Sorted_dedup_OpenChromatin_clean.bed  --outdir $OUTDIR -n out_q.0001/Combined_HCC1395/MACS -f BED -g $EFFECTIVE_GENOME_SIZE -q $FDR --keep-dup=all --nomodel --extsize 9 --verbose 0 &

macs2 callpeak -t ./bowtie1/BedFile/Combined_HCC1395BL_Sorted_dedup_OpenChromatin_clean.bed  --outdir $OUTDIR -n out_q.0001/Combined_HCC1395BL/MACS -f BED -g $EFFECTIVE_GENOME_SIZE -q $FDR --keep-dup=all --nomodel --extsize 9 --verbose 0 &




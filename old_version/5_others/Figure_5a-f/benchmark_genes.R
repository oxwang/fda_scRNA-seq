## using high, intermediate, and low expressed genes as benchmark genes ##

get_genes <- function(da)
{
	da <- da
	da_avg <- apply(da,1,mean)
	da_sd <- apply(da,1,sd)
	da_quan <- quantile(da_avg,probs=0.999)

	da <- da[da_avg > 0,]
	da_sd <- da_sd[da_avg > 0]
	da_avg <- da_avg[da_avg > 0]

	da <- da[da_sd < 1,]
	da_avg <- da_avg[da_sd < 1]

	da <- da[da_avg < da_quan,]
	da_avg <- da_avg[da_avg < da_quan]
	da_order <- order(da_avg,decreasing=T)
	da <- da[da_order,]
	da_avg <- da_avg[da_order]

	pos <- round((nrow(da)-500)/2)
	gene_list <- names(da_avg[c(1:500,pos:(pos+499),(nrow(da)-499):nrow(da))])

	return(gene_list)
}


get_summary <- function(gene,da)
{
	gene <- gene
	da <- da
	index <- match(gene,rownames(da))

	gene_summary <- rep(0,length(gene))
	da2 <- da[index[!is.na(index)],]
	gene_summary[!is.na(index)] <- apply(da2 > 0,1,sum)

	gene_summary <- cbind(gene_summary,round(gene_summary/ncol(da),4))
	colnames(gene_summary) <- c("Counts","Perct")

	return(gene_summary)
}



get_summary2 <- function(gene,da)
{
	gene <- gene
	da <- da
	index <- match(gene,da[,1])

	gene_summary <- rep(0,length(gene))
	da2 <- da[index[!is.na(index)],-c(1:6)]
	gene_summary[!is.na(index)] <- apply(da2 > 0,1,sum)

	gene_summary <- cbind(gene_summary,round(gene_summary/ncol(da),4))
	colnames(gene_summary) <- c("Counts","Perct")

	return(gene_summary)
}


## generate statistics summary of number of expressed cells and percentage per gene ##
## the summary data are available in .../Figure_5a-f/data/ ##

input <- "/genomics/1_Projects/FDA_QC/Results/RNAseq/merge/Star/NCBI_GRCh38/Gtf/"
output <- "/genomics/1_Projects/FDA_QC/Results/benchmark_gene/exprs_genes/"

load(paste(input,"txi_star_rsem_hcc1395.rdata",sep=""))
exprs_gene_tumor <- log2(txi_rsem_hcc1395$abundance+1)
exprs_gene_tumor <- get_genes(exprs_gene_tumor)

load(paste(input,"txi_star_rsem_hcc1395bl.rdata",sep=""))
exprs_gene_normal <- log2(txi_rsem_hcc1395bl$abundance+1)
exprs_gene_normal <- get_genes(exprs_gene_normal)




## 10X LLU ##

#  zumi

gene_zumi <- readRDS("/genomics/1_Projects/FDA_QC/Results/10x/merge/zUMI/HCC1395/zUMIs_output/expression/10X_LL_T.dgecounts.rds")
gene_counts <- unlist(x = gene_zumi$exons$downsampled$downsampled_100000$umicounts_downsampled,recursive = F,use.names = T)
llu_10x_zumi_tumor <- get_summary(exprs_gene_tumor,gene_counts)

gene_zumi <- readRDS("/genomics/1_Projects/FDA_QC/Results/10x/merge/zUMI/HCC1395BL/zUMIs_output/expression/10X_LL_N.dgecounts.rds")
gene_counts <- unlist(x = gene_zumi$exons$downsampled$downsampled_100000$umicounts_downsampled,recursive = F,use.names = T)
llu_10x_zumi_normal <- get_summary(exprs_gene_normal,gene_counts)

colnames(llu_10x_zumi_tumor) <- c("HCC1395 expressed cells","HCC1395 cell percentage")
colnames(llu_10x_zumi_normal) <- c("HCC1395BL expressed cells","HCC1395BL cell percentage")
rownames(llu_10x_zumi_tumor) <- exprs_gene_tumor
rownames(llu_10x_zumi_normal) <- exprs_gene_normal

write.csv(llu_10x_zumi_tumor,file=paste(output,"llu_10x_zumi_tumor.csv",sep=""),quote=F)
write.csv(llu_10x_zumi_normal,file=paste(output,"llu_10x_zumi_normal.csv",sep=""),quote=F)




## 10X NCI standard ##

#  zumi

gene_zumi <- readRDS("/genomics/1_Projects/FDA_QC/Results/10x/Third_run/NCI/zUMI/autoselect/HCC1395/zUMIs_output/expression/10X_NCI_T_S1_L001.dgecounts.rds")
gene_counts <- unlist(x = gene_zumi$exons$downsampled$downsampled_100000$umicounts_downsampled,recursive = F,use.names = T)
nci_10x_zumi_tumor_standard <- get_summary(exprs_gene_tumor,gene_counts)

gene_zumi <- readRDS("/genomics/1_Projects/FDA_QC/Results/10x/Third_run/NCI/zUMI/autoselect/HCC1395BL/zUMIs_output/expression/10X_NCI_N_S2_L001.dgecounts.rds")
gene_counts <- unlist(x = gene_zumi$exons$downsampled$downsampled_100000$umicounts_downsampled,recursive = F,use.names = T)
nci_10x_zumi_normal_standard <- get_summary(exprs_gene_normal,gene_counts)

colnames(nci_10x_zumi_tumor_standard) <- c("HCC1395 expressed cells","HCC1395 cell percentage")
colnames(nci_10x_zumi_normal_standard) <- c("HCC1395BL expressed cells","HCC1395BL cell percentage")
rownames(nci_10x_zumi_tumor_standard) <- exprs_gene_tumor
rownames(nci_10x_zumi_normal_standard) <- exprs_gene_normal

write.csv(nci_10x_zumi_tumor_standard,file=paste(output,"nci_10x_zumi_tumor_standard.csv",sep=""),quote=F)
write.csv(nci_10x_zumi_normal_standard,file=paste(output,"nci_10x_zumi_normal_standard.csv",sep=""),quote=F)



## 10X NCI short ##

#  zumi

gene_zumi <- readRDS("/genomics/1_Projects/FDA_QC/Results/10x/NCI/merge/zUMI/HCC1395/zUMIs_output/expression/10X_NC_T.dgecounts.rds")
gene_counts <- unlist(x = gene_zumi$exons$downsampled$downsampled_100000$umicounts_downsampled,recursive = F,use.names = T)
nci_10x_zumi_tumor_short <- get_summary(exprs_gene_tumor,gene_counts)

gene_zumi <- readRDS("/genomics/1_Projects/FDA_QC/Results/10x/NCI/merge/zUMI/HCC1395BL/zUMIs_output/expression/10X_NC_N.dgecounts.rds")
gene_counts <- unlist(x = gene_zumi$exons$downsampled$downsampled_100000$umicounts_downsampled,recursive = F,use.names = T)
nci_10x_zumi_normal_short <- get_summary(exprs_gene_normal,gene_counts)

colnames(nci_10x_zumi_tumor_short) <- c("HCC1395 expressed cells","HCC1395 cell percentage")
colnames(nci_10x_zumi_normal_short) <- c("HCC1395BL expressed cells","HCC1395BL cell percentage")
rownames(nci_10x_zumi_tumor_short) <- exprs_gene_tumor
rownames(nci_10x_zumi_normal_short) <- exprs_gene_normal

write.csv(nci_10x_zumi_tumor_short,file=paste(output,"nci_10x_zumi_tumor_short.csv",sep=""),quote=F)
write.csv(nci_10x_zumi_normal_short,file=paste(output,"nci_10x_zumi_normal_short.csv",sep=""),quote=F)



## C1 LLU ##

#  featureCounts

counts_hcc1395 <- read.table("/genomics/1_Projects/FDA_QC/Results/C1/Star_featureCounts/Ensembl_GRCh38/bam_subsampling/100K/featureCounts/counts_hcc1395.txt",skip=1,header=T,sep="\t")
llu_c1_featureCounts_tumor <- get_summary2(exprs_gene_tumor,counts_hcc1395)

counts_hcc1395bl <- read.table("/genomics/1_Projects/FDA_QC/Results/C1/Star_featureCounts/Ensembl_GRCh38/bam_subsampling/100K/featureCounts/counts_hcc1395bl.txt",skip=1,header=T,sep="\t")
llu_c1_featureCounts_normal <- get_summary2(exprs_gene_normal,counts_hcc1395bl)

colnames(llu_c1_featureCounts_tumor) <- c("HCC1395 expressed cells","HCC1395 cell percentage")
colnames(llu_c1_featureCounts_normal) <- c("HCC1395BL expressed cells","HCC1395BL cell percentage")
rownames(llu_c1_featureCounts_tumor) <- exprs_gene_tumor
rownames(llu_c1_featureCounts_normal) <- exprs_gene_normal

write.csv(llu_c1_featureCounts_tumor,file=paste(output,"llu_c1_featureCounts_tumor.csv",sep=""),quote=F)
write.csv(llu_c1_featureCounts_normal,file=paste(output,"llu_c1_featureCounts_normal.csv",sep=""),quote=F)




## C1 FDA ##

# featureCounts

load("/genomics/1_Projects/FDA_QC/Results/C1_ceber/Star_featureCounts/Ensembl_GRCh38/featureCounts/counts_hcc1395.rdata")
fda_c1_featureCounts_tumor <- get_summary2(exprs_gene_tumor,counts_hcc1395)
load("/genomics/1_Projects/FDA_QC/Results/C1_ceber/Star_featureCounts/Ensembl_GRCh38/featureCounts/counts_hcc1395bl.rdata")
fda_c1_featureCounts_normal <- get_summary2(exprs_gene_normal,counts_hcc1395bl)

colnames(fda_c1_featureCounts_tumor) <- c("HCC1395 expressed cells","HCC1395 cell percentage")
colnames(fda_c1_featureCounts_normal) <- c("HCC1395BL expressed cells","HCC1395BL cell percentage")
rownames(fda_c1_featureCounts_tumor) <- exprs_gene_tumor
rownames(fda_c1_featureCounts_normal) <- exprs_gene_normal

write.csv(fda_c1_featureCounts_tumor,file=paste(output,"fda_c1_featureCounts_tumor.csv",sep=""),quote=F)
write.csv(fda_c1_featureCounts_normal,file=paste(output,"fda_c1_featureCounts_normal.csv",sep=""),quote=F)



## waferGen run1 ##

# featureCounts

counts_hcc1395 <- read.table("/genomics/1_Projects/FDA_QC/Results/WaferGen/First_run/Star_featureCounts/Ensembl_GRCh38/bam_subsampling/100K/featureCounts/counts_hcc1395.txt",skip=1,header=T,sep="\t")
wafer_run1_featureCounts_tumor <- get_summary2(exprs_gene_tumor,counts_hcc1395)

counts_hcc1395bl <- read.table("/genomics/1_Projects/FDA_QC/Results/WaferGen/First_run/Star_featureCounts/Ensembl_GRCh38/bam_subsampling/100K/featureCounts/counts_hcc1395bl.txt",skip=1,header=T,sep="\t")
wafer_run1_featureCounts_normal <- get_summary2(exprs_gene_normal,counts_hcc1395bl)

colnames(wafer_run1_featureCounts_tumor) <- c("HCC1395 expressed cells","HCC1395 cell percentage")
colnames(wafer_run1_featureCounts_normal) <- c("HCC1395BL expressed cells","HCC1395BL cell percentage")
rownames(wafer_run1_featureCounts_tumor) <- exprs_gene_tumor
rownames(wafer_run1_featureCounts_normal) <- exprs_gene_normal

write.csv(wafer_run1_featureCounts_tumor,file=paste(output,"wafer_run1_featureCounts_tumor.csv",sep=""),quote=F)
write.csv(wafer_run1_featureCounts_normal,file=paste(output,"wafer_run1_featureCounts_normal.csv",sep=""),quote=F)




## waferGen run2 ##

# featureCounts

counts_hcc1395 <- read.table("/genomics/1_Projects/FDA_QC/Results/WaferGen/Second_run/Star_featureCounts/Ensembl_GRCh38/bam_subsampling/100K/featureCounts/counts_hcc1395.txt",skip=1,header=T,sep="\t")
wafer_run2_featureCounts_tumor <- get_summary2(exprs_gene_tumor,counts_hcc1395)

counts_hcc1395bl <- read.table("/genomics/1_Projects/FDA_QC/Results/WaferGen/Second_run/Star_featureCounts/Ensembl_GRCh38/bam_subsampling/100K/featureCounts/counts_hcc1395bl.txt",skip=1,header=T,sep="\t")
wafer_run2_featureCounts_normal <- get_summary2(exprs_gene_normal,counts_hcc1395bl)

colnames(wafer_run2_featureCounts_tumor) <- c("HCC1395 expressed cells","HCC1395 cell percentage")
colnames(wafer_run2_featureCounts_normal) <- c("HCC1395BL expressed cells","HCC1395BL cell percentage")
rownames(wafer_run2_featureCounts_tumor) <- exprs_gene_tumor
rownames(wafer_run2_featureCounts_normal) <- exprs_gene_normal

write.csv(wafer_run2_featureCounts_tumor,file=paste(output,"wafer_run2_featureCounts_tumor.csv",sep=""),quote=F)
write.csv(wafer_run2_featureCounts_normal,file=paste(output,"wafer_run2_featureCounts_normal.csv",sep=""),quote=F)





t_perct <- cbind(llu_10x_zumi_tumor[,2],nci_10x_zumi_tumor_standard[,2],nci_10x_zumi_tumor_short[,2],llu_c1_featureCounts_tumor[,2],fda_c1_featureCounts_tumor[,2],wafer_run1_featureCounts_tumor[,2],wafer_run2_featureCounts_tumor[,2])

colnames(t_perct) <- c("llu_10x_zumi","nci_10x_zumis_standard","nci_10x_zumis_short","llu_c1_featureCounts","fda_c1_featureCounts", "wafer_run1_featureCounts","wafer_run2_featureCounts")
rownames(t_perct) <- exprs_gene_tumor

write.csv(t_perct,file=paste(output,"tumor_perct_summary.csv",sep=""),quote=F)



n_perct <- cbind(llu_10x_zumi_normal[,2],nci_10x_zumi_normal_standard[,2],nci_10x_zumi_normal_short[,2],llu_c1_featureCounts_normal[,2],fda_c1_featureCounts_normal[,2],wafer_run1_featureCounts_normal[,2],wafer_run2_featureCounts_normal[,2])

colnames(n_perct) <- c("llu_10x_zumi","nci_10x_zumis_standard","nci_10x_zumis_short","llu_c1_featureCounts","fda_c1_featureCounts", "wafer_run1_featureCounts","wafer_run2_featureCounts")
rownames(n_perct) <- exprs_gene_normal

write.csv(n_perct,file=paste(output,"normal_perct_summary.csv",sep=""),quote=F)





## generate corrleation plot ##


library(corrplot)
library(Cairo)

output <- "../manuscript_code/4_Batch_correction/Figure_5a-f/data/"

t_perct <- read.csv(paste(output,"tumor_perct_summary.csv",sep=""),row.names=1)
n_perct <- read.csv(paste(output,"normal_perct_summary.csv",sep=""),row.names=1)

colnames(t_perct) <- colnames(n_perct) <- c("10x_LLU","10x_NCI","10x_NCI_M","C1_FDA_HT","C1_LLU","WaferGen_PE","WaferGen_SE")


CairoJPEG(filename=paste(output,"tumor_perct_correlation_all.jpeg",sep=""),width=1000, height=1000,quality=1000)

corrplot(cor(t_perct), method="number", order="hclust",addrect=2, number.cex=3,tl.cex=3)

dev.off()


CairoJPEG(filename=paste(output,"tumor_perct_correlation_high.jpeg",sep=""),width=1000, height=1000,quality=1000)

corrplot(cor(t_perct[1:500,]), method="number", order="hclust",addrect=2, number.cex=3,tl.cex=3)

dev.off()


CairoJPEG(filename=paste(output,"tumor_perct_correlation_medium.jpeg",sep=""),width=1000, height=1000,quality=1000)

corrplot(cor(t_perct[501:1000,]), method="number", order="hclust",addrect=2, number.cex=3,tl.cex=3)

dev.off()


CairoJPEG(filename=paste(output,"tumor_perct_correlation_low.jpeg",sep=""),width=1000, height=1000,quality=1000)

corrplot(cor(t_perct[1001:1500,]), method="number", order="hclust",addrect=2, number.cex=3,tl.cex=3)

dev.off()


## normal ##

CairoJPEG(filename=paste(output,"normal_perct_correlation_all.jpeg",sep=""),width=1000, height=1000,quality=1000)

corrplot(cor(n_perct), method="number", order="hclust",addrect=2, number.cex=3,tl.cex=3)

dev.off()


CairoJPEG(filename=paste(output,"normal_perct_correlation_high.jpeg",sep=""),width=1000, height=1000,quality=1000)

corrplot(cor(n_perct[1:500,]), method="number", order="hclust",addrect=2, number.cex=3,tl.cex=3)

dev.off()


CairoJPEG(filename=paste(output,"normal_perct_correlation_medium.jpeg",sep=""),width=1000, height=1000,quality=1000)

corrplot(cor(n_perct[501:1000,]), method="number", order="hclust",addrect=2, number.cex=3,tl.cex=3)

dev.off()


CairoJPEG(filename=paste(output,"normal_perct_correlation_low.jpeg",sep=""),width=1000, height=1000,quality=1000)

corrplot(cor(n_perct[1001:1500,]), method="number", order="hclust",addrect=2, number.cex=3,tl.cex=3)

dev.off()



















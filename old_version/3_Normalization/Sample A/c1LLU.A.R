##module load R/3.5.0
library(Seurat)
library(scran)
library(scater)
library(ggplot2)
library(gplots)
#library(ggbiplot)
library(cluster)
library(Linnorm)
vloc<-"/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/tumor_A/results/"
dataid <- "c1llu"
#c1llu <- read.table("/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/WaferGen/Seq2/10K/counts_hcc1395.txt",header=T,stringsAsFactors = FALSE,sep="\t", quote="")
#c1llu2 <- read.table("/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/WaferGen/Seq2/100K/counts_hcc1395.txt",header=T,stringsAsFactors = FALSE,sep="\t", quote="")
#c1llu <- read.table("/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/WaferGen/Seq1/10K/counts_hcc1395.txt",header=T,stringsAsFactors = FALSE,sep="\t", quote="")
#c1llu2 <- read.table("/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/WaferGen/Seq1/100K/counts_hcc1395.txt",header=T,stringsAsFactors = FALSE,sep="\t", quote="")

#c1llu <- read.table("/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/C1_FDA_HT/50K/counts_hcc1395.txt",header=T,stringsAsFactors = FALSE,sep="\t", quote="")
#c1llu2 <- read.table("/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/C1_FDA_HT/10K/counts_hcc1395.txt",header=T,stringsAsFactors = FALSE,sep="\t", quote="")
c1llu <- read.table("/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/C1_LLU/100K/counts_hcc1395.txt",header=T,stringsAsFactors = FALSE,sep="\t", quote="")
c1llu2 <- read.table("/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/C1_LLU/10K/counts_hcc1395.txt",header=T,stringsAsFactors = FALSE,sep="\t", quote="")

sum(sort(colnames(c1llu)) == sort(colnames(c1llu2)))

c1llu<-c1llu[,-(1:6)]
c1llu2<-c1llu2[,-(1:6)]

meta1<- colnames(c1llu)
meta1<-data.frame(names=meta1,id=unlist(strsplit(meta1,"\\."))[seq(13,length(unlist(strsplit(meta1,"\\."))),15)],stringsAsFactors = FALSE)

meta2<- colnames(c1llu2)
meta2<-data.frame(names=meta2,id=unlist(strsplit(meta2,"\\."))[seq(13,length(unlist(strsplit(meta2,"\\."))),15)],stringsAsFactors = FALSE)
ncells<-sum(meta1$id %in% meta2$id) 

sum(meta1$id == meta2$id)  ## 66

rawdata<-as.matrix(x = cbind(c1llu,c1llu2[rownames(c1llu),]))


#sce <- newSCESet(countData=zz)
rownames(rawdata) <- as.character(rownames(rawdata))

#### cpm without log
calc_cpm <-
  function (expr_mat) 
  {
    norm_factor <- colSums(expr_mat)
    return(t(t(expr_mat)/norm_factor)) * 10^6
  }
normdata <- calc_cpm(rawdata)
data1<-normdata[,meta1$names]
data2<-normdata[,meta2$names]
dist.cpmnl<-dist(t(cbind(data1,data2)))
clusterd <- rep(1:ncells,2)
#sil_cpmnl.wafer2<-silhouette(clusterd,dist.cpmnl)

assign(paste("sil_cpmnl.",dataid,sep=""),silhouette(clusterd,dist.cpmnl))
out<- paste(vloc,"sil_cpmnl.",dataid,".RData",sep="")
save(list=paste("sil_cpmnl.",dataid,sep=""),file=out)

### uq
### UQ
calc_uq <-
  function (expr_mat) 
  {
    UQ <- function(x) {
      quantile(x[x > 0], 0.75)
    }
    uq <- unlist(apply(expr_mat, 2, UQ))
    norm_factor <- uq/median(uq)
    return(t(t(expr_mat)/norm_factor))
  }

normdata <- calc_uq(rawdata)
data1<-normdata[,meta1$names]
data2<-normdata[,meta2$names]
dist.uq<-dist(t(cbind(data1,data2)))
assign(paste("sil_uq.",dataid,sep=""),silhouette(clusterd,dist.uq))
out<- paste(vloc,"sil_uq.",dataid,".RData",sep="")
save(list=paste("sil_uq.",dataid,sep=""),file=out)


### other normalization methods
## cpm

pbmc <- CreateSeuratObject(raw.data = rawdata)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
normdata<-as.matrix(x = pbmc@data)
data1<-normdata[,meta1$names]
data2<-normdata[,meta2$names]
dist.cpm<-dist(t(cbind(data1,data2)))

assign(paste("sil_cpm.",dataid,sep=""),silhouette(clusterd,dist.cpm))
out<- paste(vloc,"sil_cpm.",dataid,".RData",sep="")
save(list=paste("sil_cpm.",dataid,sep=""),file=out)


### #ngenes
pbmc2 <- ScaleData(object = pbmc, vars.to.regress = "nGene")
scaledata<-as.matrix(x = pbmc2@scale.data)
data1<-scaledata[,meta1$names]
data2<-scaledata[,meta2$names]
dist.ngene<-dist(t(cbind(data1,data2)))

assign(paste("sil_ngene.",dataid,sep=""),silhouette(clusterd,dist.ngene))
out<- paste(vloc,"sil_ngene.",dataid,".RData",sep="")
save(list=paste("sil_ngene.",dataid,sep=""),file=out)

#sil_ngene.wafer2<-silhouette(clusterd,dist.ngene)
#save(sil_ngene.wafer2,file="/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_ngene.wafer2.RData")

############# regression on seq number
seq<-rep(1:2,each=ncells)
names(seq) <- colnames(rawdata)
pbmc <- AddMetaData(object = pbmc, metadata = seq, col.name = "seqnumber")
pbmc <- ScaleData(object = pbmc, vars.to.regress = "seqnumber")
scaledata2<-as.matrix(x = pbmc@scale.data)

data1<-scaledata2[,meta1$names]
data2<-scaledata2[,meta2$names]
dist.seq<-dist(t(cbind(data1,data2)))

assign(paste("sil_seq.",dataid,sep=""),silhouette(clusterd,dist.seq))
out<- paste(vloc,"sil_seq.",dataid,".RData",sep="")
save(list=paste("sil_seq.",dataid,sep=""),file=out)

#sil_seq.wafer2<-silhouette(clusterd,dist.seq)
#save(sil_seq.wafer2,file="/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_seq.wafer2.RData")

##END

### raw data

data1<-rawdata[,meta1$names]
data2<-rawdata[,meta2$names]
dist.raw<-dist(t(cbind(data1,data2)))

assign(paste("sil_raw.",dataid,sep=""),silhouette(clusterd,dist.raw))
out<- paste(vloc,"sil_raw.",dataid,".RData",sep="")
save(list=paste("sil_raw.",dataid,sep=""),file=out)

#sil_raw.wafer2<-silhouette(clusterd,dist.raw)
#save(sil_raw.wafer2,file="/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_raw.wafer2.RData")
### DESeq

calc_sf <-
  function (expr_mat) 
  {
    geomeans <- exp(rowMeans(log(expr_mat)))
    SF <- function(cnts) {
      median((cnts/geomeans)[(is.finite(geomeans) & geomeans > 
                                0)])
    }
    norm_factor <- apply(expr_mat, 2, SF)
    return(t(t(expr_mat)/norm_factor))
  }
normdata <- calc_sf(rawdata)
data1<-normdata[,meta1$names]
data2<-normdata[,meta2$names]
dist.deseq<-dist(t(cbind(data1,data2)))

assign(paste("sil_deseq.",dataid,sep=""),silhouette(clusterd,dist.deseq))
out<- paste(vloc,"sil_deseq.",dataid,".RData",sep="")
save(list=paste("sil_deseq.",dataid,sep=""),file=out)

#sil_deseq.wafer2<-silhouette(clusterd,dist.deseq)
#save(sil_deseq.wafer2,file="/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_deseq.wafer2.RData")


###TMM
## Apply TMM normalisation taking into account all genes
sce <- SingleCellExperiment(list(counts=rawdata))
sce <- normaliseExprs(sce, method = "TMM",return_log = F)
scedata<-normcounts(sce)

data1<-scedata[,meta1$names]
data2<-scedata[,meta2$names]

dist.tmm<-dist(t(cbind(data1,data2)))
clusterd <- rep(1:ncells,2)

assign(paste("sil_tmm.",dataid,sep=""),silhouette(clusterd,dist.tmm))
out<- paste(vloc,"sil_tmm.",dataid,".RData",sep="")
save(list=paste("sil_tmm.",dataid,sep=""),file=out)
#sil_tmm.wafer2<-silhouette(clusterd,dist.tmm)
#save(sil_tmm.wafer2,file="/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_tmm.wafer2.RData")

###SCRAN
sce <- SingleCellExperiment(list(counts=rawdata))
qclust <- quickCluster(sce, min.size = 50)
sce <- computeSumFactors(sce, sizes = 50,cluster=qclust)
#sce <- computeSumFactors(sce, sizes = 1000, clusters = qclust)#,positive=T)


sce <- normalize(sce)
scedata<-exprs(sce)


data1<-scedata[,meta1$names]
data3<-scedata[,meta2$names]
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
eud_sce = NULL
for(i in 1:dim(data1)[2]) {
  eud_sce <- c(eud_sce,euc.dist(data1[,i],data3[,i]))
}
dist.sce<-dist(t(cbind(data1,data3)))
clusterd <- rep(1:ncells,2)

assign(paste("sil_scran.",dataid,sep=""),silhouette(clusterd,dist.sce))
out<- paste(vloc,"sil_scran.",dataid,".RData",sep="")
save(list=paste("sil_scran.",dataid,sep=""),file=out)

#sil_scran.wafer2<-silhouette(clusterd,dist.sce)
#save(sil_scran.wafer2,file="/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_scran.wafer2.RData")

##Linnorm
normdata <- Linnorm(rawdata)
data1<-normdata[,meta1$names]
data2<-normdata[,meta2$names]
dist.linnorm<-dist(t(cbind(data1,data2)))

assign(paste("sil_linnorm.",dataid,sep=""),silhouette(clusterd,dist.linnorm))
out<- paste(vloc,"sil_linnorm.",dataid,".RData",sep="")
save(list=paste("sil_linnorm.",dataid,sep=""),file=out)

### linnorm_norm
normdata2 <- Linnorm.Norm(rawdata)
data1<-normdata2[,meta1$names]
data2<-normdata2[,meta2$names]
dist.linnorm2<-dist(t(cbind(data1,data2)))

assign(paste("sil_linnorm2.",dataid,sep=""),silhouette(clusterd,dist.linnorm2))
out<- paste(vloc,"sil_linnorm2.",dataid,".RData",sep="")
save(list=paste("sil_linnorm2.",dataid,sep=""),file=out)


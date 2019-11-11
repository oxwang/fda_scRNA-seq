##module load R/3.5.0

#library(Seurat)
library(scran)
library(scater)
#library(ggplot2)
#library(gplots)
#library(ggbiplot)
library(cluster)

c1llu <- read.table("/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/C1_LLU/100K/counts_hcc1395bl.txt",header=T,stringsAsFactors = FALSE,sep="\t", quote="")
c1llu2 <- read.table("/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/C1_LLU/10K/counts_hcc1395bl.txt",header=T,stringsAsFactors = FALSE,sep="\t", quote="")

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


sce <- SingleCellExperiment(list(counts=rawdata))
qclust <- quickCluster(sce, min.size = 20)
#sce <- computeSumFactors(sce, cluster=qclust)
sce <- computeSumFactors(sce, sizes = 20, clusters = qclust)#,positive=T)


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
sil_scran.c1llu<-silhouette(clusterd,dist.sce)
save(sil_scran.c1llu,file="/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_scran.c1llu.RData")

### raw data

data1<-rawdata[,meta1$names]
data2<-rawdata[,meta2$names]
dist.raw<-dist(t(cbind(data1,data2)))
#save(dist.raw,file="/data/tri_data/primary/2018_03_10XGenomics/biyx/FDA/LLU_10X/results/seurat/dist_raw.RData")
sil_raw.c1llu<-silhouette(clusterd,dist.raw)
save(sil_raw.c1llu,file="/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_raw.c1llu.RData")
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
sil_deseq.c1llu<-silhouette(clusterd,dist.deseq)
save(sil_deseq.c1llu,file="/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_deseq.c1llu.RData")


###TMM
## Apply TMM normalisation taking into account all genes
sce <- SingleCellExperiment(list(counts=rawdata))
sce <- normaliseExprs(sce, method = "TMM",return_log = F)
scedata<-normcounts(sce)

data1<-scedata[,meta1$names]
data2<-scedata[,meta2$names]

dist.tmm<-dist(t(cbind(data1,data2)))
clusterd <- rep(1:ncells,2)
sil_tmm.c1llu<-silhouette(clusterd,dist.tmm)
save(sil_tmm.c1llu,file="/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_tmm.c1llu.RData")


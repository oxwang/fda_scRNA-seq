##module load R/3.5.0

#library(Seurat)
library(scran)
library(scater)
#library(ggplot2)
#library(gplots)
#library(ggbiplot)
library(cluster)

#load("/data/tri_data/primary/2018_03_10XGenomics/biyx/FDA/subsampleing/C1_LLU/100K/counts_hcc1395.txt")

#filepath <- "/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/10X_LLU/"
filepath <- "/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/10X_NCI/"

filename <- list.files(filepath,pattern = ".dgecounts.rds",recursive=T)
counts <- readRDS(paste(filepath,filename[1],sep="")) ##
counts <- counts$exons$downsampled

#Then the counts data include downsampling counts at different read depths, you need to select 10K and 100K umi counts data

c1llu <- as.matrix(counts[[2]][[2]])
c1llu2  <- as.matrix(counts[[5]][[2]])


meta1<- colnames(c1llu)
meta1<-data.frame(id=meta1,names=paste(meta1,"1",sep="_"),stringsAsFactors = FALSE)

meta2<- colnames(c1llu2)
meta2<-data.frame(id=meta2,names=paste(meta2,"2",sep="_"),stringsAsFactors = FALSE)
sum(meta1$id %in% meta2$id)  ## 2068 ?

meta1<-meta1[meta1$id %in% meta2$id,]
sum(meta1$id == meta2$id)  ## 2068

colnames(c1llu) <- paste(colnames(c1llu),"1",sep="_")
colnames(c1llu2) <- paste(colnames(c1llu2),"2",sep="_")

c1llu <- c1llu[,meta1$names]
sum(rownames(c1llu) %in% rownames(c1llu2)) ## 19000 ?

zz <- merge(c1llu, c1llu2, by = "row.names", all = TRUE)
temp<-zz[,1]
zz<-as.matrix(zz[,-1])
rownames(zz) <- temp
zz[is.na(zz)] <- 0 



#sce <- newSCESet(countData=zz)
rownames(zz) <- as.character(rownames(zz))
sce <- SingleCellExperiment(list(counts=zz))
qclust <- quickCluster(sce)#, min.size = 100)
sce <- computeSumFactors(sce, cluster=qclust)
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
#dist.sce<-dist(t(cbind(data1,data3)))
ncells<-dim(data1)[2]

clusterd <- rep(1:ncells,2)
#sil_scran.10xnciM<-silhouette(clusterd,dist.sce)
#save(sil_scran.10xnciM,file="/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_scran.10X.NCI.M.RData")

### raw data
#rawdata<-as.matrix(x = pbmc@raw.data)[,rownames(meta3)]
rawdata<-zz
data1<-rawdata[,meta1$names]
data2<-rawdata[,meta2$names]
#dist.raw<-dist(t(cbind(data1,data2)))
#sil_raw.10xnciM<-silhouette(clusterd,dist.raw)
#save(sil_raw.10xnciM,file="/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_raw.10X.NCI.M.RData")


###TMM
## Apply TMM normalisation taking into account all genes
#sce <- SingleCellExperiment(list(counts=zz))
#sce <- normaliseExprs(sce, method = "TMM",return_log = F)
#scedata<-normcounts(sce)

#data1<-scedata[,meta1$names]
#data2<-scedata[,meta2$names]

#dist.tmm<-dist(t(cbind(data1,data2)))
#sil_tmm.10xnciM<-silhouette(clusterd,dist.tmm)
#save(sil_tmm.10xnciM,file="/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_tmm.10X.NCI.M.RData")

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
sil_cpmnl.10xnciM<-silhouette(clusterd,dist.cpmnl)
save(sil_cpmnl.10xnciM,file="/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_cpmnl.10X.NCI.M.RData")

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
sil_uq.10xnciM<-silhouette(clusterd,dist.uq)
save(sil_uq.10xnciM,file="/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_uq.10X.NCI.M.RData")

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
sil_deseq.10xnciM<-silhouette(clusterd,dist.deseq)
save(sil_deseq.10xnciM,file="/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_deseq.10X.NCI.M.RData")


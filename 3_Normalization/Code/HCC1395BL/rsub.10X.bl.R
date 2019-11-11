##module load R/3.5.0

#library(Seurat)
library(scran)
library(scater)
#library(ggplot2)
#library(gplots)
#library(ggbiplot)
library(cluster)

filepath <- "/.../3_Normalization/Data/HCC1395BL/10X_LLU/"

filename <- list.files(filepath,pattern = ".dgecounts.rds",recursive=T)
counts <- readRDS(paste(filepath,filename[1],sep=""))
counts <- counts$exons$downsampled

#Then the counts data include downsampling counts at different read depths, you need to select 10K and 100K umi counts data

load(paste0(filepath,'counts_10k.Rdata'))
load(paste0(filepath,'counts_100k.Rdata'))

c1llu <- counts_10k
c1llu2 <- counts_100k

meta1<- colnames(c1llu)
meta1<-data.frame(id=meta1,names=paste(meta1,"1",sep="_"),stringsAsFactors = FALSE)

meta2<- colnames(c1llu2)
meta2<-data.frame(id=meta2,names=paste(meta2,"2",sep="_"),stringsAsFactors = FALSE)
sum(meta1$id %in% meta2$id)  ## 885

meta1<-meta1[meta1$id %in% meta2$id,]
sum(meta1$id == meta2$id)  ## 885

colnames(c1llu) <- paste(colnames(c1llu),"1",sep="_")
colnames(c1llu2) <- paste(colnames(c1llu2),"2",sep="_")

c1llu <- c1llu[,meta1$names]
sum(rownames(c1llu) %in% rownames(c1llu2)) ## 16894

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
dist.sce<-dist(t(cbind(data1,data3)))

#save(dist.sce,file="/data/tri_data/primary/2018_03_10XGenomics/biyx/FDA/LLU_10X/results/seurat/dist_sce.RData")
#meta_t <- meta3[order(meta3$cell),]
clusterd <- rep(1:885,2)
sil_scran<-silhouette(clusterd,dist.sce)
save(sil_scran,file="/.../3_Normalization/Results/HCC1395BL/sil_scran.10X_LLU.RData")

### raw data
#rawdata<-as.matrix(x = pbmc@raw.data)[,rownames(meta3)]
rawdata<-zz
data1<-rawdata[,meta1$names]
data2<-rawdata[,meta2$names]
dist.raw<-dist(t(cbind(data1,data2)))
#save(dist.raw,file="/data/tri_data/primary/2018_03_10XGenomics/biyx/FDA/LLU_10X/results/seurat/dist_raw.RData")
sil_raw<-silhouette(clusterd,dist.raw)
save(sil_raw,file="/.../3_Normalization/Results/HCC1395BL/sil_raw.10X_LLU.RData")
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
sil_deseq<-silhouette(clusterd,dist.deseq)
save(sil_deseq,file="/.../3_Normalization/Results/HCC1395BL/sil_deseq.10X_LLU.RData")


###TMM
## Apply TMM normalisation taking into account all genes
sce <- SingleCellExperiment(list(counts=zz))
sce <- normaliseExprs(sce, method = "TMM",return_log = F)
scedata<-normcounts(sce)

data1<-scedata[,meta1$names]
data2<-scedata[,meta2$names]

dist.tmm<-dist(t(cbind(data1,data2)))
sil_tmm<-silhouette(clusterd,dist.tmm)
save(sil_tmm,file="/.../3_Normalization/Results/HCC1395BL/sil_tmm.10X_LLU.RData")


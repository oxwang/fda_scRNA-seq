##module load R/3.5.0
library(Seurat)
library(scran)
library(scater)
library(ggplot2)
library(gplots)
#library(ggbiplot)
library(cluster)
library(Linnorm)
vloc<-"/.../3_Normalization/Results/HCC1395/"
dataid <- "10XLLU"
filepath <- "/.../3_Normalization/Data/HCC1395/10X_LLU/"

#Then the counts data include downsampling counts at different read depths, you need to select 10K and 100K umi counts data

load(paste0(filepath,'counts_10k.Rdata'))
load(paste0(filepath,'counts_100k.Rdata'))

c1llu <- counts_10k
c1llu2 <- counts_100k


meta1<- colnames(c1llu)
meta1<-data.frame(id=meta1,names=paste(meta1,"1",sep="_"),stringsAsFactors = FALSE)

meta2<- colnames(c1llu2)
meta2<-data.frame(id=meta2,names=paste(meta2,"2",sep="_"),stringsAsFactors = FALSE)
sum(meta1$id %in% meta2$id)
ncells<-sum(meta1$id %in% meta2$id) 

meta1<-meta1[meta1$id %in% meta2$id,]

colnames(c1llu) <- paste(colnames(c1llu),"1",sep="_")
colnames(c1llu2) <- paste(colnames(c1llu2),"2",sep="_")

c1llu <- c1llu[,meta1$names]
sum(rownames(c1llu) %in% rownames(c1llu2))

zz <- merge(c1llu, c1llu2, by = "row.names", all = TRUE)
temp<-zz[,1]
zz<-as.matrix(zz[,-1])
rownames(zz) <- temp
zz[is.na(zz)] <- 0 



#sce <- newSCESet(countData=zz)
rownames(zz) <- as.character(rownames(zz))
rawdata<-zz

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


##END

### raw data

data1<-rawdata[,meta1$names]
data2<-rawdata[,meta2$names]
dist.raw<-dist(t(cbind(data1,data2)))

assign(paste("sil_raw.",dataid,sep=""),silhouette(clusterd,dist.raw))
out<- paste(vloc,"sil_raw.",dataid,".RData",sep="")
save(list=paste("sil_raw.",dataid,sep=""),file=out)

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

###SCRAN

sce <- SingleCellExperiment(list(counts=rawdata))
qclust <- quickCluster(sce)#, min.size = 100)
sce <- computeSumFactors(sce, cluster=qclust)
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

assign(paste("sil_scran.",dataid,sep=""),silhouette(clusterd,dist.sce))
out<- paste(vloc,"sil_scran.",dataid,".RData",sep="")
save(list=paste("sil_scran.",dataid,sep=""),file=out)


##Linnorm
normdata <- Linnorm(rawdata)
data1<-normdata[,meta1$names]
data2<-normdata[,meta2$names]
dist.linnorm<-dist(t(cbind(data1,data2)))
clusterd <- rep(1:ncells,2)
assign(paste("sil_linnorm.",dataid,sep=""),silhouette(clusterd,dist.linnorm))
out<- paste(vloc,"sil_linnorm.",dataid,".RData",sep="")
save(list=paste("sil_linnorm.",dataid,sep=""),file=out)

### linnorm_norm
normdata2 <- Linnorm.Norm(rawdata)
data1<-normdata2[,meta1$names]
data2<-normdata2[,meta2$names]
dist.linnorm2<-dist(t(cbind(data1,data2)))
clusterd <- rep(1:ncells,2)
assign(paste("sil_linnorm2.",dataid,sep=""),silhouette(clusterd,dist.linnorm2))
out<- paste(vloc,"sil_linnorm2.",dataid,".RData",sep="")
save(list=paste("sil_linnorm2.",dataid,sep=""),file=out)

## sctransform
pbmc <- CreateSeuratObject(counts = rawdata)
pbmc <- SCTransform(pbmc, verbose = FALSE)
normdata<-pbmc[["SCT"]]@data

data1<-as.matrix(normdata[,meta1$names])
data2<-as.matrix(normdata[,meta2$names])
dist.sctran<-dist(t(cbind(data1,data2)))

ncells<-dim(data1)[2]
clusterd <- rep(1:ncells,2)
assign(paste("sil_sctran.",dataid,sep=""),silhouette(clusterd,dist.sctran))
out<- paste(vloc,"sil_sctran.",dataid,".RData",sep="")
save(list=paste("sil_sctran.",dataid,sep=""),file=out)


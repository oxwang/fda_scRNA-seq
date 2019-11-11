##module load R/3.5.0
library(Seurat)
library(scran)
library(scater)
library(ggplot2)
library(gplots)
#library(ggbiplot)
library(cluster)
library(Linnorm)
vloc<-"/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/"


read10X_N <- function(filename) {
  #filename <- list.files(filepath,pattern = ".dgecounts.rds",recursive=T)
  
  counts <- readRDS(filename)
  #counts <- readRDS(paste(filepath,filename[1],sep=""))
  counts <- counts$exons$downsampled
  c1llu <- as.matrix(counts[[2]][[2]])
  c1llu2  <- as.matrix(counts[[5]][[2]])
  
  
  meta1<- colnames(c1llu)
  meta1<-data.frame(id=meta1,names=paste(meta1,"1",sep="_"),stringsAsFactors = FALSE)
  
  meta2<- colnames(c1llu2)
  meta2<-data.frame(id=meta2,names=paste(meta2,"2",sep="_"),stringsAsFactors = FALSE)
  sum(meta1$id %in% meta2$id)  ## 1041
  
  meta1<-meta1[meta1$id %in% meta2$id,]
  sum(meta1$id == meta2$id)  ## 1041
  
  colnames(c1llu) <- paste(colnames(c1llu),"1",sep="_")
  colnames(c1llu2) <- paste(colnames(c1llu2),"2",sep="_")
  
  c1llu <- c1llu[,meta1$names]
  sum(rownames(c1llu) %in% rownames(c1llu2)) ## 19096
  
  zz <- merge(c1llu, c1llu2, by = "row.names", all = TRUE)
  temp<-zz[,1]
  zz<-as.matrix(zz[,-1])
  rownames(zz) <- temp
  zz[is.na(zz)] <- 0 
  
  
  
  #sce <- newSCESet(countData=zz)
  rownames(zz) <- as.character(rownames(zz))
  newList <- list("rawdata" = zz, "meta1" = meta1,"meta2"=meta2)
  return(newList)
  #rawdata<-zz
}


dataid <- "10xnciM"
filepath <- "/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/10X_NCI/"
filename<-paste(filepath,"10X_NC_N.dgecounts.rds",sep="")
newList <- read10X_N(filename)
rawdata<-newList$rawdata
meta1<-newList$meta1
meta2<-newList$meta2

##Linnorm
normdata <- Linnorm(rawdata)
data1<-normdata[,meta1$names]
data2<-normdata[,meta2$names]
dist.linnorm<-dist(t(cbind(data1,data2)))

ncells<-dim(data1)[2]
clusterd <- rep(1:ncells,2)
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
###############
dataid <- "10xnciseq4"
filename2<-paste(filepath,"10X_NCI_N_S2_L001.dgecounts.rds",sep="")
newList <- read10X_N(filename2)
rawdata<-newList$rawdata
meta1<-newList$meta1
meta2<-newList$meta2

##Linnorm
normdata <- Linnorm(rawdata)
data1<-normdata[,meta1$names]
data2<-normdata[,meta2$names]
dist.linnorm<-dist(t(cbind(data1,data2)))

ncells<-dim(data1)[2]
clusterd <- rep(1:ncells,2)
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

###########10X LLU
filepath <- "/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/10X_LLU/"

dataid <- ""
filename3<-paste(filepath,"10X_LL_N.dgecounts.rds",sep="")
newList <- read10X_N(filename3)
rawdata<-newList$rawdata
meta1<-newList$meta1
meta2<-newList$meta2

##Linnorm
normdata <- Linnorm(rawdata)
data1<-normdata[,meta1$names]
data2<-normdata[,meta2$names]
dist.linnorm<-dist(t(cbind(data1,data2)))

ncells<-dim(data1)[2]
clusterd <- rep(1:ncells,2)
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
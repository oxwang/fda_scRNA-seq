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

readc1fda<-function(filename1,filename2) {
  c1llu <- read.table(filename1,header=T,stringsAsFactors = FALSE,sep="\t", quote="")
  c1llu2 <- read.table(filename2,header=T,stringsAsFactors = FALSE,sep="\t", quote="")
  rownames(c1llu) <- c1llu[,1]
  rownames(c1llu2) <- c1llu2[,1]
  #sum(sort(colnames(c1llu)) == sort(colnames(c1llu2)))
  
  c1llu<-c1llu[,-(1:6)]
  c1llu2<-c1llu2[,-(1:6)]
  
  meta1<- colnames(c1llu)
  meta1<-data.frame(names=meta1,id=unlist(strsplit(meta1,"\\."))[seq(14,length(unlist(strsplit(meta1,"\\."))),16)],stringsAsFactors = FALSE)
  
  meta2<- colnames(c1llu2)
  meta2<-data.frame(names=meta2,id=unlist(strsplit(meta2,"\\."))[seq(14,length(unlist(strsplit(meta2,"\\."))),16)],stringsAsFactors = FALSE)
  ncells<-sum(meta1$id %in% meta2$id)  ## 187
  
  sum(meta1$id %in% meta2$id)  ## 187
  
  meta2<-meta2[meta2$id %in% meta1$id,]
  sum(meta1$id == meta2$id)  ## 187
  c1llu2 <- c1llu2[,meta2$names]
  
  
  meta1<-meta1[meta1$id %in% meta2$id,]
  
  c1llu <- c1llu[,meta1$names]
  rawdata<-as.matrix(x = cbind(c1llu,c1llu2[rownames(c1llu),]))
  
  
  #sce <- newSCESet(countData=zz)
  rownames(rawdata) <- as.character(rownames(rawdata))
  newList <- list("rawdata" = rawdata, "meta1" = meta1,"meta2"=meta2)
  return(newList)
}


readothers_N <- function(filename1,filename2) {
c1llu <- read.table(filename1,header=T,stringsAsFactors = FALSE,sep="\t", quote="")
c1llu2 <- read.table(filename2,header=T,stringsAsFactors = FALSE,sep="\t", quote="")

sum(sort(colnames(c1llu)) == sort(colnames(c1llu2)))

c1llu<-c1llu[,-(1:6)]
c1llu2<-c1llu2[,-(1:6)]

meta1<- colnames(c1llu)
meta1<-data.frame(names=meta1,id=unlist(strsplit(meta1,"\\."))[seq(13,length(unlist(strsplit(meta1,"\\."))),15)],stringsAsFactors = FALSE)

meta2<- colnames(c1llu2)
meta2<-data.frame(names=meta2,id=unlist(strsplit(meta2,"\\."))[seq(13,length(unlist(strsplit(meta2,"\\."))),15)],stringsAsFactors = FALSE)
ncells<-sum(meta1$id %in% meta2$id) 

#sum(meta1$id == meta2$id)  ## 66
meta1<-meta1[meta1$id %in% meta2$id,]

c1llu <- c1llu[,meta1$names]

#rawdata<-as.matrix(x = cbind(c1llu,c1llu2[rownames(c1llu),]))
rawdata<-as.matrix(x = cbind(c1llu,c1llu2[rownames(c1llu),]))


#sce <- newSCESet(countData=zz)
rownames(rawdata) <- as.character(rownames(rawdata))
newList <- list("rawdata" = rawdata, "meta1" = meta1,"meta2"=meta2)
return(newList)
}



#c1llu <- read.table("/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/C1_LLU/100K/counts_hcc1395bl.txt",header=T,stringsAsFactors = FALSE,sep="\t", quote="")
#c1llu2 <- read.table("/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/C1_LLU/10K/counts_hcc1395bl.txt",header=T,stringsAsFactors = FALSE,sep="\t", quote="")

##C1LLU
f1<-"/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/C1_LLU/100K/counts_hcc1395bl.txt"
f2<-"/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/C1_LLU/10K/counts_hcc1395bl.txt"

###########10X LLU
#filepath <- "/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/10X_LLU/"

dataid <- "c1llu"
#filename3<-paste(filepath,"10X_LL_N.dgecounts.rds",sep="")
newList <- readothers_N(f1,f2)
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

### C1FDA
f1<-"/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/C1_FDA_HT/50K/counts_hcc1395bl.txt"
f2 <- "/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/C1_FDA_HT/10K/counts_hcc1395bl.txt"
dataid <- "c1fda"
newList <- readc1fda(f1,f2)
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


### wafer1
f1<-"/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/WaferGen/Seq1/10K/counts_hcc1395bl.txt"
f2<-"/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/WaferGen/Seq1/100K/counts_hcc1395bl.txt"
dataid <- "wafer"
newList <- readc1fda(f1,f2)
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

### wafer2
f1<-"/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/WaferGen/Seq2/10K/counts_hcc1395bl.txt"
f2<-"/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/WaferGen/Seq2/100K/counts_hcc1395bl.txt"
dataid <- "wafer2"
newList <- readc1fda(f1,f2)
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
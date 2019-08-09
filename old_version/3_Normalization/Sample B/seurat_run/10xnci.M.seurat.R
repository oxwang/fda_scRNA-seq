library(Seurat)
library(scran)
library(ggplot2)
library(gplots)
library(ggbiplot)
library(cluster)

#load("/data/tri_data/primary/2018_03_10XGenomics/biyx/FDA/subsampleing/C1_LLU/100K/counts_hcc1395.txt")

#filepath <- "/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/10X_LLU/"
filepath <- "/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/10X_NCI/"

filename <- list.files(filepath,pattern = ".dgecounts.rds",recursive=T)
counts <- readRDS(paste(filepath,filename[1],sep=""))
counts <- counts$exons$downsampled

#Then the counts data include downsampling counts at different read depths, you need to select 10K and 100K umi counts data

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
rawdata<-zz



### other normalization methods
## cpm

pbmc <- CreateSeuratObject(raw.data = rawdata)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
normdata<-as.matrix(x = pbmc@data)
data1<-normdata[,meta1$names]
data2<-normdata[,meta2$names]
dist.cpm<-dist(t(cbind(data1,data2)))
ncells<-dim(data1)[2]
clusterd <- rep(1:ncells,2)
sil_cpm.10xnciM<-silhouette(clusterd,dist.cpm)
save(sil_cpm.10xnciM,file="/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_cpm.10X.NCI.M.RData")

### #ngenes
pbmc2 <- ScaleData(object = pbmc, vars.to.regress = "nGene")
scaledata<-as.matrix(x = pbmc2@scale.data)
data1<-scaledata[,meta1$names]
data2<-scaledata[,meta2$names]
dist.ngene<-dist(t(cbind(data1,data2)))
sil_ngene.10xnciM<-silhouette(clusterd,dist.ngene)
save(sil_ngene.10xnciM,file="/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_ngene.10X.NCI.M.RData")

############# regression on seq number
seq<-rep(1:2,each=ncells)
names(seq) <- colnames(rawdata)
pbmc <- AddMetaData(object = pbmc, metadata = seq, col.name = "seqnumber")
pbmc <- ScaleData(object = pbmc, vars.to.regress = "seqnumber")
scaledata2<-as.matrix(x = pbmc@scale.data)

data1<-scaledata2[,meta1$names]
data2<-scaledata2[,meta2$names]
dist.seq<-dist(t(cbind(data1,data2)))
sil_seq.10xnciM<-silhouette(clusterd,dist.seq)
save(sil_seq.10xnciM,file="/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_seq.10X.NCI.M.RData")


### boxplot of silhoutte
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_scran.10X.NCI.M.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_deseq.10X.NCI.M.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_raw.10X.NCI.M.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_tmm.10X.NCI.M.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_cpm.10X.NCI.M.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_ngene.10X.NCI.M.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_seq.10X.NCI.M.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_uq.10X.NCI.M.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_cpmnl.10X.NCI.M.RData")


#n2<-ncells*2
n2<-dim(sil_scran.10xnciM)[1]
dplot <- data.frame(Sil= c(sil_scran.10xnciM[,3],sil_cpm.10xnciM[,3],sil_ngene.10xnciM[,3],sil_seq.10xnciM[,3],sil_raw.10xnciM[,3],sil_deseq.10xnciM[,3],sil_tmm.10xnciM[,3],sil_uq.10xnciM[,3],sil_cpmnl.10xnciM[,3]),Methods=c(rep("scran",n2),rep("lgcpm",n2),rep("lgcpm+Regress_ngene",n2),rep("lgcpm+Regress_seq",n2),rep("raw",n2),rep("deseq",n2),rep("tmm",n2),rep("UQ",n2),rep("cpm",n2)))


p10 <- ggplot(dplot, aes(x = Methods, y = Sil,color= factor(Methods),fill = factor(Methods))) +
  geom_boxplot(fill="white",outlier.colour = NA, 
               position = position_dodge(width=0.9),width=0.5) +
  geom_point(size = 1, position = position_jitterdodge()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Silhouette Values")
ggsave("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/boxplot_10X_NCI.M.pdf", p10, width=8, height=6, units="in")


aggregate(dplot[,1,drop=F], by=list(dplot$Methods), FUN=mean, na.rm=TRUE)
#Group.1        Sil
#1                 cpm  0.1820628
#2               deseq  0.1457264
#3               lgcpm  0.1162839
#4 lgcpm+Regress_ngene  0.1119155
#5   lgcpm+Regress_seq  0.1136009
#6                 raw -0.5835620
#7               scran  0.1024825
#8                 tmm -0.3761872
#9                  UQ -0.4890488
aggregate(dplot[,1,drop=F], by=list(dplot$Methods), FUN=sd, na.rm=TRUE)
#Group.1        Sil
#1                 cpm 0.14136773
#2               deseq 0.15633930
#3               lgcpm 0.12940748
#4 lgcpm+Regress_ngene 0.08045223
#5   lgcpm+Regress_seq 0.07933639
#6                 raw 0.22339332
#7               scran 0.12121149
#8                 tmm 0.13831566
#9                  UQ 0.16889726


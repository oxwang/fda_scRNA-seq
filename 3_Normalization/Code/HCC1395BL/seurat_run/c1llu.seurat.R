library(Seurat)
library(scran)
library(ggplot2)
library(gplots)
library(ggbiplot)
library(cluster)

c1llu <- read.table("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/C1_LLU/100K/counts_hcc1395bl.txt",header=T,stringsAsFactors = FALSE,sep="\t", quote="")
c1llu2 <- read.table("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/C1_LLU/10K/counts_hcc1395bl.txt",header=T,stringsAsFactors = FALSE,sep="\t", quote="")

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
sil_cpmnl.c1llu<-silhouette(clusterd,dist.cpmnl)
save(sil_cpmnl.c1llu,file="/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_cpmnl.c1llu.RData")

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
sil_uq.c1llu<-silhouette(clusterd,dist.uq)
save(sil_uq.c1llu,file="/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_uq.c1llu.RData")


### other normalization methods
## cpm

pbmc <- CreateSeuratObject(raw.data = rawdata)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
normdata<-as.matrix(x = pbmc@data)
data1<-normdata[,meta1$names]
data2<-normdata[,meta2$names]
dist.cpm<-dist(t(cbind(data1,data2)))
clusterd <- rep(1:ncells,2)
sil_cpm.c1llu<-silhouette(clusterd,dist.cpm)
save(sil_cpm.c1llu,file="/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_cpm.c1llu.RData")

### #ngenes
pbmc2 <- ScaleData(object = pbmc, vars.to.regress = "nGene")
scaledata<-as.matrix(x = pbmc2@scale.data)
data1<-scaledata[,meta1$names]
data2<-scaledata[,meta2$names]
dist.ngene<-dist(t(cbind(data1,data2)))
sil_ngene.c1llu<-silhouette(clusterd,dist.ngene)
save(sil_ngene.c1llu,file="/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_ngene.c1llu.RData")

############# regression on seq number
seq<-rep(1:2,each=ncells)
names(seq) <- colnames(rawdata)
pbmc <- AddMetaData(object = pbmc, metadata = seq, col.name = "seqnumber")
pbmc <- ScaleData(object = pbmc, vars.to.regress = "seqnumber")
scaledata2<-as.matrix(x = pbmc@scale.data)

data1<-scaledata2[,meta1$names]
data2<-scaledata2[,meta2$names]
dist.seq<-dist(t(cbind(data1,data2)))
sil_seq.c1llu<-silhouette(clusterd,dist.seq)
save(sil_seq.c1llu,file="/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_seq.c1llu.RData")

##END
### boxplot of silhoutte
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_scran.c1llu.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_deseq.c1llu.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_raw.c1llu.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_tmm.c1llu.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_cpm.c1llu.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_cpmnl.c1llu.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_ngene.c1llu.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_seq.c1llu.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_uq.c1llu.RData")








n2<-ncells*2
dplot <- data.frame(Sil= c(sil_scran.c1llu[,3],sil_cpm.c1llu[,3],sil_ngene.c1llu[,3],sil_seq.c1llu[,3],sil_raw.c1llu[,3],sil_deseq.c1llu[,3],sil_tmm.c1llu[,3],sil_cpmnl.c1llu[,3],sil_uq.c1llu[,3]),Methods=c(rep("scran",n2),rep("lgcpm",n2),rep("lgcpm+Regress_ngene",n2),rep("lgcpm+Regress_seq",n2),rep("raw",n2),rep("deseq",n2),rep("tmm",n2),rep("cpm",n2),rep("uq",n2)))


p10 <- ggplot(dplot, aes(x = Methods, y = Sil,color= factor(Methods),fill = factor(Methods))) +
  geom_boxplot(fill="white",outlier.colour = NA, 
               position = position_dodge(width=0.9),width=0.5) +
  geom_point(size = 1, position = position_jitterdodge()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Silhouette Values")
ggsave("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/boxplot_c1llu.pdf", p10, width=8, height=6, units="in")


aggregate(dplot[,1,drop=F], by=list(dplot$Methods), FUN=mean, na.rm=TRUE)
#Group.1        Sil
#1                 cpm  0.4445600
#2               deseq  0.4022707
#3               lgcpm  0.2888680
#4 lgcpm+Regress_ngene  0.2095556
#5   lgcpm+Regress_seq  0.2132174
#6                 raw -0.4571990
#7               scran  0.2680078
#8                 tmm -0.1669875
#9                  uq -0.2923825
aggregate(dplot[,1,drop=F], by=list(dplot$Methods), FUN=sd, na.rm=TRUE)
#Group.1        Sil
#1                 cpm 0.13676004
#2               deseq 0.15238644
#3               lgcpm 0.11298727
#4 lgcpm+Regress_ngene 0.07241235
#5   lgcpm+Regress_seq 0.07224878
#6                 raw 0.14697075
#7               scran 0.10423660
#8                 tmm 0.20940953
#9                  uq 0.22717262

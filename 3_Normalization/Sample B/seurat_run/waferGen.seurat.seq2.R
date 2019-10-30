library(Seurat)
library(scran)
library(ggplot2)
library(gplots)
library(ggbiplot)
library(cluster)

c1llu <- read.table("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/WaferGen/Seq2/10K/counts_hcc1395bl.txt",header=T,stringsAsFactors = FALSE,sep="\t", quote="")
c1llu2 <- read.table("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/WaferGen/Seq2/100K/counts_hcc1395bl.txt",header=T,stringsAsFactors = FALSE,sep="\t", quote="")
rownames(c1llu) <- c1llu[,1]
rownames(c1llu2) <- c1llu2[,1]

#sum(sort(colnames(c1llu)) == sort(colnames(c1llu2)))

c1llu<-c1llu[,-(1:6)]
c1llu2<-c1llu2[,-(1:6)]

meta1<- colnames(c1llu)
meta1<-data.frame(names=meta1,id=unlist(strsplit(meta1,"\\."))[seq(14,length(unlist(strsplit(meta1,"\\."))),16)],stringsAsFactors = FALSE)

meta2<- colnames(c1llu2)
meta2<-data.frame(names=meta2,id=unlist(strsplit(meta2,"\\."))[seq(14,length(unlist(strsplit(meta2,"\\."))),16)],stringsAsFactors = FALSE)
ncells<-sum(meta1$id %in% meta2$id)  ## 539

meta1<-meta1[meta1$id %in% meta2$id,]
sum(meta1$id == meta2$id)  ## 539
c1llu <- c1llu[,meta1$names]

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
sil_cpmnl.wafer2<-silhouette(clusterd,dist.cpmnl)
save(sil_cpmnl.wafer2,file="/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_cpmnl.wafer2.RData")

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
sil_uq.wafer2<-silhouette(clusterd,dist.uq)
save(sil_uq.wafer2,file="/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_uq.wafer2.RData")

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
sil_cpm.wafer2<-silhouette(clusterd,dist.cpm)
save(sil_cpm.wafer2,file="/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_cpm.wafer2.RData")

### #ngenes
pbmc2 <- ScaleData(object = pbmc, vars.to.regress = "nGene")
scaledata<-as.matrix(x = pbmc2@scale.data)
data1<-scaledata[,meta1$names]
data2<-scaledata[,meta2$names]
dist.ngene<-dist(t(cbind(data1,data2)))
sil_ngene.wafer2<-silhouette(clusterd,dist.ngene)
save(sil_ngene.wafer2,file="/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_ngene.wafer2.RData")

############# regression on seq number
seq<-rep(1:2,each=ncells)
names(seq) <- colnames(rawdata)
pbmc <- AddMetaData(object = pbmc, metadata = seq, col.name = "seqnumber")
pbmc <- ScaleData(object = pbmc, vars.to.regress = "seqnumber")
scaledata2<-as.matrix(x = pbmc@scale.data)

data1<-scaledata2[,meta1$names]
data2<-scaledata2[,meta2$names]
dist.seq<-dist(t(cbind(data1,data2)))
sil_seq.wafer2<-silhouette(clusterd,dist.seq)
save(sil_seq.wafer2,file="/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_seq.wafer2.RData")

##END
### boxplot of silhoutte
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_scran.wafer2.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_deseq.wafer2.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_raw.wafer2.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_tmm.wafer2.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_cpm.wafer2.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_cpmnl.wafer2.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_uq.wafer2.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_seq.wafer2.RData")
load("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/sil_ngene.wafer2.RData")





n2<-dim(sil_seq.wafer2)[1]
#n2<-ncells*2
dplot <- data.frame(Sil= c(sil_scran.wafer2[,3],sil_cpm.wafer2[,3],sil_ngene.wafer2[,3],sil_seq.wafer2[,3],sil_raw.wafer2[,3],sil_deseq.wafer2[,3],sil_tmm.wafer2[,3],sil_uq.wafer2[,3],sil_cpmnl.wafer2[,3]),Methods=c(rep("scran",n2),rep("lgcpm",n2),rep("lgcpm+Regress_ngene",n2),rep("lgcpm+Regress_seq",n2),rep("raw",n2),rep("deseq",n2),rep("tmm",n2),rep("uq",n2),rep("cpm",n2)))


p10 <- ggplot(dplot, aes(x = Methods, y = Sil,color= factor(Methods),fill = factor(Methods))) +
  geom_boxplot(fill="white",outlier.colour = NA, 
               position = position_dodge(width=0.9),width=0.5) +
  geom_point(size = 1, position = position_jitterdodge()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Silhouette Values")
ggsave("/data/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/boxplot_WaferGen_Seq2.pdf", p10, width=8, height=6, units="in")


aggregate(dplot[,1,drop=F], by=list(dplot$Methods), FUN=mean, na.rm=TRUE)
#Group.1         Sil
#1                 cpm  0.46723316
#2               deseq  0.40736604
#3               lgcpm  0.36869182
#4 lgcpm+Regress_ngene  0.29865119
#5   lgcpm+Regress_seq  0.30083147
#6                 raw -0.59198782
#7               scran  0.34368002
#8                 tmm -0.03894084
#9                  uq -0.26417171
aggregate(dplot[,1,drop=F], by=list(dplot$Methods), FUN=sd, na.rm=TRUE)
#Group.1        Sil
#1                 cpm 0.09837867
#2               deseq 0.15905060
#3               lgcpm 0.08633206
#4 lgcpm+Regress_ngene 0.08463284
#5   lgcpm+Regress_seq 0.08344729
#6                 raw 0.30916676
#7               scran 0.08658451
#8                 tmm 0.17555094
#9                  uq 0.24165638


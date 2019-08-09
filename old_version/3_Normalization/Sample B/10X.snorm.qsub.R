##module load R/3.5.0
library(scater)
require("SingleCellExperiment")




filepath <- "/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/10X_LLU/"

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

#### scnorm normalization
library(SCnorm)


#rawd <- rawd[,c(meta3_1$id,meta3_3$id)]
#condition <- c(meta3_1$seq,meta3_3$seq)
#pdf("/sc/wo/tri_data/primary/2018_03_10XGenomics/biyx/FDA/LLU_10X/results/counts_depths.pdf")
#countDeptEst <- plotCountDepth(Data=rawd,Conditions=condition,FilterCellProportion=0.1,NCores=3)
#dev.off()
rownames(zz) <- as.character(rownames(zz))

sce <- SingleCellExperiment(list(counts=zz))

condition_1 <- rep(1,1770)
## one condition
#pdf("/sc/wo/tri_data/primary/2018_03_10XGenomics/biyx/FDA/LLU_10X/results/normalize_one_condition.pdf")
#par(mfrow=c(2,2))
#DataNorm1<-SCnorm(Data=rawd,Conditions=condition_1,PrintProgressPlots=TRUE,FilterCellNum=10,NCores=3)
#dev.off()

#NormalizedData1 <- results(DataNorm1)
#save(NormalizedData1,file="/sc/wo/tri_data/primary/2018_03_10XGenomics/biyx/FDA/LLU_10X/results/NormalizedData1.RData")
### two condition

##10X
DataNorm <- SCnorm(Data=sce,Conditions=condition_1,PrintProgressPlots=TRUE,FilterCellNum=10,PropToUse=0.1,Thresh=.1,ditherCounts=T)
save(DataNorm,file="/sc/wo/tri_data/primary/2018_03_10XGenomics_SLE_whole_blood_MiSeq/biyx/FDA/subsampleing/normal_BL/results/DataNorm.10X.RData")

#####

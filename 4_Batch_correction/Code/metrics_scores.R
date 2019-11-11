## calculate modified alignment score, kBET score, and silhouette score ##


#########################
## get alignment score ##
#########################

calc_align_score <- function(embedding,batch_cell,cell_perct)
{
	embedding <- embedding
	batch_cell <- batch_cell
	cell_perct <- cell_perct
	w <- batch_cell/sum(batch_cell)

	temp <- cumsum(batch_cell)
	batch_index <- cbind(c(1,temp[-length(temp)]+1),temp)

	embedding_dist <- as.matrix(dist(embedding))
	diag(embedding_dist) <- 1e5
	k <- round(nrow(embedding)*cell_perct)
	embedding_order <- t(apply(embedding_dist,1,order))

	k_mean <- c()
	for(i in 1:length(batch_cell))
	{
		temp <- embedding_order[batch_index[i,1]:batch_index[i,2],1:k]
		temp2 <- temp <= batch_index[i,2]
		temp3 <- temp >= batch_index[i,1]
		temp <- temp2+temp3
		k_mean <- c(k_mean,mean(apply(temp == 2,1,sum)))
	}
	align_score <- 1-abs(k_mean-k*w)/(k-k*w)
	align_score <- c(align_score,sum(w*align_score))

	return(align_score)
}


cell_counts <- function(da)
{
	a <- rownames(da)
	a <- sub('wta432_','wta432-',a)
	counts <- c()
	for(j in 1:5)
	{counts <- c(counts,sum(!is.na(grep(paste0('_',j),a))))}
	return(counts)
}




########################################################################
## get alighment score for each cell type in scenario 1 or scenario 4 ##
########################################################################

hvg_num <- 2000
input <- '/.../4_Batch_correction/Results/scenario_1/'
#input <- '/.../4_Batch_correction/Results/scenario_4/'

batch <- c("10x_LLU_A", "10x_NCI_A", "10x_NCI_M_A", "10x_LLU_B", "10x_NCI_B", "10x_NCI_M_B", "10x_Mix10_LLU", "10x_Mix5_NCI", "10x_Mix5_F_NCI", "10x_Mix5_NCI_M", "10x_Mix5_F_NCI_M", "10x_Mix5_F2_NCI_M","C1_FDA_HT_A", "C1_LLU_A", "ICELL8_SE_A", "ICELL8_PE_A", "C1_FDA_HT_B", "C1_LLU_B", "ICELL8_SE_B", "ICELL8_PE_B")

#batch <- c("10x_Mix10_LLU","10x_Mix5_NCI","10x_Mix5_F_NCI","10x_Mix5_F2_NCI")

correct_method <- c('Uncorrected','Seurat3.0','fastMNN','Scanorama','BBKNN','Harmony','LIMMA','ComBat')
path <- c('/unc/unc_pca.csv','/cca/cca_pca.csv','/mnn/mnn_embeddings.csv','/scanorama/scanorama_embeddings.csv','/bbknn/umap_coord_bbknn.csv','/harmony/harmony_embeddings.csv','/limma/limma_pca.csv','/combat/combat_pca.csv')

perct <- 0.01

align_score_cls1 <- align_score_cls2 <- vector('list',length=8)

cls_info <- read.csv(paste0(input,'cluster.csv'),row.names=1)

for(k in 1:length(hvg_num))
{
	temp_input <- paste0(input,'HVG_',hvg_num[k])
	
	for(j in 1:length(path))
	{
	temp <- read.csv(paste0(temp_input,path[j]),row.names=1)
	if(j ==5) {temp <- temp[,-ncol(temp)]}

	temp <- temp[rownames(cls_info),]
	temp_cls_info <- cls_info[!is.na(temp[,1]),]
	temp <- temp[!is.na(temp[,1]),]

	temp_cls1 <- temp[temp_cls_info[,2]=='HCC1395BL',]
	temp_cls1_info <- factor(temp_cls_info[temp_cls_info[,2]=='HCC1395BL',1],levels=batch)
	batch_cell1 <- table(temp_cls1_info)
	batch_cell1 <- batch_cell1[batch_cell1>0]
	align_score_cls1[[j]] <- rbind(align_score_cls1[[j]],calc_align_score(temp_cls1,batch_cell1,perct))

	temp_cls2 <- temp[temp_cls_info[,2]=='HCC1395',]
	temp_cls2_info <- factor(temp_cls_info[temp_cls_info[,2]=='HCC1395',1],levels=batch)
	batch_cell2 <- table(temp_cls2_info)
	batch_cell2 <- batch_cell2[batch_cell2>0]
	align_score_cls2[[j]] <- rbind(align_score_cls2[[j]],calc_align_score(temp_cls2,batch_cell2,perct))
	}
}

align_score_cls1_summary <- do.call(rbind,align_score_cls1)
align_score_cls2_summary <- do.call(rbind,align_score_cls2)
rownames(align_score_cls1_summary) <- rownames(align_score_cls2_summary) <- correct_method
colnames(align_score_cls1_summary)[ncol(align_score_cls1_summary)] <- colnames(align_score_cls2_summary)[ncol(align_score_cls2_summary)] <- 'Overall'

write.csv(round(align_score_cls1_summary,3),file=paste0(input,'/alignment_score_cls1.csv'),row.names=T)
write.csv(round(align_score_cls2_summary,3),file=paste0(input,'/alignment_score_cls2.csv'),row.names=T)



#######################################################
## get alighment score for scenario 2 and scenario 3 ##
#######################################################


sce <- c('scenario_2','scenario_3')
hvg_num <- 2000
batch <- c("10x_LLU","10x_NCI","C1_FDA_HT","C1_LLU","iCELL8_SE")
input <- '/.../4_Batch_correction/Results/'

correct_method <- c('Uncorrected','Seurat3.0','MNN','Scanorama','BBKNN','Harmony','LIMMA','ComBat')
path <- c('/unc/unc_pca.csv','/cca/cca_pca.csv','/mnn/mnn_embeddings.csv','/scanorama/scanorama_embeddings.csv','/bbknn/umap_coord_bbknn.csv','/harmony/harmony_embeddings.csv','/limma/limma_pca.csv','/combat/combat_pca.csv')

perct <- 0.01

for(i in 1:length(sce))
{

	align_score <- vector('list',length=8)

	for(k in 1:length(hvg_num))
	{
	temp_input <- paste0(input,sce[i],'/HVG_',hvg_num[k])

	for(j in 1:length(path))
	{
	temp <- read.csv(paste0(temp_input,path[j]),row.names=1)
	batch_cell <- cell_counts(temp)[1:length(batch)]
	if(j ==5) {temp <- temp[,-ncol(temp)]}
	align_score[[j]] <- rbind(align_score[[j]],calc_align_score(temp,batch_cell,perct))
	}

	}

	align_score_summary <- do.call(rbind,align_score)
	rownames(align_score_summary) <- paste(rep(correct_method,rep(length(hvg_num),length(correct_method))),rep(paste0('HVG_',hvg_num),length(correct_method)),sep=':')
	colnames(align_score_summary) <- c(batch[[i]],'Overall')

	write.csv(round(align_score_summary,3),file=paste0(input,sce[i],'/alignment_score.csv'),row.names=T)
}



	



###################################################################
## get kBET score for each cell type in scenario 1 or scenario 4 ##
###################################################################


library(kBET)

hvg_num <- 2000
input <- '/.../4_Batch_correction/Results/scenario_1/'
#input <- '/.../4_Batch_correction/Results/scenario_4/'

batch <- c("10x_LLU_A", "10x_NCI_A", "10x_NCI_M_A", "10x_LLU_B", "10x_NCI_B", "10x_NCI_M_B", "10x_Mix10_LLU", "10x_Mix5_NCI", "10x_Mix5_F_NCI", "10x_Mix5_NCI_M", "10x_Mix5_F_NCI_M", "10x_Mix5_F2_NCI_M","C1_FDA_HT_A", "C1_LLU_A", "ICELL8_SE_A", "ICELL8_PE_A", "C1_FDA_HT_B", "C1_LLU_B", "ICELL8_SE_B", "ICELL8_PE_B")

#batch <- c("10x_Mix10_LLU","10x_Mix5_NCI","10x_Mix5_F_NCI","10x_Mix5_F2_NCI")

correct_method <- c('Uncorrected','Seurat3.0','MNN','Scanorama','BBKNN','Harmony','LIMMA','ComBat')
path <- c('/unc/unc_pca.csv','/cca/cca_pca.csv','/mnn/mnn_embeddings.csv','/scanorama/scanorama_embeddings.csv','/bbknn/umap_coord_bbknn.csv','/harmony/harmony_embeddings.csv','/limma/limma_pca.csv','/combat/combat_pca.csv')

perct <- 0.01

kbet_score_cls1 <- kbet_score_cls2 <- vector('list',length=8)

cls_info <- read.csv(paste0(input,'cluster.csv'),row.names=1)
cls_info <- cls_info[!is.na(cls_info[,2]),]
temp <- c()
for(i in 1:length(batch))
{
	temp <- rbind(temp,cls_info[cls_info[,1]==batch[i],])
}
cls_info <- temp
cls1_info <- cls_info[cls_info[,2]=='HCC1395BL',]
cls2_info <- cls_info[cls_info[,2]=='HCC1395',]
batch_cell1 <- table(factor(cls1_info[,1],levels=batch))
batch_cell2 <- table(factor(cls2_info[,1],levels=batch))
batch_cell1 <- batch_cell1[batch_cell1>0]
batch_cell2 <- batch_cell2[batch_cell2>0]
batch_cell1 <- rep(names(batch_cell1),batch_cell1)
batch_cell2 <- rep(names(batch_cell2),batch_cell2)

for(j in 1:length(hvg_num))
{
	temp_input <- paste0(input,'/HVG_',hvg_num[j])

	for(i in 1:length(path))
	{
	temp <- read.csv(paste0(temp_input,path[i]),row.names=1)
	if(i ==5) {temp <- temp[,-ncol(temp)]}

	temp <- temp[rownames(cls_info),]
	temp_cls_info <- cls_info[!is.na(temp[,1]),]
	temp <- temp[!is.na(temp[,1]),]

	temp_cls1 <- temp[temp_cls_info[,2]=='HCC1395BL',]
	batch_cell1 <- table(temp_cls_info[temp_cls_info[,2]=='HCC1395BL',1])
	batch_cell1 <- batch_cell1[batch_cell1>0]
	batch_cell1 <- rep(names(batch_cell1),batch_cell1)
	k <- nrow(temp_cls1)*perct
#	k <- NULL
	kbet_score_cls1[[i]] <- c(kbet_score_cls1[[i]],kBET(temp_cls1,batch_cell1,do.pca=FALSE,plot=FALSE,k0=k)$summary$kBET.observed[1])

	temp_cls2 <- temp[cls_info[,2]=='HCC1395',]
	temp_cls2 <- temp[temp_cls_info[,2]=='HCC1395',]
	batch_cell2 <- table(temp_cls_info[temp_cls_info[,2]=='HCC1395',1])
	batch_cell2 <- batch_cell2[batch_cell2>0]
	batch_cell2 <- rep(names(batch_cell2),batch_cell2)
	k <- nrow(temp_cls2)*perct
#	k <- NULL
	kbet_score_cls2[[i]] <- c(kbet_score_cls2[[i]],kBET(temp_cls2,batch_cell2,do.pca=FALSE,plot=FALSE,k0=k)$summary$kBET.observed[1])
	}
}

kbet_score_cls1_summary <- do.call(rbind,kbet_score_cls1)
kbet_score_cls2_summary <- do.call(rbind,kbet_score_cls2)
rownames(kbet_score_cls1_summary) <- rownames(kbet_score_cls2_summary) <- correct_method
colnames(kbet_score_cls1_summary) <- colnames(kbet_score_cls2_summary) <- paste0('HVG_',hvg_num)
kbet_score_cls1_summary <- 1-kbet_score_cls1_summary
kbet_score_cls2_summary <- 1-kbet_score_cls2_summary
write.csv(round(kbet_score_cls1_summary,5),file=paste0(input,'/kBET_acceptence_score_cls1_1perct.csv'),row.names=T)
write.csv(round(kbet_score_cls2_summary,5),file=paste0(input,'/kBET_acceptence_score_cls2_1perct.csv'),row.names=T)



##################################################
## get kBET score for scenario 2 and scenario 3 ##
##################################################

library(kBET)

sce <- c('scenario_2','scenario_3')
input <- '/.../4_Batch_correction/Results/'
hvg_num <- 2000
batch <- c("10x_LLU","10x_NCI","C1_FDA_HT","C1_LLU","WaferGen_SE")

correct_method <- c('Uncorrected','Seurat3.0','MNN','Scanorama','BBKNN','Harmony','LIMMA','ComBat')
path <- c('/unc/unc_pca.csv','/cca/cca_pca.csv','/mnn/mnn_embeddings.csv','/scanorama/scanorama_embeddings.csv','/bbknn/umap_coord_bbknn.csv','/harmony/harmony_embeddings.csv','/limma/limma_pca.csv','/combat/combat_pca.csv')

perct <- 0.01
for(i in 1:length(sce))
{

	kbet_score <- vector('list',length=8)
	sil_score <- vector('list',length=8)

	for(j in 1:length(hvg_num))
	{
	temp_input <- paste0(input,sce[i],'/HVG_',hvg_num[j])

	for(l in 1:length(path))
	{
	temp <- read.csv(paste0(temp_input,path[l]),row.names=1)
	if(l ==5) {temp <- temp[,-ncol(temp)]}
	k <- nrow(temp)*perct
#	k <- NULL
	batch_cell <- rep(batch[[i]],cell_counts(temp)[1:length(batch[[i]])])
	kbet_score[[l]] <- c(kbet_score[[l]],kBET(temp,batch_cell,do.pca=FALSE,plot=FALSE,k0=k)$summary$kBET.observed[1])
	}

	}

	kbet_score_summary <- do.call(rbind,kbet_score)
	rownames(kbet_score_summary) <- correct_method
	colnames(kbet_score_summary) <- paste0('HVG_',hvg_num)
	kbet_score_summary <- 1-kbet_score_summary

	write.csv(round(kbet_score_summary,5),file=paste0(input,sce[i],'/kBET_acceptence_score_1perct.csv'),row.names=T)
#	write.csv(round(kbet_score_summary,5),file=paste0(input,sce[i],'/kBET_acceptence_score.csv'),row.names=T)
}
	





#######################################################
## get silhouette score for scenario 1 or scenario 4 ##
#######################################################

library(cluster)
library(ggplot2)

# read embeddings and cell identity files
input <- '/.../4_Batch_correction/Results/scenario_1/'
#input <- '/.../4_Batch_correction/Results/scenario_4/'

unc <- read.csv(paste0(input,"/HVG_2000/unc/unc_pca.csv"), row.names = 1)
cca <- read.csv(paste0(input,"/HVG_2000/cca/cca_pca.csv"), row.names = 1)
mnn <- read.csv(paste0(input,"/HVG_2000/mnn/mnn_embeddings.csv"), row.names = 1)
scanorama <- read.csv(paste0(input,"/HVG_2000/scanorama/scanorama_embeddings.csv"), row.names = 1)
bbknn <- read.csv(paste0(input,"/HVG_2000/bbknn/umap_coord_bbknn.csv"), row.names = 1)
harmony <- read.csv(paste0(input,"/HVG_2000/harmony/harmony_embeddings.csv"), row.names = 1)
limma <- read.csv(paste0(input,"/HVG_2000/limma/limma_pca.csv"), row.names = 1)
combat <- read.csv(paste0(input,"/HVG_2000/combat/combat_pca.csv"), row.names = 1)

# for already annotated cell names
cls <- read.csv(paste0(input,"cluster.csv"), row.names = 1)


# make equal of rownames of embeddings and order of cells
# different name scheme between bbknn/cca/mnn and others
unc <- unc[rownames(cls),]
cca <- cca[rownames(cls),]
mnn <- mnn[rownames(cls),]
scanorama <- scanorama[rownames(cls),]
harmony <- harmony[rownames(cls),]
bbknn <- bbknn[rownames(cls),]
limma <- limma[rownames(cls),]
combat <- combat[rownames(cls),]
cluster <- cls$cluster

# calculate dist
dist.unc <- dist(unc)
dist.cca <- dist(cca)
dist.mnn <- dist(mnn)
dist.scanorama <- dist(scanorama)
dist.bbknn <- dist(as.matrix(bbknn[,1:2])) # using UMAP_1 and UMAP_2
dist.harmony <- dist(harmony)
dist.limma <- dist(limma)
dist.combat <- dist(combat)

# Run silhouette
sil.unc <- silhouette(as.numeric(cluster), dist.unc)
sil.cca <- silhouette(as.numeric(cluster), dist.cca)
sil.mnn <- silhouette(as.numeric(cluster), dist.mnn)
sil.scanorama <- silhouette(as.numeric(cluster), dist.scanorama)
sil.bbknn <- silhouette(as.numeric(cluster), dist.bbknn)
sil.harmony <- silhouette(as.numeric(cluster), dist.harmony)
sil.limma <- silhouette(as.numeric(cluster), dist.limma)
sil.combat <- silhouette(as.numeric(cluster), dist.combat)

sil.mix <- list(unc=sil.unc,
                cca=sil.cca,
                mnn=sil.mnn,
                scanorama=sil.scanorama,
                bbknn=sil.bbknn,
                harmony=sil.harmony,
                limma=sil.limma,
                combat=sil.combat
                )

# 1. mean score per cluster (cell type)
# 2. mean score of two cell type
sil.unc.HCC1395 <- mean(sil.unc[which(cluster == "HCC1395"),3])
sil.unc.HCC1395BL <- mean(sil.unc[which(cluster == "HCC1395BL"),3])
sil.unc.total <- mean(c(sil.unc.HCC1395, sil.unc.HCC1395BL))

sil.cca.HCC1395 <- mean(sil.cca[which(cluster == "HCC1395"),3])
sil.cca.HCC1395BL <- mean(sil.cca[which(cluster == "HCC1395BL"),3])
sil.cca.total <- mean(c(sil.cca.HCC1395, sil.cca.HCC1395BL))

sil.mnn.HCC1395 <- mean(sil.mnn[which(cluster == "HCC1395"),3])
sil.mnn.HCC1395BL <- mean(sil.mnn[which(cluster == "HCC1395BL"),3])
sil.mnn.total <- mean(c(sil.mnn.HCC1395, sil.mnn.HCC1395BL))

sil.scanorama.HCC1395 <- mean(sil.scanorama[which(cluster == "HCC1395"),3])
sil.scanorama.HCC1395BL <- mean(sil.scanorama[which(cluster == "HCC1395BL"),3])
sil.scanorama.total <- mean(c(sil.scanorama.HCC1395, sil.scanorama.HCC1395BL))

sil.bbknn.HCC1395 <- mean(sil.bbknn[which(cluster == "HCC1395"),3])
sil.bbknn.HCC1395BL <- mean(sil.bbknn[which(cluster == "HCC1395BL"),3])
sil.bbknn.total <- mean(c(sil.bbknn.HCC1395, sil.bbknn.HCC1395BL))

sil.harmony.HCC1395 <- mean(sil.harmony[which(cluster == "HCC1395"),3])
sil.harmony.HCC1395BL <- mean(sil.harmony[which(cluster == "HCC1395BL"),3])
sil.harmony.total <- mean(c(sil.harmony.HCC1395, sil.harmony.HCC1395BL))

sil.limma.HCC1395 <- mean(sil.limma[which(cluster == "HCC1395"),3])
sil.limma.HCC1395BL <- mean(sil.limma[which(cluster == "HCC1395BL"),3])
sil.limma.total <- mean(c(sil.limma.HCC1395, sil.limma.HCC1395BL))

sil.combat.HCC1395 <- mean(sil.combat[which(cluster == "HCC1395"),3])
sil.combat.HCC1395BL <- mean(sil.combat[which(cluster == "HCC1395BL"),3])
sil.combat.total <- mean(c(sil.combat.HCC1395, sil.combat.HCC1395BL))


sil_matrix <- data.frame(
  method=c(rep("Uncorrected",3),rep("Seurat v3",3),rep("fastMNN",3),rep("Scanorama",3),
           rep("BBKNN",3),rep("Harmony",3),rep("limma",3),rep("ComBat",3)),
  score=c(sil.unc.HCC1395, sil.unc.HCC1395BL, sil.unc.total,
          sil.cca.HCC1395, sil.cca.HCC1395BL, sil.cca.total,
          sil.mnn.HCC1395, sil.mnn.HCC1395BL, sil.mnn.total,
          sil.scanorama.HCC1395, sil.scanorama.HCC1395BL, sil.scanorama.total,
          sil.bbknn.HCC1395, sil.bbknn.HCC1395BL, sil.bbknn.total,
          sil.harmony.HCC1395, sil.harmony.HCC1395BL, sil.harmony.total,
          sil.limma.HCC1395, sil.limma.HCC1395BL, sil.limma.total,
          sil.combat.HCC1395, sil.combat.HCC1395BL, sil.combat.total),
  type=c(rep(c("HCC1395","HCC1395BL","All"),8))
)














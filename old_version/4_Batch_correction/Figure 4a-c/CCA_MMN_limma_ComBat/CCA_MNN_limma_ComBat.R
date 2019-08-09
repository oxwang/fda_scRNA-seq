## load data, please load only one set of data (HCC1395, HCC1395BL, or mixture data) for each run ##

input <- '../manuscript_code/4_Batch_correction/Figure_4a-c/data/'

## 1. HCC1395

counts <- vector("list",5)
load(paste0(input,"10x/LLU/HCC1395/gene_counts_cellranger_filter.rdata"))
counts[[1]] <- gene_counts_filter
load(paste0(input,"10x/NCI/HCC1395/gene_counts_cellranger_filter.rdata"))
counts[[2]] <- gene_counts_filter
load(paste0(input,"C1_fda/HCC1395/gene_counts_featureCounts.rdata"))
counts[[3]] <- gene_counts
load(paste0(input,"C1/HCC1395/gene_counts_featureCounts.rdata"))
counts[[4]] <- gene_counts
load(paste0(input,"WaferGen_SE/HCC1395/gene_counts_featureCounts.rdata"))
counts[[5]] <- gene_counts
batch_cell <- c()
for(i in 1:length(counts))
{
	colnames(counts[[i]]) <- paste(colnames(counts[[i]]),i,sep="_")
	batch_cell <- c(batch_cell,ncol(counts[[i]]))
}
batch <- c("10x_LLU","10x_NCI","C1_FDA_HT","C1_LLU","WaferGen_SE")
col_type <- c("darkturquoise","darkorchid","hotpink","firebrick","seagreen4")



## 2. HCC1395BL

counts <- vector("list",5)
load(paste0(input,"10x/LLU/HCC1395BL/gene_counts_cellranger_filter.rdata"))
counts[[1]] <- gene_counts_filter
load(paste0(input,"10x/NCI/HCC1395BL/gene_counts_cellranger_filter.rdata"))
counts[[2]] <- gene_counts_filter
load(paste0(input,"C1_fda/HCC1395BL/gene_counts_featureCounts.rdata"))
counts[[3]] <- gene_counts
load(paste0(input,"C1/HCC1395BL/gene_counts_featureCounts.rdata"))
counts[[4]] <- gene_counts
load(paste0(input,"WaferGen_SE/HCC1395BL/gene_counts_featureCounts.rdata"))
counts[[5]] <- gene_counts
batch_cell <- c()
for(i in 1:length(counts))
{
	colnames(counts[[i]]) <- paste(colnames(counts[[i]]),i,sep="_")
	batch_cell <- c(batch_cell,ncol(counts[[i]]))
}
batch <- c("10x_LLU","10x_NCI","C1_FDA_HT","C1_LLU","WaferGen_SE")
col_type <- c("darkturquoise","darkorchid","hotpink","firebrick","seagreen4")



## 3. mixture data

counts <- vector("list",4)
load(paste0(input,"10x/LLU/Mix10/gene_counts_cellranger_filter.rdata"))
counts[[1]] <- gene_counts_filter
load(paste0(input,"10x/NCI/Mix5/gene_counts_cellranger_filter.rdata"))
counts[[2]] <- gene_counts_filter
load(paste0(input,"10x/NCI/Mix5_F/gene_counts_cellranger_filter.rdata"))
counts[[3]] <- gene_counts_filter
load(paste0(input,"10x/NCI/Mix5_F_2/gene_counts_cellranger_filter.rdata"))
counts[[4]] <- gene_counts_filter
batch_cell <- c()
for(i in 1:length(counts))
{
	colnames(counts[[i]]) <- paste(colnames(counts[[i]]),i,sep="_")
	batch_cell <- c(batch_cell,ncol(counts[[i]]))
}
batch <- c("10x_Mix10_LLU","10x_Mix5_NCI","10x_Mix5_F_NCI","10x_Mix5_F2_NCI")
col_type <- c("darkturquoise","darkorchid","hotpink","firebrick")



## below we used HCC1395 data as an example

######################
## CCA based method ##
######################


seurat_preprocess <- function(da)
{
	da <- da
	da <- CreateSeuratObject(raw.data=da,min.cells=3,min.genes=200)
	da <- NormalizeData(da)
	da <- ScaleData(da,display.progress=F)
	da <- FindVariableGenes(da,do.plot=F)
	return(da)
}


library(Seurat)
library(Rtsne)
library(Cairo)

hvg_num <- c(100,500,1000,2000,4000)
seurat_data <- vector("list",length(counts))
seurat_data_org <- vector("list",length(counts))

for(i in 1:length(counts))
{
	seurat_data_org[[i]] <- seurat_preprocess(counts[[i]])
}

for(k in 1:length(hvg_num))
{
	output <- paste(input,"/Eva_1/Multi_data/HCC1395/HVG_",hvg_num[k],"/",sep="")
#	output <- paste(input,"/Eva_1/Multi_data/HCC1395BL/HVG_",hvg_num[k],"/",sep="")
#	output <- paste(input,"/Eva_1/Multi_data/Mix/HVG_",hvg_num[k],"/",sep="")
	dir.create(paste(output,"CCA_figure",sep=""),recursive=T)
	
	for(i in 1:length(counts))
	{
		seurat_data[[i]] <- seurat_data_org[[i]]
		temp <- head(rownames(seurat_data[[i]]@hvg.info),hvg_num[k])
		if(i == 1)
		{genes_use <- temp}
		if(i > 1)
		{genes_use <- unique(c(genes_use,temp))}
		seurat_data[[i]]@meta.data[,"batch"] <- batch[i]
	}

	for(i in 1:length(counts))
	{
		genes_use <- intersect(genes_use, rownames(seurat_data[[i]]@scale.data))
	}

	seurat_combined <- RunMultiCCA(seurat_data, genes.use=genes_use, num.ccs=20)
	seurat_combined <- AlignSubspace(seurat_combined, reduction.type = "cca", grouping.var = "batch", dims.align = 1:15)
	t_cca <- seurat_combined@dr$cca.aligned@cell.embeddings
	set.seed(0)
	tsne_cca <- Rtsne(t_cca,perplexity=30)

	jpeg(paste(output,"CCA_figure/cca.pdf",sep=""),width=12, height=10,units="in", res=400)
	par(mar=c(4.5,5,0.5,0.5))
	plot(tsne_cca$Y[,1],tsne_cca$Y[,2],col=rep(col_type,batch_cell),cex=0.7,xlab="tSNE_1",ylab="tSNE_2",pch=16,cex.axis=2,cex.lab=2)
#	legend("bottomright",batch,col=col_type,pch=c(16,16),bty="n")
	dev.off()

	tsne_cca_coord_all <- tsne_cca$Y
	rownames(tsne_cca_coord_all) <- rownames(t_cca)
	seurat_cca_all <- seurat_combined
	save(tsne_cca_coord_all,file=paste(output,"tsne_coord_cca.rdata",sep=""))
	save(seurat_cca_all,file=paste(output,"seurat_cca.rdata",sep=""))
}






##################################
## MNN, limma and ComBat method ##
##################################


scran_preprocess <- function(da,batch)
{
	da <- da
	batch <- batch
	sample_id <- colnames(da)
	da_meta <- data.frame(Sample=sample_id,Batch=batch)
	sce <- SingleCellExperiment(list(counts=da),colData=da_meta,rowData=DataFrame(Symbol=rownames(da)))
	if(ncol(da) >= 100)
	{clusters <- quickCluster(da,min.size=ncol(da)/10,method="igraph")}
	if(ncol(da) < 100)
	{clusters <- quickCluster(da,min.size=ncol(da)/4,method="igraph")}
	sce <- computeSumFactors(sce, clusters=clusters)
	sce <- normalize(sce)
	fit <- trendVar(sce, use.spikes=F)
	dec <- decomposeVar(sce, fit=fit)
	results <- vector("list",2)
	results[[1]] <- sce
	results[[2]] <- dec
	return(results)
}




library(SingleCellExperiment)
library(scran)
library(Rtsne)
library(Matrix)
library(pracma)
library(limma)
library(sva)

hvg_num <- c(100,500,1000,2000,4000)

sce <- vector("list",length(batch))
dec <- vector("list",length(batch))
batch_id <- rep(batch,batch_cell)


temp_gene <- vector("list",length(batch))
for(i in 1:length(batch))
{
	temp <- scran_preprocess(counts[[i]],batch[i])
	sce[[i]] <- temp[[1]]
	dec[[i]] <- temp[[2]]
	temp_gene[[i]] <- rownames(dec[[i]])
}
universe <- Reduce(intersect,temp_gene)
combined_bio <- dec[[1]][universe,"bio"]
for(i in 2:length(batch))
{combined_bio <- combined_bio+dec[[i]][universe,"bio"]}

for(k in 1:length(hvg_num))
{
	output <- paste(input,"/Eva_1/Multi_data/HCC1395/HVG_",hvg_num[k],"/",sep="")
#	output <- paste(input,"/Eva_1/Multi_data/HCC1395BL/HVG_",hvg_num[k],"/",sep="")
#	output <- paste(input,"/Eva_1/Multi_data/Mix/HVG_",hvg_num[k],"/",sep="")
	dir.create(paste(output,"MNN_figure",sep=""),recursive=T)
	dir.create(paste(output,"UNC_figure",sep=""),recursive=T)
	dir.create(paste(output,"LIMMA_figure",sep=""),recursive=T)
	dir.create(paste(output,"ComBat_figure",sep=""),recursive=T)

	chosen <- universe[order(combined_bio, decreasing=TRUE)[1:hvg_num[k]]]
	sce_chosen <- vector("list",length(batch))
	for(i in 1:length(batch))
	{sce_chosen[[i]] <- sce[[i]][chosen,]}
	normed <- do.call(multiBatchNorm,sce_chosen)

## no correction

	combined <- c()
	for(i in 1:length(batch))
	{combined <- cbind(combined,logcounts(normed[[i]]))}
	t_unc <- as.matrix(t(combined))
	t_unc <- t_unc[!duplicated(t_unc),]
	set.seed(0)
	tsne_unc <- Rtsne(t_unc, perplexity = 30)
	tsne_unc_coord_all <- tsne_unc$Y
	rownames(tsne_unc_coord_all) <- rownames(t_unc)
	unc_all <- combined

	jpeg(paste(output,"UNC_figure/unc.jpeg",sep=""),width=12, height=10,units="in", res=400)
	plot(tsne_unc$Y[,1],tsne_unc$Y[,2],col=rep(col_type,batch_cell),cex=0.7,xlab="tSNE_1",ylab="tSNE_2",pch=c(16,16))
	legend("bottomright",batch,col=col_type,pch=c(16,16),bty="n")
	dev.off()


## MNN

	sce_log <- vector("list",length(batch))
	for(i in 1:length(batch))
	{sce_log[[i]] <- logcounts(normed[[i]])}
	mnn_out <- do.call(mnnCorrect,c(sce_log, list(k=20, sigma=0.1, cos.norm.in=TRUE, cos.norm.out=TRUE, compute.angle=TRUE, var.adj=TRUE)))
	t_mnn <- as.matrix(t(do.call(cbind, mnn_out$corrected)))
	rownames(t_mnn) <- colnames(combined)
	t_mnn <- t_mnn[!duplicated(t_mnn),]
	set.seed(0)
	tsne_mnn <- Rtsne(t_mnn, perplexity = 30)
	tsne_mnn_coord_all <- tsne_mnn$Y
	rownames(tsne_mnn_coord_all) <- rownames(t_mnn)
	mnn_all <- mnn_out

	jpeg(paste(output,"MNN_figure/mnn.jpeg",sep=""),width=12, height=10,units="in", res=400)
	plot(tsne_mnn$Y[,1],tsne_mnn$Y[,2],col=rep(col_type,batch_cell),cex=0.7,xlab="tSNE_1",ylab="tSNE_2",pch=c(16,16))
	legend("bottomright",batch,col=col_type,pch=c(16,16),bty="n")
	dev.off()



## limma 
	limma_all <- removeBatchEffect(as.matrix(combined),factor(batch_id))

## combat
	combat_all <- ComBat(as.matrix(combined),factor(batch_id),mod=NULL,prior.plots=F)

	t_limma <- as.matrix(t(limma_all))
	t_combat <- as.matrix(t(combat_all))
	set.seed(0)
	tsne_limma <- Rtsne(t_limma, perplexity = 30)
	tsne_combat <- Rtsne(t_combat, perplexity = 30)
	tsne_limma_coord_all <- tsne_limma$Y
	tsne_combat_coord_all <- tsne_combat$Y
	rownames(tsne_limma_coord_all) <- rownames(tsne_combat_coord_all) <- colnames(combined)
	

	jpeg(paste(output,"LIMMA_figure/limma.jpeg",sep=""),width=12, height=10,units="in", res=400)
	plot(tsne_limma$Y[,1],tsne_limma$Y[,2],col=rep(col_type,batch_cell),cex=0.7,xlab="tSNE_1",ylab="tSNE_2",pch=c(16,16))
	legend("bottomright",batch,col=col_type,pch=c(16,16),bty="n")
	dev.off()


	jpeg(paste(output,"ComBat_figure/combat.jpeg",sep=""),width=12, height=10,units="in", res=400)
	plot(tsne_combat$Y[,1],tsne_combat$Y[,2],col=rep(col_type,batch_cell),cex=0.7,xlab="tSNE_1",ylab="tSNE_2",pch=c(16,16))
	legend("bottomright",batch,col=col_type,pch=c(16,16),bty="n")
	dev.off()


	save(tsne_mnn_coord_all,file=paste(output,"tsne_coord_mnn.rdata",sep=""))
	save(tsne_unc_coord_all,file=paste(output,"tsne_coord_unc.rdata",sep=""))
	save(mnn_all,file=paste(output,"mnn.rdata",sep=""))
	save(unc_all,file=paste(output,"unc.rdata",sep=""))
	save(tsne_limma_coord_all,file=paste(output,"tsne_coord_limma.rdata",sep=""))
	save(tsne_combat_coord_all,file=paste(output,"tsne_coord_combat.rdata",sep=""))
	save(limma_all,file=paste(output,"limma.rdata",sep=""))
	save(combat_all,file=paste(output,"combat.rdata",sep=""))
}







#########################
## get alignment score ##
#########################

calc_align_score <- function(tsne_coord,batch_cell,cell_perct)
{
	tsne_coord <- tsne_coord
	batch_cell <- batch_cell
	cell_perct <- cell_perct
	w <- batch_cell/sum(batch_cell)

	temp <- cumsum(batch_cell)
	batch_index <- cbind(c(1,temp[-length(temp)]+1),temp)

	tsne_dist <- as.matrix(dist(tsne_coord))
	diag(tsne_dist) <- 1e5
	k <- round(nrow(tsne_coord)*cell_perct)
	tsne_order <- t(apply(tsne_dist,1,order))

	k_mean <- c()
	for(i in 1:length(batch_cell))
	{
		temp <- tsne_order[batch_index[i,1]:batch_index[i,2],1:k]
		temp2 <- temp <= batch_index[i,2]
		temp3 <- temp >= batch_index[i,1]
		temp <- temp2+temp3
		k_mean <- c(k_mean,mean(apply(temp == 2,1,sum)))
	}
	align_score <- 1-abs(k_mean-k*w)/(k-k*w)
	align_score <- c(sum(w*align_score),align_score)

	return(align_score)
}




## alignment score for uncorrected, mnn, and cca ##

batch <- c("10x_LLU","10x_NCI","C1_FDA_HT","C1_LLU","WaferGen_SE")
batch <- c("10x_Mix10_LLU","10x_Mix5_NCI","10x_Mix5_F_NCI","10x_Mix5_F2_NCI")
hvg_num <- c(100,500,1000,2000,4000)

align_score_unc_all <- c()
align_score_mnn_all <- c()
align_score_cca_all <- c()

for(k in 1:length(hvg_num))
{
output <- paste(input,"/Eva_1/Multi_data/HCC1395/HVG_",hvg_num[k],"/",sep="")
#output <- paste(input,"/Eva_1/Multi_data/HCC1395BL/HVG_",hvg_num[k],"/",sep="")
#output <- paste(input,"/Eva_1/Multi_data/Mix/HVG_",hvg_num[k],"/",sep="")

load(paste(output,"tsne_coord_mnn.rdata",sep=""))
load(paste(output,"tsne_coord_unc.rdata",sep=""))
load(paste(output,"tsne_coord_cca.rdata",sep=""))
load(paste(output,"mnn.rdata",sep=""))
load(paste(output,"seurat_cca.rdata",sep=""))

batch_cell_unc <- unlist(lapply(mnn_all$corrected,ncol))
batch_cell_cca <- table(seurat_cca_all@meta.data[,"batch"])
batch_cell_cca <- batch_cell_cca[match(batch,names(batch_cell_cca))]
temp_batch <- rep(1:length(batch_cell_unc),batch_cell_unc)
temp_batch <- temp_batch[!duplicated(as.matrix(t(do.call(cbind, mnn_all$corrected))))]
batch_cell_mnn <- table(temp_batch)

perct <- 0.01
align_score_unc_all <- rbind(align_score_unc_all,calc_align_score(tsne_unc_coord_all,batch_cell_unc,perct))
align_score_mnn_all <- rbind(align_score_mnn_all,calc_align_score(tsne_mnn_coord_all,batch_cell_mnn,perct))
align_score_cca_all <- rbind(align_score_cca_all,calc_align_score(tsne_cca_coord_all,batch_cell_cca,perct))

}


align_score_summary <- round(rbind(align_score_unc_all,align_score_mnn_all,align_score_cca_all),3)
colnames(align_score_summary) <- c("All",batch)
rownames(align_score_summary) <- paste(rep(c("UNC","MNN","CCA"),rep(length(hvg_num),3)),"HVG",rep(hvg_num,3),sep="_")
write.csv(align_score_summary,file=paste0(input,"/Eva_1/Multi_data/HCC1395/align_score_summary.csv"),quote=F,row.names=T)



## alignment score for limma and combat ##


batch <- c("10x_LLU","10x_NCI","C1_FDA_HT","C1_LLU","WaferGen_SE")
batch <- c("10x_Mix10_LLU","10x_Mix5_NCI","10x_Mix5_F_NCI","10x_Mix5_F2_NCI")
hvg_num <- c(100,500,1000,2000,4000)

align_score_limma_all <- c()
align_score_combat_all <- c()

for(k in 1:length(hvg_num))
{

output <- paste0(input,"/Eva_1/Multi_data/HCC1395/HVG_",hvg_num[k],"/")
load(paste(output,"tsne_coord_limma.rdata",sep=""))
load(paste(output,"tsne_coord_combat.rdata",sep=""))
load(paste(output,"mnn.rdata",sep=""))

batch_cell <- unlist(lapply(mnn_all$corrected,ncol))

perct <- 0.01
align_score_limma_all <- rbind(align_score_limma_all,calc_align_score(tsne_limma_coord_all,batch_cell,perct))
align_score_combat_all <- rbind(align_score_combat_all,calc_align_score(tsne_combat_coord_all,batch_cell,perct))

}


align_score_summary <- round(rbind(align_score_limma_all,align_score_combat_all),3)
colnames(align_score_summary) <- c("All",batch)
rownames(align_score_summary) <- paste(rep(c("LIMMA","ComBat"),rep(length(hvg_num),2)),"HVG",rep(hvg_num,2),sep="_")
write.csv(align_score_summary,file=paste0(input,"/Eva_1/Multi_data/Mix/align_score_summary2.csv"),quote=F,row.names=T)





## get cca, mnn, limma, combat results for umap plot ##

input <- paste0(input,"/Eva_1/Multi_data/")

sample <- c("HCC1395BL","HCC1395BL","Mix")
hvg <- paste("HVG",c(100,500,1000,2000,4000),sep="_")

for(i in 1:length(sample))
{
	for(j in 1:length(hvg))
	{
	input2 <- paste(input,sample[i],"/",hvg[j],sep="")
	
	da <- get(load(paste(input2,"/seurat_cca.rdata",sep="")))
	da2 <- da@dr$cca.aligned@cell.embeddings
	output <- paste(input,"/cca",sep="")
	dir.create(output,recursive=T)
	write.csv(da2,file=paste(output,"/cca_correct.csv",sep=""),row.names=T)
	
	da <- get(load(paste(input2,"/unc.rdata",sep="")))
	da2 <- as.matrix(t(da))
	cell_id <- rownames(da2)
	output <- paste(input2,"/unc",sep="")
	dir.create(output,recursive=T)
	write.csv(da2,file=paste(output,"/unc_correct.csv",sep=""),row.names=T)

	da <- get(load(paste(input2,"/limma.rdata",sep="")))
	da2 <- as.matrix(t(da))
	output <- paste(input2,"/limma",sep="")
	dir.create(output,recursive=T)
	write.csv(da2,file=paste(output,"/limma_correct.csv",sep=""),row.names=T)

	da <- get(load(paste(input2,"/combat.rdata",sep="")))
	da2 <- as.matrix(t(da))
	output <- paste(input2,"/combat",sep="")
	dir.create(output,recursive=T)
	write.csv(da2,file=paste(output,"/combat_correct.csv",sep=""),row.names=T)

	da <- get(load(paste(input2,"/mnn.rdata",sep="")))
	da2 <- as.matrix(t(do.call(cbind, da$corrected)))
	rownames(da2) <- cell_id
	da2 <- da2[!duplicated(da2),]
	output <- paste(input2,"/mnn",sep="")
	dir.create(output,recursive=T)
	write.csv(da2,file=paste(output,"/mnn_correct.csv",sep=""),row.names=T)
	}
}
















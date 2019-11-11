#############################################################
## load gene count data for four scenarios and tian's data ##
#############################################################


## scenario 1: get gene counts data for all 20 data sets

counts <- vector("list",20)

load("/.../4_Batch_correction/Data/scenario_1/HCC1395/LLU/gene_counts_cellranger_3.1.rdata")
counts[[1]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/HCC1395/NCI/standard/gene_counts_cellranger_3.1.rdata")
counts[[2]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/HCC1395/NCI/modified/gene_counts_cellranger_3.1.rdata")
counts[[3]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/HCC1395BL/LLU/gene_counts_cellranger_3.1.rdata")
counts[[4]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/HCC1395BL/NCI/standard/gene_counts_cellranger_3.1.rdata")
counts[[5]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/HCC1395BL/NCI/modified/gene_counts_cellranger_3.1.rdata")
counts[[6]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/Spike-in/Mix10/gene_counts_cellranger_3.1.rdata")
counts[[7]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/Spike-in/Mix5/standard/gene_counts_cellranger_3.1.rdata")
counts[[8]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/Spike-in/Mix5_F/standard/gene_counts_cellranger_3.1.rdata")
counts[[9]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/Spike-in/Mix5/modified/gene_counts_cellranger_3.1.rdata")
counts[[10]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/Spike-in/Mix5_F/modified/gene_counts_cellranger_3.1.rdata")
counts[[11]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/Spike-in/Mix5_F_2/gene_counts_cellranger_3.1.rdata")
counts[[12]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/HCC1395/C1_HT/gene_counts_featureCounts.rdata")
counts[[13]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/HCC1395/C1/gene_counts_featureCounts.rdata")
counts[[14]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/HCC1395/iCELL8/SE/gene_counts_featureCounts.rdata")
counts[[15]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/HCC1395/iCELL8/PE/gene_counts_featureCounts.rdata")
counts[[16]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/HCC1395BL/C1_HT/gene_counts_featureCounts.rdata")
counts[[17]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/HCC1395BL/C1/gene_counts_featureCounts.rdata")
counts[[18]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/HCC1395BL/iCELL8/SE/gene_counts_featureCounts.rdata")
counts[[19]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_1/HCC1395BL/iCELL8/PE/gene_counts_featureCounts.rdata")
counts[[20]] <- gene_counts


batch <- c("10x_LLU_A", "10x_NCI_A", "10x_NCI_M_A", "10x_LLU_B", "10x_NCI_B", "10x_NCI_M_B", "10x_Mix10_LLU", "10x_Mix5_NCI", "10x_Mix5_F_NCI", "10x_Mix5_NCI_M", "10x_Mix5_F_NCI_M", "10x_Mix5_F2_NCI_M","C1_FDA_HT_A", "C1_LLU_A", "ICELL8_SE_A", "ICELL8_PE_A", "C1_FDA_HT_B", "C1_LLU_B", "ICELL8_SE_B", "ICELL8_PE_B")
col_type <- rainbow(20)


## scenario 2: get gene counts data for HCC1395

counts <- vector("list",5)
load("/.../4_Batch_correction/Data/scenario_2/LLU/gene_counts_cellranger_3.1.rdata")
counts[[1]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_2/NCI/gene_counts_cellranger_3.1.rdata")
counts[[2]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_2/C1_HT/gene_counts_featureCounts.rdata")
counts[[3]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_2/C1/gene_counts_featureCounts.rdata")
counts[[4]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_2/iCELL8/gene_counts_featureCounts.rdata")
counts[[5]] <- gene_counts

batch <- c("10x_LLU","10x_NCI","C1_FDA_HT","C1_LLU","iCell8_SE")
col_type <- c("darkturquoise","darkorchid","hotpink","firebrick","seagreen4")



## scenario 3: get gene counts data for HCC1395BL

counts <- vector("list",5)
load("/.../4_Batch_correction/Data/scenario_3/LLU/gene_counts_cellranger_3.1.rdata")
counts[[1]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_3/NCI/gene_counts_cellranger_3.1.rdata")
counts[[2]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_3/C1_HT/gene_counts_featureCounts.rdata")
counts[[3]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_3/C1/gene_counts_featureCounts.rdata")
counts[[4]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_3/iCELL8/gene_counts_featureCounts.rdata")
counts[[5]] <- gene_counts

batch <- c("10x_LLU","10x_NCI","C1_FDA_HT","C1_LLU","iCell8_SE")
col_type <- c("darkturquoise","darkorchid","hotpink","firebrick","seagreen4")



## scenario 4: get gene counts data for spike-in data

counts <- vector("list",4)
load("/.../4_Batch_correction/Data/scenario_4/Mix10/gene_counts_cellranger_3.1.rdata")
counts[[1]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_4/Mix5/gene_counts_cellranger_3.1.rdata")
counts[[2]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_4/Mix5_F/gene_counts_cellranger_3.1.rdata")
counts[[3]] <- gene_counts
load("/.../4_Batch_correction/Data/scenario_4/Mix5_F_2/gene_counts_cellranger_3.1.rdata")
counts[[4]] <- gene_counts

batch <- c("10x_Mix10_LLU","10x_Mix5_NCI","10x_Mix5_F_NCI","10x_Mix5_F2_NCI")
col_type <- c("darkturquoise","darkorchid","hotpink","firebrick")



## supplementary Tian's data

counts <- vector("list",4)
load("/.../4_Batch_correction/Data/Tian/10x_3cl/gene_counts_cellranger_3.0.rdata")
counts[[1]] <- gene_counts
load("/.../4_Batch_correction/Data/Tian/10x_5cl/gene_counts_cellranger_3.0.rdata")
counts[[2]] <- gene_counts
load("/.../4_Batch_correction/Data/Tian/Celseq2/gene_counts_umitools_1.0.rdata")
counts[[3]] <- gene_counts
load("/.../4_Batch_correction/Data/Tian/Dropseq/gene_counts_umitools_1.0.rdata")
counts[[4]] <- gene_counts

batch <- c("sc_10x_3cl","sc_10x_5cl","sc_CEL-seq2","sc_Drop-seq")
col_type <- c("darkturquoise","darkorchid","hotpink","firebrick")




################################################################
## The following code is applied for each scenario to ##########
## perform batch correction. We used scenario 2 as an example ##
################################################################



###################################################################
## low quality cell filtering and logCPM normalization by Seurat ##
###################################################################

library(Seurat)

cca_process <- counts_norm <- vector('list',length=5)
#cca_process <- counts_norm <- vector('list',length=4)

batch_cell <- c()

for(i in 1:length(counts))
{
	colnames(counts[[i]]) <- paste(colnames(counts[[i]]),i,sep="_")
	meta_da <- data.frame(matrix(batch[i],nrow=ncol(counts[[i]]),ncol=1))
	rownames(meta_da) <- colnames(counts[[i]])
	colnames(meta_da) <- 'batch'
	tmpCells <- CreateSeuratObject(counts[[i]],min.cells = 3,min.features = 200, meta.data=meta_da)
	Total_mRNAs <- tmpCells[["nCount_RNA"]]$nCount_RNA
	mupper_bound <- 10^(mean(log10(Total_mRNAs)) + 2*sd(log10(Total_mRNAs)))
	mlower_bound <- 10^(mean(log10(Total_mRNAs)) - 2*sd(log10(Total_mRNAs)))
	Total_Genes <- tmpCells[["nFeature_RNA"]]$nFeature_RNA
	gupper_bound <- 10^(mean(log10(Total_Genes)) + 2*sd(log10(Total_Genes)))
	glower_bound <- 10^(mean(log10(Total_Genes)) - 2*sd(log10(Total_Genes)))
	tmpCells <- subset(x = tmpCells, subset = nFeature_RNA > glower_bound & nFeature_RNA < gupper_bound & nCount_RNA > mlower_bound & nCount_RNA < mupper_bound)
	tmpCells <- NormalizeData(tmpCells,verbose=F)

## generate raw counts data for scanorama and bbknn
	counts[[i]] <- GetAssayData(object=tmpCells, slot='counts')

## generate logCPM normalized data for uncorrected, limma, and ComBat
	counts_norm[[i]] <- GetAssayData(object=tmpCells, slot='data')
	batch_cell <- c(batch_cell,ncol(counts_norm[[i]]))

## generate cca_process data for seurat integration
	cca_process[[i]] <- tmpCells

## generate merged Cells data for MNN and Harmony
	if(i == 1) {Cells <- tmpCells}
	else {Cells <- merge(Cells,tmpCells)}
	
}


## generate filtered gene counts matrices for scanorama and bbknn

output <- '/.../4_Batch_correction/Results/scenario_2/'

for(i in 1:length(counts))
{
	gene_counts <- counts[[i]]
	gene_counts <- cbind.data.frame(rownames(gene_counts),as.matrix(gene_counts))
	colnames(gene_counts) <- c("gene",colnames(gene_counts)[-1])
	write.table(gene_counts,file=paste0(output,batch[i],'.txt'),sep="\t",row.names=F,quote=F)
}



######################
## batch correction ##
######################

hvg_num <- 2000


######################
## Seurat V3 method ##
######################

cca_da <- vector('list',length(cca_process))
for(i in 1:length(cca_da))
{
	cca_da[[i]] <- FindVariableFeatures(cca_process[[i]],nfeatures=hvg_num)
}
cca_anchors <- FindIntegrationAnchors(object.list=cca_da,anchor.features=hvg_num,k.filter=50)
cca_integrated <- IntegrateData(cca_anchors)
cca_integrated <- ScaleData(cca_integrated,verbose=F)
cca_integrated <- RunPCA(cca_integrated, npcs=30, verbose=F)
cca_integrated <- RunTSNE(cca_integrated,reduction='pca',dims=1:30)
cca_integrated <- RunUMAP(cca_integrated,reduction='pca',dims=1:30)

cca_pca <- Embeddings(cca_integrated[['pca']])
cca_tsne <- Embeddings(cca_integrated[['tsne']])
cca_umap <- Embeddings(cca_integrated[['umap']])

temp_output <- paste0(output,'HVG_',hvg_num,'/cca')
dir.create(temp_output,recursive=T)
write.csv(cca_pca,file=paste0(temp_output,'/cca_pca.csv'),row.names=T)
write.csv(cca_tsne,file=paste0(temp_output,'/cca_tsne.csv'),row.names=T)
write.csv(cca_umap,file=paste0(temp_output,'/cca_umap.csv'),row.names=T)

pdf(paste(temp_output,"/fig_cca_tsne.pdf",sep=""))
par(mar=c(4.5,5,0.5,0.5))
p <- DimPlot(cca_integrated, reduction = "tsne", group.by = "batch",cols=col_type)
print(p)
dev.off()

pdf(paste(temp_output,"/fig_cca_umap.pdf",sep=""))
par(mar=c(4.5,5,0.5,0.5))
p <- DimPlot(cca_integrated, reduction = "umap", group.by = "batch",cols=col_type)
print(p)
dev.off()






################
## MNN method ##
################

library(SeuratWrappers)

mnn_data <- Cells
mnn_data <- FindVariableFeatures(mnn_data,nfeatures=hvg_num)

mnn_data <- RunFastMNN(object.list=SplitObject(mnn_data, split.by = 'batch'),features=hvg_num)
mnn_data <- RunTSNE(mnn_data,reduction='mnn',dims=1:30)
mnn_data <- RunUMAP(mnn_data,reduction='mnn',dims=1:30)

mnn_embeddings <- Embeddings(mnn_data[['mnn']])
mnn_tsne <- Embeddings(mnn_data[['tsne']])
mnn_umap <- Embeddings(mnn_data[['umap']])

temp_output <- paste0(output,'HVG_',hvg_num,'/mnn')
dir.create(temp_output,recursive=T)
write.csv(mnn_embeddings,file=paste0(temp_output,'/mnn_embeddings.csv'),row.names=T)
write.csv(mnn_tsne,file=paste0(temp_output,'/mnn_tsne.csv'),row.names=T)
write.csv(mnn_umap,file=paste0(temp_output,'/mnn_umap.csv'),row.names=T)

pdf(paste(temp_output,"/fig_mnn_tsne.pdf",sep=""))
par(mar=c(4.5,5,0.5,0.5))
p <- DimPlot(mnn_data, reduction = "tsne", group.by = "batch",cols=col_type)
print(p)
dev.off()

pdf(paste(temp_output,"/fig_mnn_umap.pdf",sep=""))
par(mar=c(4.5,5,0.5,0.5))
p <- DimPlot(mnn_data, reduction = "umap", group.by = "batch",cols=col_type)
print(p)
dev.off()






####################
## Hormony method ##
####################


library(Seurat)
library(harmony)
library(Rtsne)

harmony_data <- Cells
harmony_data <- FindVariableFeatures(harmony_data,nfeatures=hvg_num)
harmony_data <- ScaleData(harmony_data,verbose=FALSE)
harmony_data <- RunPCA(harmony_data,pc.genes=harmony_data@var.genes,npcs=30,verbose=FALSE)
harmony_data <- RunHarmony(harmony_data,"batch")
harmony_data <- RunTSNE(harmony_data,reduction='harmony',dims=1:30)
harmony_data <- RunUMAP(harmony_data,reduction='harmony',dims=1:30)
harmony_embeddings <- Embeddings(harmony_data[['harmony']])
harmony_tsne <- Embeddings(harmony_data[['tsne']])
harmony_umap <- Embeddings(harmony_data[['umap']])

temp_output <- paste0(output,'HVG_',hvg_num,'/harmony')
dir.create(temp_output,recursive=T)
write.csv(harmony_embeddings,file=paste0(temp_output,'/harmony_embeddings.csv'),row.names=T)
write.csv(harmony_tsne,file=paste0(temp_output,'/harmony_tsne.csv'),row.names=T)
write.csv(harmony_umap,file=paste0(temp_output,'/harmony_umap.csv'),row.names=T)

pdf(paste(temp_output,"/fig_harmony_tsne.pdf",sep=""))
par(mar=c(4.5,5,0.5,0.5))
p <- DimPlot(harmony_data, reduction = "tsne", group.by = "batch",cols=col_type)
print(p)
dev.off()

pdf(paste(temp_output,"/fig_harmony_umap.pdf",sep=""))
par(mar=c(4.5,5,0.5,0.5))
p <- DimPlot(harmony_data, reduction = "umap", group.by = "batch",cols=col_type)
print(p)
dev.off()




######################
## scanorama method ##
######################

## please check scanorama folder for python-based batch correction

## load scanorama corrected data to seurat for TSNE and UMAP 

library(Seurat)

batch <- c("10x_LLU","10x_NCI","C1_FDA_HT","C1_LLU","iCell8_SE")
col_type <- c("darkturquoise","darkorchid","hotpink","firebrick","seagreen4")

temp_output <- paste0(output,'HVG_',hvg_num,'/scanorama')
counts <- vector('list',length(batch))
batch_cell <- c()
cell_id <- c()
for(j in 1:length(batch))
{
	counts[[j]] <- read.table(paste0(temp_output,batch[j],'.scanorama_corrected.txt'),header=T,row.names=1,sep='\t')
	batch_cell <- c(batch_cell,ncol(counts[[j]]))
	colnames(counts[[j]]) <- paste0(colnames(counts[[j]]),'_',j)
	cell_id <- c(cell_id,colnames(counts[[j]]))
}

meta_da <- data.frame(matrix(rep(batch,batch_cell),ncol=1))
rownames(meta_da) <- cell_id
colnames(meta_da) <- 'batch'

scanorama_da <- CreateSeuratObject(counts=do.call(cbind,counts),min.cells=0,min.features=0,meta.data=meta_da)
scanorama_da <- ScaleData(scanorama_da)
scanorama_da <- RunPCA(scanorama_da,features=rownames(counts[[1]]),npcs=30,verbose=FALSE)
scanorama_da <- RunTSNE(scanorama_da,reduction='pca',dims=1:30)
scanorama_da <- RunUMAP(scanorama_da,reduction='pca',dims=1:30)

scanorama_pca <- Embeddings(scanorama_da[['pca']])
scanorama_tsne <- Embeddings(scanorama_da[['tsne']])
scanorama_umap <- Embeddings(scanorama_da[['umap']])

write.csv(scanorama_pca,file=paste0(temp_output,'/scanorama_embeddings.csv'),row.names=T)
write.csv(scanorama_tsne,file=paste0(temp_output,'/scanorama_tsne.csv'),row.names=T)
write.csv(scanorama_umap,file=paste0(temp_output,'/scanorama_umap.csv'),row.names=T)

pdf(paste(temp_output,"/fig_scanorama_tsne.pdf",sep=""))
par(mar=c(4.5,5,0.5,0.5))
p <- DimPlot(scanorama_da, reduction = "tsne", group.by = "batch",cols=col_type)
print(p)
dev.off()

pdf(paste(temp_output,"/fig_scanorama_umap.pdf",sep=""))
par(mar=c(4.5,5,0.5,0.5))
p <- DimPlot(scanorama_da, reduction = "umap", group.by = "batch",cols=col_type)
print(p)
dev.off()






##################
## bbknn method ##
##################

## please check bbknn folder for python-based batch correction

## load bbknn umap coordinates for fig generation


batch <- c("10x_LLU","10x_NCI","C1_FDA_HT","C1_LLU","iCell8_SE")
col_type <- c("darkturquoise","darkorchid","hotpink","firebrick","seagreen4")

temp_output <- paste0(output,'HVG_',hvg_num,'/bbknn/')

umap <- read.csv(paste0(temp_output,'umap_coord_bbknn.csv'),row.names=1)
temp_batch <- table(umap[,3])[batch]
colnames(umap) <- c('UMAP_1','UMAP_2','sample')

pdf(paste(temp_output,"/fig_bbknn_umap.pdf",sep=""))
par(mar=c(4.5,5,0.5,0.5))
p <- ggplot(umap,aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.2) + 
scale_color_manual(values=col_type) + theme_classic() + guides(colour = guide_legend(override.aes = list(size=3))) +
theme(axis.title=element_text(size=15), axis.text=element_text(size=12,color='black'),legend.title = element_blank(),legend.text = element_text(size=12))		
print(p)
dev.off()




###########################################
## uncorrected, limma, and combat method ##
###########################################

library(limma)
library(sva)

common_genes <- rownames(counts_norm[[1]])
for(i in 2:length(counts))
{
	common_genes <- intersect(common_genes,rownames(counts_norm[[i]]))
}

common_counts_norm <- vector('list',length(counts_norm))
for(i in 1:length(counts_norm))
{	
	common_counts_norm[[i]] <- counts_norm[[i]][common_genes,]
}



## uncorrected

unc_process <- do.call(cbind,common_counts_norm)
batch_id <- rep(batch,batch_cell)
meta_da <- data.frame(batch_id)
rownames(meta_da) <- colnames(unc_process)
colnames(meta_da) <- 'batch'
unc_da <- CreateSeuratObject(unc_process,min.cells = 3,min.features = 200, meta.data=meta_da)


unc_out <- FindVariableFeatures(unc_da,nfeatures=hvg_num)
unc_out <- ScaleData(unc_out,verbose=F)
unc_out <- RunPCA(unc_out, npcs=30, verbose=F)
unc_out <- RunTSNE(unc_out,reduction='pca',dims=1:30)
unc_out <- RunUMAP(unc_out,reduction='pca',dims=1:30)

unc_pca <- Embeddings(unc_out[['pca']])
unc_tsne <- Embeddings(unc_out[['tsne']])
unc_umap <- Embeddings(unc_out[['umap']])

temp_output <- paste0(output,'HVG_',hvg_num,'/unc')
dir.create(temp_output,recursive=T)
write.csv(unc_pca,file=paste0(temp_output,'/unc_pca.csv'),row.names=T)
write.csv(unc_tsne,file=paste0(temp_output,'/unc_tsne.csv'),row.names=T)
write.csv(unc_umap,file=paste0(temp_output,'/unc_umap.csv'),row.names=T)

pdf(paste(temp_output,"/fig_unc_tsne.pdf",sep=""))
par(mar=c(4.5,5,0.5,0.5))
p <- DimPlot(unc_out, reduction = "tsne", group.by = "batch",cols=col_type)
print(p)
dev.off()

pdf(paste(temp_output,"/fig_unc_umap.pdf",sep=""))
par(mar=c(4.5,5,0.5,0.5))
p <- DimPlot(unc_out, reduction = "umap", group.by = "batch",cols=col_type)
print(p)
dev.off()




## limma ##

limma_process <- common_counts_norm
batch_id <- rep(batch,batch_cell)
limma_combined <- do.call(cbind,limma_process)
limma_all <- removeBatchEffect(as.matrix(limma_combined),factor(batch_id))
meta_da <- data.frame(batch_id)
rownames(meta_da) <- colnames(limma_combined)
colnames(meta_da) <- 'batch'
limma_da <- CreateSeuratObject(limma_all,min.cells = 3,min.features = 200, meta.data=meta_da)

limma_out <- FindVariableFeatures(limma_da,nfeatures=hvg_num)
limma_out <- ScaleData(limma_out,verbose=F)
limma_out <- RunPCA(limma_out, npcs=30, verbose=F)
limma_out <- RunTSNE(limma_out,reduction='pca',dims=1:30,check_duplicates = FALSE)
limma_out <- RunUMAP(limma_out,reduction='pca',dims=1:30)

limma_pca <- Embeddings(limma_out[['pca']])
limma_tsne <- Embeddings(limma_out[['tsne']])
limma_umap <- Embeddings(limma_out[['umap']])

temp_output <- paste0(output,'HVG_',hvg_num,'/limma')
dir.create(temp_output,recursive=T)
write.csv(limma_pca,file=paste0(temp_output,'/limma_pca.csv'),row.names=T)
write.csv(limma_tsne,file=paste0(temp_output,'/limma_tsne.csv'),row.names=T)
write.csv(limma_umap,file=paste0(temp_output,'/limma_umap.csv'),row.names=T)

pdf(paste(temp_output,"/fig_limma_tsne.pdf",sep=""))
par(mar=c(4.5,5,0.5,0.5))
p <- DimPlot(limma_out, reduction = "tsne", group.by = "batch",cols=col_type)
print(p)
dev.off()

pdf(paste(temp_output,"/fig_limma_umap.pdf",sep=""))
par(mar=c(4.5,5,0.5,0.5))
p <- DimPlot(limma_out, reduction = "umap", group.by = "batch",cols=col_type)
print(p)
dev.off()



## combat ##

combat_process <- common_counts_norm
batch_id <- rep(batch,batch_cell)
combat_combined <- do.call(cbind,combat_process)
combat_all <- ComBat(as.matrix(combat_combined),factor(batch_id),mod=NULL,prior.plots=F)
meta_da <- data.frame(batch_id)
rownames(meta_da) <- colnames(combat_combined)
colnames(meta_da) <- 'batch'
combat_da <- CreateSeuratObject(combat_all,min.cells = 3,min.features = 200, meta.data=meta_da)

combat_out <- FindVariableFeatures(combat_da,nfeatures=hvg_num)
combat_out <- ScaleData(combat_out,verbose=F)
combat_out <- RunPCA(combat_out, npcs=30, verbose=F,features=common_genes)
combat_out <- RunTSNE(combat_out,reduction='pca',dims=1:30)
combat_out <- RunUMAP(combat_out,reduction='pca',dims=1:30)

combat_pca <- Embeddings(combat_out[['pca']])
combat_tsne <- Embeddings(combat_out[['tsne']])
combat_umap <- Embeddings(combat_out[['umap']])

temp_output <- paste0(output,'HVG_',hvg_num[i],'/combat')
dir.create(temp_output,recursive=T)
write.csv(combat_pca,file=paste0(temp_output,'/combat_pca.csv'),row.names=T)
write.csv(combat_tsne,file=paste0(temp_output,'/combat_tsne.csv'),row.names=T)
write.csv(combat_umap,file=paste0(temp_output,'/combat_umap.csv'),row.names=T)

pdf(paste(temp_output,"/fig_combat_tsne.pdf",sep=""))
par(mar=c(4.5,5,0.5,0.5))
p <- DimPlot(combat_out, reduction = "tsne", group.by = "batch",cols=col_type)
print(p)
dev.off()

pdf(paste(temp_output,"/fig_combat_umap.pdf",sep=""))
par(mar=c(4.5,5,0.5,0.5))
p <- DimPlot(combat_out, reduction = "umap", group.by = "batch",cols=col_type)
print(p)
dev.off()



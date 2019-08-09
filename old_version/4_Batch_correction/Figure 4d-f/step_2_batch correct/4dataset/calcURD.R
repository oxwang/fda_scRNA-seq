library(Seurat)
library(scran)
library(DropletUtils)
library(scater)

library(URD)

sessionInfo()

seurat_preprocess <- function(da)
{
  da <- da
  da <- CreateSeuratObject(raw.data=da,min.cells=3,min.genes=200)
  da <- NormalizeData(da)
  da <- ScaleData(da,display.progress=F)
  da <- FindVariableGenes(da,do.plot=F)
  return(da)
}

setwd('FilteredNew/4dataset/')
load('./FilteredSamples/LLU-LLU/HCC1395(A)/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Merge_LLU_HCC1395(A)", gene_counts@Dimnames[[2]], sep="_")
sample1 <- seurat_preprocess(gene_counts)
sample1@meta.data$sample <- "10X_LLU_A"

load('./FilteredSamples/LLU-NCI/HCC1395(A)/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Seq4_LLU_NCI_HCC1395(A)", gene_counts@Dimnames[[2]], sep="_")
sample2 <- seurat_preprocess(gene_counts)
sample2@meta.data$sample <- "10X_NCI_A"

load('./FilteredSamples/LLU-LLU/HCC1395BL(B)/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Merge_LLU_HCC1395BL(B)", gene_counts@Dimnames[[2]], sep="_")
sample3 <- seurat_preprocess(gene_counts)
sample3@meta.data$sample <- "10X_LLU_B"

load('./FilteredSamples/LLU-NCI/HCC1395BL(B)/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Seq4_LLU_NCI_HCC1395BL(B)", gene_counts@Dimnames[[2]], sep="_")
sample4 <- seurat_preprocess(gene_counts)
sample4@meta.data$sample <- "10X_NCI_B"

merged <- MergeSeurat(sample1, sample2)
merged <- MergeSeurat(merged, sample3)
merged <- MergeSeurat(merged, sample4)
merged <- NormalizeData(merged)
merged <- ScaleData(merged)
test <- createURD(count.data = merged@data)
test <- calcPCA(test, mp.factor = 2)
sum(test@pca.sig)
write.table(sum(test@pca.sig), "URD-seurat.txt")



load('./FilteredSamples/LLU-LLU/HCC1395(A)/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Merge_LLU_HCC1395(A)", gene_counts@Dimnames[[2]], sep="_")
sample1 <- CreateSeuratObject(gene_counts)
sample1@meta.data$sample <- "10X_LLU_A"

load('./FilteredSamples/LLU-NCI/HCC1395(A)/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Seq4_LLU_NCI_HCC1395(A)", gene_counts@Dimnames[[2]], sep="_")
sample2 <- CreateSeuratObject(gene_counts)
sample2@meta.data$sample <- "10X_NCI_A"

load('./FilteredSamples/LLU-LLU/HCC1395BL(B)/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Merge_LLU_HCC1395BL(B)", gene_counts@Dimnames[[2]], sep="_")
sample3 <- CreateSeuratObject(gene_counts)
sample3@meta.data$sample <- "10X_LLU_B"

load('./FilteredSamples/LLU-NCI/HCC1395BL(B)/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Seq4_LLU_NCI_HCC1395BL(B)", gene_counts@Dimnames[[2]], sep="_")
sample4 <- CreateSeuratObject(gene_counts)
sample4@meta.data$sample <- "10X_NCI_B"

merged <- MergeSeurat(sample1, sample2)
merged <- MergeSeurat(merged, sample3)
merged <- MergeSeurat(merged, sample4)
test <- createURD(count.data = merged@raw.data)
test <- calcPCA(test, mp.factor = 2)
sum(test@pca.sig)
write.table(sum(test@pca.sig), "URD-raw.txt")



setwd('FilteredNew/4dataset/mnnCorrect/')
load("combined_mnnCorrectAdjustedObject.Robj")
test <- createURD(count.data = assay(combined, "mnnCounts"))
test <- calcPCA(test, mp.factor = 2)
sum(test@pca.sig)
write.table(sum(test@pca.sig), "URD.txt")

setwd('../ComBat')
result <- readRDS('combat.rds')
test <- createURD(count.data = assay(result, "combatCorrect"))
test <- calcPCA(test, mp.factor = 2)
sum(test@pca.sig)
write.table(sum(test@pca.sig), "URD.txt")

setwd('../limma')
result <- readRDS('limma.rds')
test <- createURD(count.data = assay(result, "limmaCorrect"))
test <- calcPCA(test, mp.factor = 2)
sum(test@pca.sig)
write.table(sum(test@pca.sig), "URD.txt")

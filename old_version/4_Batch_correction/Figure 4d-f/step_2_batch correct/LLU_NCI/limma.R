library(Seurat)
library(scater)
library(scran)
library(limma)

sessionInfo()

numVarGenes = 1000
pcs <- read.table('URD.txt')

scran_preprocess <- function(da,batch)
{
  da <- da
  batch <- batch
  sample_id <- colnames(da)
  da_meta <- data.frame(Sample=sample_id,Batch=batch)
  #sce <- da
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

batch <- c("10X_LLU_B", "10X_LLU_Spikein_10%A", "10X_LLU_A", "10X_NCI_A", "10X_NCI_B", "10X_NCI_Spikein_5%A", "10X_NCI_Spikein_5%A_Fix1")
counts <- vector("list",length(batch))

load('FilteredSamples/LLU-LLU/HCC1395BL(B)/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Merge_LLU_HCC1395BL(B)", gene_counts@Dimnames[[2]], sep="_")
LLU_Tumor <- CreateSeuratObject(gene_counts)
LLU_Tumor@meta.data$sample <- "10X_LLU_B"

load('FilteredSamples/LLU-LLU/Spikein_10/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Merge_LLU_Spikein", gene_counts@Dimnames[[2]], sep="_")
LLU_Mix <- CreateSeuratObject(gene_counts)
LLU_Mix@meta.data$sample <- "10X_LLU_Spikein_10%A"

load('FilteredSamples/LLU-LLU/HCC1395(A)/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Merge_LLU_HCC1395(A)", gene_counts@Dimnames[[2]], sep="_")
LLU_Normal <- CreateSeuratObject(gene_counts)
LLU_Normal@meta.data$sample <- "10X_LLU_A"

load('FilteredSamples/LLU-NCI/HCC1395(A)/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("MergeNCI_HCC1395(A)", gene_counts@Dimnames[[2]], sep="_")
NCI_Normal <- CreateSeuratObject(gene_counts)
NCI_Normal@meta.data$sample <- "10X_NCI_A"

load('FilteredSamples/LLU-NCI/HCC1395BL(B)/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("MergeNCI_HCC1395BL(B)", gene_counts@Dimnames[[2]], sep="_")
NCI_Tumor <- CreateSeuratObject(gene_counts)
NCI_Tumor@meta.data$sample <- "10X_NCI_B"

load('FilteredSamples/LLU-NCI/Spikein_5/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Merge_NCI_Spikein", gene_counts@Dimnames[[2]], sep="_")
NCI_Mix <- CreateSeuratObject(gene_counts)
NCI_Mix@meta.data$sample <- "10X_NCI_Spikein_5%A"

load('FilteredSamples/LLU-NCI/Spikein_5_Fix1/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Merge_NCI_Spikein_Fix1", gene_counts@Dimnames[[2]], sep="_")
NCI_Fix <- CreateSeuratObject(gene_counts)
NCI_Fix@meta.data$sample <- "10X_NCI_Spikein_5%A_Fix1"

setwd('FilteredNew/LLU_NCI/limma/')

merged <- MergeSeurat(LLU_Tumor, LLU_Mix)
merged <- MergeSeurat(merged, LLU_Normal)
merged <- MergeSeurat(merged, NCI_Normal)
merged <- MergeSeurat(merged, NCI_Tumor)
merged <- MergeSeurat(merged, NCI_Mix)
merged <- MergeSeurat(merged, NCI_Fix)

#Steps to take
#Load in data, create Seurat objects, merge using Seurat, then create Single Cell Experiment objects

batch <- c("10X_LLU_B", "10X_LLU_Spikein_10%A", "10X_LLU_A", "10X_NCI_A", "10X_NCI_B", "10X_NCI_Spikein_5%A", "10X_NCI_Spikein_5%A_Fix1")

#How to create single cell experiment
index <- grep("^10X_LLU_B$", merged@meta.data$sample)
counts[[1]] <- merged@raw.data[,index]

index <- grep("^10X_LLU_Spikein_10%A$", merged@meta.data$sample)
counts[[2]] <- merged@raw.data[,index]

index <- grep("^10X_LLU_A$", merged@meta.data$sample)
counts[[3]] <- merged@raw.data[,index]

index <- grep("^10X_NCI_A$", merged@meta.data$sample)
counts[[4]] <- merged@raw.data[,index]

index <- grep("^10X_NCI_B$", merged@meta.data$sample)
counts[[5]] <- merged@raw.data[,index]

index <- grep("^10X_NCI_Spikein_5%A$", merged@meta.data$sample)
counts[[6]] <- merged@raw.data[,index]

index <- grep("^10X_NCI_Spikein_5%A_Fix1$", merged@meta.data$sample)
counts[[7]] <- merged@raw.data[,index]

dim(counts[[1]])
dim(counts[[2]])
dim(counts[[3]])
dim(counts[[4]])
dim(counts[[5]])
dim(counts[[6]])
dim(counts[[7]])

sce <- vector("list",length(batch))
dec <- vector("list",length(batch))
temp_gene <- vector("list",length(batch))
for(i in 1:length(batch))
{
  temp <- scran_preprocess(as(counts[[i]], "dgCMatrix"),batch[i])
  sce[[i]] <- temp[[1]]
  dec[[i]] <- temp[[2]]
  temp_gene[[i]] <- rownames(dec[[i]])
}

decomp <- list()
hvgs <- list()
genes.use <- c()
for (i in 1:length(sce)) {
  fit <- trendVar(sce[[i]], use.spikes=FALSE)
  decomp[[i]] <- decomposeVar(sce[[i]], fit)
  hvgs[[i]] <- order(decomp[[i]]$bio, decreasing=TRUE)
  genes.use <- union(genes.use, rownames(decomp[[i]][hvgs[[i]][1:numVarGenes],]))
}
print(length(genes.use))

save(decomp, hvgs, genes.use, file="hvgs_data.Robj")

temp <- vector("list",length(batch))
for(i in 1:length(batch)) {
  da <- assay(sce[[i]], "counts")[genes.use,]
  da2 <- assay(sce[[i]], "logcounts")[genes.use,]
  da_meta <- data.frame(Sample=colnames(da),Batch=batch[i])
  temp[[i]] <- SingleCellExperiment(list(counts=da, logcounts=da2),colData=da_meta,rowData=DataFrame(Symbol=rownames(da)))
}

dim(temp[[1]])
dim(temp[[2]])
dim(temp[[3]])
dim(temp[[4]])
dim(temp[[5]])
dim(temp[[6]])
dim(temp[[7]])

combined <- cbind(temp[[1]],
                  temp[[2]],
                  temp[[3]],
                  temp[[4]],
                  temp[[5]],
                  temp[[6]],
                  temp[[7]])

limma_result <- removeBatchEffect(as.matrix(assay(combined, "logcounts")), factor(combined@colData$Batch))

assay(combined, "limmaCorrect") <- as(limma_result, "dgCMatrix")
saveRDS(combined, 'limma.rds')

combined <- runTSNE(combined, exprs_values="counts", n_dimred=pcs)
pdf('./counts.pdf')
plotTSNE(combined, colour_by="Batch")
dev.off()

combined <- runTSNE(combined, exprs_values="limmaCorrect", n_dimred=pcs)
pdf('./limma_factor-sample.pdf')
plotTSNE(combined, colour_by="Batch")
dev.off()

saveRDS(combined, 'limma.rds')

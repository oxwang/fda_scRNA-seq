library(Seurat)
library(scran)
library(DropletUtils)
library(scater)

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

batch <- c("10X_LLU_Spikein_10%A", "10X_LLU_A")
counts <- vector("list",length(batch))

load('FilteredSamples/LLU-LLU/Spikein_10/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Merge_LLU_Spikein", gene_counts@Dimnames[[2]], sep="_")
LLU_Mix <- CreateSeuratObject(gene_counts)
LLU_Mix@meta.data$sample <- "10X_LLU_Spikein_10%A"

load('FilteredSamples/LLU-LLU/HCC1395(A)/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Merge_LLU_HCC1395(A)", gene_counts@Dimnames[[2]], sep="_")
LLU_Normal <- CreateSeuratObject(gene_counts)
LLU_Normal@meta.data$sample <- "10X_LLU_A"

setwd('FilteredNew/2dataset-tumor/mnnCorrect/')

merged <- MergeSeurat(LLU_Mix, LLU_Normal)

#Steps to take
#Load in data, create Seurat objects, merge using Seurat, then create Single Cell Experiment objects

#How to create single cell experiment
index <- grep("10X_LLU_Spikein_10%A", merged@meta.data$sample)
counts[[1]] <- merged@raw.data[,index]

index <- grep("10X_LLU_A", merged@meta.data$sample)
counts[[2]] <- merged@raw.data[,index]

sce <- vector("list",length(batch))
dec <- vector("list",length(batch))
temp_gene <- vector("list",length(batch))
for(i in 1:length(batch))
{
  temp <- scran_preprocess(as(counts[[i]], "dgCMatrix"),batch[i])
  #temp <- scran_preprocess(batch[[i]])
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

out <- mnnCorrect(assay(sce[[1]], "logcounts"),
                  assay(sce[[2]], "logcounts"),
                  subset.row = genes.use)

save(out, file="out_mnnCorrectObject.Robj")

assay(sce[[1]], "mnnCounts") <- out$corrected[[1]]
assay(sce[[2]], "mnnCounts") <- out$corrected[[2]]

combined <- cbind(sce[[1]],
                  sce[[2]])
save(combined, file = "combined_mnnCorrectAdjustedObject.Robj")

pdf('PCA_mnnCount.pdf')
combined <- runPCA(combined, exprs_values = "mnnCounts")
print(plotReducedDim(combined, use_dimred = "PCA",
                     colour_by = "Batch", shape_by = "Batch"))
dev.off()

pdf('PCA_logCounts.pdf')
combined <- runPCA(combined, exprs_values = "logcounts")
print(plotReducedDim(combined, use_dimred = "PCA",
                     colour_by = "Batch", shape_by = "Batch"))
dev.off()

pdf('TSNE_logCounts.pdf')
combined <- runTSNE(combined, exprs_values = "logcounts", n_dimred=pcs)
print(plotTSNE(combined, colour_by = "Batch", shape_by = "Batch"))
dev.off()

pdf('TSNE_mnnCounts.pdf')
combined <- runTSNE(combined, exprs_values = "mnnCounts", n_dimred=pcs)
print(plotTSNE(combined, colour_by = "Batch", shape_by = "Batch"))
dev.off()

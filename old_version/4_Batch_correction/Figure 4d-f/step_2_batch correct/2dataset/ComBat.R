library(Seurat)
library(scater)
library(scran)
library(sva)

sessionInfo()

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

batch <- c("10X_LLU_Spikein_10%A", "10X_LLU_B")
counts <- vector("list",length(batch))

load('FilteredSamples/LLU-LLU/Spikein_10/gene_counts_cellranger.rdata')
counts[[1]] <- gene_counts
gene_counts@Dimnames[[2]] <- paste("Merge_LLU_Spikein", gene_counts@Dimnames[[2]], sep="_")
LLU_Mix <- CreateSeuratObject(gene_counts)
LLU_Mix@meta.data$sample <- "10X_LLU_Spikein_10%A"

load('FilteredSamples/LLU-LLU/HCC1395BL(B)/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Merge_LLU_HCC1395BL(B)", gene_counts@Dimnames[[2]], sep="_")
LLU_Normal <- CreateSeuratObject(gene_counts)
LLU_Normal@meta.data$sample <- "10X_LLU_B"

setwd('FilteredNew/2dataset/ComBat/')

merged <- MergeSeurat(LLU_Mix, LLU_Normal)

#Steps to take
#Load in data, create Seurat objects, merge using Seurat, then create Single Cell Experiment objects

#How to create single cell experiment
index <- grep("10X_LLU_Spikein_10%A", merged@meta.data$sample)
counts[[1]] <- merged@raw.data[,index]

index <- grep("10X_LLU_B", merged@meta.data$sample)
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

temp <- vector("list",length(batch))
for(i in 1:length(batch)) {
  da <- assay(sce[[i]], "counts")[genes.use,]
  da2 <- assay(sce[[i]], "logcounts")[genes.use,]
  da_meta <- data.frame(Sample=colnames(da),Batch=batch[i])
  temp[[i]] <- SingleCellExperiment(list(counts=da, logcounts=da2),colData=da_meta,rowData=DataFrame(Symbol=rownames(da)))
}

dim(temp[[1]])
dim(temp[[2]])

combined <- cbind(temp[[1]],
                  temp[[2]])

combined <- runTSNE(combined, exprs_values="counts", n_dimred=pcs)
pdf('./counts.pdf')
plotTSNE(combined, colour_by="Batch", shape_by="Batch")
dev.off()

combat_result <- ComBat(as.matrix(assay(combined, "logcounts")), factor(combined@colData$Batch), mod=NULL, prior.plots=F)

assay(combined, "combatCorrect") <- as(combat_result, "dgCMatrix")

saveRDS(combined, 'combat.rds')

combined <- runTSNE(combined, exprs_values="combatCorrect", n_dimred=pcs)
pdf('./combat_batch-sample_mod-none.pdf')
plotTSNE(combined, colour_by="Batch", shape_by="Batch")
dev.off()

saveRDS(combined, 'combat.rds')

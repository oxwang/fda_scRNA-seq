library(Seurat)
library(grid)
sessionInfo()
#setwd('/Volumes/ccr-sf/static/illumina/NGS_SHARE/BSG/SEQC2/SGC/NC/10x_Genomics_merged/')

library(RColorBrewer)
n <- 7
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

numVarGenes = 1000

seurat_preprocess <- function(da)
{
  da <- da
  da <- CreateSeuratObject(raw.data=da,min.cells=3,min.genes=200)
  da <- NormalizeData(da)
  da <- ScaleData(da,display.progress=F)
  da <- FindVariableGenes(da,do.plot=F)
  return(da)
}

setwd('/is2/projects/CCR-SF/static/illumina/NGS_SHARE/BSG/SEQC2/SGC/Biogenlink/correction/')

load('./FilteredSamples/LLU-LLU/HCC1395BL(B)/gene_counts_cellranger.rdata')
sample1 <- CreateSeuratObject(gene_counts, project="LLU_LLU_HCC1395BL(B)")
sample1 <- seurat_preprocess(gene_counts)
sample1@meta.data$sample <- "10X_LLU_B"

load('./FilteredSamples/LLU-LLU/Spikein_10/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Merge_LLU_Spikein", gene_counts@Dimnames[[2]], sep="_")
sample2 <- seurat_preprocess(gene_counts)
sample2@meta.data$sample <- "10X_LLU_Spikein_10%A"

load('./FilteredSamples/LLU-LLU/HCC1395(A)/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Merge_LLU_HCC1395(A)", gene_counts@Dimnames[[2]], sep="_")
sample3 <- seurat_preprocess(gene_counts)
sample3@meta.data$sample <- "10X_LLU_A"

load('./FilteredSamples/LLU-NCI/HCC1395(A)/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Seq4_LLU_NCI_HCC1395(A)", gene_counts@Dimnames[[2]], sep="_")
sample4 <- seurat_preprocess(gene_counts)
sample4@meta.data$sample <- "10X_NCI_A"

load('./FilteredSamples/LLU-NCI/HCC1395BL(B)/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Seq4_LLU_NCI_HCC1395BL(B)", gene_counts@Dimnames[[2]], sep="_")
sample5 <- seurat_preprocess(gene_counts)
sample5@meta.data$sample <- "10X_NCI_B"

load('./FilteredSamples/LLU-NCI/Spikein_5/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Merge_NCI_Spikein", gene_counts@Dimnames[[2]], sep="_")
sample6 <- seurat_preprocess(gene_counts)
sample6@meta.data$sample <- "10X_NCI_Spikein_5%A"

load('./FilteredSamples/LLU-NCI/Spikein_5_Fix1/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Merge_NCI_Spikein_Fix1", gene_counts@Dimnames[[2]], sep="_")
sample7 <- seurat_preprocess(gene_counts)
sample7@meta.data$sample <- "10X_NCI_Spikein_5%A_Fix1"

setwd('./FilteredNew/LLU_NCI/CCA_sean-pcs/')

pcs <- read.table('URD.txt')

# lastly, we set the 'identity' in each dataset for easy identification
# later it will be transferred to the merged object in RunCCA
# sample1@meta.data[, "sample"] <- "Merge_LLU_HCC1395(A)"
# sample2@meta.data[, "sample"] <- "Seq4_LLU_NCI_HCC1395(A)"

ob.list <- list(sample1, sample2, sample3, sample4, sample5, sample6, sample7)
genes.use <- c()
for (i in 1:length(ob.list)) {
	genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), numVarGenes))
}
genes.use <- unique(genes.use)
print(length(genes.use))
for (i in 1:length(ob.list)) {
	genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
}
print(length(genes.use))

#Run multi-set CCA
samples.integrated <- RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = pcs[[1]])

p1 <- DimPlot(object = samples.integrated, reduction.use = "cca", group.by = "sample", pt.size = 0.5,
              do.return = TRUE)
p2 <- VlnPlot(object = samples.integrated, features.plot = "CC1", group.by = "sample", do.return = TRUE)
pdf("Prior_CCA_and_Violin.pdf")
plot_grid(p1, p2)
dev.off()

pdf("Dim_Heatmap.pdf")
DimHeatmap(object = samples.integrated, reduction.type = "cca", dim.use = 1:9,
           do.balanced = TRUE)
DimHeatmap(object = samples.integrated, reduction.type = "cca", dim.use = 10:18,
           do.balanced = TRUE)
dev.off()


#CC Selection
MetageneBicorPlot(samples.integrated, grouping.var = "sample", dims.eval = 1:pcs[[1]])

#Run rare non-overlapping filtering
samples.integrated <- CalcVarExpRatio(object = samples.integrated, reduction.type = "pca",
                                      grouping.var = "sample", dims.use = 1:pcs[[1]])
samples.all <- samples.integrated
samples.integrated <- SubsetData(samples.integrated, subset.name = "var.ratio.pca",
                                 accept.low = 0.5)

#Alignment
samples.integrated <- AlignSubspace(samples.integrated, reduction.type = "cca",
                                    grouping.var = "sample", dims.align = 1:pcs[[1]])

#Want to compare with not discarding? How many cells do we lose in this dataset?
p1 <- VlnPlot(object = samples.integrated, features.plot = "ACC1", group.by = "sample",
              do.return = TRUE)
p2 <- VlnPlot(object = samples.integrated, features.plot = "ACC2", group.by = "sample",
              do.return = TRUE)
pdf("Post_alignment_CCA.pdf")
plot_grid(p1, p2)
dev.off()

#t-SNE and Clustering
#Cluster labels saved in @ident
samples.integrated <- FindClusters(samples.integrated, reduction.type = "cca.aligned",
                                   dims.use = 1:pcs[[1]], save.SNN = T, resolution = 0.4)
samples.integrated <- RunTSNE(samples.integrated, reduction.use = "cca.aligned",
                              dims.use = 1:pcs[[1]])

save(samples.integrated, file ='./seurat_CCA_object.Robj')

#Visualization
pdf(paste('Post_alignment_TSNEPlotwith', pcs[[1]], 'PCs_Samples.pdf', sep=""))
TSNEPlot(samples.integrated, group.by='sample')
dev.off()
pdf(paste('Post_alignment_TSNEPlotwith', pcs[[1]], 'PCs.pdf', sep=""))
TSNEPlot(samples.integrated)
dev.off()

#No Discard
samples.all <- AlignSubspace(object = samples.all, reduction.type = "cca", grouping.var = "sample",
                             dims.align = 1:pcs[[1]])
p1 <- VlnPlot(object = samples.all, features.plot = "ACC1", group.by = "sample",
              do.return = TRUE)
p2 <- VlnPlot(object = samples.all, features.plot = "ACC2", group.by = "sample",
              do.return = TRUE)
pdf("No_discard_Post_alignment_CCA.pdf")
plot_grid(p1, p2)
dev.off()

#TSNE
samples.all <- RunTSNE(object = samples.all, reduction.use = "cca.aligned", dims.use = 1:pcs[[1]],
                       do.fast = TRUE)
samples.all <- FindClusters(object = samples.all, reduction.type = "cca.aligned", dims.use = 1:pcs[[1]],
                            save.SNN = TRUE)
pdf(paste('No_discard_Post_alignment_TSNEPlotwith', pcs[[1]], 'PCs_Samples.pdf', sep=""))
TSNEPlot(samples.all, group.by='sample')
dev.off()
pdf(paste('No_discard_Post_alignment_TSNEPlotwith', pcs[[1]], 'PCs.pdf', sep=""))
TSNEPlot(samples.all)
dev.off()

save(samples.all, file ='./seurat_CCA_no_discard_object.Robj')

temp <- TSNEPlot(samples.all, group.by="sample") + scale_colour_manual(values = col_vector)
pdf(paste('No_discard_Post_alignment_TSNEPlotwith', pcs[[1]], 'PCs_Samples.pdf', sep=""))
mylegend <- g_legend(temp)
grid.draw(mylegend)
dev.off()
temp <- TSNEPlot(samples.all, group.by="sample") + scale_colour_manual(values = col_vector) + theme(legend.position="none")
pdf(paste('No_discard_Post_alignment_TSNEPlotwith', pcs[[1]], 'PCs.pdf', sep=""))
temp
dev.off()

temp <- TSNEPlot(samples.all, group.by="sample") + scale_colour_manual(values = col_vector)
pdf(paste('No_discard_Post_alignment_TSNEPlotwith', pcs[[1]], 'PCs_Samples_Legend.pdf', sep=""))
mylegend <- g_legend(temp)
grid.draw(mylegend)
dev.off()
temp <- TSNEPlot(samples.all, group.by="sample") + scale_colour_manual(values = col_vector) + theme(legend.position="none")
pdf(paste('No_discard_Post_alignment_TSNEPlotwith', pcs[[1]], 'PCs_Samples_NoLegend.pdf', sep=""))
temp
dev.off()

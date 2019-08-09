library(Seurat)
library(grid)
sessionInfo()

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

load('./FilteredSamples/LLU-LLU/Spikein_10/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Merge_LLU_Spikein", gene_counts@Dimnames[[2]], sep="_")
sample1 <- seurat_preprocess(gene_counts)
sample1@meta.data$sample <- "10X_LLU_Spikein_10%A"

load('./FilteredSamples/LLU-LLU/HCC1395(A)/gene_counts_cellranger.rdata')
gene_counts@Dimnames[[2]] <- paste("Merge_LLU_HCC1395(A)", gene_counts@Dimnames[[2]], sep="_")
sample2 <- seurat_preprocess(gene_counts)
sample2@meta.data$sample <- "10X_LLU_A"

setwd('./FilteredNew/2dataset-tumor/CCA_sean-pcs/')

pcs <- read.table('URD.txt')

# lastly, we set the 'identity' in each dataset for easy identification
# later it will be transferred to the merged object in RunCCA
# sample1@meta.data[, "sample"] <- "Merge_LLU_HCC1395(A)"
# sample2@meta.data[, "sample"] <- "Seq4_LLU_NCI_HCC1395(A)"

ob.list <- list(sample1, sample2)
genes.use <- c()
for (i in 1:length(ob.list)) {
	genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), numVarGenes))
}
genes.use <- unique(genes.use)
print(length(genes.use))
for (i in 1:length(ob.list)) {
	genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
}
hvg.llu <- rownames(x = head(x = sample1@hvg.info, n = numVarGenes))
hvg.nci <- rownames(x = head(x = sample2@hvg.info, n = numVarGenes))
hvg.union <- union(x = hvg.llu, y = hvg.nci)
print(length(hvg.union))
print(length(genes.use))

#Run multi-set CCA
samples.integrated <- RunCCA(sample1, sample2, genes.use = hvg.union, num.cc = pcs[[1]])

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

temp <- TSNEPlot(samples.all, group.by="sample") + scale_colour_manual(values = c("darkturquoise", "darkorchid", "hotpink", "firebrick", "seagreen4"))
pdf(paste('No_discard_Post_alignment_TSNEPlotwith', pcs[[1]], 'PCs_Samples_Legend.pdf', sep=""))
mylegend <- g_legend(temp)
grid.draw(mylegend)
dev.off()
temp <- TSNEPlot(samples.all, group.by="sample") + scale_colour_manual(values = c("darkturquoise", "darkorchid", "hotpink", "firebrick", "seagreen4")) + theme(legend.position="none")
pdf(paste('No_discard_Post_alignment_TSNEPlotwith', pcs[[1]], 'PCs_Samples_NoLegend.pdf', sep=""))
temp
dev.off()

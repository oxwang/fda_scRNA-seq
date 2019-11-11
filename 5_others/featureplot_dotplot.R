## generate featureplot and dotplot based on mnn methods ##

# 1. load data

counts <- vector("list",20)

batch <- c("10x_Mix10_LLU", "10x_Mix5_NCI", "10x_Mix5_F_NCI", "10x_Mix5_NCI_M", "10x_Mix5_F_NCI_M", "10x_Mix5_F2_NCI_M", "10x_LLU_A", "10x_NCI_A", "10x_NCI_M_A", "10x_LLU_B", "10x_NCI_B", "10x_NCI_M_B", "C1_FDA_HT_A", "C1_LLU_A", "ICELL8_SE_A", "ICELL8_PE_A", "C1_FDA_HT_B", "C1_LLU_B", "ICELL8_SE_B", "ICELL8_PE_B")
col_type <- rainbow(20)

input <- output <- '/.../4_Batch_correction/Results/scenario_1/'

for(i in 1:20)
{
	counts[[i]] <- read.table(paste0(input,batch[i],'.txt.gz'),header=T,row.names=1)
}


# 2. create seurat object

library(Seurat)

for(i in 1:length(counts))
{
	meta_da <- data.frame(matrix(batch[i],nrow=ncol(counts[[i]]),ncol=1))
	rownames(meta_da) <- colnames(counts[[i]])
	colnames(meta_da) <- 'batch'
	tmpCells <- CreateSeuratObject(counts[[i]],meta.data=meta_da)
	tmpCells <- NormalizeData(tmpCells,verbose=F)

## generate merged Cells data for MNN
	if(i == 1) {Cells <- tmpCells}
	else {Cells <- merge(Cells,tmpCells)}
	
}

hvg_num <- 2000
mnn_data <- Cells
mnn_data <- FindVariableFeatures(mnn_data,nfeatures=hvg_num)


# 3. fastMNN correction

library(SeuratWrappers)

start_time <- Sys.time()
mnn_data <- RunFastMNN(object.list=SplitObject(mnn_data, split.by = 'batch'),features=hvg_num)
end_time <- Sys.time()		## performance time: 408 seconds
mnn_data <- RunTSNE(mnn_data,reduction='mnn',dims=1:30)
mnn_data <- RunUMAP(mnn_data,reduction='mnn',dims=1:30)

# 4. generate featureplot

Idents(mnn_data) <- mnn_data$batch
DimPlot(mnn_data)
DimPlot(mnn_data, reduction="tsne")

library(org.Hs.eg.db)
genesA1 <-c("TGFBI","AKR1C2","BGN","KRT81","THBS1","APP","NUPR1","LY6K","TM4SF1","SERPINE2")
genesA2 <- mapIds(org.Hs.eg.db, keys=genesA1,keytype="SYMBOL",column="ENSEMBL")
genesB1 <-c("MS4A1","CCR7","IRF4","LSP1","CD40","CD79A","CD53","CCL3")
genesB2 <- mapIds(org.Hs.eg.db, keys=genesB1,keytype="SYMBOL",column="ENSEMBL")

features.plot <- c("ENSG00000101017","ENSG00000026508","ENSG00000141736","ENSG00000146648","ENSG00000091831","ENSG00000136997","ENSG00000100030","ENSG00000142208")

feature_plot_A_tsne <- vector("list", length(genesA1))
feature_plot_A_umap <- vector("list", length(genesA1))
for (i in 1:length(genesA1)) {
  feature_plot_A_tsne[[i]] <- FeaturePlot(object = mnn_data, features = genesA2[i],
                                     cols= c("lightgrey","red"),reduction="tsne",
                                     pt.size = .3,ncol = 4,order =T)
  feature_plot_A_umap[[i]] <- FeaturePlot(object = mnn_data, features = genesA2[i],
                                          cols= c("lightgrey","red"),reduction="umap",
                                          pt.size = .3,ncol = 4,order =T)
}

feature_plot_B_tsne <- vector("list", length(genesB1))
feature_plot_B_umap <- vector("list", length(genesB1))
for (i in 1:length(genesB1)) {
  feature_plot_B_tsne[[i]] <- FeaturePlot(object = mnn_data, features = genesB2[i],
                                          cols= c("lightgrey","blue"),reduction="tsne",
                                          pt.size = .3,ncol = 4,order =T)
  feature_plot_B_umap[[i]] <- FeaturePlot(object = mnn_data, features = genesB2[i],
                                          cols= c("lightgrey","blue"),reduction="umap",
                                          pt.size = .3,ncol = 4,order =T)
}
saveRDS(genesA1, file = "genesA1.rds")
saveRDS(genesB1, file = "genesB1.rds")

saveRDS(feature_plot_A_tsne, file = "feature_plot_A_tsne.rds")
saveRDS(feature_plot_A_umap, file = "feature_plot_A_umap.rds")
saveRDS(feature_plot_B_tsne, file = "feature_plot_B_tsne.rds")
saveRDS(feature_plot_B_umap, file = "feature_plot_B_umap.rds")




# 5. mnnCorrect correction

library(SummarizedExperiment)
library(batchelor)
library(ggplot2)
library(org.Hs.eg.db)

mnn_data <- Cells
mnn_data <- FindVariableFeatures(mnn_data,nfeatures=hvg_num)
temp_hvgs <- VariableFeatures(mnn_data)

batch_cell <- c()
for(i in 1:length(counts))
{batch_cell <- c(batch_cell,ncol(counts[[i]]))}
batch_id <- rep(batch,batch_cell)

mnn_process <- GetAssayData(mnn_data,'data')
start_time <- Sys.time()
mnn_data <- mnnCorrect(as.matrix(mnn_process), batch=factor(batch_id), subset.row=temp_hvgs)
end_time <- Sys.time()		## performance time: 2843 seconds

mnn_corrected <- assay(mnn_data,'corrected')




## 6. generate dot plot for uncorrected data

library(ggplot2)

Idents(Cells) <- 'batch'
Cells <- RenameIdents(object = Cells, '10x_LLU_A'= 'A','10x_NCI_A'= 'A','10x_NCI_M_A'= 'A','C1_FDA_HT_A'= 'A','C1_LLU_A'= 'A','ICELL8_SE_A'= 'A','ICELL8_PE_A'= 'A','10x_LLU_B'= 'B','10x_NCI_B'= 'B','10x_NCI_M_B'= 'B','C1_FDA_HT_B'= 'B','C1_LLU_B'= 'B','ICELL8_SE_B'= 'B','ICELL8_PE_B'= 'B','10x_Mix10_LLU'= 'B','10x_Mix5_NCI'= 'B','10x_Mix5_F_NCI'= 'B','10x_Mix5_NCI_M'= 'B','10x_Mix5_F_NCI_M'= 'B','10x_Mix5_F2_NCI_M'= 'B')
Cells[["CellType"]] <- Idents(object = Cells)
Idents(object = Cells) <- 'batch'

Cells <- RenameIdents(object = Cells, '10x_LLU_A'= '10x-LLU-A','10x_NCI_A'= '10x-NCI-A','10x_NCI_M_A'= '10x-NCI-M-A','C1_FDA_HT_A'= 'C1-FDA-HT-A','C1_LLU_A'= 'C1-LLU-A','ICELL8_SE_A'= 'ICELL8-SE-A','ICELL8_PE_A'= 'ICELL8-PE-A','10x_LLU_B'= '10x-LLU-B','10x_NCI_B'= '10x-NCI-B','10x_NCI_M_B'= '10x-NCI-M-B','C1_FDA_HT_B'= 'C1-FDA-HT-B','C1_LLU_B'= 'C1-LLU-B','ICELL8_SE_B'= 'ICELL8-SE-B','ICELL8_PE_B'= 'ICELL8-PE-B','10x_Mix10_LLU'= '10x-Mix10-LLU','10x_Mix5_NCI'= '10x-Mix5-NCI','10x_Mix5_F_NCI'= '10x-Mix5-F-NCI','10x_Mix5_NCI_M'= '10x-Mix5-NCI-M','10x_Mix5_F_NCI_M'= '10x-Mix5-F-NCI-M','10x_Mix5_F2_NCI_M'= '10x-Mix5-F2-NCI-M')

Idents(Cells) <- factor(Idents(Cells), levels = c( 'C1-FDA-HT-A','C1-LLU-A','ICELL8-SE-A','ICELL8-PE-A','10x-LLU-A','10x-NCI-A','10x-NCI-M-A','10x-Mix10-LLU','10x-Mix5-NCI','10x-Mix5-F-NCI','10x-Mix5-NCI-M','10x-Mix5-F-NCI-M','10x-Mix5-F2-NCI-M','C1-FDA-HT-B','C1-LLU-B','ICELL8-SE-B','ICELL8-PE-B','10x-LLU-B','10x-NCI-B','10x-NCI-M-B'))

ensembl_ID <- row.names(Cells)
Gene_symbl <-mapIds(org.Hs.eg.db, keys=ensembl_ID,keytype="ENSEMBL",column="SYMBOL");
Cells@assays$RNA@data@Dimnames[[1]] <-Gene_symbl
Cells@assays$RNA@counts <-Cells@assays$RNA@data


markers.to.plot <-c("TGFBI","AKR1C2","BGN","KRT81","THBS1","APP","NUPR1","LY6K","TM4SF1","SERPINE2","MS4A1","CCR7","IRF4","LSP1","CD40","CD79A","CD53","CCL3","CD74","CD19")
markers.to.plot2 <-markers.to.plot[markers.to.plot %in% rownames(Cells)]

unc_Dotplot <-DotPlot(Cells, features =markers.to.plot2, cols= c("Blue","Red"), split.by ='CellType', col.max = 2, dot.scale = 12) + theme(axis.text=element_text(size=16,face="bold"), axis.title=element_text(size=14,face="bold"))+ RotatedAxis()

unc_data <- unc_Dotplot$data
save(unc_data,file="/.../5_others/Results/unc_dotplot.Rdata")



## 7. generate dot plot for mnn corrected data

library(ggplot2)

Cells <- CreateSeuratObject(mnn_corrected)
Cells[['batch']] <- batch_id

Idents(Cells) <- 'batch'
Cells <- RenameIdents(object = Cells, '10x_LLU_A'= 'A','10x_NCI_A'= 'A','10x_NCI_M_A'= 'A','C1_FDA_HT_A'= 'A','C1_LLU_A'= 'A','ICELL8_SE_A'= 'A','ICELL8_PE_A'= 'A','10x_LLU_B'= 'B','10x_NCI_B'= 'B','10x_NCI_M_B'= 'B','C1_FDA_HT_B'= 'B','C1_LLU_B'= 'B','ICELL8_SE_B'= 'B','ICELL8_PE_B'= 'B','10x_Mix10_LLU'= 'B','10x_Mix5_NCI'= 'B','10x_Mix5_F_NCI'= 'B','10x_Mix5_NCI_M'= 'B','10x_Mix5_F_NCI_M'= 'B','10x_Mix5_F2_NCI_M'= 'B')
Cells[["CellType"]] <- Idents(object = Cells)
Idents(object = Cells) <- 'batch'

Cells <- RenameIdents(object = Cells, '10x_LLU_A'= '10x-LLU-A','10x_NCI_A'= '10x-NCI-A','10x_NCI_M_A'= '10x-NCI-M-A','C1_FDA_HT_A'= 'C1-FDA-HT-A','C1_LLU_A'= 'C1-LLU-A','ICELL8_SE_A'= 'ICELL8-SE-A','ICELL8_PE_A'= 'ICELL8-PE-A','10x_LLU_B'= '10x-LLU-B','10x_NCI_B'= '10x-NCI-B','10x_NCI_M_B'= '10x-NCI-M-B','C1_FDA_HT_B'= 'C1-FDA-HT-B','C1_LLU_B'= 'C1-LLU-B','ICELL8_SE_B'= 'ICELL8-SE-B','ICELL8_PE_B'= 'ICELL8-PE-B','10x_Mix10_LLU'= '10x-Mix10-LLU','10x_Mix5_NCI'= '10x-Mix5-NCI','10x_Mix5_F_NCI'= '10x-Mix5-F-NCI','10x_Mix5_NCI_M'= '10x-Mix5-NCI-M','10x_Mix5_F_NCI_M'= '10x-Mix5-F-NCI-M','10x_Mix5_F2_NCI_M'= '10x-Mix5-F2-NCI-M')

Idents(Cells) <- factor(Idents(Cells), levels = c( 'C1-FDA-HT-A','C1-LLU-A','ICELL8-SE-A','ICELL8-PE-A','10x-LLU-A','10x-NCI-A','10x-NCI-M-A','10x-Mix10-LLU','10x-Mix5-NCI','10x-Mix5-F-NCI','10x-Mix5-NCI-M','10x-Mix5-F-NCI-M','10x-Mix5-F2-NCI-M','C1-FDA-HT-B','C1-LLU-B','ICELL8-SE-B','ICELL8-PE-B','10x-LLU-B','10x-NCI-B','10x-NCI-M-B'))


ensembl_ID <- row.names(Cells)
Gene_symbl <-mapIds(org.Hs.eg.db, keys=ensembl_ID,keytype="ENSEMBL",column="SYMBOL");
Cells@assays$RNA@data@Dimnames[[1]] <-Gene_symbl
Cells@assays$RNA@counts <-Cells@assays$RNA@data


markers.to.plot <-c("TGFBI","AKR1C2","BGN","KRT81","THBS1","APP","NUPR1","LY6K","TM4SF1","SERPINE2","MS4A1","CCR7","IRF4","LSP1","CD40","CD79A","CD53","CCL3","CD74","CD19")
markers.to.plot2 <-markers.to.plot[markers.to.plot %in% rownames(Cells)]

mnn_Dotplot <-DotPlot(Cells, features =markers.to.plot2, cols= c("Blue","Red"), split.by ='CellType', col.max = 2, dot.scale = 12) + theme(axis.text=element_text(size=16,face="bold"), axis.title=element_text(size=14,face="bold"))+ RotatedAxis()

mnn_data <- mnn_Dotplot$data
save(mnn_data,file="/.../5_others/Results/mnn_dotplot.Rdata")






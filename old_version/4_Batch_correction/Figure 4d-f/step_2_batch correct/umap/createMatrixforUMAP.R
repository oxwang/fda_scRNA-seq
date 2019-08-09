setwd('FilteredNew/2dataset/')

load('FilteredNew/2dataset/CCA_sean/seurat_CCA_no_discard_object.Robj')
write.csv(samples.all@dr$cca.aligned@cell.embeddings, 'FilteredNew/2dataset/CCA_sean/cca_aligned.csv')

load('CCA/seurat_CCA_no_discard_object.Robj')
write.csv(samples.all@dr$cca.aligned@cell.embeddings, 'CCA/cca_aligned.csv')

load('mnnCorrect/combined_mnnCorrectAdjustedObject.Robj')
write.csv(as.matrix(assay(combined, "mnnCounts")), 'mnnCorrect/mnnCounts.csv')

result <- readRDS('ComBat/combat.rds')
write.csv(as.matrix(assay(result, "combatCorrect")), 'combat/combat.csv')

result <- readRDS('limma/limma.rds')
write.csv(as.matrix(assay(result, "limmaCorrect")), 'limma/limma.csv')


setwd('FilteredNew/2dataset-tumor/')
load('CCA_sean/seurat_CCA_no_discard_object.Robj')
write.csv(samples.all@dr$cca.aligned@cell.embeddings, 'CCA_sean/cca_aligned.csv')

load('mnnCorrect/combined_mnnCorrectAdjustedObject.Robj')
write.csv(as.matrix(assay(combined, "mnnCounts")), 'mnnCorrect/mnnCounts.csv')

result <- readRDS('ComBat/combat.rds')
write.csv(as.matrix(assay(result, "combatCorrect")), 'ComBat/combat.csv')

result <- readRDS('limma/limma.rds')
write.csv(as.matrix(assay(result, "limmaCorrect")), 'limma/limma.csv')


setwd('FilteredNew/4dataset/')
load('CCA_sean/seurat_CCA_no_discard_object.Robj')
write.csv(samples.all@dr$cca.aligned@cell.embeddings, 'CCA_sean/cca_aligned.csv')

load('mnnCorrect/combined_mnnCorrectAdjustedObject.Robj')
write.csv(as.matrix(assay(combined, "mnnCounts")), 'mnnCorrect/mnnCounts.csv')

result <- readRDS('ComBat/combat.rds')
write.csv(as.matrix(assay(result, "combatCorrect")), 'combat/combat.csv')

result <- readRDS('limma/limma.rds')
write.csv(as.matrix(assay(result, "limmaCorrect")), 'limma/limma.csv')

setwd('FilteredNew/')
load('2dataset/CCA_sean-pcs/seurat_CCA_no_discard_object.Robj')
write.csv(samples.all@dr$cca.aligned@cell.embeddings, '2dataset/CCA_sean-pcs/cca_aligned.csv')

load('2dataset-tumor/CCA_sean-pcs/seurat_CCA_no_discard_object.Robj')
write.csv(samples.all@dr$cca.aligned@cell.embeddings, '2dataset-tumor/CCA_sean-pcs/cca_aligned.csv')

load('4dataset/CCA_sean-pcs/seurat_CCA_no_discard_object.Robj')
write.csv(samples.all@dr$cca.aligned@cell.embeddings, '4dataset/CCA_sean-pcs/cca_aligned.csv')

load('LLU_NCI/CCA_sean-pcs/seurat_CCA_no_discard_object.Robj')
write.csv(samples.all@dr$cca.aligned@cell.embeddings, 'LLU_NCI/CCA_sean-pcs/cca_aligned.csv')

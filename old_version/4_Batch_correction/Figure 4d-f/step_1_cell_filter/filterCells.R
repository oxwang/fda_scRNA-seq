library(Seurat)
sessionInfo()

numVarGenes = 1000

cellFilter <- function(gene_counts) {
  par(mfrow=c(1,2))
  print(dim(gene_counts))
  total_mrna <- colSums(gene_counts)
  hist(total_mrna, breaks = 20)
  print(sum(total_mrna < 1e6))
  
  upper_bound <- 10^(mean(log10(total_mrna)) +
                       2*sd(log10(total_mrna)))
  lower_bound <- 10^(mean(log10(total_mrna)) -
                       2*sd(log10(total_mrna)))
  print(c(lower_bound, upper_bound))
  filter_mrna <- (total_mrna > lower_bound) & (total_mrna < upper_bound)
  print(sum(filter_mrna))
  hist(total_mrna[filter_mrna], breaks=20)
  return(gene_counts[,filter_mrna])
}

load('../LLU-LLU/Spikein_10/merge/gene_counts_cellranger.rdata')
gene_counts <- cellFilter(gene_counts)
save(gene_counts, file='./FilteredSamples/LLU-LLU/Spikein_10/gene_counts_cellranger.rdata')

load('../LLU-LLU/HCC1395(A)/merge/gene_counts_cellranger.rdata')
gene_counts <- cellFilter(gene_counts)
save(gene_counts, file='./FilteredSamples/LLU-LLU/HCC1395(A)/gene_counts_cellranger.rdata')

load('../LLU-LLU/HCC1395BL(B)/merge/gene_counts_cellranger.rdata')
gene_counts <- cellFilter(gene_counts)
save(gene_counts, file='./FilteredSamples/LLU-LLU/HCC1395BL(B)/gene_counts_cellranger.rdata')

load('../LLU-NCI/HCC1395(A)/Seq4/gene_counts_cellranger.rdata')
gene_counts <- cellFilter(gene_counts)
save(gene_counts, file='./FilteredSamples/LLU-NCI/HCC1395(A)/gene_counts_cellranger.rdata')

load('../LLU-NCI/HCC1395BL(B)/Seq4/gene_counts_cellranger.rdata')
gene_counts <- cellFilter(gene_counts)
save(gene_counts, file='./FilteredSamples/LLU-NCI/HCC1395BL(B)/gene_counts_cellranger.rdata')

load('../LLU-NCI/Spikein_5/Seq4/gene_counts_cellranger.rdata')
gene_counts <- cellFilter(gene_counts)
save(gene_counts, file='./FilteredSamples/LLU-NCI/Spikein_5/gene_counts_cellranger.rdata')

load('../LLU-NCI/Spikein_5_Fix1/Seq4/gene_counts_cellranger.rdata')
gene_counts <- cellFilter(gene_counts)
save(gene_counts, file='./FilteredSamples/LLU-NCI/Spikein_5_Fix1/gene_counts_cellranger.rdata')


load('../NCI-NCI/HCC1395(A)/merge/gene_counts_cellranger.rdata')
gene_counts <- cellFilter(gene_counts)
save(gene_counts, file='./FilteredSamples/NCI-NCI/HCC1395(A)/gene_counts_cellranger.rdata')

load('../NCI-NCI/HCC1395BL(B)/merge/gene_counts_cellranger.rdata')
gene_counts <- cellFilter(gene_counts)
save(gene_counts, file='./FilteredSamples/NCI-NCI/HCC1395BL(B)/gene_counts_cellranger.rdata')

load('../NCI-NCI/Spikein_5/merge/gene_counts_cellranger.rdata')
gene_counts <- cellFilter(gene_counts)
save(gene_counts, file='./FilteredSamples/NCI-NCI/Spikein_5/gene_counts_cellranger.rdata')

load('../NCI-NCI/Spikein_5_Fix1/merge/gene_counts_cellranger.rdata')
gene_counts <- cellFilter(gene_counts)
save(gene_counts, file='./FilteredSamples/NCI-NCI/Spikein_5_Fix1/gene_counts_cellranger.rdata')

load('../NCI-NCI/Spikein_5_Fix2/merge/gene_counts_cellranger.rdata')
gene_counts <- cellFilter(gene_counts)
save(gene_counts, file='./FilteredSamples/NCI-NCI/Spikein_5_Fix2/gene_counts_cellranger.rdata')
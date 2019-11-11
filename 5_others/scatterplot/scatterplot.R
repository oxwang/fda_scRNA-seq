#######################################
# R script for plotting scatterplots #
#######################################

library(psych)
library("ggplot2")

## plot A samples scatter plot

A<-read.table("/.../5_others/scatterplot/Data/A_CPM_merged_all.new.csv", header=TRUE, sep="," )
colnames(A) <-c("Gene","BK_RNA-seq","10X_LLU","10X_NCI","10X_NCI_M","C1_FDA","C1_LLU","ICELL8_PE","ICELL8_SE")
png(file="scatterplot_Sample-A.png",width=12,height=12, units = "in", res=300)
pairs.panels(A[,2:9],
             method = "pearson",
             hist.col = "#00AFBB",
             density = TRUE,
             ellipses = TRUE)
dev.off()


## plot B samples scatter plot

B<-read.table("/.../5_others/scatterplot/Data/B_CPM_merged_all.new.csv", header=TRUE, sep="," )
colnames(B) <-c("gene","BK_RNA-seq","10X_LLU","10X_NCI","10X_NCI_M","C1_FDA","C1_LLU","ICELL8_PE","ICELL8_SE")
png(file="scatterplot_Sample-B.png",width=12,height=12, units = "in", res=300)
pairs.panels(B[,2:9],
             method = "pearson",
             hist.col = "#00AFBB",
             density = TRUE,
             ellipses = TRUE)
 dev.off()

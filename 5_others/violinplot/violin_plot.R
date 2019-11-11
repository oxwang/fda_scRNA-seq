#############################
# R script for violin plot  #
#############################

library(dplyr)
library("ggplot2")

## read data tables
data<-read.table("/.../5_others/violinplot/Data/proteinCodingRNA_all_A_mod.txt", header=FALSE, sep="\t", stringsAsFactors = FALSE)  %>% mutate(row = row_number())
data<-read.table("/.../5_others/violinplot/Data/antisenseRNA_all_B_mod.txt", header=FALSE, sep="\t", stringsAsFactors = FALSE)  %>% mutate(row = row_number())
data<-read.table("/.../5_others/violinplot/Data/miscRNA_all_B_mod.txt", header=FALSE, sep="\t", stringsAsFactors = FALSE)  %>% mutate(row = row_number())
data<-read.table("/.../5_others/violinplot/Data/lincRNA_all_B_mod.txt", header=FALSE, sep="\t", stringsAsFactors = FALSE)  %>% mutate(row = row_number())

## plot proteinCodingRNA violin plot
png("proteinCodingRNA_all_A_plot.png", width=8, height=4, units="in", res=300)
ggplot(data, aes(x = reorder(V3,row), y = data$V2)) + geom_violin(aes(fill =data$V3)) + geom_boxplot(width=0.05) +
  labs(x=NULL, y="log(CPM)") + guides(fill=FALSE) + scale_y_log10(breaks=c(0,0.001,0.01,0.05,0.1,0.5,1,5,10,100,500,1000,5000),
  labels=c(0,0.001,0.01,0.05,0.1,0.5,1,5,10,100,500,1000,5000)) +
  theme(panel.grid.major = element_line(color = "grey")) + theme_bw() +
  theme(text = element_text(size=10)) + theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major.x  = element_blank()) +
  labs(title="Protein Coding RNA")
dev.off()

## plot antisenseRNA violin plot
png("antisenseRNA_all_B_plot.png", width=8, height=4, units="in", res=300)
ggplot(data, aes(x = reorder(V3,row), y = data$V2)) + geom_violin(aes(fill =data$V3)) + geom_boxplot(width=0.05) +
  labs(x=NULL, y="log(CPM)") + guides(fill=FALSE) + scale_y_log10(breaks=c(0,0.001,0.01,0.05,0.1,0.5,1,5,10,100,500,1000,5000),
  labels=c(0,0.001,0.01,0.05,0.1,0.5,1,5,10,100,500,1000,5000)) +
  theme(panel.grid.major = element_line(color = "grey")) + theme_bw() +
  theme(text = element_text(size=10)) + theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major.x  = element_blank()) +
  labs(title="Antisense RNA")
dev.off()

## plot miscRNA  violin plot
png("miscRNA_all_B_plot.png", width=8, height=4, units="in", res=300)
ggplot(data, aes(x = reorder(V3,row), y = data$V2)) + geom_violin(aes(fill =data$V3)) + geom_boxplot(width=0.05) +
  labs(x=NULL, y="log(CPM)") + guides(fill=FALSE) + scale_y_log10(breaks=c(0,0.001,0.01,0.05,0.1,0.5,1,5,10,100,500,1000,5000),
labels=c(0,0.001,0.01,0.05,0.1,0.5,1,5,10,100,500,1000,5000)) +
  theme(panel.grid.major = element_line(color = "grey")) + theme_bw() +
  theme(text = element_text(size=10)) + theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major.x  = element_blank()) +
  labs(title="miscRNA")
dev.off()

## plot lincRNAs violin plot
png("lincRNA_all_B_plot.png", width=8, height=4, units="in", res=300)
ggplot(data, aes(x = reorder(V3,row), y = data$V2)) + geom_violin(aes(fill =data$V3)) + geom_boxplot(width=0.05) +
  labs(x=NULL, y="log(CPM)") + guides(fill=FALSE) + scale_y_log10(breaks=c(0,0.001,0.01,0.05,0.1,0.5,1,5,10,100,500,1000,5000),
  labels=c(0,0.001,0.01,0.05,0.1,0.5,1,5,10,100,500,1000,5000)) +
  theme(panel.grid.major = element_line(color = "grey")) + theme_bw() +
  theme(text = element_text(size=10)) + theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major.x  = element_blank()) +
  labs(title="lincRNA")
dev.off()

sessionInfo()

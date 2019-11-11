
## generate corrleation plot ##


library(corrplot)
library(Cairo)

output <- "/.../5_others/Figure_5a-f/Data/"

t_perct <- read.csv(paste(output,"tumor_perct_summary.csv",sep=""),row.names=1)
n_perct <- read.csv(paste(output,"normal_perct_summary.csv",sep=""),row.names=1)

colnames(t_perct) <- colnames(n_perct) <- c("10x_LLU","10x_NCI","10x_NCI_M","C1_FDA_HT","C1_LLU","WaferGen_PE","WaferGen_SE")


CairoJPEG(filename=paste(output,"tumor_perct_correlation_all.jpeg",sep=""),width=1000, height=1000,quality=1000)

corrplot(cor(t_perct), method="number", order="hclust",addrect=2, number.cex=3,tl.cex=3)

dev.off()


CairoJPEG(filename=paste(output,"tumor_perct_correlation_high.jpeg",sep=""),width=1000, height=1000,quality=1000)

corrplot(cor(t_perct[1:500,]), method="number", order="hclust",addrect=2, number.cex=3,tl.cex=3)

dev.off()


CairoJPEG(filename=paste(output,"tumor_perct_correlation_medium.jpeg",sep=""),width=1000, height=1000,quality=1000)

corrplot(cor(t_perct[501:1000,]), method="number", order="hclust",addrect=2, number.cex=3,tl.cex=3)

dev.off()


CairoJPEG(filename=paste(output,"tumor_perct_correlation_low.jpeg",sep=""),width=1000, height=1000,quality=1000)

corrplot(cor(t_perct[1001:1500,]), method="number", order="hclust",addrect=2, number.cex=3,tl.cex=3)

dev.off()


## normal ##

CairoJPEG(filename=paste(output,"normal_perct_correlation_all.jpeg",sep=""),width=1000, height=1000,quality=1000)

corrplot(cor(n_perct), method="number", order="hclust",addrect=2, number.cex=3,tl.cex=3)

dev.off()


CairoJPEG(filename=paste(output,"normal_perct_correlation_high.jpeg",sep=""),width=1000, height=1000,quality=1000)

corrplot(cor(n_perct[1:500,]), method="number", order="hclust",addrect=2, number.cex=3,tl.cex=3)

dev.off()


CairoJPEG(filename=paste(output,"normal_perct_correlation_medium.jpeg",sep=""),width=1000, height=1000,quality=1000)

corrplot(cor(n_perct[501:1000,]), method="number", order="hclust",addrect=2, number.cex=3,tl.cex=3)

dev.off()


CairoJPEG(filename=paste(output,"normal_perct_correlation_low.jpeg",sep=""),width=1000, height=1000,quality=1000)

corrplot(cor(n_perct[1001:1500,]), method="number", order="hclust",addrect=2, number.cex=3,tl.cex=3)

dev.off()



















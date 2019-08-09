## generate pdf figure ##

input <- "../manuscript_code/4_Batch_correction/Figure_4a-c/data//Eva_1/Multi_data/"
sample <- c("HCC1395","HCC1395BL","Mix")
hvg <- paste("HVG",c(100,500,1000,2000,4000),sep="_")
method_type <- c("CCA","ComBat","LIMMA","MNN","UNC")
method_type2 <- tolower(method_type)
col_type <- c("darkturquoise","darkorchid","hotpink","firebrick","seagreen4")

for(i in 1:length(sample))
{
	for(j in 1:length(hvg))
	{
	input2 <- paste(input,sample[i],"/",hvg[j],sep="")
	for(k in 1:length(method_type2))
	{
	da <- read.csv(paste(input2,"/",method_type2[k],"/",method_type2[k],"_umap.csv",sep=""),row.names=1)
	cell_names <- rownames(da)
	batch_cell <- c(length(grep("_1",cell_names)),length(grep("_2",cell_names)),length(grep("_3",cell_names)),length(grep("_4",cell_names)),length(grep("_5",cell_names)))
	col_type2 <- col_type[batch_cell!=0]
	batch_cell <- batch_cell[batch_cell!=0]

	pdf(paste(input2,"/",method_type[k],"_figure/",sample[i],"_",method_type2[k],"_umap.pdf",sep=""))
	par(mar=c(4.5,5,0.5,0.5))
	plot(da[,1],da[,2],col=rep(col_type2,batch_cell),cex=0.7,xlab="UMAP_1",ylab="UMAP_2",pch=16,cex.axis=1.8,cex.lab=2)
	dev.off()
	}
	
	}
}











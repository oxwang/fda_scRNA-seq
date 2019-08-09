##########################
## generate umap figure ##
##########################

batch <- c("10x_LLU","10x_NCI","C1_FDA","C1_LLU","WaferGen")
batch2 <- c("10x_LLU","10x_NCI","C1_FDA_HT","C1_LLU","WaferGen_SE")
col_type <- c("darkturquoise","darkorchid","hotpink","firebrick","seagreen4")

batch <- c("10x_Mix10_LLU","10x_Mix5_NCI","10x_Mix5_F_NCI","10x_Mix5_F2_NCI")
col_type <- c("darkturquoise","darkorchid","hotpink","firebrick")


## the code use mix samples as an example

hvg_num <- c(100,500,1000,2000,4000)

for(k in 1:length(hvg_num))
{

	input <- paste("../manuscript_code/4_Batch_correction/Figure_4a-c/bbknn/Mix/HVG_",hvg_num[k],"/",sep="")
	umap_coord <- read.csv(paste(input,"umap_coord_bbknn.csv",sep=""),row.names=1)

	batch_cell <- table(umap_coord[,3])
	batch_cell <- batch_cell[match(batch,names(batch_cell))]
	dir.create(paste(input,"bbknn_figure",sep=""),recursive=T)

	jpeg(paste(input,"bbknn_figure/Mix_bbknn.jpeg",sep=""),width=12, height=10,units="in", res=400)
	plot(umap_coord[,1],umap_coord[,2],col=rep(col_type,batch_cell),cex=0.7,xlab="UMAP_1",ylab="UMAP_2",pch=c(16,16))
	legend("bottomright",batch2,col=col_type,pch=c(16,16),bty="n")
	dev.off()

}





#########################
## get alignment score ##
#########################

calc_align_score <- function(tsne_coord,batch_cell,cell_perct)
{
	tsne_coord <- tsne_coord
	batch_cell <- batch_cell
	cell_perct <- cell_perct
	w <- batch_cell/sum(batch_cell)

	temp <- cumsum(batch_cell)
	batch_index <- cbind(c(1,temp[-length(temp)]+1),temp)

	tsne_dist <- as.matrix(dist(tsne_coord))
	diag(tsne_dist) <- 1e5
	k <- round(nrow(tsne_coord)*cell_perct)
	tsne_order <- t(apply(tsne_dist,1,order))

	k_mean <- c()
	for(i in 1:length(batch_cell))
	{
		temp <- tsne_order[batch_index[i,1]:batch_index[i,2],1:k]
		temp2 <- temp <= batch_index[i,2]
		temp3 <- temp >= batch_index[i,1]
		temp <- temp2+temp3
		k_mean <- c(k_mean,mean(apply(temp == 2,1,sum)))
	}
	align_score <- 1-abs(k_mean-k*w)/(k-k*w)
	align_score <- c(sum(w*align_score),align_score)

	return(align_score)
}






batch <- c("10x_LLU","10x_NCI","C1_FDA_HT","C1_LLU","WaferGen_SE")
batch2 <- c("10x_LLU","10x_NCI","C1_FDA","C1_LLU","WaferGen")

batch <- c("10x_Mix10_LLU","10x_Mix5_NCI","10x_Mix5_F_NCI","10x_Mix5_F2_NCI")
batch2 <- c("10x_Mix10_LLU","10x_Mix5_NCI","10x_Mix5_F_NCI","10x_Mix5_F2_NCI")

hvg_num <- c(100,500,1000,2000,4000)

align_score_bbknn_all <- c()

for(k in 1:length(hvg_num))
{

output <- paste("../manuscript_code/4_Batch_correction/Figure_4a-c/bbknn/Mix/HVG_",hvg_num[k],"/",sep="")
umap_coord_bbknn <- read.csv(paste(output,"umap_coord_bbknn.csv",sep=""),row.names=1)

batch_cell_bbknn <- table(umap_coord_bbknn[,3])
batch_cell_bbknn <- as.numeric(batch_cell_bbknn[match(batch2,names(batch_cell_bbknn))])

perct <- 0.01
align_score_bbknn_all <- rbind(align_score_bbknn_all,calc_align_score(umap_coord_bbknn[,-3],batch_cell_bbknn,perct))

}


align_score_summary <- round(align_score_bbknn_all,3)
colnames(align_score_summary) <- c("All",batch)
rownames(align_score_summary) <- paste("bbknn","HVG",hvg_num,sep="_")
write.csv(align_score_summary,file="../manuscript_code/4_Batch_correction/Figure_4a-c/bbknn/Mix/align_score_summary.csv",quote=F,row.names=T)










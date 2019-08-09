
## save count data for scanorama and bbknn (data are stored in .../Figure 4a-c/scanorama)

## HCC1395 ##

output <- "/media/genomicslab/Data/1_Projects/FDA_QC/Results/batch/Eva_1/Multi_data/scanorama/HCC1395/"
dir.create(output,recursive=T)

load("/media/genomicslab/Data/1_Projects/FDA_QC/Results/10x/merge/gene_counts/HCC1395/gene_counts_cellranger_filter.rdata")
gene_counts <- gene_counts_filter
gene_counts <- cbind.data.frame(rownames(gene_counts),as.matrix(gene_counts))
colnames(gene_counts) <- c("gene",colnames(gene_counts)[-1])
write.table(gene_counts,file=paste(output,"10x_LLU.txt",sep=""),sep="\t",row.names=F,quote=F)

load("/media/genomicslab/Data/1_Projects/FDA_QC/Results/10x/Third_run/NCI/gene_counts/HCC1395/gene_counts_cellranger_filter.rdata")
gene_counts <- gene_counts_filter
gene_counts <- cbind.data.frame(rownames(gene_counts),as.matrix(gene_counts))
colnames(gene_counts) <- c("gene",colnames(gene_counts)[-1])
write.table(gene_counts,file=paste(output,"10x_NCI.txt",sep=""),sep="\t",row.names=F,quote=F)

load("/media/genomicslab/Data/1_Projects/FDA_QC/Results/C1_ceber/gene_counts/HCC1395/gene_counts_featureCounts.rdata")
gene_counts <- cbind.data.frame(rownames(gene_counts),as.matrix(gene_counts))
colnames(gene_counts) <- c("gene",colnames(gene_counts)[-1])
write.table(gene_counts,file=paste(output,"C1_FDA.txt",sep=""),sep="\t",row.names=F,quote=F)

load("/media/genomicslab/Data/1_Projects/FDA_QC/Results/C1/gene_counts/HCC1395/gene_counts_featureCounts.rdata")
gene_counts <- cbind.data.frame(rownames(gene_counts),as.matrix(gene_counts))
colnames(gene_counts) <- c("gene",colnames(gene_counts)[-1])
write.table(gene_counts,file=paste(output,"C1_LLU.txt",sep=""),sep="\t",row.names=F,quote=F)

load("/media/genomicslab/Data/1_Projects/FDA_QC/Results/WaferGen/Second_run/gene_counts/HCC1395/gene_counts_featureCounts.rdata")
gene_counts <- cbind.data.frame(rownames(gene_counts),as.matrix(gene_counts))
colnames(gene_counts) <- c("gene",colnames(gene_counts)[-1])
write.table(gene_counts,file=paste(output,"WaferGen.txt",sep=""),sep="\t",row.names=F,quote=F)



## HCC1395BL ##

output <- "/media/genomicslab/Data/1_Projects/FDA_QC/Results/batch/Eva_1/Multi_data/scanorama/HCC1395BL/"
dir.create(output,recursive=T)

load("/media/genomicslab/Data/1_Projects/FDA_QC/Results/10x/merge/gene_counts/HCC1395BL/gene_counts_cellranger_filter.rdata")
gene_counts <- gene_counts_filter
gene_counts <- cbind.data.frame(rownames(gene_counts),as.matrix(gene_counts))
colnames(gene_counts) <- c("gene",colnames(gene_counts)[-1])
write.table(gene_counts,file=paste(output,"10x_LLU.txt",sep=""),sep="\t",row.names=F,quote=F)

load("/media/genomicslab/Data/1_Projects/FDA_QC/Results/10x/Third_run/NCI/gene_counts/HCC1395BL/gene_counts_cellranger_filter.rdata")
gene_counts <- gene_counts_filter
gene_counts <- cbind.data.frame(rownames(gene_counts),as.matrix(gene_counts))
colnames(gene_counts) <- c("gene",colnames(gene_counts)[-1])
write.table(gene_counts,file=paste(output,"10x_NCI.txt",sep=""),sep="\t",row.names=F,quote=F)

load("/media/genomicslab/Data/1_Projects/FDA_QC/Results/C1_ceber/gene_counts/HCC1395BL/gene_counts_featureCounts.rdata")
gene_counts <- cbind.data.frame(rownames(gene_counts),as.matrix(gene_counts))
colnames(gene_counts) <- c("gene",colnames(gene_counts)[-1])
write.table(gene_counts,file=paste(output,"C1_FDA.txt",sep=""),sep="\t",row.names=F,quote=F)

load("/media/genomicslab/Data/1_Projects/FDA_QC/Results/C1/gene_counts/HCC1395BL/gene_counts_featureCounts.rdata")
gene_counts <- cbind.data.frame(rownames(gene_counts),as.matrix(gene_counts))
colnames(gene_counts) <- c("gene",colnames(gene_counts)[-1])
write.table(gene_counts,file=paste(output,"C1_LLU.txt",sep=""),sep="\t",row.names=F,quote=F)

load("/media/genomicslab/Data/1_Projects/FDA_QC/Results/WaferGen/Second_run/gene_counts/HCC1395BL/gene_counts_featureCounts.rdata")
gene_counts <- cbind.data.frame(rownames(gene_counts),as.matrix(gene_counts))
colnames(gene_counts) <- c("gene",colnames(gene_counts)[-1])
write.table(gene_counts,file=paste(output,"WaferGen.txt",sep=""),sep="\t",row.names=F,quote=F)



## Mix ##

output <- "/media/genomicslab/Data/1_Projects/FDA_QC/Results/batch/Eva_1/Multi_data/scanorama/Mix/"
dir.create(output,recursive=T)

load("/media/genomicslab/Data/1_Projects/FDA_QC/Results/10x/merge/gene_counts/Mix10/gene_counts_cellranger_filter.rdata")
gene_counts <- gene_counts_filter
gene_counts <- cbind.data.frame(rownames(gene_counts),as.matrix(gene_counts))
colnames(gene_counts) <- c("gene",colnames(gene_counts)[-1])
write.table(gene_counts,file=paste(output,"10x_Mix10_LLU.txt",sep=""),sep="\t",row.names=F,quote=F)

load("/media/genomicslab/Data/1_Projects/FDA_QC/Results/10x/NCI/merge/gene_counts/Mix5/gene_counts_cellranger_filter.rdata")
gene_counts <- gene_counts_filter
gene_counts <- cbind.data.frame(rownames(gene_counts),as.matrix(gene_counts))
colnames(gene_counts) <- c("gene",colnames(gene_counts)[-1])
write.table(gene_counts,file=paste(output,"10x_Mix5_NCI.txt",sep=""),sep="\t",row.names=F,quote=F)

load("/media/genomicslab/Data/1_Projects/FDA_QC/Results/10x/NCI/merge/gene_counts/Mix5_F/gene_counts_cellranger_filter.rdata")
gene_counts <- gene_counts_filter
gene_counts <- cbind.data.frame(rownames(gene_counts),as.matrix(gene_counts))
colnames(gene_counts) <- c("gene",colnames(gene_counts)[-1])
write.table(gene_counts,file=paste(output,"10x_Mix5_F_NCI.txt",sep=""),sep="\t",row.names=F,quote=F)

load("/media/genomicslab/Data/1_Projects/FDA_QC/Results/10x/NCI/merge/gene_counts/Mix5_F_2/gene_counts_cellranger_filter.rdata")
gene_counts <- gene_counts_filter
gene_counts <- cbind.data.frame(rownames(gene_counts),as.matrix(gene_counts))
colnames(gene_counts) <- c("gene",colnames(gene_counts)[-1])
write.table(gene_counts,file=paste(output,"10x_Mix5_F2_NCI.txt",sep=""),sep="\t",row.names=F,quote=F)





## scanorama process

## HCC1395 ##

python_path=/home/genomics/bioinformatics/app/miniconda/miniconda3/bin/python
cd /home/genomics/bioinformatics/app/scanorama
$python_path bin/process.py conf/hcc1395.txt
$python_path bin/hcc1395.py


## HCC1395BL ##

python_path=/home/genomics/bioinformatics/app/miniconda/miniconda3/bin/python
cd /home/genomics/bioinformatics/app/scanorama
$python_path bin/process.py conf/hcc1395bl.txt
$python_path bin/hcc1395bl.py


## Mix ##

python_path=/home/genomics/bioinformatics/app/miniconda/miniconda3/bin/python
cd /home/genomics/bioinformatics/app/scanorama
$python_path bin/process.py conf/mix.txt
$python_path bin/mix.py







## generate tsne by R

library(Rtsne)

batch <- c("10x_LLU","10x_NCI","C1_FDA_HT","C1_LLU","WaferGen_SE")
col_type <- c("darkturquoise","darkorchid","hotpink","firebrick","seagreen4")

batch <- c("10x_Mix10_LLU","10x_Mix5_NCI","10x_Mix5_F_NCI","10x_Mix5_F2_NCI")
col_type <- c("darkturquoise","darkorchid","hotpink","firebrick")


hvg_num <- c(100,500,1000,2000,4000)

for(k in 1:length(hvg_num))
{

	input <- paste("../manuscript_code/4_Batch_correction/Figure_4a-c/scanorama/data/HCC1395/HVG_",hvg_num[k],"/",sep="")
	file <- list.files(input,pattern = ".txt")
	combined <- c()
	batch_cell <- c()

	for(i in 1:length(file))
	{
	temp <- read.table(paste(input,file[i],sep=""),sep="\t",header=T,row.names=1)
	colnames(temp) <- paste(colnames(temp),i,sep="_")
	combined <- cbind(combined,as.matrix(temp))
	batch_cell <- c(batch_cell,ncol(temp))
	}

	names(batch_cell) <- batch
	write.csv(batch_cell,file=paste(input,"batch_cell.csv",sep=""))

	t_scanorama <- as.matrix(t(combined))
	set.seed(0)
	tsne_scanorama <- Rtsne(t_scanorama, perplexity = 30)
	tsne_scanorama_coord_all <- tsne_scanorama$Y

	dir.create(paste(input,"Scanorama_figure",sep=""),recursive=T)
	jpeg(paste(input,"Scanorama_figure/scanorama.jpeg",sep=""),width=12, height=10,units="in", res=400)
	plot(tsne_scanorama$Y[,1],tsne_scanorama$Y[,2],col=rep(col_type,batch_cell),cex=0.7,xlab="tSNE_1",ylab="tSNE_2",pch=c(16,16))
	legend("bottomright",batch,col=col_type,pch=c(16,16),bty="n")
	dev.off()

	pdf(paste(input,"Scanorama_figure/scanorama.pdf",sep=""))
	par(mar=c(4.5,5,0.5,0.5))
	plot(tsne_scanorama$Y[,1],tsne_scanorama$Y[,2],col=rep(col_type,batch_cell),cex=0.7,xlab="tSNE_1",ylab="tSNE_2",pch=16,cex.axis=2,cex.lab=2)
	dev.off()

	save(tsne_scanorama_coord_all,file=paste(input,"tsne_coord_scanorama.rdata",sep=""))
}



## generate umap coordinates by python ##

samples=['HCC1395','HCC1395BL','Mix']
datasets=['10x_LLU','10x_NCI','C1_FDA','C1_LLU','WaferGen']
datasets2=['10x_Mix10_LLU','10x_Mix5_NCI','10x_Mix5_F_NCI','10x_Mix5_F2_NCI']
hvgs=[100,500,1000,2000,4000]


for sample in samples:
	for hvg in hvgs:
		dataset_chosen=datasets
		holder=[]
		if sample=='Mix':
		 dataset_chosen=datasets2
		for dataset in dataset_chosen:
		 input='../manuscript_code/4_Batch_correction/Figure_4a-c/scanorama/data/'+sample+'/HVG_'+str(hvg)+'/'+dataset+'.scanorama_corrected.txt'
		 counts=pd.read_table(input,index_col=0)
		 holder.append(anndata.AnnData(X=counts.values).T)
		 holder[-1].var_names=counts.index
		 holder[-1].obs_names=counts.columns
		 holder[-1].obs['sample']=dataset
		adata = holder[0].concatenate(holder[1:], join='outer')
		sc.tl.pca(adata)
		adata.obsm['X_pca']*=-1
		sc.pp.neighbors(adata,n_pcs=20,n_neighbors=20)
		sc.tl.umap(adata)
		umap_coord=pd.DataFrame(adata.obsm['X_umap'])
		umap_coord.index=adata.obs['sample']
		output='../manuscript_code/4_Batch_correction/Figure_4a-c/scanorama/data/'+sample+'/HVG_'+str(hvg)+'/scanorama_umap.csv'
		umap_coord.to_csv(output)



## generate umap pdf figure ##

input <- "../manuscript_code/4_Batch_correction/Figure_4a-c/scanorama/data/"
sample <- c("HCC1395","HCC1395BL","Mix")
hvg <- paste("HVG",c(100,500,1000,2000,4000),sep="_")
col_type <- c("darkturquoise","darkorchid","hotpink","firebrick","seagreen4")
batch <- c("10x_LLU","10x_NCI","C1_FDA","C1_LLU","WaferGen")
batch2 <- c("10x_Mix10_LLU","10x_Mix5_NCI","10x_Mix5_F_NCI","10x_Mix5_F2_NCI")


for(i in 1:length(sample))
{
	for(j in 1:length(hvg))
	{
	input2 <- paste(input,sample[i],"/",hvg[j],sep="")
	da <- read.csv(paste(input2,"/scanorama_umap.csv",sep=""))
	if(i < 3) {temp_batch <- batch}
	if(i == 3) {temp_batch <- batch2}
	batch_cell <- table(da[,1])
	batch_cell <- batch_cell[match(names(batch_cell),temp_batch)]
	col_type2 <- col_type[1:length(batch_cell)]

	pdf(paste(input2,"/Scanorama_figure/",sample[i],"_scanorama_umap.pdf",sep=""))
	par(mar=c(4.5,5,1,0.5))
	plot(da[,2],da[,3],col=rep(col_type2,batch_cell),cex=0.7,xlab="UMAP_1",ylab="UMAP_2",pch=16,cex.axis=1.8,cex.lab=2)
	dev.off()	
	}
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
batch <- c("10x_Mix10_LLU","10x_Mix5_NCI","10x_Mix5_F_NCI","10x_Mix5_F2_NCI")
hvg_num <- c(100,500,1000,2000,4000)

align_score_scanorama_all <- c()

for(k in 1:length(hvg_num))
{

output <- paste("../manuscript_code/4_Batch_correction/Figure_4a-c/scanorama/data/HCC1395/HVG_",hvg_num[k],"/",sep="")
load(paste(output,"tsne_coord_scanorama.rdata",sep=""))

batch_cell_scanorama <- as.numeric(t(read.csv(paste(output,"batch_cell.csv",sep=""),row.names=1)))

perct <- 0.01
align_score_scanorama_all <- rbind(align_score_scanorama_all,calc_align_score(tsne_scanorama_coord_all,batch_cell_scanorama,perct))

}


align_score_summary <- round(align_score_scanorama_all,3)
colnames(align_score_summary) <- c("All",batch)
rownames(align_score_summary) <- paste("Scanorama","HVG",hvg_num,sep="_")
write.csv(align_score_summary,file="../manuscript_code/4_Batch_correction/Figure_4a-c/scanorama/data/HCC1395/align_score_summary.csv",quote=F,row.names=T)










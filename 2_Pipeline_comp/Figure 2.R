
## barplot of total number of expressed genes,
## boxplot of number of expressed genes per cell, and
## violin plot of gene expression corrleation and IoU

input <- '/.../2_Pipeline_comp/Data/'
output <- '/.../2_Pipeline_comp/Results/'

## 10x, figure 2a-c ##

filepath <- vector("list",3)

filepath[[1]] <- paste0(input,'/10x/LLU/')
filepath[[2]] <- paste0(input,'/10x/NCI/')
filepath[[3]] <- paste0(input,'/10x/NCI_M/')

seq_type <- c("10x_LLU","10x_NCI","10x_NCI_M")
samples <- c("HCC1395","HCC1395BL")
samples2 <- c("A","B")

da <- c()
gene_exprs_num <- c()

for(i in 1:length(filepath))
{
	temp2 <- c()
	temp4 <- c()
	for(j in 1:length(samples))
	{
	temp <- read.csv(paste(filepath[[i]],samples[j],"/ranger_zumi_umitools_gene.csv",sep=""))
	temp2 <- rbind(temp2,temp)

	temp3 <- read.csv(paste(filepath[[i]],samples[j],"/gene_exprs_num_per_cell.csv",sep=""))
	temp3 <- cbind(temp3,samples2[j])
	temp4 <- rbind(temp4,temp3)
	}
	temp2 <- cbind(temp2,seq_type[i])
	temp4 <- cbind(temp4,seq_type[i])
	da <- rbind(da,temp2)
	gene_exprs_num <- rbind(gene_exprs_num,temp4)
}

levels(da[,9])[levels(da[,9])=="HCC1395"] <- "A"
levels(da[,9])[levels(da[,9])=="HCC1395BL"] <- "B"

da <- cbind(da,paste(da[,10],da[,9],sep="_"))
temp <- colnames(da)
temp[c(10,11)] <- c("Seq_site","Cell_Seq")
colnames(da) <- temp

gene_exprs_num <- cbind(gene_exprs_num,paste(gene_exprs_num[,5],gene_exprs_num[,4],sep="_"))
colnames(gene_exprs_num) <- c("Cell_id","Expressed_gene","Pipeline","Cell_type","Seq_site","Cell_seq")

Cell_seq <- names(table(gene_exprs_num$Cell_seq))
cell_num <- c()
gene_exprs_num2 <- c()

for(i in 1:length(Cell_seq))
{
	temp <- gene_exprs_num[gene_exprs_num$Cell_seq==Cell_seq[i],]
	cell_num <- rbind(cell_num,table(temp[,3]))
	temp_1 <- temp[temp[,3]=="cellrange",]
	temp_2 <- temp[temp[,3]=="umi_tools",]
	temp_3 <- temp[temp[,3]=="zUMIs",]

	temp_2[,1] <- substr(temp_2[,1],1,16)
	temp_3[,1] <- substr(temp_3[,1],1,16)

	overlap_cell <- intersect(temp_1[,1],temp_2[,1])
	overlap_cell <- intersect(overlap_cell,temp_3[,1])

	temp_1 <- temp_1[match(overlap_cell,temp_1[,1]),]
	temp_2 <- temp_2[match(overlap_cell,temp_2[,1]),]
	temp_3 <- temp_3[match(overlap_cell,temp_3[,1]),]

	gene_exprs_num2 <- rbind(gene_exprs_num2,temp_1,temp_2,temp_3)
}

rownames(cell_num) <- Cell_seq





library(ggplot2)
library(reshape2)


## plot cor and IoU together


da_2 <- melt(da[,6:8], id="Label")
da_2 <- cbind(da_2,da[,11])
colnames(da_2) <- c(colnames(da_2)[-4],"Cell_Seq")

da_2[,4] <- factor(da_2[,4],levels=paste(rep(seq_type,2),rep(samples2,rep(3,2)),sep="_"))

jpeg(paste0(output,"/10x/summary/IoU_cor_10x.jpeg"),width=12,height=6,units="in",res=600)

p <- ggplot(da_2,aes(x=Cell_Seq,y=value,fill=Label))+geom_violin(position=position_dodge())+facet_grid(.~variable)+ 
labs(title="Gene expression (counts) correlation and fraction of consensus genes per cell",x="", y="") + 
theme(plot.title = element_text(hjust = 0.5))

dev.off()




da_3 <- gene_exprs_num2
da_3[,6] <- factor(da_3[,6],levels=paste(rep(seq_type,2),rep(samples2,rep(3,2)),sep="_"))

jpeg(paste0(output,"/10x/summary/gene_exprs_num_per_cell_10x.jpeg"),width=8,height=6,units="in",res=600)

p <- ggplot(da_3, aes(x=Cell_seq,y=Expressed_gene,fill=Pipeline)) + geom_boxplot(position=position_dodge()) + 
labs(title="Number of expressed genes per cell",x="", y="") + theme(plot.title = element_text(hjust = 0.5))

dev.off()


da_4 <- cbind.data.frame(matrix(t(cell_num),ncol=1),rep(colnames(cell_num),nrow(cell_num)),rep(rownames(cell_num),rep(ncol(cell_num),nrow(cell_num))))
colnames(da_4) <- c("Cell_num","Pipeline","Cell_seq")
da_4[,3] <- factor(da_4[,3],levels=paste(rep(seq_type,2),rep(samples2,rep(3,2)),sep="_"))


jpeg(paste0(output,"/10x/summary/num_cell_10x.jpeg"),width=8,height=6,units="in",res=600)

p <- ggplot(da_4, aes(x=Cell_seq,y=Cell_num,fill=Pipeline)) + geom_bar(stat="identity",position=position_dodge()) + 
labs(title="Number of selected cells",x="", y="") + theme(plot.title = element_text(hjust = 0.5))

dev.off()








## other data, figure 2d-e  ##


## barplot of total number of expressed genes for 4 cell types and 
## boxplot of number of expressed genes per cell for 4 cell types

filepath <- vector("list",4)

filepath[[1]] <- paste0(input,"/C1/")
filepath[[2]] <- paste0(input,"/C1_HT/")
filepath[[3]] <- paste0(input,"/iCELL8/PE/")
filepath[[4]] <- paste0(input,"/iCELL8/SE/")

seq_type <- c("C1_LLU","C1_FDA_HT","iCELL8_PE","iCELL8_SE")
samples <- c("HCC1395","HCC1395BL")
samples2 <- c("A","B")

da <- c()
gene_exprs_num <- c()

for(i in 1:length(filepath))
{
	temp2 <- c()
	temp4 <- c()
	for(j in 1:length(samples))
	{
	temp <- read.csv(paste(filepath[[i]],samples[j],"/rsem_featureCounts_kallisto_gene.csv",sep=""))
	temp2 <- rbind(temp2,temp)

	temp3 <- read.csv(paste(filepath[[i]],samples[j],"/gene_exprs_num_per_cell.csv",sep=""))
	temp3 <- cbind(temp3,samples2[j])
	temp4 <- rbind(temp4,temp3)
	}
	temp2 <- cbind(temp2,seq_type[i])
	temp4 <- cbind(temp4,seq_type[i])
	da <- rbind(da,temp2)
	gene_exprs_num <- rbind(gene_exprs_num,temp4)
}

levels(da[,9])[levels(da[,9])=="HCC1395"] <- "A"
levels(da[,9])[levels(da[,9])=="HCC1395BL"] <- "B"

da <- cbind(da,paste(da[,10],da[,9],sep="_"))
temp <- colnames(da)
temp[c(10,11)] <- c("Seq_site","Cell_Seq")
colnames(da) <- temp

gene_exprs_num <- cbind(gene_exprs_num,paste(gene_exprs_num[,5],gene_exprs_num[,4],sep="_"))
colnames(gene_exprs_num) <- c("Cell_id","Expressed_gene","Pipeline","Cell_type","Seq_site","Cell_seq")







library(ggplot2)
library(reshape2)


## plot cor and IoU together


da_2 <- melt(da[,6:8], id="Label")
da_2 <- cbind(da_2,da[,11])
colnames(da_2) <- c(colnames(da_2)[-4],"Cell_Seq")

da_2[,4] <- factor(da_2[,4],levels=paste(rep(seq_type,2),rep(samples2,rep(4,2)),sep="_"))

jpeg(paste0(output,"/other/summary/IoU_cor_others.jpeg",width=12,height=6,units="in",res=600)

p <- ggplot(da_2,aes(x=Cell_Seq,y=value,fill=Label))+geom_violin(position=position_dodge())+facet_grid(.~variable)+ 
labs(title="Gene expression (counts) correlation and fraction of consensus genes per cell",x="", y="") + 
theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=5))


dev.off()




da_3 <- gene_exprs_num
da_3[,6] <- factor(da_3[,6],levels=paste(rep(seq_type,2),rep(samples2,rep(4,2)),sep="_"))

jpeg(paste0(output,"other/summary/gene_exprs_num_per_cell_others.jpeg",sep=""),width=10,height=6,units="in",res=600)

p <- ggplot(da_3, aes(x=Cell_seq,y=Expressed_gene,fill=Pipeline)) + geom_boxplot(position=position_dodge()) + 
labs(title="Number of expressed genes per cell",x="", y="") + theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=8))

dev.off()





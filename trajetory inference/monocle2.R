if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("monocle")
library("monocle")
library(monocle)
YDL<- readRDS("YDL_NONPRO.rds")
#准备monocle分析需要的文件
monocle.matrix=as.matrix(YDL@assays$RNA@data)
monocle.matrix=cbind(id=row.names(monocle.matrix),monocle.matrix)
write.table(monocle.matrix,file="monocleMatrix.txt",quote=F,sep="\t",row.names=F)
monocle.sample=as.matrix(YDL@meta.data)
monocle.sample=cbind(id=row.names(monocle.sample),monocle.sample)
write.table(monocle.sample,file="monocleSample.txt",quote=F,sep="\t",row.names=F)
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.geneAnn=cbind(id=row.names(monocle.geneAnn),monocle.geneAnn)
write.table(monocle.geneAnn,file="monocleGene.txt",quote=F,sep="\t",row.names=F)
write.table(singler$other,file="monocleClusterAnn.txt",quote=F,sep="\t",col.names=F)
write.table(YDL.markers,file="monocleMarkers.txt",sep="\t",row.names=F,quote=F)


#设置工作目录
monocle.matrix=read.table("monocleMatrix.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.sample=read.table("monocleSample.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.geneAnn=read.table("monocleGene.txt",sep="\t",header=T,row.names=1,check.names=F)
marker=read.table("monocleMarkers.txt",sep="\t",header=T,check.names=F)

#将Seurat结果转换为monocle需要的细胞矩阵，细胞注释表和基因注释表表
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)

#给其中一列数据重命名
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])
#伪时间分析流程
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, marker$gene)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,reduction_method = 'DDRTree')
cds <- orderCells(cds)
pdf(file="cluster.trajectory.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, color_by = "Cluster")
plot_cell_trajectory(cds,color_by="Cluster")+facet_wrap(~Cluster,nrow=3,ncol = 3)
plot_cell_trajectory(cds, color_by = "YDL@active.ident")
plot_cell_trajectory(cds, color_by="Pseudotime", show_backbone=FALSE)
# 可以很明显看到细胞的发育轨迹 
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "State") +facet_wrap(~State, nrow = 1)
dev.off()
pdf(file="cellType.trajectory.pdf",width=6.5,height=6)]
plot_cell_trajectory(cds,color_by = "cell_type2")
dev.off()
saveRDS(cds,"cds.rds")

library(monocle)
rm(list=ls())
cds<-readRDS("cds.rds")
hm<-plot_pseudotime_heatmap(cds,num_clusters=6,show_rownames=F,return_heatmap = T)
#hm<-plot_pseudotime_heatmap(cds[ordering_genes,],num_clusters=4,show_rownames=F,return_heatmap = T)
save(hm,file="heatmap.Rdata")
c1 <- as.data.frame(cutree(hm$tree_row, k=4)) 
colnames(c1) <- "Cluster"
c1$Gene <- rownames(c1)
write.table(c1,"genes_in_heatmap_clusters.txt",sep="\t",row.names = F)

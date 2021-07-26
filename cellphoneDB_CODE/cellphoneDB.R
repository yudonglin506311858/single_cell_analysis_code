Sys.setenv(LANGUAGE = "en") #显示英文报错信息
gc()
memory.limit(9999999999)
set.seed(123)
rm(list = ls())  
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)
setwd("D:/analysis/bone_marrow_niche")

load("NicheData10x.rda")
YDL<-NicheData10x

library(Seurat)
library(tidyverse)

YDL_counts <- as.matrix(YDL@assays$RNA@data)
YDL_counts <- data.frame(gene=rownames(YDL_counts), YDL_counts)
YDL_meta <- data.frame(Cell=rownames(YDL@meta.data), cell_type=YDL@meta.data$celltype)
#YDL_meta <- data.frame(Cell=rownames(YDL@meta.data), cell_type=YDL@meta.data$seurat_clusters)
#write.table(YDL_counts, "test_counts.txt", row.names=F, sep='\t')#如果是人类数据直接读出文件
write.table(YDL_meta, "test_meta.txt", row.names=F, sep='\t', quote=F)
head(YDL_meta)
head(YDL_counts[1:4,1:4])
#可能需要将test_meta.txt的.改为-


library(biomaRt)
# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"),
                   filters = "mgi_symbol",
                   values = x , mart = mouse,
                   attributesL = c("hgnc_symbol"),
                   martL = human, uniqueRows=T)
  #humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  #print(head(genes))
  return(genesV2)
}
a<-convertMouseGeneList(rownames(YDL_counts))
head(a)
colnames(a)<-c("gene","Gene")
a<-a[!duplicated(a$Gene),]#去除人类基因重复值
test_counts<-merge(a,YDL_counts,by=c("gene"))
test_counts[1:4,1:4]
test_count<-test_counts[,-1]#删除第一列
test_count[1:4,1:4]
dim(test_count)
dim(YDL_meta)
NAME<-c("Gene",meta_data$cell_type)
head(NAME)
NAME<-c(NAME)
head(meta_data$Cell)
head(colnames(test_count))
colnames(test_count)<-NAME
#或者更加简单
#colnames(test_count)<-c("Gene",meta_data$Cell)
write.table(test_count, "test_counts.txt", row.names=F, sep='\t', quote=F)

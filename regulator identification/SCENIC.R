if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::version()
# If your bioconductor version is previous to 3.9, see the section bellow

## Required
BiocManager::install(c("AUCell", "RcisTarget"))
BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost

## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "pheatmap", "R2HTML", "Rtsne"))
# To support paralell execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"))
# To export/visualize in http://scope.aertslab.org
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)

library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(foreach)
library(Seurat)
library(GENIE3)
library(AUCell)
library(RcisTarget)
library(loomR) 
##==分析准备==##
set.seed(123)
dir.create("SCENIC")
dir.create("SCENIC/int")
YDL<-readRDS("ydl.rds")
scRNA <- YDL
setwd("./SCENIC") 
##准备细胞meta信息
##准备细胞meta信息
cellInfo <- data.frame(scRNA@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
#colnames(cellInfo)[which(colnames(cellInfo)=="celltype_Monaco")] <- "celltype"
cellInfo <- cellInfo[,c("sample","cluster")]
saveRDS(cellInfo, file="int/cellInfo.Rds")
##准备表达矩阵
#为了节省计算资源，随机抽取1000个细胞的数据子集
subcell <- sample(colnames(scRNA),1000)
scRNAsub <- scRNA[,subcell]
saveRDS(scRNAsub, "scRNAsub.rds")
exprMat <- as.matrix(scRNAsub@assays$RNA@counts)
dim(exprMat)
# [1] 32285  1000
head(cellInfo)
mydbDIR <- "/data/yudonglin/scenic/cisTarget_databases/cisTarget_databases/"
mydbs <- c("mm9-500bp-upstream-7species.mc9nr.feather",
           "mm9-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")
scenicOptions <- initializeScenic(org="mgi", 
                                  nCores=40,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "SCENIC")
saveRDS(scenicOptions, "int/scenicOptions.rds")
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
##计算相关性矩阵
runCorrelation(exprMat_filtered, scenicOptions)
##TF-Targets相关性回归分析
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 20)
##推断共表达模块
runSCENIC_1_coexNetwork2modules(scenicOptions)
##推断转录调控网络（regulon）
runSCENIC_2_createRegulons(scenicOptions)
#以上代码可增加参数coexMethod=c("w001", "w005", "top50", "top5perTarget", "top10perTarget", "top50perTarget"))
#默认6种方法的共表达网络都计算，可以少选几种方法以减少计算量
##==regulon活性评分与可视化==##
##regulons计算AUC值并进行下游分析
exprMat_all <- as.matrix(scRNA@assays$RNA@counts)
exprMat_all <- log2(exprMat_all+1)
saveRDS(exprMat_all,"exprMat_all.rds")
# rm(list=ls())
# scenicOptions<-readRDS("scenicOptions.rds")
# exprMat_all<-readRDS("exprMat_all.rds")
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)
# #使用shiny互动调整阈值
# aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_all)
# savedSelections <- shiny::runApp(aucellApp)
# #保存调整后的阈值
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
# saveRDS(scenicOptions, file="int/scenicOptions.Rds")
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all)

#REFERENCE:https://mp.weixin.qq.com/s/QBBQ2TXQzNrjPILNi4EHdA

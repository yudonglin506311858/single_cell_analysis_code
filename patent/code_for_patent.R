
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

setwd("D:/analysis/forpublication/summary/nonpro1_own_nonpro2")

YDL<-readRDS("D:/analysis/forpublication/summary/nonpro1_own_nonpro2/nonpro1_nonpro2_SINGLET.RDS")

setwd("D:/课题总结/CD44-专利")


current.cluster.ids <- c(0, 1, 2, 3, 4, 5,6)
new.cluster.ids <- c(
  "BasoE2",
  "Ery/Mono",
  "BasoE1",
  "Ery/B",
  "CD44hiOrthoE1",
  "CD44hiOrthoE2","ProE")
names(new.cluster.ids) <- levels(YDL)
YDL <- RenameIdents(YDL, new.cluster.ids)
DimPlot(YDL, reduction = "tsne", label = TRUE, pt.size = 1.5)
DimPlot(YDL, reduction = "umap", label = TRUE, pt.size = 1.5)

YDL$celltype<-Idents(YDL)
YDL$celltype<-factor(YDL$celltype,levels = c("ProE","BasoE1","BasoE2","CD44hiOrthoE1","CD44hiOrthoE2","Ery/B","Ery/Mono"))
Idents(YDL)<-YDL$celltype

cols = c("#FB61D7","#53B400","#F8766D","#00B6EB","#A58AFF","#C49A00","#00C094")
pdf("专利_figure2B.pdf")
DimPlot(YDL, reduction = "umap", label = TRUE, pt.size = 1.5,cols = c("#FB61D7","#53B400","#F8766D","#00B6EB","#A58AFF","#00C094","#C49A00"))
DimPlot(YDL, reduction = "umap", label = F, pt.size = 1.5,cols = c("#FB61D7","#53B400","#F8766D","#00B6EB","#A58AFF","#00C094","#C49A00"))
DimPlot(YDL, reduction = "umap", label = TRUE, pt.size = 1.5,cols = c("#FB61D7","#53B400","#F8766D","#00B6EB","#A58AFF","#00C094","#C49A00"),split.by = "orig.ident")
DimPlot(YDL, reduction = "umap", label = F, pt.size = 1.5,cols = c("#FB61D7","#53B400","#F8766D","#00B6EB","#A58AFF","#00C094","#C49A00"),split.by = "orig.ident")

# Color palette
#colors <- c("#F8766D","#E68613","#CD9600","#ABA300","#7CAE00","#0CB702","#00BE67","#00C19A","#00BFC4","#00B8E7","#00A9FF","#8494FF","#C77CFF","#ED68ED","#FF61CC","#FF68A1")
colors <- c("#FB61D7","#53B400","#F8766D","#00B6EB","#A58AFF","#C49A00","#00C094")

DimPlot(YDL,cols =colors,reduction="umap")
a<-as.data.frame(table(YDL@meta.data$celltype))
# 数据准备
info = a$Freq
# 命名
names = as.character(a$Var1)
#names = c("0","1","2","3","4","5","6","7")
# 涂色（可选）
#cols = c("#ED1C24","#22B14C","#FFC90E","#3f48CC","#3f90CC","#22B17C","#FFC93E")
# 计算百分比
piepercent = paste(round(100*info/sum(info),2), "%")
# 绘图
pie(info, labels=piepercent, main = "total cluster ratio", col=colors, family='GB1')
# 添加颜色样本标注
legend("topright", names, cex=0.7, fill=colors)


library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","nCount_RNA"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = log(nCount_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","nFeature_RNA"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = log(nFeature_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

dev.off()





pdf("专利_figure2C.pdf")
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Hbb-bs","Gypa","Gata1","Nfe2"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Tfrc","Cd44","Xpo7","Malat1"),cols = c("gray", "red"))#actin
dev.off()

pdf("专利_figure4A.pdf")
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Ptprc"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Cd44"),cols = c("gray", "red"))#actin
dev.off()

pdf("专利_figure8AB.pdf")
DimPlot(YDL, reduction = "umap", label = TRUE, pt.size = 1.5,cols = c("#FB61D7","#53B400","#F8766D","#00B6EB","#A58AFF","#00C094","#C49A00"))
DimPlot(YDL, reduction = "umap", label = F, pt.size = 1.5,cols = c("#FB61D7","#53B400","#F8766D","#00B6EB","#A58AFF","#00C094","#C49A00"))
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("S100a8"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("S100a9"),cols = c("gray", "red"))#actin
YDL<-subset(YDL,idents = c("Ery/Mono"))
ElbowPlot(YDL)#选择top20个PC
pcSelect=20
YDL <- FindNeighbors(object = YDL, dims = 1:pcSelect)                #计算邻接距离
##接着优化模型,resolution参数决定下游聚类分析得到的分群数,对于3K左右的细胞,设为0.4-1.2 能得到较好的结果(官方说明)；如果数据量增大,该参数也应该适当增大。
YDL <- FindClusters(object = YDL, resolution = 0.2)                  #对细胞分组,优化标准模块化
##使用Idents（）函数可查看不同细胞的分群；
head(Idents(YDL), 5)
DimPlot(YDL, reduction = "umap", label = TRUE, pt.size = 1.5)
DimPlot(YDL, reduction = "umap", label = F, pt.size = 1.5)

VlnPlot(YDL, "S100a8",slot = "data",pt.size =0)
VlnPlot(YDL, "S100a9",slot = "data",pt.size =0)
VlnPlot(YDL, "Cd44",slot = "data",pt.size =0)
VlnPlot(YDL, "Ptprc",slot = "data",pt.size =0)

dev.off()






pdf("专利_figure8C生物过程打分.pdf")

nonpro1<-YDL
colors<-c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")
response_to_calcium_ion_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/response to calcium ion.csv",header = F)
response_to_calcium_ion_score<-response_to_calcium_ion_score[,1]
head(response_to_calcium_ion_score)
response_to_calcium_ion_score<-as.data.frame(response_to_calcium_ion_score)
colnames(response_to_calcium_ion_score)<-c("response_to_calcium_ion_score")
response_to_calcium_ion_score<-list(response_to_calcium_ion_score$response_to_calcium_ion_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = response_to_calcium_ion_score,
  name = "response_to_calcium_ion_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[19] <- 'response_to_calcium_ion_score' 
# VlnPlot(object = nonpro1 ,cols =colors,features = c("response_to_calcium_ion_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("response_to_calcium_ion_score"),pt.size = 0)
# VlnPlot(object = nonpro1 ,cols =colors,features = c("response_to_calcium_ion_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1 ,cols =colors,features = c("response_to_calcium_ion_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)


#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","response_to_calcium_ion_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","response_to_calcium_ion_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


cellular_response_to_calcium_ion_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/cellular response to calcium ion.csv",header = F)
cellular_response_to_calcium_ion_score<-toupper(cellular_response_to_calcium_ion_score[,1])
head(cellular_response_to_calcium_ion_score)
cellular_response_to_calcium_ion_score<-as.data.frame(cellular_response_to_calcium_ion_score)
colnames(cellular_response_to_calcium_ion_score)<-c("cellular_response_to_calcium_ion_score")
cellular_response_to_calcium_ion_score<-list(cellular_response_to_calcium_ion_score$cellular_response_to_calcium_ion_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = cellular_response_to_calcium_ion_score,
  name = "cellular_response_to_calcium_ion_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[20] <- 'cellular_response_to_calcium_ion_score' 
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

# 
# 
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","cellular_response_to_calcium_ion_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = cellular_response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","cellular_response_to_calcium_ion_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = cellular_response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




MAPK_cascade_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/MAPK cascade.csv",header = F)
MAPK_cascade_score<-MAPK_cascade_score[,1]
head(MAPK_cascade_score)
MAPK_cascade_score<-as.data.frame(MAPK_cascade_score)
colnames(MAPK_cascade_score)<-c("MAPK_cascade_score")
MAPK_cascade_score<-list(MAPK_cascade_score$MAPK_cascade_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = MAPK_cascade_score,
  name = "MAPK_cascade_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[21] <- 'MAPK_cascade_score' 

# VlnPlot(object = nonpro1 ,cols =colors,features = c("MAPK_cascade_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("MAPK_cascade_score"),pt.size = 0)
# VlnPlot(object = nonpro1 ,cols =colors,features = c("MAPK_cascade_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("MAPK_cascade_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1 ,cols =colors,features = c("MAPK_cascade_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("MAPK_cascade_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("MAPK_cascade_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","MAPK_cascade_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = MAPK_cascade_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","MAPK_cascade_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = MAPK_cascade_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



vesicle_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/vesicle.csv",header = F)
vesicle_score<-vesicle_score[,1]
head(vesicle_score)
vesicle_score<-as.data.frame(vesicle_score)
colnames(vesicle_score)<-c("vesicle_score")
vesicle_score<-list(vesicle_score$vesicle_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = vesicle_score,
  name = "vesicle_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[22] <- 'vesicle_score' 
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vesicle_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("vesicle_score"),pt.size = 0)
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vesicle_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("vesicle_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vesicle_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("vesicle_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("vesicle_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","vesicle_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = vesicle_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","vesicle_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = vesicle_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

actin_binding_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/actin_binding.csv",header = F)
actin_binding_score<-actin_binding_score[,1]
head(actin_binding_score)
actin_binding_score<-as.data.frame(actin_binding_score)
colnames(actin_binding_score)<-c("actin_binding_score")
actin_binding_score<-list(actin_binding_score$actin_binding_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = actin_binding_score,
  name = "actin_binding_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[23] <- 'actin_binding_score' 
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"),pt.size = 0)
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("actin_binding_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

# 
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","actin_binding_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = actin_binding_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","actin_binding_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = actin_binding_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




regulation_of_cell_shape_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/regulation_of_cell_shape.csv",header = F)
regulation_of_cell_shape_score<-regulation_of_cell_shape_score[,1]
head(regulation_of_cell_shape_score)
regulation_of_cell_shape_score<-as.data.frame(regulation_of_cell_shape_score)
colnames(regulation_of_cell_shape_score)<-c("regulation_of_cell_shape_score")
regulation_of_cell_shape_score<-list(regulation_of_cell_shape_score$regulation_of_cell_shape_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = regulation_of_cell_shape_score,
  name = "regulation_of_cell_shape_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[24] <- 'regulation_of_cell_shape_score' 
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"),pt.size = 0)
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("regulation_of_cell_shape_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","regulation_of_cell_shape_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = regulation_of_cell_shape_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","regulation_of_cell_shape_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = regulation_of_cell_shape_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



vacuole_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/vacuole.csv",header = F)
vacuole_score<-vacuole_score[,1]
head(vacuole_score)
vacuole_score<-as.data.frame(vacuole_score)
colnames(vacuole_score)<-c("vacuole_score")
vacuole_score<-list(vacuole_score$vacuole_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = vacuole_score,
  name = "vacuole_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[25] <- 'vacuole_score' 
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vacuole_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("vacuole_score"),pt.size = 0)
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vacuole_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("vacuole_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vacuole_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("vacuole_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("vacuole_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","vacuole_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = vacuole_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","vacuole_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = vacuole_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

YDL.markers <- FindAllMarkers(YDL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'roc')
# install.packages("magrittr") # package installations are only needed the first time you use it
# install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run evYDL time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
YDL.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
##存储marker

write.csv(YDL.markers,file="allmarker_celltype_cluster1.csv")
YDL.markers<-read.csv("allmarker_celltype_cluster1.csv")
#绘制分cluster的热图
top10 <- YDL.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#绘制marker在各个cluster的热图
pdf(file="专利_figure8D热图.pdf",width=8,height=13)
DoHeatmap(object = YDL, features = top10$gene) + NoLegend()
DoHeatmap(subset(YDL, downsample = 100), features = top10$gene, size = 3)+ NoLegend()
DoHeatmap(object = YDL, features = YDL.markers$gene) + NoLegend()
DoHeatmap(subset(YDL, downsample = 100), features = top10$gene, size = 3)+ NoLegend()
dev.off()


library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)



#合并四个时期的BP，以热图展示
a <- read.csv("allmarker_celltype_cluster1.csv")
b<-a[a$cluster=="0","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
E0 <- eg[,2]
E0<-as.data.frame(E0)
colnames(E0)<-c("cluster0")
head(E0)
#E0<-c(E0)

a <- read.csv("allmarker_celltype_cluster1.csv")
b<-a[a$cluster=="1","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
E1 <- eg[,2]
E1<-as.data.frame(E1)
colnames(E1)<-c("cluster1")
head(E1)
#E1<-c(E1)

a <- read.csv("allmarker_celltype_cluster1.csv")
b<-a[a$cluster=="2","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
E2 <- eg[,2]
E2<-as.data.frame(E2)
colnames(E2)<-c("cluster2")
head(E2)
#E2<-c(E2)

a <- read.csv("allmarker_celltype_cluster1.csv")
b<-a[a$cluster=="3","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
E3 <- eg[,2]
E3<-as.data.frame(E3)
colnames(E3)<-c("cluster3")
head(E3)
#E3<-c(E3)

a <- read.csv("allmarker_celltype_cluster1.csv")
b<-a[a$cluster=="4","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
E4 <- eg[,2]
E4<-as.data.frame(E4)
colnames(E4)<-c("cluster4")
head(E4)
#E4<-c(E4)

a <- read.csv("allmarker_celltype_cluster1.csv")
b<-a[a$cluster=="5","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
E5 <- eg[,2]
E5<-as.data.frame(E5)
colnames(E5)<-c("cluster5")
head(E5)
#E5<-c(E5)



data<-list(cluster0 =E0$cluster0,cluster1 =E1$cluster1, cluster2 = E2$cluster2, cluster3 =E3$cluster3,cluster4 =E4$cluster4,cluster5 =E5$cluster5)

lapply(data, head)


head(as.data.frame(ck))

ck <- compareCluster(geneCluster = data,OrgDb = org.Mm.eg.db, fun = "enrichGO", pvalueCutoff=0.05)
dotplot(ck, showCategory =10)
ego <- simplify(ck,cutoff=0.7,by="p.adjust",select_fun=min)
dotplot(ego, showCategory =15)

write.csv(ck,"ck_enrichGO.csv")


ck_enrichKEGG <- compareCluster(geneCluster = data, fun = "enrichKEGG",organism="mmu")
dotplot(ck_enrichKEGG, showCategory =25)
write.csv(ck_enrichKEGG,"ck_enrichKEGG.csv")

#visualize the result using dotplot method.

pdf("专利_figure8E富集分析.pdf",width=15,height=15)
dotplot(ego, showCategory =10)
dotplot(ck_enrichKEGG, showCategory =10)
dev.off()







#热图，全部差异基因都展示
#热图的bar颜色需要修改

#dotplot修改，加上脱核相关的基因
Idents(YDL) <- YDL@meta.data$celltype
YDL.markers <- FindAllMarkers(YDL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'roc')
# install.packages("magrittr") # package installations are only needed the first time you use it
# install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run evYDL time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
YDL.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
##存储marker

write.csv(YDL.markers,file="allmarker_celltype_nonpro12.csv")
YDL.markers<-read.csv("allmarker_celltype_nonpro12.csv")
#绘制分cluster的热图
top10 <- YDL.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#绘制marker在各个cluster的热图
pdf(file="专利_figure2D.pdf",width=8,height=13)
# DoHeatmap(object = YDL, features = top10$gene) + NoLegend()
# DoHeatmap(subset(YDL, downsample = 100), features = top10$gene, size = 3)+ NoLegend()
# DoHeatmap(object = YDL, features = YDL.markers$gene) + NoLegend()
# DoHeatmap(subset(YDL, downsample = 100), features = top10$gene, size = 3)+ NoLegend()


library(pheatmap)
library(ggplot2)
library(colorRamps)
library(RColorBrewer)
library(viridis)
library(cowplot)

DoHeatmap(object = YDL, features = top10$gene,group.colors= c("#FB61D7","#53B400","#F8766D","#00B6EB","#A58AFF","#00C094","#C49A00")) + NoLegend()+scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
DoHeatmap(object = YDL, features = top10$gene,draw.lines = F,group.colors= c("#FB61D7","#53B400","#F8766D","#00B6EB","#A58AFF","#00C094","#C49A00")) + NoLegend()+scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

DoHeatmap(object = YDL, features = YDL.markers$gene,draw.lines = T,group.colors= c("#FB61D7","#53B400","#F8766D","#00B6EB","#A58AFF","#00C094","#C49A00")) + NoLegend()+scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
DoHeatmap(object = YDL, features = YDL.markers$gene,draw.lines = F,group.colors= c("#FB61D7","#53B400","#F8766D","#00B6EB","#A58AFF","#00C094","#C49A00")) + NoLegend()+scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

dev.off()
table(YDL.markers$cluster)

feature<-c("Hbb-bs","Gata1","Nfe2","Xpo7","Tubb5","Tubb4b","Actb","Tmsb4x","Tmsb10","Vim","Jun","Sox4","Jund","Cebpb","Cd79a","Igkc","Ighm","Ly6d","Lyz2","Lgals3","Elane","S100a8")

pdf("专利_figure2E.pdf")

DotPlot(object = YDL, features =feature)+ RotatedAxis()
DotPlot(object = YDL, features =feature)+ RotatedAxis()+ coord_flip()

#加边框
DotPlot(object = YDL, features =feature)+ RotatedAxis()+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

dev.off()



setwd("D:/analysis/forpublication/summary/nonpro1_own_nonpro2")

marker0<-read.csv("GO_BP_marker0.csv")
marker0_data<-marker0[,c("ID","Description","GeneRatio")]

marker0_data$ID<-paste(marker0_data$ID,marker0_data$Description,sep="_")
marker0_data<-marker0_data[,c("ID","GeneRatio")]
marker0_data$GeneRatio<-sapply(marker0_data$GeneRatio, function(x) eval(parse(text=x)))
head(marker0_data)
dim(marker0_data)

marker1<-read.csv("GO_BP_marker1.csv")
marker1_data<-marker1[,c("ID","Description","GeneRatio")]

marker1_data$ID<-paste(marker1_data$ID,marker1_data$Description,sep="_")
marker1_data<-marker1_data[,c("ID","GeneRatio")]
marker1_data$GeneRatio<-sapply(marker1_data$GeneRatio, function(x) eval(parse(text=x)))
head(marker1_data)
dim(marker1_data)

marker2<-read.csv("GO_BP_marker2.csv")
marker2_data<-marker2[,c("ID","Description","GeneRatio")]

marker2_data$ID<-paste(marker2_data$ID,marker2_data$Description,sep="_")
marker2_data<-marker2_data[,c("ID","GeneRatio")]
marker2_data$GeneRatio<-sapply(marker2_data$GeneRatio, function(x) eval(parse(text=x)))
head(marker2_data)
dim(marker2_data)

marker3<-read.csv("GO_BP_marker3.csv")
marker3_data<-marker3[,c("ID","Description","GeneRatio")]

marker3_data$ID<-paste(marker3_data$ID,marker3_data$Description,sep="_")
marker3_data<-marker3_data[,c("ID","GeneRatio")]
marker3_data$GeneRatio<-sapply(marker3_data$GeneRatio, function(x) eval(parse(text=x)))
head(marker3_data)
dim(marker3_data)


marker4<-read.csv("GO_BP_marker4.csv")
marker4_data<-marker4[,c("ID","Description","GeneRatio")]

marker4_data$ID<-paste(marker4_data$ID,marker4_data$Description,sep="_")
marker4_data<-marker4_data[,c("ID","GeneRatio")]
marker4_data$GeneRatio<-sapply(marker4_data$GeneRatio, function(x) eval(parse(text=x)))
head(marker4_data)
dim(marker4_data)


marker5<-read.csv("GO_BP_marker5.csv")
marker5_data<-marker5[,c("ID","Description","GeneRatio")]

marker5_data$ID<-paste(marker5_data$ID,marker5_data$Description,sep="_")
marker5_data<-marker5_data[,c("ID","GeneRatio")]
marker5_data$GeneRatio<-sapply(marker5_data$GeneRatio, function(x) eval(parse(text=x)))
head(marker5_data)
dim(marker5_data)


marker6<-read.csv("GO_BP_marker6.csv")
marker6_data<-marker6[,c("ID","Description","GeneRatio")]

marker6_data$ID<-paste(marker6_data$ID,marker6_data$Description,sep="_")
marker6_data<-marker6_data[,c("ID","GeneRatio")]
marker6_data$GeneRatio<-sapply(marker6_data$GeneRatio, function(x) eval(parse(text=x)))
head(marker6_data)
dim(marker6_data)





#合并不同的矩阵，列名弄成一样长，然后补上0值。

ID<-c(marker0_data$ID,marker1_data$ID,marker2_data$ID,marker3_data$ID,marker4_data$ID,marker5_data$ID,marker6_data$ID)
name<-as.data.frame(ID)
name<-name[!duplicated(name$ID),] #删掉所有列上都重复的
name<-as.data.frame(name)
colnames(name)<-c("ID")
head(name)



DATA_cluster0<-dplyr::bind_rows(name,marker0_data)
DATA_cluster0<-DATA_cluster0[order(DATA_cluster0[,2],decreasing=T),]#以第一列降序排列
DATA_cluster0<-DATA_cluster0[!duplicated(DATA_cluster0$ID),] #删掉所有列上都重复的
DATA_cluster0[is.na(DATA_cluster0)] <- 0#把NA值全部替换为0
colnames(DATA_cluster0)<-c("ID","cluster0")
head(DATA_cluster0)


DATA_cluster1<-dplyr::bind_rows(name,marker1_data)
DATA_cluster1<-DATA_cluster1[order(DATA_cluster1[,2],decreasing=T),]#以第一列降序排列
DATA_cluster1<-DATA_cluster1[!duplicated(DATA_cluster1$ID),] #删掉所有列上都重复的
DATA_cluster1[is.na(DATA_cluster1)] <- 0#把NA值全部替换为0
colnames(DATA_cluster1)<-c("ID","cluster1")
head(DATA_cluster1)


DATA_cluster2<-dplyr::bind_rows(name,marker2_data)
DATA_cluster2<-DATA_cluster2[order(DATA_cluster2[,2],decreasing=T),]#以第一列降序排列
DATA_cluster2<-DATA_cluster2[!duplicated(DATA_cluster2$ID),] #删掉所有列上都重复的
DATA_cluster2[is.na(DATA_cluster2)] <- 0#把NA值全部替换为0
colnames(DATA_cluster2)<-c("ID","cluster2")
head(DATA_cluster2)


DATA_cluster3<-dplyr::bind_rows(name,marker3_data)
DATA_cluster3<-DATA_cluster3[order(DATA_cluster3[,2],decreasing=T),]#以第一列降序排列
DATA_cluster3<-DATA_cluster3[!duplicated(DATA_cluster3$ID),] #删掉所有列上都重复的
DATA_cluster3[is.na(DATA_cluster3)] <- 0#把NA值全部替换为0
colnames(DATA_cluster3)<-c("ID","cluster3")
head(DATA_cluster3)


DATA_cluster4<-dplyr::bind_rows(name,marker4_data)
DATA_cluster4<-DATA_cluster4[order(DATA_cluster4[,2],decreasing=T),]#以第一列降序排列
DATA_cluster4<-DATA_cluster4[!duplicated(DATA_cluster4$ID),] #删掉所有列上都重复的
DATA_cluster4[is.na(DATA_cluster4)] <- 0#把NA值全部替换为0
colnames(DATA_cluster4)<-c("ID","cluster4")
head(DATA_cluster4)


DATA_cluster5<-dplyr::bind_rows(name,marker5_data)
DATA_cluster5<-DATA_cluster5[order(DATA_cluster5[,2],decreasing=T),]#以第一列降序排列
DATA_cluster5<-DATA_cluster5[!duplicated(DATA_cluster5$ID),] #删掉所有列上都重复的
DATA_cluster5[is.na(DATA_cluster5)] <- 0#把NA值全部替换为0
colnames(DATA_cluster5)<-c("ID","cluster5")
head(DATA_cluster5)


DATA_cluster6<-dplyr::bind_rows(name,marker6_data)
DATA_cluster6<-DATA_cluster6[order(DATA_cluster6[,2],decreasing=T),]#以第一列降序排列
DATA_cluster6<-DATA_cluster6[!duplicated(DATA_cluster6$ID),] #删掉所有列上都重复的
DATA_cluster6[is.na(DATA_cluster6)] <- 0#把NA值全部替换为0
colnames(DATA_cluster6)<-c("ID","cluster6")
head(DATA_cluster6)





cluster01<-merge(DATA_cluster0,DATA_cluster1,by=c("ID"))
cluster23<-merge(DATA_cluster2,DATA_cluster3,by=c("ID"))
cluster45<-merge(DATA_cluster4,DATA_cluster5,by=c("ID"))

cluster0123<-merge(cluster01,cluster23,by=c("ID"))
cluster456<-merge(cluster45,DATA_cluster6,by=c("ID"))

DATA<-merge(cluster0123,cluster456,by=c("ID"))
# DATA<-DATA[,c("ID","cluster0","cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","cluster7")]
# 
# dim(DATA)
# head(DATA)
# DATA<-DATA[!duplicated(DATA$ID),] #删掉所有列上都重复的
rownames(DATA)<-DATA[,1]
DATA<-DATA[,-1]
DATA<-as.data.frame(DATA)
# https://blog.csdn.net/c1z2w3456789/article/details/79467095



DATA<-read.csv("D:/analysis/forpublication/summary/nonpro1_own_nonpro2/results/GO_7cluster.csv",row.names = 1)

head(DATA)
colnames(DATA)<- c(
  "BasoE2",
  "Ery/Mono",
  "BasoE1",
  "Ery/B",
  "CD44hiOrthoE1",
  "CD44hiOrthoE2","ProE")
#data<-DATA[c(3,0,2,5,4,6,1,7)]
data<-DATA[c("ProE","BasoE1","BasoE2","CD44hiOrthoE1","CD44hiOrthoE2","Ery/B","Ery/Mono")]
write.csv(data,"GO_nonpro12.csv")

#选择特定的GO term作图
term<-c("GO:0034101_erythrocyte homeostasis","GO:0030218_erythrocyte differentiation",
        "GO:0000280_nuclear division",
        "GO:0044770_cell cycle phase transition",
        "GO:0006281_DNA repair",
        "GO:0061572_actin filament bundle organization",
        "GO:0030041_actin filament polymerization",
        "GO:0045010_actin nucleation","GO:0017038_protein import","GO:0050657_nucleic acid transport",
        "GO:0050658_RNA transport",
        "GO:0006611_protein export from nucleus",
        "GO:0022613_ribonucleoprotein complex biogenesis",
        "GO:0006364_rRNA processing",
        "GO:0006457_protein folding",
        "GO:0002181_cytoplasmic translation",
        "GO:0016032_viral process",
        "GO:0045088_regulation of innate immune response",
        "GO:0035455_response to interferon-alpha",
        "GO:0035456_response to interferon-beta",
        "GO:0034341_response to interferon-gamma",
        "GO:0006900_vesicle budding from membrane",
        "GO:0048193_Golgi vesicle transport",
        "GO:0030098_lymphocyte differentiation",
        "GO:0022409_positive regulation of cell-cell adhesion",
        "GO:0007264_small GTPase mediated signal transduction",
        "GO:0008360_regulation of cell shape",
        "GO:0002274_myeloid leukocyte activation",
        "GO:0051235_maintenance of location",
        "GO:0050900_leukocyte migration",
        "GO:0006909_phagocytosis")

term<-c("GO:0034101_erythrocyte homeostasis","GO:0030218_erythrocyte differentiation",
        "GO:0000280_nuclear division",
        "GO:0044770_cell cycle phase transition",
        "GO:0006281_DNA repair",
        "GO:0061572_actin filament bundle organization",
        "GO:0030041_actin filament polymerization",
        "GO:0045010_actin nucleation","GO:0017038_protein import","GO:0050657_nucleic acid transport",
        "GO:0050658_RNA transport",
        "GO:0006611_protein export from nucleus",
        "GO:0022613_ribonucleoprotein complex biogenesis",
        "GO:0006364_rRNA processing",
        "GO:0006457_protein folding",
        "GO:0002181_cytoplasmic translation",
        "GO:0016032_viral process",
        "GO:0045088_regulation of innate immune response",
        "GO:0035455_response to interferon-alpha",
        "GO:0035456_response to interferon-beta",
        "GO:0034341_response to interferon-gamma",
        "GO:0006900_vesicle budding from membrane",
        "GO:0048193_Golgi vesicle transport",
        "GO:0030098_lymphocyte differentiation",
        "GO:0022409_positive regulation of cell-cell adhesion",
        "GO:0007264_small GTPase mediated signal transduction",
        "GO:0008360_regulation of cell shape",
        "GO:0002274_myeloid leukocyte activation",
        "GO:0051235_maintenance of location",
        "GO:0050900_leukocyte migration",
        "GO:0006909_phagocytosis")



#删掉所有列上都重复的
newdata<-data[c(term),]
newdata<-na.omit(newdata)
#低值为蓝色，高值为红色，中间值为白色：
#pdf("fig1g.pdf")
pheatmap(newdata,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))
pheatmap(newdata,filename = "fig1g_new.pdf",
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))
pheatmap(newdata,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(newdata,scale = "none",filename = "fig1g_new.pdf",
         cluster rows = F,border color = NA,
         cluster cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))



newdata<-read.csv("GO_7cluster_BP.csv",row.names = 1)
newdata_MF<-read.csv("GO_7cluster_MF.csv",row.names = 1)

#合并BP和MF
DATA_BP_MF<-rbind(newdata,newdata_MF)
head(DATA_BP_MF)
colnames(DATA_BP_MF)<- c(
  "BasoE2",
  "Ery/Mono",
  "BasoE1",
  "Ery/B",
  "CD44hiOrthoE1",
  "CD44hiOrthoE2","ProE")
#data<-DATA[c(3,0,2,5,4,6,1,7)]
data<-DATA_BP_MF[c("ProE","BasoE1","BasoE2","CD44hiOrthoE1","CD44hiOrthoE2","Ery/B","Ery/Mono")]

#选择特定的GO term作图
term<-c("GO:0044389_ubiquitin-like protein ligase binding","GO:0030218_erythrocyte differentiation",
        "GO:0016209_antioxidant activity",
        "GO:0000280_nuclear division",
        "GO:0044770_cell cycle phase transition",
        "GO:0006281_DNA repair",
        "GO:0050658_RNA transport",
        "GO:0003712_transcription coregulator activity",
        "GO:0035326_enhancer binding",
        "GO:0006611_protein export from nucleus",
        "GO:0003735_structural constituent of ribosome",
        "GO:0006457_protein folding",
        "GO:0016032_viral process",
        "GO:0045088_regulation of innate immune response",
        "GO:0035455_response to interferon-alpha",
        "GO:0034341_response to interferon-gamma",
        "GO:0006900_vesicle budding from membrane",
        "GO:0030041_actin filament polymerization",
        "GO:0045010_actin nucleation",
        "GO:0030098_lymphocyte differentiation",
        "GO:0022409_positive regulation of cell-cell adhesion",
        "GO:0007264_small GTPase mediated signal transduction",
        "GO:0008360_regulation of cell shape",
        "GO:0051235_maintenance of location",
        "GO:0050900_leukocyte migration",
        "GO:0006909_phagocytosis",
        "GO:0044548_S100 protein binding",
        "GO:0035325_Toll-like receptor binding")


term<-c("GO:0044389_ubiquitin-like protein ligase binding",
        "GO:0030218_erythrocyte differentiation",
        "GO:0016209_antioxidant activity",
        "GO:0000280_nuclear division",
        "GO:0044770_cell cycle phase transition",
        "GO:0006281_DNA repair",
        "GO:0050658_RNA transport",
        "GO:0003712_transcription coregulator activity",
        "GO:0006611_protein export from nucleus",
        "GO:0003735_structural constituent of ribosome",
        "GO:0006457_protein folding",
        "GO:0016032_viral process",
        "GO:0045088_regulation of innate immune response",
        "GO:0030041_actin filament polymerization",
        "GO:0008360_regulation of cell shape",
        "GO:0030098_lymphocyte differentiation",
        "GO:0022409_positive regulation of cell−cell adhesion",
        "GO:0050900_leukocyte migration",
        "GO:0019882_antigen processing and presentation",
        "GO:0006909_phagocytosis",
        "GO:0050727_regulation of inflammatory response",
        "GO:0007264_small GTPase mediated signal transduction")

#删掉所有列上都重复的
DATA_BP_MF<-data[c(term),]
DATA_BP_MF<-na.omit(DATA_BP_MF)
pheatmap(DATA_BP_MF,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))
pheatmap(DATA_BP_MF, legend_labels = c("0","0.1", "0.2", "0.3"),
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(DATA_BP_MF,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))
pheatmap(DATA_BP_MF,filename = "专利_figure2F.pdf",width = 7,height = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))


pheatmap(newdata,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(newdata,scale = "none",filename = "fig1g_new.pdf",
         cluster rows = F,border color = NA,
         cluster cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))




pheatmap(DATA_BP_MF,filename = "fig1g_bp_mf1.pdf",
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))
pheatmap(DATA_BP_MF,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(DATA_BP_MF,scale = "none",filename = "fig1g_new.pdf",
         cluster rows = F,border color = NA,
         cluster cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))
dev.off()















load("D:/analysis/bone_marrow_niche/NicheData10x.rda")

NicheData10x[["celltype"]]<-NicheData10x@active.ident
NicheData10x[["orig.ident"]]<-NicheData10x@meta.data$metadata....experiment..

immunecells<-subset(NicheData10x,ident = c("B cell","Dendritic cells","Monocytes","NK cells","Neutrophils","T cells"))
DimPlot(immunecells,reduction = "tsne",label = TRUE,pt.size = 1.5)
immunecells.markers <- FindAllMarkers(immunecells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'roc')
# install.packages("magrittr") # package installations are only needed the first time you use it
# install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run evYDL time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
immunecells.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
##存储marker

write.csv(immunecells.markers,file="allmarker_immunecells.csv")
#绘制分cluster的热图
top10 <- immunecells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#绘制marker在各个cluster的热图
pdf(file="tsneHeatmap_immunecells.pdf",width=12,height=9)
DoHeatmap(object = immunecells, features = top10$gene) + NoLegend()
DoHeatmap(subset(immunecells, downsample = 100), features = top10$gene, size = 3)+ NoLegend()
dev.off()




YDL.AVERAGE<-AverageExpression(object =YDL,return.seurat=F)
YDL.AVERAGE<-as.data.frame(YDL.AVERAGE)
write.csv(YDL.AVERAGE,"nonpro12.AVERAGE.csv")
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(immunecells,reduction = "tsne",label = TRUE,pt.size = 1.5)


immunecells.AVERAGE<-AverageExpression(object =immunecells,return.seurat=F)
immunecells.AVERAGE<-as.data.frame(immunecells.AVERAGE)



head(immunecells.AVERAGE)
head(YDL.AVERAGE)
dim(immunecells.AVERAGE)
dim(YDL.AVERAGE)
head(c(rownames(immunecells.AVERAGE)))
#让矩阵变成一样大小，但是损失部分红细胞信息
#重新运行代码，保留所有的基因。
YDL.AVERAGE<-YDL.AVERAGE[c(rownames(immunecells.AVERAGE)),]
# data<-cbind(immunecells.AVERAGE,YDL.AVERAGE)

data<-cbind(immunecells.AVERAGE,YDL.AVERAGE)
data<-as.data.frame(data)


#data<-na.omit(data)
data[is.na(data)] <- 0#把NA值全部替换为0；
head(data)
dim(data)


library(pheatmap)
library(ggplot2)
library(colorRamps)
library(RColorBrewer)
library(viridis)
library(cowplot)

mat<-data
#3，0，2，5，4，6，1，7
colnames(mat)<-c("B cell","Dendritic cells","Monocytes","NK cells","Neutrophils","T cells","ProE","BasoE1",
                 "BasoE2","CD44hiOrthoE1","CD44hiOrthoE2","Ery/B","Ery/Mono")
mat<-mat[c("B cell","Dendritic cells","Monocytes","NK cells","Neutrophils","T cells","Ery/B","Ery/Mono")]
mat<-mat[c(unique(top10$gene),"Hbb-a1","Hbb-a2","Hbb-bs","Hbb-bt","Gata1","Klf1","Tal1","Nfe2","Xpo7"),]
mat<-na.omit(mat)

pheatmap(
  mat   = log10(mat+1),
  scale = "row",color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  cluster_rows = F,
  cluster_cols = F,
  border_color = NA,
  #annotation_col = annotation_col,
  #annotation_row = annotation_row,
  show_colnames = TRUE,
  show_rownames = TRUE,
  drop_levels   = TRUE,
  fontsize  = 8,filename = "专利_figure3A.pdf",width = 7,height = 9
)
pheatmap(
  mat   = mat,
  scale = "row",color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  cluster_rows = F,
  cluster_cols = F,
  border_color = NA,
  #annotation_col = annotation_col,
  #annotation_row = annotation_row,
  show_colnames = TRUE,
  show_rownames = TRUE,
  drop_levels   = TRUE,
  fontsize  = 8#,filename = "fig2h.pdf"
)

dev.off()
pheatmap(
  mat   = log10(mat+1),
  scale = "row",
  cluster_rows = F,
  cluster_cols = T,
  border_color = NA,
  #annotation_col = annotation_col,
  #annotation_row = annotation_row,
  show_colnames = TRUE,
  show_rownames = TRUE,
  drop_levels   = TRUE,
  fontsize  = 8,
  #main = "Pheatmap",
  filename = "fig2h_聚类.pdf"
)

pheatmap(
  mat   = log10(mat+1),
  scale = "row",
  cluster_rows = F,
  cluster_cols = T,
  border_color = NA,
  #annotation_col = annotation_col,
  #annotation_row = annotation_row,
  show_colnames = TRUE,
  show_rownames = TRUE,
  drop_levels   = TRUE,
  fontsize  = 8,
  main = "Pheatmap",filename = "heatmap_1.pdf"
)


pheatmap(
  mat   = log10(mat+1),
  scale = "row",
  cluster_rows = F,
  cluster_cols = F,
  border_color = NA,
  #annotation_col = annotation_col,
  #annotation_row = annotation_row,
  show_colnames = TRUE,
  show_rownames = TRUE,
  drop_levels   = TRUE,
  fontsize  = 8,
  #main = "Pheatmap",
  filename = "heatmap_top10.pdf"
)
ggsave("heatmap.pdf" , device ="pdf" ,width = 8, height = 9,units = c( "cm"))

dev.off()


pdf("专利_figure3B.pdf")


nonpro<-read.csv("D:/analysis/yanxingcell/nonpro细胞的富集分析.csv")
head(nonpro)
immune<-read.csv("D:/analysis/yanxingcell/炎性细胞的富集分析.csv")
head(immune)


merge_tf0<-nonpro$ID
merge_tf1<-DATA2$ID


library (ggvenn)
x<-list("Immune cells"=immune$X,"CD45-Nonpro"=nonpro$X
        )
ggvenn(x)
pdf("专利_figure3B.pdf")
ggvenn(x)
dev.off()


pdf("专利_figure3C.pdf")
all<-read.csv("D:/analysis/yanxingcell/结合免疫细胞的炎性红细胞的富集分析结果.csv",row.names=1)
head(all)
colnames(all)<-c("B cells","T cells","NK cells","Dendritic cells","Monocytes","Neutrophils","BasoE2",
                 "Ery/Mono",
                 "BasoE1",
                 "Ery/B",
                 "CD44hiOrthoE1",
                 "CD44hiOrthoE2","ProE")

data<-all
data<-data[which(rowSums(data13) > 0),] 
pheatmap::pheatmap(data,show_rownames = F,#filename = "结合炎性细胞的炎性红细胞.pdf",width = 20,height = 300,
                   cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("white","red","darkred"))(100))
pheatmap::pheatmap(data,filename = "专利_figure3结合炎性细胞的炎性红细胞的富集分析.pdf",width = 20,height = 300,
                   cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("white","red","darkred"))(100))
pheatmap::pheatmap(data,filename = "结合炎性细胞的炎性红细胞的富集分析.jpg",width = 20,height = 300,
                   cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("white","red","darkred"))(100))

write.csv(data,"专利_figure3BC对应的富集分析结果.csv")
pdf("专利_figure4A分.pdf")


rm(nonpro1)
nonpro1<-YDL

pdf("专利_figure3D生物过程打分.pdf")


response_to_calcium_ion_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/response to calcium ion.csv",header = F)
response_to_calcium_ion_score<-response_to_calcium_ion_score[,1]
head(response_to_calcium_ion_score)
response_to_calcium_ion_score<-as.data.frame(response_to_calcium_ion_score)
colnames(response_to_calcium_ion_score)<-c("response_to_calcium_ion_score")
response_to_calcium_ion_score<-list(response_to_calcium_ion_score$response_to_calcium_ion_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = response_to_calcium_ion_score,
  name = "response_to_calcium_ion_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[19] <- 'response_to_calcium_ion_score' 
# VlnPlot(object = nonpro1 ,cols =colors,features = c("response_to_calcium_ion_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("response_to_calcium_ion_score"),pt.size = 0)
# VlnPlot(object = nonpro1 ,cols =colors,features = c("response_to_calcium_ion_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1 ,cols =colors,features = c("response_to_calcium_ion_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)


#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","response_to_calcium_ion_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","response_to_calcium_ion_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


cellular_response_to_calcium_ion_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/cellular response to calcium ion.csv",header = F)
cellular_response_to_calcium_ion_score<-toupper(cellular_response_to_calcium_ion_score[,1])
head(cellular_response_to_calcium_ion_score)
cellular_response_to_calcium_ion_score<-as.data.frame(cellular_response_to_calcium_ion_score)
colnames(cellular_response_to_calcium_ion_score)<-c("cellular_response_to_calcium_ion_score")
cellular_response_to_calcium_ion_score<-list(cellular_response_to_calcium_ion_score$cellular_response_to_calcium_ion_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = cellular_response_to_calcium_ion_score,
  name = "cellular_response_to_calcium_ion_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[20] <- 'cellular_response_to_calcium_ion_score' 
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

# 
# 
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","cellular_response_to_calcium_ion_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = cellular_response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","cellular_response_to_calcium_ion_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = cellular_response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




MAPK_cascade_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/MAPK cascade.csv",header = F)
MAPK_cascade_score<-MAPK_cascade_score[,1]
head(MAPK_cascade_score)
MAPK_cascade_score<-as.data.frame(MAPK_cascade_score)
colnames(MAPK_cascade_score)<-c("MAPK_cascade_score")
MAPK_cascade_score<-list(MAPK_cascade_score$MAPK_cascade_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = MAPK_cascade_score,
  name = "MAPK_cascade_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[21] <- 'MAPK_cascade_score' 

# VlnPlot(object = nonpro1 ,cols =colors,features = c("MAPK_cascade_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("MAPK_cascade_score"),pt.size = 0)
# VlnPlot(object = nonpro1 ,cols =colors,features = c("MAPK_cascade_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("MAPK_cascade_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1 ,cols =colors,features = c("MAPK_cascade_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("MAPK_cascade_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("MAPK_cascade_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","MAPK_cascade_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = MAPK_cascade_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","MAPK_cascade_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = MAPK_cascade_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



vesicle_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/vesicle.csv",header = F)
vesicle_score<-vesicle_score[,1]
head(vesicle_score)
vesicle_score<-as.data.frame(vesicle_score)
colnames(vesicle_score)<-c("vesicle_score")
vesicle_score<-list(vesicle_score$vesicle_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = vesicle_score,
  name = "vesicle_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[22] <- 'vesicle_score' 
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vesicle_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("vesicle_score"),pt.size = 0)
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vesicle_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("vesicle_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vesicle_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("vesicle_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("vesicle_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","vesicle_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = vesicle_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","vesicle_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = vesicle_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

actin_binding_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/actin_binding.csv",header = F)
actin_binding_score<-actin_binding_score[,1]
head(actin_binding_score)
actin_binding_score<-as.data.frame(actin_binding_score)
colnames(actin_binding_score)<-c("actin_binding_score")
actin_binding_score<-list(actin_binding_score$actin_binding_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = actin_binding_score,
  name = "actin_binding_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[23] <- 'actin_binding_score' 
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"),pt.size = 0)
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("actin_binding_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

# 
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","actin_binding_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = actin_binding_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","actin_binding_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = actin_binding_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




regulation_of_cell_shape_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/regulation_of_cell_shape.csv",header = F)
regulation_of_cell_shape_score<-regulation_of_cell_shape_score[,1]
head(regulation_of_cell_shape_score)
regulation_of_cell_shape_score<-as.data.frame(regulation_of_cell_shape_score)
colnames(regulation_of_cell_shape_score)<-c("regulation_of_cell_shape_score")
regulation_of_cell_shape_score<-list(regulation_of_cell_shape_score$regulation_of_cell_shape_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = regulation_of_cell_shape_score,
  name = "regulation_of_cell_shape_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[24] <- 'regulation_of_cell_shape_score' 
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"),pt.size = 0)
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("regulation_of_cell_shape_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","regulation_of_cell_shape_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = regulation_of_cell_shape_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","regulation_of_cell_shape_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = regulation_of_cell_shape_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



vacuole_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/vacuole.csv",header = F)
vacuole_score<-vacuole_score[,1]
head(vacuole_score)
vacuole_score<-as.data.frame(vacuole_score)
colnames(vacuole_score)<-c("vacuole_score")
vacuole_score<-list(vacuole_score$vacuole_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = vacuole_score,
  name = "vacuole_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[25] <- 'vacuole_score' 
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vacuole_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("vacuole_score"),pt.size = 0)
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vacuole_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("vacuole_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vacuole_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("vacuole_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("vacuole_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","vacuole_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = vacuole_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","vacuole_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = vacuole_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()




YDL<-readRDS("D:/analysis/forpublication/summary/nonpro1_own_nonpro2_nonpro3_ortho/nonpro1_nonpro2_nonpro3_SINGLET.RDS")
current.cluster.ids <- c(0, 1, 2, 3, 4, 5,6,7,8)
new.cluster.ids <- c(
  "Ortho",
  "EB4",
  "Ery/Mono1",
  "EB3",
  "EB2","Ery/B","EB1",
  "Ery/ApoE+","Ery/Mono2")
names(new.cluster.ids) <- levels(YDL)
YDL <- RenameIdents(YDL, new.cluster.ids)
DimPlot(YDL, reduction = "tsne", label = TRUE, pt.size = 1.5)
DimPlot(YDL, reduction = "umap", label = TRUE, pt.size = 1.5)


library("Seurat")
library("ggplot2")
gg <- TSNEPlot(YDL)
col<-ggplot_build(gg)$data
col<-as.data.frame(col)
table(col$colour)
table(YDL$seurat_clusters)

a<-as.data.frame(table(col$colour))
b<-as.data.frame(table(YDL$seurat_clusters))
c<-merge(a,b,by = "Freq")
c<-c[order(c$Freq,decreasing = T),]

YDL$celltype<-Idents(YDL)
YDL$celltype<-factor(YDL$celltype,levels = c("EB1","EB2","EB3","EB4","Ortho","Ery/Mono1","Ery/Mono2","Ery/B","Ery/ApoE+"))
Idents(YDL)<-YDL$celltype

pdf("专利_figure5B.pdf")
cols = c("#619CFF","#00C19F","#00BA38","#D39200","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB")

DimPlot(YDL, reduction = "umap", label = TRUE, pt.size = 1.5,cols = cols)
DimPlot(YDL, reduction = "umap", label = F, pt.size = 1.5,cols = c("#619CFF","#00C19F","#00BA38","#D39200","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"))
DimPlot(YDL, reduction = "umap", label = T, pt.size = 1.5,cols = c("#619CFF","#00C19F","#00BA38","#D39200","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"),split.by = "orig.ident")
DimPlot(YDL, reduction = "umap", label = F, pt.size = 1.5,cols = c("#619CFF","#00C19F","#00BA38","#D39200","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"),split.by = "orig.ident")


cols = c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB")

DimPlot(YDL, reduction = "umap", label = TRUE, pt.size = 1.5,cols = cols)
DimPlot(YDL, reduction = "umap", label = F, pt.size = 1.5,cols = c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"))
DimPlot(YDL, reduction = "umap", label = T, pt.size = 1.5,cols = c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"),split.by = "orig.ident")
DimPlot(YDL, reduction = "umap", label = F, pt.size = 1.5,cols = c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"),split.by = "orig.ident")


# Color palette
#colors <- c("#F8766D","#E68613","#CD9600","#ABA300","#7CAE00","#0CB702","#00BE67","#00C19A","#00BFC4","#00B8E7","#00A9FF","#8494FF","#C77CFF","#ED68ED","#FF61CC","#FF68A1")
colors <- c("#619CFF","#00C19F","#00BA38","#D39200","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB")

DimPlot(YDL,cols =colors,reduction="umap")
a<-as.data.frame(table(YDL@meta.data$celltype))
# 数据准备
info = a$Freq
# 命名
names = as.character(a$Var1)
#names = c("0","1","2","3","4","5","6","7")
# 涂色（可选）
#cols = c("#ED1C24","#22B14C","#FFC90E","#3f48CC","#3f90CC","#22B17C","#FFC93E")
# 计算百分比
piepercent = paste(round(100*info/sum(info),2), "%")
# 绘图
pie(info, labels=piepercent, main = "total cluster ratio", col=colors, family='GB1')
# 添加颜色样本标注
legend("topright", names, cex=0.6, fill=colors)


library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","nCount_RNA"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = log(nCount_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","nFeature_RNA"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = log(nFeature_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

dev.off()

pdf("专利_figure5C.pdf",width = 14,height = 5)
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Hbb-bs","Xpo7","Gata1"),cols = c("gray", "red"),ncol = 3)#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Cd79a","S100a8","Ptprc"),cols = c("gray", "red"),ncol = 3)#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("C1qa","Apoe","Vcam1"),cols = c("gray", "red"),ncol = 3)#actin
dev.off()



pdf("专利_figure5D-1.pdf")
# Color palette
colors <- c("#619CFF","#00C19F","#00BA38","#D39200","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB")

DimPlot(YDL, pt.size = 1.5,cols =colors,reduction="umap")
table(YDL@meta.data$celltype)
a<-as.data.frame(table(YDL@meta.data$celltype))
# 数据准备
info = a$Freq
# 命名
names = as.character(a$Var1)
#names = c("0","1","2","3","4","5","6","7")
# 涂色（可选）
#cols = c("#ED1C24","#22B14C","#FFC90E","#3f48CC","#3f90CC","#22B17C","#FFC93E")
# 计算百分比
piepercent = paste(round(100*info/sum(info),2), "%")
# 绘图
pie(info, labels=piepercent, main = "total cluster ratio", col=colors, family='GB1')
# 添加颜色样本标注
legend("topright", names, cex=0.6, fill=colors)



table(YDL@meta.data$orig.ident)

nonpro1<-subset(YDL, cells= rownames(YDL@meta.data[YDL@meta.data$orig.ident=="nonpro1",]))
DimPlot(nonpro1, reduction = "umap", label = T, pt.size = 1.5,cols = c("#619CFF","#00C19F","#00BA38","#D39200","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"),split.by = "orig.ident")
DimPlot(nonpro1, reduction = "umap", label = F, pt.size = 1.5,cols = c("#619CFF","#00C19F","#00BA38","#D39200","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"),split.by = "orig.ident")
table(nonpro1@meta.data$celltype)
a<-as.data.frame(table(nonpro1@meta.data$celltype))
rm(nonpro1)
# 数据准备
info = a$Freq
# 命名
names = as.character(a$Var1)
#names = c("0","1","2","3","4","5","6","7")
# 涂色（可选）
#cols = c("#ED1C24","#22B14C","#FFC90E","#3f48CC","#3f90CC","#22B17C","#FFC93E")
# 计算百分比
piepercent = paste(round(100*info/sum(info)), "%")
# 绘图
pie(info, labels=piepercent, main = "nonpro1 cluster ratio", col=colors, family='GB1')
# 添加颜色样本标注
legend("topright", names, cex=0.6, fill=colors)

nonpro2<-subset(YDL, cells= rownames(YDL@meta.data[YDL@meta.data$orig.ident=="nonpro2",]))
DimPlot(nonpro2, reduction = "umap", label = T, pt.size = 1.5,cols = c("#619CFF","#00C19F","#00BA38","#D39200","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"),split.by = "orig.ident")
DimPlot(nonpro2, reduction = "umap", label = F, pt.size = 1.5,cols = c("#619CFF","#00C19F","#00BA38","#D39200","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"),split.by = "orig.ident")
table(nonpro2@meta.data$celltype)
a<-as.data.frame(table(nonpro2@meta.data$celltype))
rm(nonpro2)
# 数据准备
info = a$Freq
# 命名
names = as.character(a$Var1)
#names = c("0","1","2","3","4","5","6","7")
# 涂色（可选）
#cols = c("#ED1C24","#22B14C","#FFC90E","#3f48CC","#3f90CC","#22B17C","#FFC93E")
# 计算百分比
piepercent = paste(round(100*info/sum(info)), "%")
# 绘图
pie(info, labels=piepercent, main = "nonpro2 cluster ratio", col=colors, family='GB1')
# 添加颜色样本标注
legend("topright", names, cex=0.6, fill=colors)

nonpro3<-subset(YDL, cells= rownames(YDL@meta.data[YDL@meta.data$orig.ident=="nonpro3",]))
DimPlot(nonpro3, reduction = "umap", label = T, pt.size = 1.5,cols = c("#619CFF","#00C19F","#00BA38","#D39200","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"),split.by = "orig.ident")
DimPlot(nonpro3, reduction = "umap", label = F, pt.size = 1.5,cols = c("#619CFF","#00C19F","#00BA38","#D39200","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"),split.by = "orig.ident")
table(nonpro3@meta.data$celltype)
a<-as.data.frame(table(nonpro3@meta.data$celltype))
rm(nonpro3)
# 数据准备
info = a$Freq
# 命名
names = as.character(a$Var1)
#names = c("0","1","2","3","4","5","6","7")
# 涂色（可选）
#cols = c("#ED1C24","#22B14C","#FFC90E","#3f48CC","#3f90CC","#22B17C","#FFC93E")
# 计算百分比
piepercent = paste(round(100*info/sum(info)), "%")
# 绘图
pie(info, labels=piepercent, main = "nonpro3 cluster ratio", col=colors, family='GB1')
# 添加颜色样本标注
legend("topright", names, cex=0.6, fill=colors)

orth<-subset(YDL, cells= rownames(YDL@meta.data[YDL@meta.data$orig.ident=="orth",]))
DimPlot(orth, reduction = "umap", label = T, pt.size = 1.5,cols = c("#619CFF","#00C19F","#00BA38","#D39200","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"),split.by = "orig.ident")
DimPlot(orth, reduction = "umap", label = F, pt.size = 1.5,cols = c("#619CFF","#00C19F","#00BA38","#D39200","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"),split.by = "orig.ident")
table(orth@meta.data$celltype)
a<-as.data.frame(table(orth@meta.data$celltype))
rm(orth)
# 数据准备
info = a$Freq
# 命名
names = as.character(a$Var1)
#names = c("0","1","2","3","4","5","6","7")
# 涂色（可选）
#cols = c("#ED1C24","#22B14C","#FFC90E","#3f48CC","#3f90CC","#22B17C","#FFC93E")
# 计算百分比
piepercent = paste(round(100*info/sum(info)), "%")
# 绘图
pie(info, labels=piepercent, main = "ortho cluster ratio", col=colors, family='GB1')
# 添加颜色样本标注
legend("topright", names, cex=0.6, fill=colors)


dev.off()


rm(nonpro1)
nonpro1<-YDL

pdf("专利_figure5F生物过程打分.pdf")


response_to_calcium_ion_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/response to calcium ion.csv",header = F)
response_to_calcium_ion_score<-response_to_calcium_ion_score[,1]
head(response_to_calcium_ion_score)
response_to_calcium_ion_score<-as.data.frame(response_to_calcium_ion_score)
colnames(response_to_calcium_ion_score)<-c("response_to_calcium_ion_score")
response_to_calcium_ion_score<-list(response_to_calcium_ion_score$response_to_calcium_ion_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = response_to_calcium_ion_score,
  name = "response_to_calcium_ion_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[26] <- 'response_to_calcium_ion_score' 
# VlnPlot(object = nonpro1 ,cols =colors,features = c("response_to_calcium_ion_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("response_to_calcium_ion_score"),pt.size = 0)
# VlnPlot(object = nonpro1 ,cols =colors,features = c("response_to_calcium_ion_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1 ,cols =colors,features = c("response_to_calcium_ion_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)


#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","response_to_calcium_ion_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","response_to_calcium_ion_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


cellular_response_to_calcium_ion_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/cellular response to calcium ion.csv",header = F)
cellular_response_to_calcium_ion_score<-toupper(cellular_response_to_calcium_ion_score[,1])
head(cellular_response_to_calcium_ion_score)
cellular_response_to_calcium_ion_score<-as.data.frame(cellular_response_to_calcium_ion_score)
colnames(cellular_response_to_calcium_ion_score)<-c("cellular_response_to_calcium_ion_score")
cellular_response_to_calcium_ion_score<-list(cellular_response_to_calcium_ion_score$cellular_response_to_calcium_ion_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = cellular_response_to_calcium_ion_score,
  name = "cellular_response_to_calcium_ion_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[27] <- 'cellular_response_to_calcium_ion_score' 
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

# 
# 
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","cellular_response_to_calcium_ion_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = cellular_response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","cellular_response_to_calcium_ion_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = cellular_response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




MAPK_cascade_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/MAPK cascade.csv",header = F)
MAPK_cascade_score<-MAPK_cascade_score[,1]
head(MAPK_cascade_score)
MAPK_cascade_score<-as.data.frame(MAPK_cascade_score)
colnames(MAPK_cascade_score)<-c("MAPK_cascade_score")
MAPK_cascade_score<-list(MAPK_cascade_score$MAPK_cascade_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = MAPK_cascade_score,
  name = "MAPK_cascade_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[28] <- 'MAPK_cascade_score' 

# VlnPlot(object = nonpro1 ,cols =colors,features = c("MAPK_cascade_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("MAPK_cascade_score"),pt.size = 0)
# VlnPlot(object = nonpro1 ,cols =colors,features = c("MAPK_cascade_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("MAPK_cascade_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1 ,cols =colors,features = c("MAPK_cascade_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("MAPK_cascade_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("MAPK_cascade_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","MAPK_cascade_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = MAPK_cascade_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","MAPK_cascade_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = MAPK_cascade_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



vesicle_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/vesicle.csv",header = F)
vesicle_score<-vesicle_score[,1]
head(vesicle_score)
vesicle_score<-as.data.frame(vesicle_score)
colnames(vesicle_score)<-c("vesicle_score")
vesicle_score<-list(vesicle_score$vesicle_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = vesicle_score,
  name = "vesicle_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[29] <- 'vesicle_score' 
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vesicle_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("vesicle_score"),pt.size = 0)
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vesicle_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("vesicle_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vesicle_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("vesicle_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("vesicle_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","vesicle_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = vesicle_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","vesicle_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = vesicle_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

actin_binding_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/actin_binding.csv",header = F)
actin_binding_score<-actin_binding_score[,1]
head(actin_binding_score)
actin_binding_score<-as.data.frame(actin_binding_score)
colnames(actin_binding_score)<-c("actin_binding_score")
actin_binding_score<-list(actin_binding_score$actin_binding_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = actin_binding_score,
  name = "actin_binding_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[30] <- 'actin_binding_score' 
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"),pt.size = 0)
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("actin_binding_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

# 
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","actin_binding_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = actin_binding_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","actin_binding_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = actin_binding_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




regulation_of_cell_shape_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/regulation_of_cell_shape.csv",header = F)
regulation_of_cell_shape_score<-regulation_of_cell_shape_score[,1]
head(regulation_of_cell_shape_score)
regulation_of_cell_shape_score<-as.data.frame(regulation_of_cell_shape_score)
colnames(regulation_of_cell_shape_score)<-c("regulation_of_cell_shape_score")
regulation_of_cell_shape_score<-list(regulation_of_cell_shape_score$regulation_of_cell_shape_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = regulation_of_cell_shape_score,
  name = "regulation_of_cell_shape_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[31] <- 'regulation_of_cell_shape_score' 
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"),pt.size = 0)
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("regulation_of_cell_shape_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","regulation_of_cell_shape_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = regulation_of_cell_shape_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","regulation_of_cell_shape_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = regulation_of_cell_shape_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



vacuole_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/vacuole.csv",header = F)
vacuole_score<-vacuole_score[,1]
head(vacuole_score)
vacuole_score<-as.data.frame(vacuole_score)
colnames(vacuole_score)<-c("vacuole_score")
vacuole_score<-list(vacuole_score$vacuole_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = vacuole_score,
  name = "vacuole_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[32] <- 'vacuole_score' 
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vacuole_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("vacuole_score"),pt.size = 0)
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vacuole_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("vacuole_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vacuole_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =colors, features = c("vacuole_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#F8766D","#93AA00","#FF61C3","#00B9E3","#DB72FB"), features = c("vacuole_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","vacuole_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = vacuole_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","vacuole_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = vacuole_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()








YDL<-readRDS("D:/analysis/forpublication/human_bm_ubc/human_bm_all.RDS")

setwd("D:/课题总结/CD44-专利")

library("Seurat")
library("ggplot2")
gg <- TSNEPlot(YDL)
col<-ggplot_build(gg)$data
col<-as.data.frame(col)
table(col$colour)
table(YDL$seurat_clusters)

a<-as.data.frame(table(col$colour))
b<-as.data.frame(table(YDL$seurat_clusters))
c<-merge(a,b,by = "Freq")
c<-c[order(c$Freq,decreasing = T),]

current.cluster.ids <- c(0, 1, 2, 3, 4, 5,6,7,8)
new.cluster.ids <- c(
  "Late-Ortho","Early-Ortho",
  "Pro/Baso","Early-Poly",
  "T-like","Late-Poly",
  "NK/T-like",
  "Mono/N-like","B-like")
names(new.cluster.ids) <- levels(YDL)
YDL <- RenameIdents(YDL, new.cluster.ids)
DimPlot(YDL, reduction = "tsne", label = TRUE, pt.size = 1.5)
DimPlot(YDL, reduction = "umap", label = TRUE, pt.size = 1.5)

YDL$celltype<-Idents(YDL)
YDL$celltype<-factor(YDL$celltype,levels = c("Pro/Baso","Early-Poly","Late-Poly","Early-Ortho","Late-Ortho","T-like",
                                             "NK/T-like",
                                             "Mono/N-like","B-like"))
Idents(YDL)<-YDL$celltype

cols = c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3")


pdf("专利_figure10AB.pdf")
DimPlot(YDL, reduction = "umap", label = TRUE, pt.size = 1.5,cols = c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"))
DimPlot(YDL, reduction = "umap", label = F, pt.size = 1.5,cols = c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"))
DimPlot(YDL, reduction = "umap", label = TRUE, pt.size = 1.5,cols = c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"),split.by = "orig.ident")
DimPlot(YDL, reduction = "umap", label = F, pt.size = 1.5,cols = c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"),split.by = "orig.ident")
DimPlot(YDL, reduction = "umap", label = TRUE, pt.size = 1.5,group.by = "Phase")
# Color palette
#colors <- c("#F8766D","#E68613","#CD9600","#ABA300","#7CAE00","#0CB702","#00BE67","#00C19A","#00BFC4","#00B8E7","#00A9FF","#8494FF","#C77CFF","#ED68ED","#FF61CC","#FF68A1")
colors <- c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3")

DimPlot(YDL,cols =colors,reduction="umap")
a<-as.data.frame(table(YDL@meta.data$celltype))
# 数据准备
info = a$Freq
# 命名
names = as.character(a$Var1)
#names = c("0","1","2","3","4","5","6","7")
# 涂色（可选）
#cols = c("#ED1C24","#22B14C","#FFC90E","#3f48CC","#3f90CC","#22B17C","#FFC93E")
# 计算百分比
piepercent = paste(round(100*info/sum(info),2), "%")
# 绘图
pie(info, labels=piepercent, main = "total cluster ratio", col=colors, family='GB1')
# 添加颜色样本标注
legend("topright", names, cex=0.6, fill=colors)


library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","nCount_RNA"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = log(nCount_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","nFeature_RNA"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = log(nFeature_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

dev.off()


pdf("专利_figure10B.pdf")
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("HBB","GYPA","CD44","PTPRC"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("CD74","TMSB4X","CD79A","S100A8"),cols = c("gray", "red"))#actin
dev.off()


rm(nonpro1)
nonpro1<-YDL

pdf("专利_figure10C生物过程打分.pdf")


response_to_calcium_ion_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/response to calcium ion.csv",header = F)
response_to_calcium_ion_score<-response_to_calcium_ion_score[,1]
head(response_to_calcium_ion_score)
response_to_calcium_ion_score<-as.data.frame(response_to_calcium_ion_score)
colnames(response_to_calcium_ion_score)<-c("response_to_calcium_ion_score")
response_to_calcium_ion_score<-list(response_to_calcium_ion_score$response_to_calcium_ion_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = response_to_calcium_ion_score,
  name = "response_to_calcium_ion_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[12] <- 'response_to_calcium_ion_score' 
# VlnPlot(object = nonpro1 ,cols =colors,features = c("response_to_calcium_ion_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("response_to_calcium_ion_score"),pt.size = 0)
# VlnPlot(object = nonpro1 ,cols =colors,features = c("response_to_calcium_ion_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1 ,cols =colors,features = c("response_to_calcium_ion_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)


#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","response_to_calcium_ion_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","response_to_calcium_ion_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


cellular_response_to_calcium_ion_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/cellular response to calcium ion.csv",header = F)
cellular_response_to_calcium_ion_score<-toupper(cellular_response_to_calcium_ion_score[,1])
head(cellular_response_to_calcium_ion_score)
cellular_response_to_calcium_ion_score<-as.data.frame(cellular_response_to_calcium_ion_score)
colnames(cellular_response_to_calcium_ion_score)<-c("cellular_response_to_calcium_ion_score")
cellular_response_to_calcium_ion_score<-list(cellular_response_to_calcium_ion_score$cellular_response_to_calcium_ion_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = cellular_response_to_calcium_ion_score,
  name = "cellular_response_to_calcium_ion_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[13] <- 'cellular_response_to_calcium_ion_score' 
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("cellular_response_to_calcium_ion_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

# 
# 
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","cellular_response_to_calcium_ion_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = cellular_response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","cellular_response_to_calcium_ion_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = cellular_response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




MAPK_cascade_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/MAPK cascade.csv",header = F)
MAPK_cascade_score<-MAPK_cascade_score[,1]
head(MAPK_cascade_score)
MAPK_cascade_score<-as.data.frame(MAPK_cascade_score)
colnames(MAPK_cascade_score)<-c("MAPK_cascade_score")
MAPK_cascade_score<-list(MAPK_cascade_score$MAPK_cascade_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = MAPK_cascade_score,
  name = "MAPK_cascade_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[14] <- 'MAPK_cascade_score' 

# VlnPlot(object = nonpro1 ,cols =colors,features = c("MAPK_cascade_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("MAPK_cascade_score"),pt.size = 0)
# VlnPlot(object = nonpro1 ,cols =colors,features = c("MAPK_cascade_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("MAPK_cascade_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1 ,cols =colors,features = c("MAPK_cascade_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("MAPK_cascade_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("MAPK_cascade_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","MAPK_cascade_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = MAPK_cascade_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","MAPK_cascade_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = MAPK_cascade_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



vesicle_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/vesicle.csv",header = F)
vesicle_score<-vesicle_score[,1]
head(vesicle_score)
vesicle_score<-as.data.frame(vesicle_score)
colnames(vesicle_score)<-c("vesicle_score")
vesicle_score<-list(vesicle_score$vesicle_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = vesicle_score,
  name = "vesicle_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[15] <- 'vesicle_score' 
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vesicle_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("vesicle_score"),pt.size = 0)
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vesicle_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("vesicle_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vesicle_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("vesicle_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("vesicle_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)



#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","vesicle_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = vesicle_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","vesicle_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = vesicle_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

actin_binding_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/actin_binding.csv",header = F)
actin_binding_score<-actin_binding_score[,1]
head(actin_binding_score)
actin_binding_score<-as.data.frame(actin_binding_score)
colnames(actin_binding_score)<-c("actin_binding_score")
actin_binding_score<-list(actin_binding_score$actin_binding_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = actin_binding_score,
  name = "actin_binding_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[16] <- 'actin_binding_score' 
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"),pt.size = 0)
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("actin_binding_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("actin_binding_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("actin_binding_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

# 
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","actin_binding_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = actin_binding_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","actin_binding_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = actin_binding_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




regulation_of_cell_shape_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/regulation_of_cell_shape.csv",header = F)
regulation_of_cell_shape_score<-regulation_of_cell_shape_score[,1]
head(regulation_of_cell_shape_score)
regulation_of_cell_shape_score<-as.data.frame(regulation_of_cell_shape_score)
colnames(regulation_of_cell_shape_score)<-c("regulation_of_cell_shape_score")
regulation_of_cell_shape_score<-list(regulation_of_cell_shape_score$regulation_of_cell_shape_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = regulation_of_cell_shape_score,
  name = "regulation_of_cell_shape_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[17] <- 'regulation_of_cell_shape_score' 
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"),pt.size = 0)
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("regulation_of_cell_shape_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("regulation_of_cell_shape_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("regulation_of_cell_shape_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","regulation_of_cell_shape_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = regulation_of_cell_shape_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","regulation_of_cell_shape_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = regulation_of_cell_shape_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



vacuole_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/vacuole.csv",header = F)
vacuole_score<-vacuole_score[,1]
head(vacuole_score)
vacuole_score<-as.data.frame(vacuole_score)
colnames(vacuole_score)<-c("vacuole_score")
vacuole_score<-list(vacuole_score$vacuole_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = vacuole_score,
  name = "vacuole_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[18] <- 'vacuole_score' 
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vacuole_score"))
# VlnPlot(object = nonpro1,cols =colors, features = c("vacuole_score"),pt.size = 0)
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vacuole_score"))+ coord_flip()
# VlnPlot(object = nonpro1,cols =colors, features = c("vacuole_score"),pt.size = 0)+ coord_flip()
# VlnPlot(object = nonpro1 ,cols =colors,features = c("vacuole_score"))+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("vacuole_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("vacuole_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","vacuole_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = vacuole_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","vacuole_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = vacuole_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))






dev.off()







rm(nonpro1)
nonpro1<-YDL
colors<-c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3")

pdf("生物过程打分_人类骨髓——重命名.pdf")
#自噬过程的打分

AUTOPHAGE_score<-read.table("D:/analysis/forpublication/immunecell/李程_各个时期的中性粒细胞marker/AUTOPHAGE_MOUSE.txt",header = T)
#AUTOPHAGE_score<-AUTOPHAGE_score[,-1]
head(AUTOPHAGE_score)
AUTOPHAGE_score<-as.data.frame(AUTOPHAGE_score)
colnames(AUTOPHAGE_score)<-c("AUTOPHAGE_score")
AUTOPHAGE_score<-list(AUTOPHAGE_score$AUTOPHAGE_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = AUTOPHAGE_score,
  name = "AUTOPHAGE_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[12] <- 'AUTOPHAGE_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("AUTOPHAGE_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("AUTOPHAGE_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)


#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","AUTOPHAGE_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = AUTOPHAGE_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","AUTOPHAGE_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = AUTOPHAGE_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))





cellular_response_to_calcium_ion_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/cellular response to calcium ion.csv",header = F)
cellular_response_to_calcium_ion_score<-cellular_response_to_calcium_ion_score[,1]
head(cellular_response_to_calcium_ion_score,20)
#cellular_response_to_calcium_ion_score<-toupper(cellular_response_to_calcium_ion_score$V1)
cellular_response_to_calcium_ion_score<-as.data.frame(cellular_response_to_calcium_ion_score)
colnames(cellular_response_to_calcium_ion_score)<-c("cellular_response_to_calcium_ion_score")
cellular_response_to_calcium_ion_score<-list(cellular_response_to_calcium_ion_score$cellular_response_to_calcium_ion_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = cellular_response_to_calcium_ion_score,
  name = "cellular_response_to_calcium_ion_score"
)
head(nonpro1@meta.data,20)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[13] <- 'cellular_response_to_calcium_ion_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","cellular_response_to_calcium_ion_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = cellular_response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","cellular_response_to_calcium_ion_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = cellular_response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



response_to_calcium_ion_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/response to calcium ion.csv",header = F)
response_to_calcium_ion_score<-response_to_calcium_ion_score[,1]
head(response_to_calcium_ion_score)
response_to_calcium_ion_score<-as.data.frame(response_to_calcium_ion_score)
colnames(response_to_calcium_ion_score)<-c("response_to_calcium_ion_score")
response_to_calcium_ion_score<-list(response_to_calcium_ion_score$response_to_calcium_ion_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = response_to_calcium_ion_score,
  name = "response_to_calcium_ion_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[14] <- 'response_to_calcium_ion_score' 

VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)



#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","response_to_calcium_ion_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","response_to_calcium_ion_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



calcium_channel_activity_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/calcium_channel_activity.csv",header = F)
calcium_channel_activity_score<-calcium_channel_activity_score[,1]
head(calcium_channel_activity_score)
calcium_channel_activity_score<-as.data.frame(calcium_channel_activity_score)
colnames(calcium_channel_activity_score)<-c("calcium_channel_activity_score")
calcium_channel_activity_score<-list(calcium_channel_activity_score$calcium_channel_activity_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = calcium_channel_activity_score,
  name = "calcium_channel_activity_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[15] <- 'calcium_channel_activity_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("calcium_channel_activity_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("calcium_channel_activity_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","calcium_channel_activity_score"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = calcium_channel_activity_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","calcium_channel_activity_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = calcium_channel_activity_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



MAPK_cascade_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/MAPK cascade.csv",header = F)
MAPK_cascade_score<-MAPK_cascade_score[,1]
head(MAPK_cascade_score)
MAPK_cascade_score<-as.data.frame(MAPK_cascade_score)
colnames(MAPK_cascade_score)<-c("MAPK_cascade_score")
MAPK_cascade_score<-list(MAPK_cascade_score$MAPK_cascade_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = MAPK_cascade_score,
  name = "MAPK_cascade_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[16] <- 'MAPK_cascade_score' 

VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("MAPK_cascade_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("MAPK_cascade_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","MAPK_cascade_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = MAPK_cascade_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","MAPK_cascade_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = MAPK_cascade_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



activation_of_GTPase_activity_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/activation of GTPase activity.csv",header = F)
activation_of_GTPase_activity_score<-activation_of_GTPase_activity_score[,1]
head(activation_of_GTPase_activity_score)
activation_of_GTPase_activity_score<-as.data.frame(activation_of_GTPase_activity_score)
colnames(activation_of_GTPase_activity_score)<-c("activation_of_GTPase_activity_score")
activation_of_GTPase_activity_score<-list(activation_of_GTPase_activity_score$activation_of_GTPase_activity_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = activation_of_GTPase_activity_score,
  name = "activation_of_GTPase_activity_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[17] <- 'activation_of_GTPase_activity_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)


#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","activation_of_GTPase_activity_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = activation_of_GTPase_activity_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","activation_of_GTPase_activity_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = activation_of_GTPase_activity_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



GTPase_activity_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/GTPase activity.csv",header = F)
GTPase_activity_score<-GTPase_activity_score[,1]
head(GTPase_activity_score)
GTPase_activity_score<-as.data.frame(GTPase_activity_score)
colnames(GTPase_activity_score)<-c("GTPase_activity_score")
GTPase_activity_score<-list(GTPase_activity_score$GTPase_activity_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = GTPase_activity_score,
  name = "GTPase_activity_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[18] <- 'GTPase_activity_score' 

VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","GTPase_activity_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = GTPase_activity_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","GTPase_activity_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = GTPase_activity_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



vesicle_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/vesicle.csv",header = F)
vesicle_score<-vesicle_score[,1]
head(vesicle_score)
vesicle_score<-as.data.frame(vesicle_score)
colnames(vesicle_score)<-c("vesicle_score")
vesicle_score<-list(vesicle_score$vesicle_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = vesicle_score,
  name = "vesicle_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[19] <- 'vesicle_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("vesicle_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("vesicle_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","vesicle_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = vesicle_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","vesicle_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = vesicle_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



autophagy_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/autophagy.csv",header = F)
autophagy_score<-autophagy_score[,1]
head(autophagy_score)
autophagy_score<-as.data.frame(autophagy_score)
colnames(autophagy_score)<-c("autophagy_score")
autophagy_score<-list(autophagy_score$autophagy_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = autophagy_score,
  name = "autophagy_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[20] <- 'autophagy_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("autophagy_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("autophagy_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","autophagy_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = autophagy_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","autophagy_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = autophagy_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



vacuole_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/vacuole.csv",header = F)
vacuole_score<-vacuole_score[,1]
head(vacuole_score)
vacuole_score<-as.data.frame(vacuole_score)
colnames(vacuole_score)<-c("vacuole_score")
vacuole_score<-list(vacuole_score$vacuole_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = vacuole_score,
  name = "vacuole_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[21] <- 'vacuole_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("vacuole_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("vacuole_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","vacuole_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = vacuole_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","vacuole_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = vacuole_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



regulation_of_actin_filament_polymerization_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/regulation of actin filament polymerization.csv",header = F)
regulation_of_actin_filament_polymerization_score<-regulation_of_actin_filament_polymerization_score[,1]
head(regulation_of_actin_filament_polymerization_score)
regulation_of_actin_filament_polymerization_score<-as.data.frame(regulation_of_actin_filament_polymerization_score)
colnames(regulation_of_actin_filament_polymerization_score)<-c("regulation_of_actin_filament_polymerization_score")
regulation_of_actin_filament_polymerization_score<-list(regulation_of_actin_filament_polymerization_score$regulation_of_actin_filament_polymerization_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = regulation_of_actin_filament_polymerization_score,
  name = "regulation_of_actin_filament_polymerization_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[22] <- 'regulation_of_actin_filament_polymerization_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("regulation_of_actin_filament_polymerization_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("regulation_of_actin_filament_polymerization_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","regulation_of_actin_filament_polymerization_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = regulation_of_actin_filament_polymerization_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","regulation_of_actin_filament_polymerization_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = regulation_of_actin_filament_polymerization_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



microtubule_based_process_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/microtubule-based process.csv",header = F)
microtubule_based_process_score<-microtubule_based_process_score[,1]
head(microtubule_based_process_score)
microtubule_based_process_score<-as.data.frame(microtubule_based_process_score)
colnames(microtubule_based_process_score)<-c("microtubule_based_process_score")
microtubule_based_process_score<-list(microtubule_based_process_score$microtubule_based_process_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = microtubule_based_process_score,
  name = "microtubule_based_process_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[23] <- 'microtubule_based_process_score' 

VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("microtubule_based_process_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("microtubule_based_process_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","microtubule_based_process_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = microtubule_based_process_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","microtubule_based_process_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = microtubule_based_process_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))







cytokinesis_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/cytokinesis.csv",header = F)
cytokinesis_score<-cytokinesis_score[,1]
head(cytokinesis_score)
cytokinesis_score<-as.data.frame(cytokinesis_score)
colnames(cytokinesis_score)<-c("cytokinesis_score")
cytokinesis_score<-list(cytokinesis_score$cytokinesis_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = cytokinesis_score,
  name = "cytokinesis_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[24] <- 'cytokinesis_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("cytokinesis_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("cytokinesis_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","cytokinesis_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = cytokinesis_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","cytokinesis_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = cytokinesis_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))





microtubule_based_process_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/microtubule-based process.csv",header = F)
microtubule_based_process_score<-microtubule_based_process_score[,1]
head(microtubule_based_process_score)
microtubule_based_process_score<-as.data.frame(microtubule_based_process_score)
colnames(microtubule_based_process_score)<-c("microtubule_based_process_score")
microtubule_based_process_score<-list(microtubule_based_process_score$microtubule_based_process_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = microtubule_based_process_score,
  name = "microtubule_based_process_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[25] <- 'microtubule_based_process_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("microtubule_based_process_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("microtubule_based_process_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","microtubule_based_process_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = microtubule_based_process_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","microtubule_based_process_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = microtubule_based_process_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))







cytokinesis_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/cytokinesis.csv",header = F)
cytokinesis_score<-cytokinesis_score[,1]
head(cytokinesis_score)
cytokinesis_score<-as.data.frame(cytokinesis_score)
colnames(cytokinesis_score)<-c("cytokinesis_score")
cytokinesis_score<-list(cytokinesis_score$cytokinesis_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = cytokinesis_score,
  name = "cytokinesis_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[26] <- 'cytokinesis_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("cytokinesis_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("cytokinesis_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","cytokinesis_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = cytokinesis_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","cytokinesis_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = cytokinesis_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# 
# dev.off()
# 
# 
# pdf("炎性红细胞特异的GO term打分.pdf")
actin_binding_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/actin_binding.csv",header = F)
actin_binding_score<-actin_binding_score[,1]
head(actin_binding_score)
actin_binding_score<-as.data.frame(actin_binding_score)
colnames(actin_binding_score)<-c("actin_binding_score")
actin_binding_score<-list(actin_binding_score$actin_binding_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = actin_binding_score,
  name = "actin_binding_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[27] <- 'actin_binding_score' 

VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("actin_binding_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("actin_binding_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","actin_binding_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = actin_binding_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","actin_binding_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = actin_binding_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



protein_maturation_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/protein_maturation.csv",header = F)
protein_maturation_score<-protein_maturation_score[,1]
head(protein_maturation_score)
protein_maturation_score<-as.data.frame(protein_maturation_score)
colnames(protein_maturation_score)<-c("protein_maturation_score")
protein_maturation_score<-list(protein_maturation_score$protein_maturation_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = protein_maturation_score,
  name = "protein_maturation_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[28] <- 'protein_maturation_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("protein_maturation_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("protein_maturation_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","protein_maturation_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = protein_maturation_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","protein_maturation_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = protein_maturation_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



positive_regulation_of_protein_transport_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/positive_regulation_of_protein_transport.csv",header = F)
positive_regulation_of_protein_transport_score<-positive_regulation_of_protein_transport_score[,1]
head(positive_regulation_of_protein_transport_score)
positive_regulation_of_protein_transport_score<-as.data.frame(positive_regulation_of_protein_transport_score)
colnames(positive_regulation_of_protein_transport_score)<-c("positive_regulation_of_protein_transport_score")
positive_regulation_of_protein_transport_score<-list(positive_regulation_of_protein_transport_score$positive_regulation_of_protein_transport_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = positive_regulation_of_protein_transport_score,
  name = "positive_regulation_of_protein_transport_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[29] <- 'positive_regulation_of_protein_transport_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("positive_regulation_of_protein_transport_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("positive_regulation_of_protein_transport_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","positive_regulation_of_protein_transport_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = positive_regulation_of_protein_transport_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","positive_regulation_of_protein_transport_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = positive_regulation_of_protein_transport_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

calcium_dependent_protein_binding_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/calcium_dependent_protein_binding.csv",header = F)
calcium_dependent_protein_binding_score<-calcium_dependent_protein_binding_score[,1]
head(calcium_dependent_protein_binding_score)
calcium_dependent_protein_binding_score<-as.data.frame(calcium_dependent_protein_binding_score)
colnames(calcium_dependent_protein_binding_score)<-c("calcium_dependent_protein_binding_score")
calcium_dependent_protein_binding_score<-list(calcium_dependent_protein_binding_score$calcium_dependent_protein_binding_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = calcium_dependent_protein_binding_score,
  name = "calcium_dependent_protein_binding_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[30] <- 'calcium_dependent_protein_binding_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("calcium_dependent_protein_binding_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("calcium_dependent_protein_binding_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","calcium_dependent_protein_binding_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = calcium_dependent_protein_binding_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","calcium_dependent_protein_binding_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = calcium_dependent_protein_binding_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ATP_biosynthetic_process_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/ATP_biosynthetic_process.csv",header = F)
ATP_biosynthetic_process_score<-ATP_biosynthetic_process_score[,1]
head(ATP_biosynthetic_process_score)
ATP_biosynthetic_process_score<-as.data.frame(ATP_biosynthetic_process_score)
colnames(ATP_biosynthetic_process_score)<-c("ATP_biosynthetic_process_score")
ATP_biosynthetic_process_score<-list(ATP_biosynthetic_process_score$ATP_biosynthetic_process_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = ATP_biosynthetic_process_score,
  name = "ATP_biosynthetic_process_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[31] <- 'ATP_biosynthetic_process_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("ATP_biosynthetic_process_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("ATP_biosynthetic_process_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","ATP_biosynthetic_process_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = ATP_biosynthetic_process_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","ATP_biosynthetic_process_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = ATP_biosynthetic_process_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




regulation_of_endothelial_cell_migration_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/regulation_of_endothelial_cell_migration.csv",header = F)
regulation_of_endothelial_cell_migration_score<-regulation_of_endothelial_cell_migration_score[,1]
head(regulation_of_endothelial_cell_migration_score)
regulation_of_endothelial_cell_migration_score<-as.data.frame(regulation_of_endothelial_cell_migration_score)
colnames(regulation_of_endothelial_cell_migration_score)<-c("regulation_of_endothelial_cell_migration_score")
regulation_of_endothelial_cell_migration_score<-list(regulation_of_endothelial_cell_migration_score$regulation_of_endothelial_cell_migration_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = regulation_of_endothelial_cell_migration_score,
  name = "regulation_of_endothelial_cell_migration_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[32] <- 'regulation_of_endothelial_cell_migration_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("regulation_of_endothelial_cell_migration_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("regulation_of_endothelial_cell_migration_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","regulation_of_endothelial_cell_migration_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = regulation_of_endothelial_cell_migration_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","regulation_of_endothelial_cell_migration_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = regulation_of_endothelial_cell_migration_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



endothelial_cell_migration_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/endothelial_cell_migration.csv",header = F)
endothelial_cell_migration_score<-endothelial_cell_migration_score[,1]
head(endothelial_cell_migration_score)
endothelial_cell_migration_score<-as.data.frame(endothelial_cell_migration_score)
colnames(endothelial_cell_migration_score)<-c("endothelial_cell_migration_score")
endothelial_cell_migration_score<-list(endothelial_cell_migration_score$endothelial_cell_migration_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = endothelial_cell_migration_score,
  name = "endothelial_cell_migration_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[33] <- 'endothelial_cell_migration_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("endothelial_cell_migration_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("endothelial_cell_migration_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","endothelial_cell_migration_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = endothelial_cell_migration_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","endothelial_cell_migration_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = endothelial_cell_migration_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



regulation_of_cell_shape_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/regulation_of_cell_shape.csv",header = F)
regulation_of_cell_shape_score<-regulation_of_cell_shape_score[,1]
head(regulation_of_cell_shape_score)
regulation_of_cell_shape_score<-as.data.frame(regulation_of_cell_shape_score)
colnames(regulation_of_cell_shape_score)<-c("regulation_of_cell_shape_score")
regulation_of_cell_shape_score<-list(regulation_of_cell_shape_score$regulation_of_cell_shape_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = regulation_of_cell_shape_score,
  name = "regulation_of_cell_shape_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[34] <- 'regulation_of_cell_shape_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("regulation_of_cell_shape_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("regulation_of_cell_shape_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","regulation_of_cell_shape_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = regulation_of_cell_shape_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","regulation_of_cell_shape_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = regulation_of_cell_shape_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


structural_constituent_of_cytoskeleton_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/structural_constituent_of_cytoskeleton.csv",header = F)
structural_constituent_of_cytoskeleton_score<-structural_constituent_of_cytoskeleton_score[,1]
head(structural_constituent_of_cytoskeleton_score)
structural_constituent_of_cytoskeleton_score<-as.data.frame(structural_constituent_of_cytoskeleton_score)
colnames(structural_constituent_of_cytoskeleton_score)<-c("structural_constituent_of_cytoskeleton_score")
structural_constituent_of_cytoskeleton_score<-list(structural_constituent_of_cytoskeleton_score$structural_constituent_of_cytoskeleton_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = structural_constituent_of_cytoskeleton_score,
  name = "structural_constituent_of_cytoskeleton_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[35] <- 'structural_constituent_of_cytoskeleton_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("structural_constituent_of_cytoskeleton_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("structural_constituent_of_cytoskeleton_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","structural_constituent_of_cytoskeleton_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = structural_constituent_of_cytoskeleton_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","structural_constituent_of_cytoskeleton_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = structural_constituent_of_cytoskeleton_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



actin_filament_bundle_assembly_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/actin_filament_bundle_assembly.csv",header = F)
actin_filament_bundle_assembly_score<-actin_filament_bundle_assembly_score[,1]
head(actin_filament_bundle_assembly_score)
actin_filament_bundle_assembly_score<-as.data.frame(actin_filament_bundle_assembly_score)
colnames(actin_filament_bundle_assembly_score)<-c("actin_filament_bundle_assembly_score")
actin_filament_bundle_assembly_score<-list(actin_filament_bundle_assembly_score$actin_filament_bundle_assembly_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = actin_filament_bundle_assembly_score,
  name = "actin_filament_bundle_assembly_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[36] <- 'actin_filament_bundle_assembly_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("actin_filament_bundle_assembly_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("actin_filament_bundle_assembly_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","actin_filament_bundle_assembly_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = actin_filament_bundle_assembly_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","actin_filament_bundle_assembly_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = actin_filament_bundle_assembly_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



regulation_of_vesicle_fusion_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/regulation_of_vesicle_fusion.csv",header = F)
regulation_of_vesicle_fusion_score<-regulation_of_vesicle_fusion_score[,1]
head(regulation_of_vesicle_fusion_score)
regulation_of_vesicle_fusion_score<-as.data.frame(regulation_of_vesicle_fusion_score)
colnames(regulation_of_vesicle_fusion_score)<-c("regulation_of_vesicle_fusion_score")
regulation_of_vesicle_fusion_score<-list(regulation_of_vesicle_fusion_score$regulation_of_vesicle_fusion_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = regulation_of_vesicle_fusion_score,
  name = "regulation_of_vesicle_fusion_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[37] <- 'regulation_of_vesicle_fusion_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("regulation_of_vesicle_fusion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("regulation_of_vesicle_fusion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","regulation_of_vesicle_fusion_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = regulation_of_vesicle_fusion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","regulation_of_vesicle_fusion_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = regulation_of_vesicle_fusion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



lytic_vacuole_organization_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/lytic_vacuole_organization.csv",header = F)
lytic_vacuole_organization_score<-lytic_vacuole_organization_score[,1]
head(lytic_vacuole_organization_score)
lytic_vacuole_organization_score<-as.data.frame(lytic_vacuole_organization_score)
colnames(lytic_vacuole_organization_score)<-c("lytic_vacuole_organization_score")
lytic_vacuole_organization_score<-list(lytic_vacuole_organization_score$lytic_vacuole_organization_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = lytic_vacuole_organization_score,
  name = "lytic_vacuole_organization_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[38] <- 'lytic_vacuole_organization_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("lytic_vacuole_organization_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("lytic_vacuole_organization_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","lytic_vacuole_organization_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = lytic_vacuole_organization_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","lytic_vacuole_organization_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = lytic_vacuole_organization_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



regulation_of_histone_H3_K9_methylation_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/regulation_of_histone_H3_K9_methylation.csv",header = F)
regulation_of_histone_H3_K9_methylation_score<-regulation_of_histone_H3_K9_methylation_score[,1]
head(regulation_of_histone_H3_K9_methylation_score)
regulation_of_histone_H3_K9_methylation_score<-as.data.frame(regulation_of_histone_H3_K9_methylation_score)
colnames(regulation_of_histone_H3_K9_methylation_score)<-c("regulation_of_histone_H3_K9_methylation_score")
regulation_of_histone_H3_K9_methylation_score<-list(regulation_of_histone_H3_K9_methylation_score$regulation_of_histone_H3_K9_methylation_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = regulation_of_histone_H3_K9_methylation_score,
  name = "regulation_of_histone_H3_K9_methylation_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[39] <- 'regulation_of_histone_H3_K9_methylation_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("regulation_of_histone_H3_K9_methylation_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("regulation_of_histone_H3_K9_methylation_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","regulation_of_histone_H3_K9_methylation_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = regulation_of_histone_H3_K9_methylation_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","regulation_of_histone_H3_K9_methylation_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = regulation_of_histone_H3_K9_methylation_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




Toll_like_receptor_binding_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/Toll_like_receptor_binding.csv",header = F)
Toll_like_receptor_binding_score<-Toll_like_receptor_binding_score[,1]
head(Toll_like_receptor_binding_score)
Toll_like_receptor_binding_score<-as.data.frame(Toll_like_receptor_binding_score)
colnames(Toll_like_receptor_binding_score)<-c("Toll_like_receptor_binding_score")
Toll_like_receptor_binding_score<-list(Toll_like_receptor_binding_score$Toll_like_receptor_binding_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = Toll_like_receptor_binding_score,
  name = "Toll_like_receptor_binding_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[40] <- 'Toll_like_receptor_binding_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("Toll_like_receptor_binding_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("Toll_like_receptor_binding_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","Toll_like_receptor_binding_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = Toll_like_receptor_binding_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","Toll_like_receptor_binding_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = Toll_like_receptor_binding_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


cell_fate_commitment_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/cell_fate_commitment.csv",header = F)
cell_fate_commitment_score<-cell_fate_commitment_score[,1]
head(cell_fate_commitment_score)
cell_fate_commitment_score<-as.data.frame(cell_fate_commitment_score)
colnames(cell_fate_commitment_score)<-c("cell_fate_commitment_score")
cell_fate_commitment_score<-list(cell_fate_commitment_score$cell_fate_commitment_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = cell_fate_commitment_score,
  name = "cell_fate_commitment_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[41] <- 'cell_fate_commitment_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("cell_fate_commitment_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("cell_fate_commitment_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","cell_fate_commitment_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = cell_fate_commitment_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","cell_fate_commitment_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = cell_fate_commitment_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



positive_regulation_of_MAP_kinase_activity_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/positive_regulation_of_MAP_kinase_activity.csv",header = F)
positive_regulation_of_MAP_kinase_activity_score<-positive_regulation_of_MAP_kinase_activity_score[,1]
head(positive_regulation_of_MAP_kinase_activity_score)
positive_regulation_of_MAP_kinase_activity_score<-as.data.frame(positive_regulation_of_MAP_kinase_activity_score)
colnames(positive_regulation_of_MAP_kinase_activity_score)<-c("positive_regulation_of_MAP_kinase_activity_score")
positive_regulation_of_MAP_kinase_activity_score<-list(positive_regulation_of_MAP_kinase_activity_score$positive_regulation_of_MAP_kinase_activity_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = positive_regulation_of_MAP_kinase_activity_score,
  name = "positive_regulation_of_MAP_kinase_activity_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[42] <- 'positive_regulation_of_MAP_kinase_activity_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("positive_regulation_of_MAP_kinase_activity_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("positive_regulation_of_MAP_kinase_activity_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","positive_regulation_of_MAP_kinase_activity_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = positive_regulation_of_MAP_kinase_activity_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","positive_regulation_of_MAP_kinase_activity_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = positive_regulation_of_MAP_kinase_activity_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



chromatin_remodeling_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/chromatin_remodeling.csv",header = F)
chromatin_remodeling_score<-chromatin_remodeling_score[,1]
head(chromatin_remodeling_score)
chromatin_remodeling_score<-as.data.frame(chromatin_remodeling_score)
colnames(chromatin_remodeling_score)<-c("chromatin_remodeling_score")
chromatin_remodeling_score<-list(chromatin_remodeling_score$chromatin_remodeling_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = chromatin_remodeling_score,
  name = "chromatin_remodeling_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[43] <- 'chromatin_remodeling_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("chromatin_remodeling_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("chromatin_remodeling_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","chromatin_remodeling_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = chromatin_remodeling_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","chromatin_remodeling_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = chromatin_remodeling_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


activation_of_innate_immune_response_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/activation_of_innate_immune_response.csv",header = F)
activation_of_innate_immune_response_score<-activation_of_innate_immune_response_score[,1]
head(activation_of_innate_immune_response_score)
activation_of_innate_immune_response_score<-as.data.frame(activation_of_innate_immune_response_score)
colnames(activation_of_innate_immune_response_score)<-c("activation_of_innate_immune_response_score")
activation_of_innate_immune_response_score<-list(activation_of_innate_immune_response_score$activation_of_innate_immune_response_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = activation_of_innate_immune_response_score,
  name = "activation_of_innate_immune_response_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[44] <- 'activation_of_innate_immune_response_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("activation_of_innate_immune_response_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("activation_of_innate_immune_response_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#nonpro1@meta.data$activation_of_innate_immune_response_score<-scale(nonpro1@meta.data$activation_of_innate_immune_response_score)
#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","activation_of_innate_immune_response_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = activation_of_innate_immune_response_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","activation_of_innate_immune_response_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = activation_of_innate_immune_response_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



cellular_response_to_calcium_ion_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/cellular_response_to_calcium_ion.csv",header = F)
cellular_response_to_calcium_ion_score<-cellular_response_to_calcium_ion_score[,1]
head(cellular_response_to_calcium_ion_score)
cellular_response_to_calcium_ion_score<-as.data.frame(cellular_response_to_calcium_ion_score)
colnames(cellular_response_to_calcium_ion_score)<-c("cellular_response_to_calcium_ion_score")
cellular_response_to_calcium_ion_score<-list(cellular_response_to_calcium_ion_score$cellular_response_to_calcium_ion_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = cellular_response_to_calcium_ion_score,
  name = "cellular_response_to_calcium_ion_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[45] <- 'cellular_response_to_calcium_ion_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("cellular_response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#nonpro1@meta.data$cellular_response_to_calcium_ion_score<-scale(nonpro1@meta.data$cellular_response_to_calcium_ion_score)
#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","cellular_response_to_calcium_ion_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = cellular_response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","cellular_response_to_calcium_ion_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = cellular_response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



response_to_calcium_ion_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/response_to_calcium_ion.csv",header = F)
response_to_calcium_ion_score<-response_to_calcium_ion_score[,1]
head(response_to_calcium_ion_score)
response_to_calcium_ion_score<-as.data.frame(response_to_calcium_ion_score)
colnames(response_to_calcium_ion_score)<-c("response_to_calcium_ion_score")
response_to_calcium_ion_score<-list(response_to_calcium_ion_score$response_to_calcium_ion_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = response_to_calcium_ion_score,
  name = "response_to_calcium_ion_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[46] <- 'response_to_calcium_ion_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("response_to_calcium_ion_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)


summary(nonpro1@meta.data$response_to_calcium_ion_score_Feature)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.259960 -0.004752  0.097868  0.084557  0.152367  2.377477 
#boxplot(nonpro1@meta.data$response_to_calcium_ion_score_Feature)

#nonpro1@meta.data$response_to_calcium_ion_score<-scale(nonpro1@meta.data$response_to_calcium_ion_score)
#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","response_to_calcium_ion_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","response_to_calcium_ion_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = response_to_calcium_ion_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



regulation_of_MAP_kinase_activity_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/regulation_of_MAP_kinase_activity.csv",header = F)
regulation_of_MAP_kinase_activity_score<-regulation_of_MAP_kinase_activity_score[,1]
head(regulation_of_MAP_kinase_activity_score)
regulation_of_MAP_kinase_activity_score<-as.data.frame(regulation_of_MAP_kinase_activity_score)
colnames(regulation_of_MAP_kinase_activity_score)<-c("regulation_of_MAP_kinase_activity_score")
regulation_of_MAP_kinase_activity_score<-list(regulation_of_MAP_kinase_activity_score$regulation_of_MAP_kinase_activity_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = regulation_of_MAP_kinase_activity_score,
  name = "regulation_of_MAP_kinase_activity_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[47] <- 'regulation_of_MAP_kinase_activity_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("regulation_of_MAP_kinase_activity_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("regulation_of_MAP_kinase_activity_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#nonpro1@meta.data$regulation_of_MAP_kinase_activity_score<-scale(nonpro1@meta.data$regulation_of_MAP_kinase_activity_score)
#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","regulation_of_MAP_kinase_activity_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = regulation_of_MAP_kinase_activity_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","regulation_of_MAP_kinase_activity_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = regulation_of_MAP_kinase_activity_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



actin_filament_binding_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/actin_filament_binding.csv",header = F)
actin_filament_binding_score<-actin_filament_binding_score[,1]
head(actin_filament_binding_score)
actin_filament_binding_score<-as.data.frame(actin_filament_binding_score)
colnames(actin_filament_binding_score)<-c("actin_filament_binding_score")
actin_filament_binding_score<-list(actin_filament_binding_score$actin_filament_binding_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = actin_filament_binding_score,
  name = "actin_filament_binding_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[48] <- 'actin_filament_binding_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("actin_filament_binding_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("actin_filament_binding_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#nonpro1@meta.data$actin_filament_binding_score<-scale(nonpro1@meta.data$actin_filament_binding_score)
#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","actin_filament_binding_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = actin_filament_binding_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","actin_filament_binding_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = actin_filament_binding_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




actin_monomer_binding_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/actin_monomer_binding.csv",header = F)
actin_monomer_binding_score<-actin_monomer_binding_score[,1]
head(actin_monomer_binding_score)
actin_monomer_binding_score<-as.data.frame(actin_monomer_binding_score)
colnames(actin_monomer_binding_score)<-c("actin_monomer_binding_score")
actin_monomer_binding_score<-list(actin_monomer_binding_score$actin_monomer_binding_score)
nonpro1 <- AddModuleScore(
  object = nonpro1,
  features = actin_monomer_binding_score,
  name = "actin_monomer_binding_score"
)
head(nonpro1@meta.data)
dim(nonpro1@meta.data)
colnames(nonpro1@meta.data)[49] <- 'actin_monomer_binding_score' 
VlnPlot(object = nonpro1,cols =c("#93AA00","#00BA38","#00B9E3","#D39200","#F8766D","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("actin_monomer_binding_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
VlnPlot(object = nonpro1,cols =c("#989898","#989898","#989898","#989898","#989898","#00C19F","#619CFF","#DB72FB","#FF61C3"), features = c("actin_monomer_binding_score"),pt.size = 0)+ coord_flip()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

#nonpro1@meta.data$actin_monomer_binding_score<-scale(nonpro1@meta.data$actin_monomer_binding_score)
#绘制分数的分布图
# library(ggplot2)
# mydata<- FetchData(nonpro1,vars = c("tSNE_1","tSNE_2","actin_monomer_binding_score"))
# a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = actin_monomer_binding_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))
# 
# a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro1,vars = c("UMAP_1","UMAP_2","actin_monomer_binding_score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = actin_monomer_binding_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



dev.off()






# grammar
library(tidyverse)
library(magrittr)
library(glue)
library(data.table)

# analysis
library(DESeq2)
library(org.Mm.eg.db)
library(clusterProfiler)

# graphics
library(ggplot2)
library(ggsci)
library(latex2exp)
library(patchwork)

# * 2. Load data ----------------------------------------------------------

GOfile <- list.files(".", "GO.csv")

# * 3. Plot ---------------------------------------------------------------

GOdata <- map(GOfile, ~ fread(.x)) %>% set_names(str_remove(GOfile, ".GO.*"))
plotData <- GOdata %>% imap(~ {.x[pvalue < 0.05 & Count >= 3][1:8][, type := .y]})
plotData %<>% map(~ {.x[, p := -log10(pvalue)]})
head(plotData$GO.csv,20)
GOplot <- imap(plotData, ~ {
  ggplot(.x, aes(x = p, y = fct_reorder(Description, p,))) +
    geom_col(aes(alpha = p), fill = "black", width = .8, show.legend = F) +
    scale_alpha_continuous(range = c(.5, 1)) +
    scale_x_continuous(expand = expansion(c(0, 0.05))) +
    labs(x = TeX("$-log_{10}(\\textit{P}\\,value)$"), y = "", title = .y) +
    theme(
      aspect.ratio = 0.75,
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line())
})

purrr::reduce(GOplot, `+`) + plot_layout(ncol = 1)

ggsave("GO.pdf", width = 10, height = 4 * length(GOfile))



YDL<-readRDS("D:/analysis/forpublication/human_abc/GSE149938/humanbloodcell.rds")
ery<-subset(YDL, cells= rownames(YDL@meta.data[YDL@meta.data$orig.ident=="ery",]))
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(ery,reduction = "tsne",label = TRUE,pt.size = 1.5)

##鉴定表达高变基因(2000个）,用于下游分析,如PCA；
ery <- FindVariableFeatures(object = ery, selection.method = "vst", nfeatures = 2000)
ery <- ScaleData(object = ery, features = rownames(ery))
#线性降维（PCA）,默认用高变基因集,但也可通过features参数自己指定；
ery=RunPCA(object= ery,npcs = 20,pc.genes=VariableFeatures(object = ery))     #PCA分析
ElbowPlot(ery)#选择top20个PC
pcSelect=15
ery <- FindNeighbors(object = ery, dims = 1:pcSelect)     
##接着优化模型,resolution参数决定下游聚类分析得到的分群数,对于3K左右的细胞,设为0.4-1.2 能得到较好的结果(官方说明)；如果数据量增大,该参数也应该适当增大。
ery<- FindClusters(object = ery, resolution = 0.3)                  #对细胞分组,优化标准模块化
##使用Idents（）函数可查看不同细胞的分群；
head(Idents(ery), 5)

##这里采用基于TSNE的聚类方法。
ery<- RunTSNE(object = ery, dims = 1:pcSelect,check_duplicates = FALSE)                      #TSNE聚类
pdf(file="TSNE_rty.pdf",width=6.5,height=6)
TSNEPlot(object = ery, pt.size = 2, label = TRUE)    #TSNE可视化

#另一个可视化的方法
DimPlot(ery,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(ery,reduction = "tsne",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
dev.off()
write.table(ery$seurat_clusters,file="tsneCluster_ery.txt",quote=F,sep="\t",col.names=F)


##这里采用基于umap的聚类方法
#这里采用基于图论的聚类方法
ery <- RunUMAP(object = ery, dims = 1:pcSelect)                      #umap聚类
pdf(file="umap__ery.pdf",width=6.5,height=6)
UMAPPlot(object = ery, pt.size = 1.5, label = TRUE)    #umap可视化
#另一个可视化的方法
DimPlot(object=ery,label = TRUE,reduction="umap")
##用DimPlot()函数绘制散点图,reduction = "tsne",指定绘制类型；如果不指定,默认先从搜索 umap,然后 tsne, 再然后 pca；也可以直接使用这3个函数PCAPlot()、TSNEPlot()、UMAPPlot()； cols,pt.size分别调整分组颜色和点的大小；
dev.off()

##细胞周期归类

ery<- CellCycleScoring(object = ery, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

head(x = ery@meta.data)
pdf(file="CellCycle_ery.pdf",width=6.5,height=6)
DimPlot(ery,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5)
DimPlot(ery,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5,split.by = "orig.ident")
DimPlot(ery,reduction = "umap",label = TRUE,group.by="Phase",pt.size = 1.5)
dev.off()


ery.markers <- FindAllMarkers(ery, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'roc')
# install.packages("magrittr") # package installations are only needed the first time you use it
# install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run evYDL time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
ery.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
##存储marker

write.csv(ery.markers,file="allmarker_celltpye_ery.csv")

#绘制分cluster的热图
top10 <- ery.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#绘制marker在各个cluster的热图
pdf(file="tsneHeatmap_ery.pdf",width=12,height=36)
DoHeatmap(object = ery, features = top10$gene) + NoLegend()
DoHeatmap(subset(ery, downsample = 100), features = top10$gene, size = 3)+ NoLegend()
dev.off()




ery<-readRDS("D:/analysis/forpublication/human_abc/GSE149938/ery.rds")
library("Seurat")
library("ggplot2")
gg <- TSNEPlot(ery)
col<-ggplot_build(gg)$data
col<-as.data.frame(col)
table(col$colour)
table(ery$seurat_clusters)

a<-as.data.frame(table(col$colour))
b<-as.data.frame(table(ery$seurat_clusters))
c<-merge(a,b,by = "Freq")
c<-c[order(c$Freq,decreasing = T),]

pdf("专利_figure9D.pdf")
DimPlot(ery,reduction = "tsne",label = TRUE,pt.size = 1.5)
colors <- c("#F8766D","#C49A00","#53B400","#00C094","#00B6EB","#A58AFF","#FB61D7")

#DimPlot(ery,cols =colors,reduction="tsne")
a<-as.data.frame(table(ery@meta.data$seurat_clusters))
# 数据准备
info = a$Freq
# 命名
names = as.character(a$Var1)
#names = c("0","1","2","3","4","5","6","7")
# 涂色（可选）
#cols = c("#ED1C24","#22B14C","#FFC90E","#3f48CC","#3f90CC","#22B17C","#FFC93E")
# 计算百分比
piepercent = paste(round(100*info/sum(info),2), "%")
# 绘图
pie(info, labels=piepercent, main = "total cluster ratio", col=colors, family='GB1')
# 添加颜色样本标注
legend("topright", names, cex=0.7, fill=colors)
dev.off()


pdf("专利_figure9D热图.pdf",height = 10,width = 14)
FeaturePlot(object = ery, reduction = "tsne",pt.size = 1.5,features = c("HBB","GATA1","S100A8","CD44","PTPRC","CD79A"),cols = c("gray", "red"),ncol = 3)#actin
dev.off()

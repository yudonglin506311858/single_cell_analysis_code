set.seed(123)
library("Seurat")
nonpro<-readRDS("D:/analysis/nonpro/nonpro_singlet.rds")
DimPlot(nonpro,reduction = "tsne",label = TRUE,pt.size = 1.5)
#平均表达值



#pro的打分

pro_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/pro_symbol.csv",header = T)
pro_score<-pro_score[,-1]
head(pro_score)
pro_score<-as.data.frame(pro_score)
colnames(pro_score)<-c("pro_score")
pro_score<-list(pro_score$pro_score)
nonpro <- AddModuleScore(
  object = nonpro,
  features = pro_score,
  name = "pro_score"
)
head(nonpro@meta.data)
colnames(nonpro@meta.data)[9] <- 'pro_score' 
VlnPlot(object = nonpro, features = c("pro_score"))
summary(nonpro@meta.data$pro_score_Feature)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.259960 -0.004752  0.097868  0.084557  0.152367  2.377477 
boxplot(nonpro@meta.data$pro_score_Feature)


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro,vars = c("tSNE_1","tSNE_2","pro_score"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = pro_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




#baso的打分

baso_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/baso_symbol_1.csv",header = T)
baso_score<-baso_score[,-1]
head(baso_score)
baso_score<-as.data.frame(baso_score)
colnames(baso_score)<-c("baso_score")
baso_score<-list(baso_score$baso_score)
nonpro <- AddModuleScore(
  object = nonpro,
  features = baso_score,
  name = "baso_score"
)
head(nonpro@meta.data)
colnames(nonpro@meta.data)[10] <- 'baso_score' 
VlnPlot(object = nonpro, features = c("baso_score"))
summary(nonpro@meta.data$baso_score)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.259960 -0.004752  0.097868  0.084557  0.152367  2.377477 
boxplot(nonpro@meta.data$baso_score)

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro,vars = c("tSNE_1","tSNE_2","baso_score"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = baso_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))





#ortho的打分

ortho_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/ortho_symbol.csv",header = T)
ortho_score<-ortho_score[,-1]
head(ortho_score)
ortho_score<-as.data.frame(ortho_score)
colnames(ortho_score)<-c("ortho_score")
ortho_score<-list(ortho_score$ortho_score)
nonpro <- AddModuleScore(
  object = nonpro,
  features = ortho_score,
  name = "ortho_score"
)
head(nonpro@meta.data)
colnames(nonpro@meta.data)[11] <- 'ortho_score' 
VlnPlot(object = nonpro, features = c("ortho_score"))
summary(nonpro@meta.data$ortho_score_Feature)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.259960 -0.004752  0.097868  0.084557  0.152367  2.377477 
boxplot(nonpro@meta.data$ortho_score_Feature)


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro,vars = c("tSNE_1","tSNE_2","ortho_score"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = ortho_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




#poly的打分

poly_score<-read.csv("D:/analysis/anxiuli_bulkdata/symbolgene/poly_symbol_1.csv",header = T)
poly_score<-poly_score[,-1]
head(poly_score)
poly_score<-as.data.frame(poly_score)
colnames(poly_score)<-c("poly_score")
poly_score<-list(poly_score$poly_score)
nonpro <- AddModuleScore(
  object = nonpro,
  features = poly_score,
  name = "poly_score"
)
head(nonpro@meta.data)
colnames(nonpro@meta.data)[12] <- 'poly_score' 
VlnPlot(object = nonpro, features = c("poly_score"))
summary(nonpro@meta.data$poly_score_Feature)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.259960 -0.004752  0.097868  0.084557  0.152367  2.377477 
boxplot(nonpro@meta.data$poly_score_Feature)

#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro,vars = c("tSNE_1","tSNE_2","poly_score"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = poly_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

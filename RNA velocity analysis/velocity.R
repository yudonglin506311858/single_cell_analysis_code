conda create -n velocyto
conda activate velocyto
conda install R
conda install -c conda-forge openmp
conda install -c anaconda boost

#R environment
install.packages("devtools")

library(devtools)
install_github("velocyto-team/velocyto.R")
library(velocyto.R)
#导入包
library(velocyto.R)
library(pagoda2)
YDL<- readRDS("YDL_NONPRO.rds")
ldat <- read.loom.matrices(file = "/data/yudonglin/singlecell/ABFC20190816-04-分析结果/NoPro/velocyto/NoPro.loom")
emat <- ldat$spliced
nmat <- ldat$unspliced
seurat.object<-YDL
emb <- seurat.object@reductions$umap@cell.embeddings
# Estimate the cell-cell distances 
cell.dist <- as.dist(1-armaCor(t(seurat.object@reductions$umap@cell.embeddings)))
head(colnames(emat))
head(colnames(nmat))
colnames(emat) <- paste(substring(colnames(emat),7,22),"-1",sep="")
colnames(nmat) <- paste(substring(colnames(nmat),7,22),"-1",sep="")
head(colnames(emat))
head(colnames(nmat))
# I'm not sure what this parameter does to be honest. 0.02 default
# perform gamma fit on a top/bottom quantiles of expression magnitudes
fit.quantile <- 0.02
# Main velocity estimation

rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=2,
                                            kCells=10,
                                            cell.dist=cell.dist,
                                            fit.quantile=fit.quantile,
                                            n.cores=24)
# This section gets the colors out of the seurat tSNE object so that my seurat and velocyto plots use the same color scheme.
library("Seurat")
library("ggplot2")

gg <- UMAPPlot(seurat.object)
ggplot_build(gg)$data
colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
names(colors) <- rownames(emb)
p1 <- show.velocity.on.embedding.cor(emb,rvel.cd,n=30,scale='sqrt',
                                      cell.colors=ac(colors,alpha=0.5),
                                      cex=0.8,arrow.scale=2,show.grid.flow=T,
                                      min.grid.cell.mass=1.0,grid.n=50,arrow.lwd=1,
                                      do.par=F,cell.border.alpha = 0.1,
                                      n.cores=24,main="RNA Velocity")
DimPlot(YDL, reduction = 'umap', label=T)
dev.off()
#导入包
library(velocyto.R)
library(pagoda2)
YDL<- readRDS("YDL_NONPRO.rds")
ldat <- read.loom.matrices(file = "/data/yudonglin/singlecell/ABFC20190816-04-分析结果/NoPro/velocyto/NoPro.loom")
emat <- ldat$spliced
nmat <- ldat$unspliced
seurat.object<-YDL
emb <- seurat.object@reductions$tsne@cell.embeddings
# Estimate the cell-cell distances 
cell.dist <- as.dist(1-armaCor(t(seurat.object@reductions$tsne@cell.embeddings)))
head(colnames(emat))
head(colnames(nmat))
colnames(emat) <- paste(substring(colnames(emat),7,22),"-1",sep="")
colnames(nmat) <- paste(substring(colnames(nmat),7,22),"-1",sep="")
head(colnames(emat))
head(colnames(nmat))
# I'm not sure what this parameter does to be honest. 0.02 default
# perform gamma fit on a top/bottom quantiles of expression magnitudes
fit.quantile <- 0.02
# Main velocity estimation

rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=2,
                                            kCells=10,
                                            cell.dist=cell.dist,
                                            fit.quantile=fit.quantile,
                                            n.cores=48)
# This section gets the colors out of the seurat tSNE object so that my seurat and velocyto plots use the same color scheme.
library("Seurat")
library("ggplot2")

gg <- TSNEPlot(seurat.object)
ggplot_build(gg)$data
colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
names(colors) <- rownames(emb)
p1 <- show.velocity.on.embedding.cor(emb,rvel.cd,n=30,scale='sqrt',
                                      cell.colors=ac(colors,alpha=0.5),
                                      cex=0.8,arrow.scale=2,show.grid.flow=T,
                                      min.grid.cell.mass=1.0,grid.n=50,arrow.lwd=1,
                                      do.par=F,cell.border.alpha = 0.1,
                                      n.cores=48,main="RNA Velocity")
DimPlot(YDL, reduction = 'tsne', label=T)
dev.off()

#REFERNCE:https://cloud.tencent.com/developer/article/1620345

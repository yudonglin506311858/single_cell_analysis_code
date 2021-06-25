# conda activate velocity
#export OPENBLAS_NUM_THREADS=1
#导入包
library(velocyto.R)
library(pagoda2)
YDL<- readRDS("ortho_nonpro1_SINGLET.RDS")

head(colnames(YDL))
x1 <-velocyto.R::read.loom.matrices(file = "/data/yudonglin/singlecell/ABFC20190816-04-分析结果/NoPro/velocyto/NoPro.loom", engine = "hdf5r")
splice<-x1$spliced
unsplice<-x1$unspliced
head(colnames(splice))
head(colnames(unsplice))
colnames(splice) <- paste(substring(colnames(splice),7,22),"-1_1",sep="")
colnames(unsplice) <- paste(substring(colnames(unsplice),7,22),"-1_1",sep="")
head(colnames(splice))
head(colnames(unsplice))
x1$spliced<-splice
x1$unspliced<-unsplice

x2 <-velocyto.R::read.loom.matrices(file = "/data/yudonglin/software/cellranger-4.0.0/ortho/velocyto/ortho.loom", engine = "hdf5r")
splice<-x2$spliced
unsplice<-x2$unspliced
head(colnames(splice))
head(colnames(unsplice))
colnames(splice) <- paste(substring(colnames(splice),7,22),"-1_2",sep="")
colnames(unsplice) <- paste(substring(colnames(unsplice),7,22),"-1_2",sep="")
head(colnames(splice))
head(colnames(unsplice))
x2$spliced<-splice
x2$unspliced<-unsplice

spliced <- cbind(x1[["spliced"]], x2[["spliced"]])
unspliced <- cbind(x1[["unspliced"]], x2[["unspliced"]])

emat <- spliced
nmat <- unspliced
seurat.object<-YDL
emb <- seurat.object@reductions$tsne@cell.embeddings
# Estimate the cell-cell distances 
cell.dist <- as.dist(1-armaCor(t(seurat.object@reductions$tsne@cell.embeddings)))


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
p2<-show.velocity.on.embedding.cor(emb, rvel.cd, n = 200, scale = 'sqrt',
                                   cell.colors = ac(colors, alpha = 0.5),
                                   cex = 0.8, arrow.scale = 20, show.grid.flow = T,
                                   min.grid.cell.mass = 0.5, grid.n = 40,
                                   arrow.lwd = 1,do.par = T, cell.border.alpha = 0.1,n.cores=48,main="RNA Velocity")

dev.off()
savehistory("ortho_velocity.txt")


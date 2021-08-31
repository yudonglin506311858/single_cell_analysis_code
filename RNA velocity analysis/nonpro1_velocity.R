#export OPENBLAS_NUM_THREADS=1
## conda activate velocity
#导入包
set.seed(123)
library(Seurat)
library(velocyto.R)
library(pagoda2)
YDL<- readRDS("nonpro1_SINGLET.RDS")
ldat <- read.loom.matrices("/data/yudonglin/software/cellranger-4.0.0/erythropoiesis_nonpro/velocyto/nonpro1.loom")
emat <- ldat$spliced
nmat <- ldat$unspliced
seurat.object<-YDL
emb <- seurat.object@reductions$tsne@cell.embeddings
# Estimate the cell-cell distances 
cell.dist <- as.dist(1-armaCor(t(seurat.object@reductions$tsne@cell.embeddings)))
head(colnames(emat))
head(colnames(nmat))
colnames(emat) <- paste(substring(colnames(emat),9,24),"-1",sep="")
colnames(nmat) <- paste(substring(colnames(nmat),9,24),"-1",sep="")
head(colnames(emat))
head(colnames(nmat))
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
pdf("nonpro_rna_velocity.pdf")
p1 <- show.velocity.on.embedding.cor(emb,rvel.cd,n=30,scale='sqrt',
                                      cell.colors=ac(colors,alpha=0.5),
                                      cex=0.8,arrow.scale=60,show.grid.flow=T,
                                      min.grid.cell.mass=1.0,grid.n=50,arrow.lwd=1,
                                      do.par=F,cell.border.alpha = 0.1,
                                      n.cores=48,main="RNA Velocity")
DimPlot(YDL, reduction = 'tsne', label=T)
p2<-show.velocity.on.embedding.cor(emb, rvel.cd, n = 200, scale = 'sqrt',
                               cell.colors = ac(colors, alpha = 0.5),
                               cex = 0.8, arrow.scale = 10, show.grid.flow = T,
                               min.grid.cell.mass = 0.5, grid.n = 40,
                               arrow.lwd = 1,do.par = T, cell.border.alpha = 0.1,n.cores=48,main="RNA Velocity")
p3<-show.velocity.on.embedding.cor(emb, rvel.cd, n = 200, scale = 'sqrt',
                               cell.colors = ac(colors, alpha = 0.5),
                               cex = 0.8, arrow.scale = 50, show.grid.flow = T,
                               min.grid.cell.mass = 0.5, grid.n = 40,
                               arrow.lwd = 1,do.par = T, cell.border.alpha = 0.1,n.cores=48,main="RNA Velocity")
DimPlot(YDL, reduction = 'tsne')

dev.off()
savehistory("code.txt")

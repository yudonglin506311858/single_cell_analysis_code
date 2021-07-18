library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)

data(gcSample)
lapply(gcSample, head)
ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG")
head(as.data.frame(ck))

#visualize the result using dotplot method.
dotplot(ck)

ck <- compareCluster(geneCluster = gcSample, fun = "enrichGO",OrgDb = org.Mm.eg.db)
dotplot(ck, showCategory=10)

#合并四个时期的BP，以热图展示
a <- read.csv("allmarker_PURE.csv")
b<-a[a$cluster=="proe","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene_proe <- eg[,2]
gene_proe<-as.data.frame(gene_proe)
colnames(gene_proe)<-c("pro")
head(gene_proe)
gene_proe<-c(gene_proe)


a <- read.csv("allmarker_PURE.csv")
b<-a[a$cluster=="baso","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene_baso <- eg[,2]
gene_baso<-as.data.frame(gene_baso)
colnames(gene_baso)<-c("baso")
head(gene_baso)
gene_baso<-c(gene_baso)

a <- read.csv("allmarker_PURE.csv")
b<-a[a$cluster=="poly","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene_poly <- eg[,2]
gene_poly<-as.data.frame(gene_poly)
colnames(gene_poly)<-c("poly")
head(gene_poly)
gene_poly<-c(gene_poly)

a <- read.csv("allmarker_PURE.csv")
b<-a[a$cluster=="orth","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene_ortho <- eg[,2]
gene_ortho<-as.data.frame(gene_ortho)
colnames(gene_ortho)<-c("ortho")
head(gene_ortho)
class(gene_ortho)
gene_ortho<-c(gene_ortho)

data<-structure(list(pro =gene_proe, baso = gene_baso, poly =gene_poly, ortho=gene_ortho))
data<-list(pro =gene_proe$pro, baso = gene_baso$baso, poly =gene_poly$poly, ortho=gene_ortho$ortho)

lapply(data, head)


head(as.data.frame(ck))

ck <- compareCluster(geneCluster = data,OrgDb = org.Mm.eg.db, fun = "enrichGO", pvalueCutoff=0.05)
ck_1 <- compareCluster(geneCluster = data, fun = "enrichKEGG")



ck_enrichKEGG <- compareCluster(geneCluster = data, fun = "enrichKEGG")
ck_enrichDO <- compareCluster(geneCluster = data, fun = "enrichGO",OrgDb = org.Mm.eg.db)
dotplot(ck_enrichKEGG)
dotplot(ck_enrichDO, showCategory =15)
#visualize the result using dotplot method.

pdf("多样本富集分析.pdf",width=10,height=8)
dotplot(ck, showCategory =20)
dotplot(ck, showCategory =20,split="ONTOLOGY") 
dev.off()


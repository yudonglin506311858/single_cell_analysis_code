if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install('clusterProfiler')

library(clusterProfiler)
library("org.Mm.eg.db")

a <- read.table("genes_in_heatmap_clusters.txt",colClasses = "character",header = T)
b<-a[a$Cluster=="1","Gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)

ego <- enrichGO(
  
  gene          = gene,
  
  keyType = "ENTREZID",
  
  OrgDb         = org.Mm.eg.db,
  
  ont           = "MF",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff  = 0.05,
  
  qvalueCutoff  = 0.05,
  
  readable      = TRUE)

barplot(ego)



dotplot(ego)



pdf("GO_MF_cluster1.pdf")

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_cluster1.csv")



#生物过程(biologicalprocess)

ego <- enrichGO(
  
  gene          = gene,
  
  keyType = "ENTREZID",
  
  OrgDb         = org.Mm.eg.db,
  
  ont           = "BP",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff  = 0.05,
  
  qvalueCutoff  = 0.05,
  
  readable      = TRUE)

barplot(ego)



dotplot(ego)



pdf("GO_BP_cluster1.pdf")

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_cluster1.csv")



#细胞组成(cellularcomponent)

ego <- enrichGO(
  
  gene          = gene,
  
  keyType = "ENTREZID",
  
  OrgDb         = org.Mm.eg.db,
  
  ont           = "CC",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff  = 0.05,
  
  qvalueCutoff  = 0.05,
  
  readable      = TRUE)

barplot(ego)



dotplot(ego)



pdf("GO_CC_cluster1.pdf")

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_cluster1.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_cluster1.csv")

pdf("KEGG_cluster1.pdf")

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_cluster1.pdf")
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#模式2
a <- read.table("genes_in_heatmap_clusters.txt",colClasses = "character",header = T)
b<-a[a$Cluster=="2","Gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)

ego <- enrichGO(
  
  gene          = gene,
  
  keyType = "ENTREZID",
  
  OrgDb         = org.Mm.eg.db,
  
  ont           = "MF",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff  = 0.05,
  
  qvalueCutoff  = 0.05,
  
  readable      = TRUE)

barplot(ego)



dotplot(ego)



pdf("GO_MF_cluster2.pdf")

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_cluster2.csv")



#生物过程(biologicalprocess)

ego <- enrichGO(
  
  gene          = gene,
  
  keyType = "ENTREZID",
  
  OrgDb         = org.Mm.eg.db,
  
  ont           = "BP",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff  = 0.05,
  
  qvalueCutoff  = 0.05,
  
  readable      = TRUE)

barplot(ego)



dotplot(ego)



pdf("GO_BP_cluster2.pdf")

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_cluster2.csv")



#细胞组成(cellularcomponent)

ego <- enrichGO(
  
  gene          = gene,
  
  keyType = "ENTREZID",
  
  OrgDb         = org.Mm.eg.db,
  
  ont           = "CC",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff  = 0.05,
  
  qvalueCutoff  = 0.05,
  
  readable      = TRUE)

barplot(ego)



dotplot(ego)



pdf("GO_CC_cluster2.pdf")

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_cluster2.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_cluster2.csv")

pdf("KEGG_cluster2.pdf")

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_cluster2.pdf")
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()



#模式三：
a <- read.table("genes_in_heatmap_clusters.txt",colClasses = "character",header = T)
b<-a[a$Cluster=="3","Gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)

ego <- enrichGO(
  
  gene          = gene,
  
  keyType = "ENTREZID",
  
  OrgDb         = org.Mm.eg.db,
  
  ont           = "MF",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff  = 0.05,
  
  qvalueCutoff  = 0.05,
  
  readable      = TRUE)

barplot(ego)



dotplot(ego)



pdf("GO_MF_cluster3.pdf")

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_cluster3.csv")



#生物过程(biologicalprocess)

ego <- enrichGO(
  
  gene          = gene,
  
  keyType = "ENTREZID",
  
  OrgDb         = org.Mm.eg.db,
  
  ont           = "BP",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff  = 0.05,
  
  qvalueCutoff  = 0.05,
  
  readable      = TRUE)

barplot(ego)



dotplot(ego)



pdf("GO_BP_cluster3.pdf")

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_cluster3.csv")



#细胞组成(cellularcomponent)

ego <- enrichGO(
  
  gene          = gene,
  
  keyType = "ENTREZID",
  
  OrgDb         = org.Mm.eg.db,
  
  ont           = "CC",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff  = 0.05,
  
  qvalueCutoff  = 0.05,
  
  readable      = TRUE)

barplot(ego)



dotplot(ego)



pdf("GO_CC_cluster3.pdf")

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_cluster3.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_cluster3.csv")

pdf("KEGG_cluster3.pdf")

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_cluster3.pdf")
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()




#模式4：
a <- read.table("genes_in_heatmap_clusters.txt",colClasses = "character",header = T)
b<-a[a$Cluster=="4","Gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)

ego <- enrichGO(
  
  gene          = gene,
  
  keyType = "ENTREZID",
  
  OrgDb         = org.Mm.eg.db,
  
  ont           = "MF",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff  = 0.05,
  
  qvalueCutoff  = 0.05,
  
  readable      = TRUE)

barplot(ego)



dotplot(ego)



pdf("GO_MF_cluster4.pdf")

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_cluster4.csv")



#生物过程(biologicalprocess)

ego <- enrichGO(
  
  gene          = gene,
  
  keyType = "ENTREZID",
  
  OrgDb         = org.Mm.eg.db,
  
  ont           = "BP",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff  = 0.05,
  
  qvalueCutoff  = 0.05,
  
  readable      = TRUE)

barplot(ego)



dotplot(ego)



pdf("GO_BP_cluster4.pdf")

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_cluster4.csv")



#细胞组成(cellularcomponent)

ego <- enrichGO(
  
  gene          = gene,
  
  keyType = "ENTREZID",
  
  OrgDb         = org.Mm.eg.db,
  
  ont           = "CC",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff  = 0.05,
  
  qvalueCutoff  = 0.05,
  
  readable      = TRUE)

barplot(ego)



dotplot(ego)



pdf("GO_CC_cluster4.pdf")

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_cluster4.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_cluster4.csv")

pdf("KEGG_cluster4.pdf")

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_cluster4.pdf")
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

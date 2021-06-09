#regulon的富集分析

library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)


#cluster0的regulon
data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "cell_0"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "cell_0"], decreasing = T)), fixed = T)
)
a <- data$label
a<-strsplit(a, "_")
a<-as.data.frame(unlist(a))
b<-a$`unlist(a)`
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



pdf("GO_MF_regulon0.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_regulon0.csv")



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



pdf("GO_BP_regulon0.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_regulon0.csv")



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



pdf("GO_CC_regulon0.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_regulon0.csv")



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

write.csv(ekegg,file="kegg_regulon0.csv")

pdf("KEGG_regulon0.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_regulon0.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()


#cluster1的regulon
data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "cell_1"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "cell_1"], decreasing = T)), fixed = T)
)
a <- data$label
a<-strsplit(a, "_")
a<-as.data.frame(unlist(a))
b<-a$`unlist(a)`
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



pdf("GO_MF_regulon1.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_regulon1.csv")



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



pdf("GO_BP_regulon1.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_regulon1.csv")



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



pdf("GO_CC_regulon1.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_regulon1.csv")



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

write.csv(ekegg,file="kegg_regulon1.csv")

pdf("KEGG_regulon1.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_regulon1.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()


#cluster2的regulon
data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "cell_2"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "cell_2"], decreasing = T)), fixed = T)
)
a <- data$label
a<-strsplit(a, "_")
a<-as.data.frame(unlist(a))
b<-a$`unlist(a)`
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



pdf("GO_MF_regulon2.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_regulon2.csv")



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



pdf("GO_BP_regulon2.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_regulon2.csv")



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



pdf("GO_CC_regulon2.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_regulon2.csv")



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

write.csv(ekegg,file="kegg_regulon2.csv")

pdf("KEGG_regulon2.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_regulon2.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()



#cluster3的regulon
data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "cell_3"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "cell_3"], decreasing = T)), fixed = T)
)
a <- data$label
a<-strsplit(a, "_")
a<-as.data.frame(unlist(a))
b<-a$`unlist(a)`
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



pdf("GO_MF_regulon3.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_regulon3.csv")



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



pdf("GO_BP_regulon3.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_regulon3.csv")



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



pdf("GO_CC_regulon3.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_regulon3.csv")



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

write.csv(ekegg,file="kegg_regulon3.csv")

pdf("KEGG_regulon3.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_regulon3.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()




#cluster4的regulon
data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "cell_4"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "cell_4"], decreasing = T)), fixed = T)
)
a <- data$label
a<-strsplit(a, "_")
a<-as.data.frame(unlist(a))
b<-a$`unlist(a)`
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



pdf("GO_MF_regulon4.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_regulon4.csv")



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



pdf("GO_BP_regulon4.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_regulon4.csv")



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



pdf("GO_CC_regulon4.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_regulon4.csv")



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

write.csv(ekegg,file="kegg_regulon4.csv")

pdf("KEGG_regulon4.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_regulon4.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()



#cluster5的regulon
data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "cell_5"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "cell_5"], decreasing = T)), fixed = T)
)
a <- data$label
a<-strsplit(a, "_")
a<-as.data.frame(unlist(a))
b<-a$`unlist(a)`
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



pdf("GO_MF_regulon5.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_regulon5.csv")



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



pdf("GO_BP_regulon5.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_regulon5.csv")



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



pdf("GO_CC_regulon5.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_regulon5.csv")



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

write.csv(ekegg,file="kegg_regulon5.csv")

pdf("KEGG_regulon5.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_regulon5.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()



#cluster6的regulon
data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "cell_6"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "cell_6"], decreasing = T)), fixed = T)
)
a <- data$label
a<-strsplit(a, "_")
a<-as.data.frame(unlist(a))
b<-a$`unlist(a)`
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



pdf("GO_MF_regulon6.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_regulon6.csv")



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



pdf("GO_BP_regulon6.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_regulon6.csv")



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



pdf("GO_CC_regulon6.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_regulon6.csv")



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

write.csv(ekegg,file="kegg_regulon6.csv")

pdf("KEGG_regulon6.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_regulon6.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()



#cluster7的regulon
data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "cell_7"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "cell_7"], decreasing = T)), fixed = T)
)
a <- data$label
a<-strsplit(a, "_")
a<-as.data.frame(unlist(a))
b<-a$`unlist(a)`
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



pdf("GO_MF_regulon7.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_regulon7.csv")



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



pdf("GO_BP_regulon7.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_regulon7.csv")



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



pdf("GO_CC_regulon7.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_regulon7.csv")



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

write.csv(ekegg,file="kegg_regulon7.csv")

pdf("KEGG_regulon7.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_regulon7.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#cluster8的regulon
data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "cell_8"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "cell_8"], decreasing = T)), fixed = T)
)
a <- data$label
a<-strsplit(a, "_")
a<-as.data.frame(unlist(a))
b<-a$`unlist(a)`
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



pdf("GO_MF_regulon8.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_regulon8.csv")



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



pdf("GO_BP_regulon8.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_regulon8.csv")



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



pdf("GO_CC_regulon8.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_regulon8.csv")



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

write.csv(ekegg,file="kegg_regulon8.csv")

pdf("KEGG_regulon8.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_regulon8.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#cluster9的regulon
data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "cell_9"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "cell_9"], decreasing = T)), fixed = T)
)
a <- data$label
a<-strsplit(a, "_")
a<-as.data.frame(unlist(a))
b<-a$`unlist(a)`
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



pdf("GO_MF_regulon9.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_regulon9.csv")



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



pdf("GO_BP_regulon9.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_regulon9.csv")



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



pdf("GO_CC_regulon9.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_regulon9.csv")



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

write.csv(ekegg,file="kegg_regulon9.csv")

pdf("KEGG_regulon9.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_regulon9.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#cluster10的regulon
data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "cell_10"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "cell_10"], decreasing = T)), fixed = T)
)
a <- data$label
a<-strsplit(a, "_")
a<-as.data.frame(unlist(a))
b<-a$`unlist(a)`
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



pdf("GO_MF_regulon10.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_regulon10.csv")



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



pdf("GO_BP_regulon10.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_regulon10.csv")



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



pdf("GO_CC_regulon10.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_regulon10.csv")



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

write.csv(ekegg,file="kegg_regulon10.csv")

pdf("KEGG_regulon10.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_regulon10.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#cluster11的regulon
data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "cell_11"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "cell_11"], decreasing = T)), fixed = T)
)
a <- data$label
a<-strsplit(a, "_")
a<-as.data.frame(unlist(a))
b<-a$`unlist(a)`
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



pdf("GO_MF_regulon11.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_regulon11.csv")



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



pdf("GO_BP_regulon11.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_regulon11.csv")



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



pdf("GO_CC_regulon11.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_regulon11.csv")



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

write.csv(ekegg,file="kegg_regulon11.csv")

pdf("KEGG_regulon11.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_regulon11.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#cluster12的regulon
data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "cell_12"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "cell_12"], decreasing = T)), fixed = T)
)
a <- data$label
a<-strsplit(a, "_")
a<-as.data.frame(unlist(a))
b<-a$`unlist(a)`
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



pdf("GO_MF_regulon12.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_regulon12.csv")



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



pdf("GO_BP_regulon12.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_regulon12.csv")



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



pdf("GO_CC_regulon12.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_regulon12.csv")



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

write.csv(ekegg,file="kegg_regulon12.csv")

pdf("KEGG_regulon12.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_regulon12.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#cluster13的regulon
data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "cell_13"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "cell_13"], decreasing = T)), fixed = T)
)
a <- data$label
a<-strsplit(a, "_")
a<-as.data.frame(unlist(a))
b<-a$`unlist(a)`
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



pdf("GO_MF_regulon13.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_regulon13.csv")



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



pdf("GO_BP_regulon13.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_regulon13.csv")



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



pdf("GO_CC_regulon13.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_regulon13.csv")



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

write.csv(ekegg,file="kegg_regulon13.csv")

pdf("KEGG_regulon13.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_regulon13.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

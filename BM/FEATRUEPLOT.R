pdf("血红蛋白_Erythroblasts .pdf")
FeaturePlot(object = Erythroblasts, reduction = "tsne",pt.size = 1.5,features = c("Hbb-bs","Hbb-bt","Hba-a1","Hba-a2"),cols = c("gray", "red"))#actin
dev.off()

pdf("转录因子_Erythroblasts .pdf")
FeaturePlot(object = Erythroblasts, reduction = "tsne",pt.size = 1.5,features = c("Gata1","Klf1","Tal1","Nfe2"),cols = c("gray", "red"))#actin
dev.off()

pdf("晚期基因_Erythroblasts .pdf")
FeaturePlot(object = Erythroblasts, reduction = "tsne",pt.size = 1.5,features = c("Xpo7","Fam46c","Cd36","Trim58"),cols = c("gray", "red"))#actin
dev.off()

pdf("B细胞1_Erythroblasts .pdf")
FeaturePlot(object = Erythroblasts, reduction = "tsne",pt.size = 1.5,features = c("Igkc","Ly6d","Cd79a","Cd79b"),cols = c("gray", "red"))#actin
dev.off()

pdf("中性粒细胞_Erythroblasts .pdf")
FeaturePlot(object = Erythroblasts, reduction = "tsne",pt.size = 1.5,features = c("Elane","Mpo","Lyz2","Prtn3"),cols = c("gray", "red"))#actin
FeaturePlot(object = Erythroblasts, reduction = "tsne",pt.size = 1.5,features = c("S100a8","S100a9","Lgals1","Lgals3"),cols = c("gray", "red"))#actin
dev.off()

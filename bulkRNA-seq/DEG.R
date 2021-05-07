#D:\analysis\anxiuli_bulkdata
#https://www-ncbi-nlm-nih-gov-s.webvpn.cams.cn/geo/query/acc.cgi?acc=GSE53983
setwd("D:/analysis/anxiuli_bulkdata")
GSE53983_All_mm_countData<-read.table("GSE53983_All_mm_countData.txt",header = T,row.names = 1)
GSE53983_All_mm_rpkmData<-read.table("GSE53983_All_mm_rpkmData.txt",header = T,row.names = 1)



library("biomaRt")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("DESeq2")
library(edgeR)
library(limma)
library("Rgraphviz")
#rm(list=ls())
getwd()
#count文件存储
#寻找pro相对于其他各个时期的
cts<-GSE53983_All_mm_countData

#PCA
condition<-factor(c("mm_proerythroblast","mm_proerythroblast","mm_proerythroblast","mm_basophilic","mm_basophilic","mm_basophilic","mm_polychromatic","mm_polychromatic","mm_polychromatic","mm_orthochromatic","mm_orthochromatic","mm_orthochromatic"),levels = c("mm_proerythroblast","mm_basophilic","mm_polychromatic","mm_orthochromatic"))
head(cts)
tmp<-cts[,c(1:12)]
head(tmp)
colData <- data.frame(row.names=colnames(cts[,c(1:12)]), condition)
head(colData,10)
dds_all<- DESeqDataSetFromMatrix(countData = cts[,c(1:12)],colData = colData,design= ~condition)
head(dds_all)
dds_all<- DESeq(dds_all)
vsd_all<-vst(dds_all,blind=FALSE)
head(vsd_all)
#dist(t(assay(vsd_all)))
plotPCA(vsd_all,intgroup="condition")



#样本的聚类图
sampleDists <- dist(t(assay(vsd_all)))
library("RColorBrewer")
library(pheatmap)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_all$condition, vsd_all$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)



#做一下差异表达分析:pro的差异表达基因（相对于其他所有细胞）
mycounts<-cts[,c(1:12)]
head(mycounts)

# 这一步很关键，要明白condition这里是因子，不是样本名称；小鼠数据有对照组和处理组，各两个重复
condition <- factor(c(rep("mm_proerythroblast",3),rep("control",9)), levels = c("control","mm_proerythroblast"))
condition

#colData也可以自己在excel做好另存为.csv格式，再导入即可
colData <- data.frame(row.names=colnames(mycounts), condition)
colData
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)

res = results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)
#所有结果先进行输出
write.csv(res,file="All_results_pro.csv")
table(res$padj<0.05)#一共3351个genes上调，3365个genes下调，没有离群值。padj小于0.05的共有7234个。


diff_gene <- subset(res, padj < 0.05 & (log2FoldChange < -1|log2FoldChange > 1))
head(diff_gene)
dim(diff_gene)
head(res)
write.csv(diff_gene,file = "diff_gene_all_pro.csv",row.names = T)
diff_gene_up <- subset(res, padj < 0.05 & (log2FoldChange > 1))
write.csv(diff_gene_up,file = "diff_gene_up_pro.csv",row.names = T)
diff_gene_down <- subset(res, padj < 0.05 & (log2FoldChange < -1))
write.csv(diff_gene_down,file = "diff_gene_down_pro.csv",row.names = T)
dim(diff_gene)
dim(diff_gene_up)
dim(diff_gene_down)
head(diff_gene)
head(diff_gene_up)
head(diff_gene_down)


pro_up <- subset(res, padj < 0.01 & (log2FoldChange > 2))
dim(pro_up)
head(pro_up)
result<-rownames(pro_up)
write.csv(result,"pro_symbol.csv")

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file = "DEG_all_pro.csv",row.names = F)
head(resdata)
summary(res[order(resdata$pvalue),])
table(resdata$padj<0.05)#number of true 小于0.05 的基因个数


pdf("vocanoplot_pro.pdf")
library(ggplot2)
dataset <-read.csv('DEG_all_pro.csv',header = TRUE)
dim(dataset)
head(dataset)
dataset <-na.omit(dataset)
cut_off_pvalue = 0.0000001
cut_off_log2FoldChange = 1
dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange>cut_off_log2FoldChange ,'Up','Down'), 'Stable')
p<- ggplot(dataset,aes(x = log2FoldChange, y = -log10(pvalue),colour=change))+labs(title= "vocano plot") +geom_point(alpha=0.4, size=3.5) +scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +geom_hline(yintercept =-log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +labs(x="log2(fold change)", y="-log10 (p-value)")+ theme_bw()+theme(plot.title = element_text(hjust = 0.5), legend.position="right",        legend.title = element_blank() )

library(ggrepel)
dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange> cut_off_log2FoldChange ,'Up','Down'), 'Stable')

dataset$label = ifelse(dataset$pvalue <cut_off_pvalue & dataset$log2FoldChange <= -5| dataset$log2FoldChange >= 2,as.character(dataset$Row.names),"")

p+geom_text_repel(data = dataset, aes(x =log2FoldChange, y =-log10(pvalue), label =label), size = 3,box.padding =unit(0.5, "lines"),point.padding = unit(0.8,"lines"),segment.color ="black",show.legend = FALSE)



# grammar
library(tidyverse)
library(magrittr)
library(glue)
library(data.table)

# analysis
library(DESeq2)

# graphics
library(ggplot2)
library(ggrepel)
library(ggsci)
library(scales)
library(latex2exp)
diffData <- fread('DEG_all_pro.csv')

colnames(diffData)[1] <- "gene"

diffData[is.na(padj), padj := 1][]
diffData[, p := -log10(padj)][]


diffData[, type := "ns"][]
diffData[log2FoldChange > 1 & padj < 0.05, type := "up"][log2FoldChange < -1 & padj < 0.05, type := "down"][]

labelGene <- diffData[order(p, decreasing = T)][type == "up"][1:10]

pal_nejm()(8) %>% show_col()
typeColor <- structure(
  c(pal_nejm()(2), "gray80"),
  names = c("up", "down", "ns")
)

ggplot(diffData, aes(x = log2FoldChange, y = p)) +
  geom_point(aes(color = type, size = p), show.legend = F) +
  geom_hline(yintercept = -log10(0.05), color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = 1, color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = -1, color = "gray60", linetype = "dashed") +
  geom_text_repel(
    data = labelGene, aes(label = gene),
    size = 3, fontface = 3,
    nudge_x = .5, nudge_y = .5) +
  scale_radius(range = c(.1, 2)) +
  scale_color_manual(values = typeColor) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  labs(
    x = TeX("$log_{2}(Fold\\,Change)$"),
    y = TeX("$-log_{10}(\\textit{P}\\,value)$")) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line())
dev.off()







#做一下差异表达分析:baso的差异表达基因（相对于其他所有细胞）
mycounts<-cts[,c(1:12)]
head(mycounts)

# 这一步很关键，要明白condition这里是因子，不是样本名称；小鼠数据有对照组和处理组，各两个重复
condition <- factor(c(rep("control",3),rep("mm_basophilic",3),rep("control",6)), levels = c("control","mm_basophilic"))
condition

#colData也可以自己在excel做好另存为.csv格式，再导入即可
colData <- data.frame(row.names=colnames(mycounts), condition)
colData
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)

res = results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)
#所有结果先进行输出
write.csv(res,file="All_results_baso.csv")
table(res$padj<0.05)#一共3351个genes上调，3365个genes下调，没有离群值。padj小于0.05的共有7234个。


diff_gene <- subset(res, padj < 0.05 & (log2FoldChange < -1|log2FoldChange > 1))
head(diff_gene)
dim(diff_gene)
head(res)
write.csv(diff_gene,file = "diff_gene_all_baso.csv",row.names = T)
diff_gene_up <- subset(res, padj < 0.05 & (log2FoldChange > 1))
write.csv(diff_gene_up,file = "diff_gene_up_baso.csv",row.names = T)
diff_gene_down <- subset(res, padj < 0.05 & (log2FoldChange < -1))
write.csv(diff_gene_down,file = "diff_gene_down_baso.csv",row.names = T)
dim(diff_gene)
dim(diff_gene_up)
dim(diff_gene_down)
head(diff_gene)
head(diff_gene_up)
head(diff_gene_down)


baso_up <- subset(res, padj < 0.05)
dim(baso_up)
head(baso_up)
result<-rownames(baso_up)
write.csv(result,"baso_symbol.csv")

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file = "DEG_all_baso.csv",row.names = F)
head(resdata)
summary(res[order(resdata$pvalue),])
table(resdata$padj<0.05)#number of true 小于0.05 的基因个数


pdf("vocanoplot_baso.pdf")
library(ggplot2)
dataset <-read.csv('DEG_all_baso.csv',header = TRUE)
dim(dataset)
head(dataset)
dataset <-na.omit(dataset)
cut_off_pvalue = 0.0000001
cut_off_log2FoldChange = 1
dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange>cut_off_log2FoldChange ,'Up','Down'), 'Stable')
p<- ggplot(dataset,aes(x = log2FoldChange, y = -log10(pvalue),colour=change))+labs(title= "vocano plot") +geom_point(alpha=0.4, size=3.5) +scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +geom_hline(yintercept =-log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +labs(x="log2(fold change)", y="-log10 (p-value)")+ theme_bw()+theme(plot.title = element_text(hjust = 0.5), legend.position="right",        legend.title = element_blank() )

library(ggrepel)
dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange> cut_off_log2FoldChange ,'Up','Down'), 'Stable')

dataset$label = ifelse(dataset$pvalue <cut_off_pvalue & dataset$log2FoldChange <= -5| dataset$log2FoldChange >= 2,as.character(dataset$Row.names),"")

p+geom_text_repel(data = dataset, aes(x =log2FoldChange, y =-log10(pvalue), label =label), size = 3,box.padding =unit(0.5, "lines"),point.padding = unit(0.8,"lines"),segment.color ="black",show.legend = FALSE)



# grammar
library(tidyverse)
library(magrittr)
library(glue)
library(data.table)

# analysis
library(DESeq2)

# graphics
library(ggplot2)
library(ggrepel)
library(ggsci)
library(scales)
library(latex2exp)
diffData <- fread('DEG_all_baso.csv')

colnames(diffData)[1] <- "gene"

diffData[is.na(padj), padj := 1][]
diffData[, p := -log10(padj)][]


diffData[, type := "ns"][]
diffData[log2FoldChange > 1 & padj < 0.05, type := "up"][log2FoldChange < -1 & padj < 0.05, type := "down"][]

labelGene <- diffData[order(p, decreasing = T)][type == "up"][1:10]

pal_nejm()(8) %>% show_col()
typeColor <- structure(
  c(pal_nejm()(2), "gray80"),
  names = c("up", "down", "ns")
)

ggplot(diffData, aes(x = log2FoldChange, y = p)) +
  geom_point(aes(color = type, size = p), show.legend = F) +
  geom_hline(yintercept = -log10(0.05), color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = 1, color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = -1, color = "gray60", linetype = "dashed") +
  geom_text_repel(
    data = labelGene, aes(label = gene),
    size = 3, fontface = 3,
    nudge_x = .5, nudge_y = .5) +
  scale_radius(range = c(.1, 2)) +
  scale_color_manual(values = typeColor) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  labs(
    x = TeX("$log_{2}(Fold\\,Change)$"),
    y = TeX("$-log_{10}(\\textit{P}\\,value)$")) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line())
dev.off()













#做一下差异表达分析:baso的差异表达基因（相对于其他所有细胞）
mycounts<-cts[,c(1:12)]
head(mycounts)

# 这一步很关键，要明白condition这里是因子，不是样本名称；小鼠数据有对照组和处理组，各两个重复
condition <- factor(c(rep("control",3),rep("mm_basophilic",3),rep("control",6)), levels = c("control","mm_basophilic"))
condition

#colData也可以自己在excel做好另存为.csv格式，再导入即可
colData <- data.frame(row.names=colnames(mycounts), condition)
colData
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)

res = results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)
#所有结果先进行输出
write.csv(res,file="All_results_baso.csv")
table(res$padj<0.05)#一共3351个genes上调，3365个genes下调，没有离群值。padj小于0.05的共有7234个。


diff_gene <- subset(res, padj < 0.05 & (log2FoldChange < -1|log2FoldChange > 1))
head(diff_gene)
dim(diff_gene)
head(res)
write.csv(diff_gene,file = "diff_gene_all_baso.csv",row.names = T)
diff_gene_up <- subset(res, padj < 0.05 & (log2FoldChange > 1))
write.csv(diff_gene_up,file = "diff_gene_up_baso.csv",row.names = T)
diff_gene_down <- subset(res, padj < 0.05 & (log2FoldChange < -1))
write.csv(diff_gene_down,file = "diff_gene_down_baso.csv",row.names = T)
dim(diff_gene)
dim(diff_gene_up)
dim(diff_gene_down)
head(diff_gene)
head(diff_gene_up)
head(diff_gene_down)


baso_up <- subset(res, padj < 0.05)
dim(baso_up)
head(baso_up)
result<-rownames(baso_up)
write.csv(result,"baso_symbol.csv")

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file = "DEG_all_baso.csv",row.names = F)
head(resdata)
summary(res[order(resdata$pvalue),])
table(resdata$padj<0.05)#number of true 小于0.05 的基因个数


pdf("vocanoplot_baso.pdf")
library(ggplot2)
dataset <-read.csv('DEG_all_baso.csv',header = TRUE)
dim(dataset)
head(dataset)
dataset <-na.omit(dataset)
cut_off_pvalue = 0.0000001
cut_off_log2FoldChange = 1
dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange>cut_off_log2FoldChange ,'Up','Down'), 'Stable')
p<- ggplot(dataset,aes(x = log2FoldChange, y = -log10(pvalue),colour=change))+labs(title= "vocano plot") +geom_point(alpha=0.4, size=3.5) +scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +geom_hline(yintercept =-log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +labs(x="log2(fold change)", y="-log10 (p-value)")+ theme_bw()+theme(plot.title = element_text(hjust = 0.5), legend.position="right",        legend.title = element_blank() )

library(ggrepel)
dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange> cut_off_log2FoldChange ,'Up','Down'), 'Stable')

dataset$label = ifelse(dataset$pvalue <cut_off_pvalue & dataset$log2FoldChange <= -5| dataset$log2FoldChange >= 2,as.character(dataset$Row.names),"")

p+geom_text_repel(data = dataset, aes(x =log2FoldChange, y =-log10(pvalue), label =label), size = 3,box.padding =unit(0.5, "lines"),point.padding = unit(0.8,"lines"),segment.color ="black",show.legend = FALSE)



# grammar
library(tidyverse)
library(magrittr)
library(glue)
library(data.table)

# analysis
library(DESeq2)

# graphics
library(ggplot2)
library(ggrepel)
library(ggsci)
library(scales)
library(latex2exp)
diffData <- fread('DEG_all_baso.csv')

colnames(diffData)[1] <- "gene"

diffData[is.na(padj), padj := 1][]
diffData[, p := -log10(padj)][]


diffData[, type := "ns"][]
diffData[log2FoldChange > 1 & padj < 0.05, type := "up"][log2FoldChange < -1 & padj < 0.05, type := "down"][]

labelGene <- diffData[order(p, decreasing = T)][type == "up"][1:10]

pal_nejm()(8) %>% show_col()
typeColor <- structure(
  c(pal_nejm()(2), "gray80"),
  names = c("up", "down", "ns")
)

ggplot(diffData, aes(x = log2FoldChange, y = p)) +
  geom_point(aes(color = type, size = p), show.legend = F) +
  geom_hline(yintercept = -log10(0.05), color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = 1, color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = -1, color = "gray60", linetype = "dashed") +
  geom_text_repel(
    data = labelGene, aes(label = gene),
    size = 3, fontface = 3,
    nudge_x = .5, nudge_y = .5) +
  scale_radius(range = c(.1, 2)) +
  scale_color_manual(values = typeColor) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  labs(
    x = TeX("$log_{2}(Fold\\,Change)$"),
    y = TeX("$-log_{10}(\\textit{P}\\,value)$")) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line())
dev.off()













#做一下差异表达分析:poly的差异表达基因（相对于其他所有细胞）
mycounts<-cts[,c(1:12)]
head(mycounts)

# 这一步很关键，要明白condition这里是因子，不是样本名称；小鼠数据有对照组和处理组，各两个重复
condition <- factor(c(rep("control",6),rep("mm_polychromatic",3),rep("control",3)), levels = c("control","mm_polychromatic"))
condition

#colData也可以自己在excel做好另存为.csv格式，再导入即可
colData <- data.frame(row.names=colnames(mycounts), condition)
colData
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)

res = results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)
#所有结果先进行输出
write.csv(res,file="All_results_poly.csv")
table(res$padj<0.05)#一共3351个genes上调，3365个genes下调，没有离群值。padj小于0.05的共有7234个。


diff_gene <- subset(res, padj < 0.05)
head(diff_gene)
dim(diff_gene)
head(res)
write.csv(diff_gene,file = "diff_gene_all_poly.csv",row.names = T)
diff_gene_up <- subset(res, padj < 0.05 & (log2FoldChange > 0))
write.csv(diff_gene_up,file = "diff_gene_up_poly.csv",row.names = T)
diff_gene_down <- subset(res, padj < 0.05 & (log2FoldChange < 0))
write.csv(diff_gene_down,file = "diff_gene_down_poly.csv",row.names = T)
dim(diff_gene)
dim(diff_gene_up)
dim(diff_gene_down)
head(diff_gene)
head(diff_gene_up)
head(diff_gene_down)


poly_up <- subset(res, padj < 0.05)
dim(poly_up)
head(poly_up)
result<-rownames(poly_up)
write.csv(result,"poly_symbol.csv")

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file = "DEG_all_poly.csv",row.names = F)
head(resdata)
summary(res[order(resdata$pvalue),])
table(resdata$padj<0.05)#number of true 小于0.05 的基因个数


pdf("vocanoplot_poly.pdf")
library(ggplot2)
dataset <-read.csv('DEG_all_poly.csv',header = TRUE)
dim(dataset)
head(dataset)
dataset <-na.omit(dataset)
cut_off_pvalue = 0.0000001
cut_off_log2FoldChange = 1
dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange>cut_off_log2FoldChange ,'Up','Down'), 'Stable')
p<- ggplot(dataset,aes(x = log2FoldChange, y = -log10(pvalue),colour=change))+labs(title= "vocano plot") +geom_point(alpha=0.4, size=3.5) +scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +geom_hline(yintercept =-log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +labs(x="log2(fold change)", y="-log10 (p-value)")+ theme_bw()+theme(plot.title = element_text(hjust = 0.5), legend.position="right",        legend.title = element_blank() )

library(ggrepel)
dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange> cut_off_log2FoldChange ,'Up','Down'), 'Stable')

dataset$label = ifelse(dataset$pvalue <cut_off_pvalue & dataset$log2FoldChange <= -5| dataset$log2FoldChange >= 2,as.character(dataset$Row.names),"")

p+geom_text_repel(data = dataset, aes(x =log2FoldChange, y =-log10(pvalue), label =label), size = 3,box.padding =unit(0.5, "lines"),point.padding = unit(0.8,"lines"),segment.color ="black",show.legend = FALSE)



# grammar
library(tidyverse)
library(magrittr)
library(glue)
library(data.table)

# analysis
library(DESeq2)

# graphics
library(ggplot2)
library(ggrepel)
library(ggsci)
library(scales)
library(latex2exp)
diffData <- fread('DEG_all_poly.csv')

colnames(diffData)[1] <- "gene"

diffData[is.na(padj), padj := 1][]
diffData[, p := -log10(padj)][]


diffData[, type := "ns"][]
diffData[log2FoldChange > 1 & padj < 0.05, type := "up"][log2FoldChange < -1 & padj < 0.05, type := "down"][]

labelGene <- diffData[order(p, decreasing = T)][type == "up"][1:10]

pal_nejm()(8) %>% show_col()
typeColor <- structure(
  c(pal_nejm()(2), "gray80"),
  names = c("up", "down", "ns")
)

ggplot(diffData, aes(x = log2FoldChange, y = p)) +
  geom_point(aes(color = type, size = p), show.legend = F) +
  geom_hline(yintercept = -log10(0.05), color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = 1, color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = -1, color = "gray60", linetype = "dashed") +
  geom_text_repel(
    data = labelGene, aes(label = gene),
    size = 3, fontface = 3,
    nudge_x = .5, nudge_y = .5) +
  scale_radius(range = c(.1, 2)) +
  scale_color_manual(values = typeColor) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  labs(
    x = TeX("$log_{2}(Fold\\,Change)$"),
    y = TeX("$-log_{10}(\\textit{P}\\,value)$")) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line())
dev.off()




#做一下差异表达分析:ortho的差异表达基因（相对于其他所有细胞）
mycounts<-cts[,c(1:12)]
head(mycounts)

# 这一步很关键，要明白condition这里是因子，不是样本名称；小鼠数据有对照组和处理组，各两个重复
condition <- factor(c(rep("control",9),rep("mm_orthochromatic",3)), levels = c("control","mm_orthochromatic"))
condition

#colData也可以自己在excel做好另存为.csv格式，再导入即可
colData <- data.frame(row.names=colnames(mycounts), condition)
colData
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)

res = results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)
#所有结果先进行输出
write.csv(res,file="All_results_ortho.csv")
table(res$padj<0.05)#一共3351个genes上调，3365个genes下调，没有离群值。padj小于0.05的共有7234个。


diff_gene <- subset(res, padj < 0.05 & (log2FoldChange < -1|log2FoldChange > 1))
head(diff_gene)
dim(diff_gene)
head(res)
write.csv(diff_gene,file = "diff_gene_all_ortho.csv",row.names = T)
diff_gene_up <- subset(res, padj < 0.05 & (log2FoldChange > 1))
write.csv(diff_gene_up,file = "diff_gene_up_ortho.csv",row.names = T)
diff_gene_down <- subset(res, padj < 0.05 & (log2FoldChange < -1))
write.csv(diff_gene_down,file = "diff_gene_down_ortho.csv",row.names = T)
dim(diff_gene)
dim(diff_gene_up)
dim(diff_gene_down)
head(diff_gene)
head(diff_gene_up)
head(diff_gene_down)


ortho_up <- subset(res, padj < 0.05 & (log2FoldChange > 1))

dim(ortho_up)
head(ortho_up)
result<-rownames(ortho_up)
write.csv(result,"ortho_symbol.csv")

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file = "DEG_all_ortho.csv",row.names = F)
head(resdata)
summary(res[order(resdata$pvalue),])
table(resdata$padj<0.05)#number of true 小于0.05 的基因个数


pdf("vocanoplot_ortho.pdf")
library(ggplot2)
dataset <-read.csv('DEG_all_ortho.csv',header = TRUE)
dim(dataset)
head(dataset)
dataset <-na.omit(dataset)
cut_off_pvalue = 0.0000001
cut_off_log2FoldChange = 1
dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange>cut_off_log2FoldChange ,'Up','Down'), 'Stable')
p<- ggplot(dataset,aes(x = log2FoldChange, y = -log10(pvalue),colour=change))+labs(title= "vocano plot") +geom_point(alpha=0.4, size=3.5) +scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +geom_hline(yintercept =-log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +labs(x="log2(fold change)", y="-log10 (p-value)")+ theme_bw()+theme(plot.title = element_text(hjust = 0.5), legend.position="right",        legend.title = element_blank() )

library(ggrepel)
dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange> cut_off_log2FoldChange ,'Up','Down'), 'Stable')

dataset$label = ifelse(dataset$pvalue <cut_off_pvalue & dataset$log2FoldChange <= -5| dataset$log2FoldChange >= 2,as.character(dataset$Row.names),"")

p+geom_text_repel(data = dataset, aes(x =log2FoldChange, y =-log10(pvalue), label =label), size = 3,box.padding =unit(0.5, "lines"),point.padding = unit(0.8,"lines"),segment.color ="black",show.legend = FALSE)



# grammar
library(tidyverse)
library(magrittr)
library(glue)
library(data.table)

# analysis
library(DESeq2)

# graphics
library(ggplot2)
library(ggrepel)
library(ggsci)
library(scales)
library(latex2exp)
diffData <- fread('DEG_all_ortho.csv')

colnames(diffData)[1] <- "gene"

diffData[is.na(padj), padj := 1][]
diffData[, p := -log10(padj)][]


diffData[, type := "ns"][]
diffData[log2FoldChange > 1 & padj < 0.05, type := "up"][log2FoldChange < -1 & padj < 0.05, type := "down"][]

labelGene <- diffData[order(p, decreasing = T)][type == "up"][1:10]

pal_nejm()(8) %>% show_col()
typeColor <- structure(
  c(pal_nejm()(2), "gray80"),
  names = c("up", "down", "ns")
)

ggplot(diffData, aes(x = log2FoldChange, y = p)) +
  geom_point(aes(color = type, size = p), show.legend = F) +
  geom_hline(yintercept = -log10(0.05), color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = 1, color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = -1, color = "gray60", linetype = "dashed") +
  geom_text_repel(
    data = labelGene, aes(label = gene),
    size = 3, fontface = 3,
    nudge_x = .5, nudge_y = .5) +
  scale_radius(range = c(.1, 2)) +
  scale_color_manual(values = typeColor) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  labs(
    x = TeX("$log_{2}(Fold\\,Change)$"),
    y = TeX("$-log_{10}(\\textit{P}\\,value)$")) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line())
dev.off()







#baso和poly的差异表达基因太少，因此去除太相近的细胞进行对比

#做一下差异表达分析:baso的差异表达基因（相对于其他所有细胞（除poly））
mycounts<-cts[,c(1,2,3,4,5,6,10,11,12)]
head(mycounts)

# 这一步很关键，要明白condition这里是因子，不是样本名称；小鼠数据有对照组和处理组，各两个重复
condition <- factor(c(rep("control",3),rep("mm_basophilic",3),rep("control",3)), levels = c("control","mm_basophilic"))
condition

#colData也可以自己在excel做好另存为.csv格式，再导入即可
colData <- data.frame(row.names=colnames(mycounts), condition)
colData
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)

res = results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)
#所有结果先进行输出
write.csv(res,file="All_results_baso_1.csv")
table(res$padj<0.05)#一共3351个genes上调，3365个genes下调，没有离群值。padj小于0.05的共有7234个。


diff_gene <- subset(res, padj < 0.05 & (log2FoldChange < -1|log2FoldChange > 1))
head(diff_gene)
dim(diff_gene)
head(res)
write.csv(diff_gene,file = "diff_gene_all_baso_1.csv",row.names = T)
diff_gene_up <- subset(res, padj < 0.05 & (log2FoldChange > 1))
write.csv(diff_gene_up,file = "diff_gene_up_baso_1.csv",row.names = T)
diff_gene_down <- subset(res, padj < 0.05 & (log2FoldChange < -1))
write.csv(diff_gene_down,file = "diff_gene_down_baso_1.csv",row.names = T)
dim(diff_gene)
dim(diff_gene_up)
dim(diff_gene_down)
head(diff_gene)
head(diff_gene_up)
head(diff_gene_down)


baso_up <- subset(res, padj < 0.05& (log2FoldChange > 1))
dim(baso_up)
head(baso_up)
result<-rownames(baso_up)
write.csv(result,"baso_symbol_1.csv")

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file = "DEG_all_baso_1.csv",row.names = F)
head(resdata)
summary(res[order(resdata$pvalue),])
table(resdata$padj<0.05)#number of true 小于0.05 的基因个数


pdf("vocanoplot_baso_1.pdf")
library(ggplot2)
dataset <-read.csv('DEG_all_baso_1.csv',header = TRUE)
dim(dataset)
head(dataset)
dataset <-na.omit(dataset)
cut_off_pvalue = 0.0000001
cut_off_log2FoldChange = 1
dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange>cut_off_log2FoldChange ,'Up','Down'), 'Stable')
p<- ggplot(dataset,aes(x = log2FoldChange, y = -log10(pvalue),colour=change))+labs(title= "vocano plot") +geom_point(alpha=0.4, size=3.5) +scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +geom_hline(yintercept =-log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +labs(x="log2(fold change)", y="-log10 (p-value)")+ theme_bw()+theme(plot.title = element_text(hjust = 0.5), legend.position="right",        legend.title = element_blank() )

library(ggrepel)
dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange> cut_off_log2FoldChange ,'Up','Down'), 'Stable')

dataset$label = ifelse(dataset$pvalue <cut_off_pvalue & dataset$log2FoldChange <= -5| dataset$log2FoldChange >= 2,as.character(dataset$Row.names),"")

p+geom_text_repel(data = dataset, aes(x =log2FoldChange, y =-log10(pvalue), label =label), size = 3,box.padding =unit(0.5, "lines"),point.padding = unit(0.8,"lines"),segment.color ="black",show.legend = FALSE)



# grammar
library(tidyverse)
library(magrittr)
library(glue)
library(data.table)

# analysis
library(DESeq2)

# graphics
library(ggplot2)
library(ggrepel)
library(ggsci)
library(scales)
library(latex2exp)
diffData <- fread('DEG_all_baso_1.csv')

colnames(diffData)[1] <- "gene"

diffData[is.na(padj), padj := 1][]
diffData[, p := -log10(padj)][]


diffData[, type := "ns"][]
diffData[log2FoldChange > 1 & padj < 0.05, type := "up"][log2FoldChange < -1 & padj < 0.05, type := "down"][]

labelGene <- diffData[order(p, decreasing = T)][type == "up"][1:10]

pal_nejm()(8) %>% show_col()
typeColor <- structure(
  c(pal_nejm()(2), "gray80"),
  names = c("up", "down", "ns")
)

ggplot(diffData, aes(x = log2FoldChange, y = p)) +
  geom_point(aes(color = type, size = p), show.legend = F) +
  geom_hline(yintercept = -log10(0.05), color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = 1, color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = -1, color = "gray60", linetype = "dashed") +
  geom_text_repel(
    data = labelGene, aes(label = gene),
    size = 3, fontface = 3,
    nudge_x = .5, nudge_y = .5) +
  scale_radius(range = c(.1, 2)) +
  scale_color_manual(values = typeColor) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  labs(
    x = TeX("$log_{2}(Fold\\,Change)$"),
    y = TeX("$-log_{10}(\\textit{P}\\,value)$")) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line())
dev.off()











#做一下差异表达分析:poly的差异表达基因（相对于其他所有细胞）
mycounts<-cts[,c(1,2,3,7,8,9,10,11,12)]
head(mycounts)

# 这一步很关键，要明白condition这里是因子，不是样本名称；小鼠数据有对照组和处理组，各两个重复
condition <- factor(c(rep("control",3),rep("mm_polychromatic",3),rep("control",3)), levels = c("control","mm_polychromatic"))
condition

#colData也可以自己在excel做好另存为.csv格式，再导入即可
colData <- data.frame(row.names=colnames(mycounts), condition)
colData
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)

res = results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)
#所有结果先进行输出
write.csv(res,file="All_results_poly_1.csv")
table(res$padj<0.05)#一共3351个genes上调，3365个genes下调，没有离群值。padj小于0.05的共有7234个。


diff_gene <- subset(res, padj < 0.05)
head(diff_gene)
dim(diff_gene)
head(res)
write.csv(diff_gene,file = "diff_gene_all_poly_1.csv",row.names = T)
diff_gene_up <- subset(res, padj < 0.05 & (log2FoldChange > 0))
write.csv(diff_gene_up,file = "diff_gene_up_poly_1.csv",row.names = T)
diff_gene_down <- subset(res, padj < 0.05 & (log2FoldChange < 0))
write.csv(diff_gene_down,file = "diff_gene_down_poly_1.csv",row.names = T)
dim(diff_gene)
dim(diff_gene_up)
dim(diff_gene_down)
head(diff_gene)
head(diff_gene_up)
head(diff_gene_down)


poly_up <- subset(res, padj < 0.05& (log2FoldChange > 0))
dim(poly_up)
head(poly_up)
result<-rownames(poly_up)
write.csv(result,"poly_symbol_1.csv")

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file = "DEG_all_poly_1.csv",row.names = F)
head(resdata)
summary(res[order(resdata$pvalue),])
table(resdata$padj<0.05)#number of true 小于0.05 的基因个数


pdf("vocanoplot_poly_1.pdf")
library(ggplot2)
dataset <-read.csv('DEG_all_poly_1.csv',header = TRUE)
dim(dataset)
head(dataset)
dataset <-na.omit(dataset)
cut_off_pvalue = 0.0000001
cut_off_log2FoldChange = 1
dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange>cut_off_log2FoldChange ,'Up','Down'), 'Stable')
p<- ggplot(dataset,aes(x = log2FoldChange, y = -log10(pvalue),colour=change))+labs(title= "vocano plot") +geom_point(alpha=0.4, size=3.5) +scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +geom_hline(yintercept =-log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +labs(x="log2(fold change)", y="-log10 (p-value)")+ theme_bw()+theme(plot.title = element_text(hjust = 0.5), legend.position="right",        legend.title = element_blank() )

library(ggrepel)
dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange> cut_off_log2FoldChange ,'Up','Down'), 'Stable')

dataset$label = ifelse(dataset$pvalue <cut_off_pvalue & dataset$log2FoldChange <= -5| dataset$log2FoldChange >= 2,as.character(dataset$Row.names),"")

p+geom_text_repel(data = dataset, aes(x =log2FoldChange, y =-log10(pvalue), label =label), size = 3,box.padding =unit(0.5, "lines"),point.padding = unit(0.8,"lines"),segment.color ="black",show.legend = FALSE)



# grammar
library(tidyverse)
library(magrittr)
library(glue)
library(data.table)

# analysis
library(DESeq2)

# graphics
library(ggplot2)
library(ggrepel)
library(ggsci)
library(scales)
library(latex2exp)
diffData <- fread('DEG_all_poly_1.csv')

colnames(diffData)[1] <- "gene"

diffData[is.na(padj), padj := 1][]
diffData[, p := -log10(padj)][]


diffData[, type := "ns"][]
diffData[log2FoldChange > 1 & padj < 0.05, type := "up"][log2FoldChange < -1 & padj < 0.05, type := "down"][]

labelGene <- diffData[order(p, decreasing = T)][type == "up"][1:10]

pal_nejm()(8) %>% show_col()
typeColor <- structure(
  c(pal_nejm()(2), "gray80"),
  names = c("up", "down", "ns")
)

ggplot(diffData, aes(x = log2FoldChange, y = p)) +
  geom_point(aes(color = type, size = p), show.legend = F) +
  geom_hline(yintercept = -log10(0.05), color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = 1, color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = -1, color = "gray60", linetype = "dashed") +
  geom_text_repel(
    data = labelGene, aes(label = gene),
    size = 3, fontface = 3,
    nudge_x = .5, nudge_y = .5) +
  scale_radius(range = c(.1, 2)) +
  scale_color_manual(values = typeColor) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  labs(
    x = TeX("$log_{2}(Fold\\,Change)$"),
    y = TeX("$-log_{10}(\\textit{P}\\,value)$")) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line())
dev.off()

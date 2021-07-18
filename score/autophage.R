AUTOPHAGE<-read.table("D:/analysis/forpublication/immunecell/李程_各个时期的中性粒细胞marker/AUTOPHAGE_HUMAN.txt",header = T)
library(biomaRt)
# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"),
                   filters = "mgi_symbol",
                   values = x , mart = mouse,
                   attributesL = c("hgnc_symbol"),
                   martL = human, uniqueRows=T)
  #humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  #print(head(genes))
  return(genesV2)
}
# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"),
                   filters = "hgnc_symbol",
                   values = x , mart = human,
                   attributesL = c("mgi_symbol"),
                   martL = mouse, uniqueRows=T)
  #mousex <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  #print(head(genes))
  return(genesV2)
}

a<-convertMouseGeneList(AUTOPHAGE$Symbol)
AUTOPHAGE<-a$MGI.symbol
AUTOPHAGE<-as.data.frame(AUTOPHAGE)
write.table(AUTOPHAGE,"AUTOPHAGE_MOUSE.txt",row.names = F)



YDL<-readRDS("YDL.RDS")
nonpro<-YDL

#autophage的打分

AUTOPHAGE_score<-read.table("D:/analysis/forpublication/immunecell/李程_各个时期的中性粒细胞marker/AUTOPHAGE_MOUSE.txt",header = T)
AUTOPHAGE_score<-AUTOPHAGE_score[,-1]
head(AUTOPHAGE_score)
AUTOPHAGE_score<-as.data.frame(AUTOPHAGE_score)
colnames(AUTOPHAGE_score)<-c("AUTOPHAGE_score")
AUTOPHAGE_score<-list(AUTOPHAGE_score$AUTOPHAGE_score)
nonpro <- AddModuleScore(
  object = nonpro,
  features = AUTOPHAGE_score,
  name = "AUTOPHAGE_score"
)
head(nonpro@meta.data)
colnames(nonpro@meta.data)[13] <- 'AUTOPHAGE_score' 
VlnPlot(object = nonpro, features = c("AUTOPHAGE_score"))
summary(nonpro@meta.data$AUTOPHAGE_score_Feature)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.259960 -0.004752  0.097868  0.084557  0.152367  2.377477 
boxplot(nonpro@meta.data$pro_score_Feature)


#绘制分数的分布图
library(ggplot2)
mydata<- FetchData(nonpro,vars = c("tSNE_1","tSNE_2","AUTOPHAGE_score"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = AUTOPHAGE_score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


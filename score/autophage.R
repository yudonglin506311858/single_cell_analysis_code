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

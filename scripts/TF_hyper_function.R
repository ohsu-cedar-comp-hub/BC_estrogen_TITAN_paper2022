library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(parallel)
factors <- read.delim("human_factor.txt", header=T, sep = "\t")
MCF7 <- factors %>% filter(Cell_line == "MCF-7") 
T47D <- factors %>% filter(Cell_line == "T47D")
ZR751 <- factors %>% filter(Cell_line == "ZR-75-1")
TFs <- rbind(ZR751, T47D, MCF7)
Genes <- readRDS("Genes_top2000.RDS") #This is the entire geneset inputted into the topic model (e.g. top 2000 variable genes in the dataset)
Genes <- unlist(Genes)
TFs <- TFs[-1*which(Genes == "error"),]
Genes <- Genes[-1*which(Genes == "error")]
ModelList <- read.table("Model_MCF7_E2_20T_CLR_5000Variable_M10_top50_genes_topics.txt", #This is a dataframe of the top 50 genes for each topic, can be extracted from a model using TopTopicGenes function in TITAN
                        header = T, sep = "\t")
TFPhyper <- function(GeneListTopic, Gene) {
  TopicGenes <- GeneListTopic
  ChipGenes <- unique(unlist(strsplit(Gene, ",")))
  TotalGenes <- length(as.data.frame(org.Hs.egSYMBOL2EG)[,2]) -50
  overlapLegnth <- length(intersect(ChipGenes, TopicGenes)) -1
  pValue <- phyper(q = overlapLegnth, 
                   m = length(TopicGenes),
                   n = TotalGenes,
                   k =length(ChipGenes), lower.tail = F)
  return(pValue)
}
for (i in 1:length(colnames(ModelList))) {
  print(colnames(ModelList)[i])
  PValues <-  mclapply(TFs$Genes, FUN = function(x) {TFPhyper(GeneListTopic = ModelList[,colnames(ModelList)[i]], Gene = x)}, mc.cores = 4)
  TFs$Pvalue.Tmp <-unlist(PValues)
  colnames(TFs)[ncol(TFs)] <- paste0(colnames(ModelList)[i], "PValue", sep = "_")
}
TFs <- as.data.frame(TFs)[,-1 *which(colnames(TFs) == "Genes")]
write.table(TFs, "TranscriptionFactors.txt",
            quote = F, sep = "\t", row.names = F, col.names = T)
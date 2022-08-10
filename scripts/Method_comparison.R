#Author: Alex Chitsazan
library(Seurat)
library(tidyverse)
library(lda)
library(RColorBrewer)
## Initialize the Values
seed <- 8
iterations <- 500
burnin <- 250
alpha <- 50
alphaByTopic <- TRUE
beta <- 0.1
returnType <- 'allModels'
addModels <- TRUE
set.seed(8)

###################################
###        MCF7             #######
###################################
## PEPE data
MCF7_E2 <- readRDS("E2_MCF7.rds") #load in MCF7 seurat object
DefaultAssay(MCF7_E2) <- "RNA"

## Cluster 1
Cluster1.Markers <- FindMarkers(MCF7_E2, ident.1 = 1)
GeneMarkers <- head(rownames(Cluster1.Markers), 50)
MCF7_E2<- AddModuleScore(MCF7_E2, features = list(GeneMarkers), name = "Cluster")
MCF7_E2$Cluster1 <- scale(MCF7_E2$Cluster1, center=T, scale = T)

## Titan

library(TITAN)
Model_MCF7 <- readRDS("")
Model_MCF7_E2_20T_CLR_5000Variable_M10 <- readRDS("Model_newMCF7_E2_20T_CLR_5000Variable_M10.rds")
MCF7_E2 <- TITAN::addTopicsToSeuratObject(Object = MCF7_E2, model = Model_MCF7_E2_20T_CLR_5000Variable_M10)

# topicClus
load("MCF_E2.rda") #load countClust object
MCF7.topicClus <- Topic_clus
MCF7.topicClus$D
omega <- MCF7.topicClus$omega
colnames(omega) <- paste0("countClustCluster", 1:20)
omegaToTitan <- as.data.frame(omega)
omegaToTitan <- t(scale(t(omegaToTitan), center=TRUE, scale = T))
MCF7_E2@meta.data <- cbind(MCF7_E2@meta.data, omegaToTitan)
HeatmapTopic2(Object = MCF7_E2,
              topics =  omegaToTitan,
              AnnoVector = MCF7_E2@meta.data$hash.ID,
              AnnoName = "Time Point", clusterTopics = T)
dev.off()
library(pheatmap)
library(CountClust)
MCF7.gene_names <- rownames(MCF7_E2)
theta_mat <- MCF7.topicClus$theta;
top_features <- ExtractTopFeatures(theta_mat, top_features=50,
                                   method="poisson", options="min");
gene_list <- t(do.call(rbind, lapply(1:dim(top_features$indices)[1],
                                     function(x) MCF7.gene_names[top_features$indices[x,]])))
colnames(gene_list) <- paste0("Topic_", 1:20)
write.table(gene_list, 
            "MCF7_geneList.txt",
            sep = "\t",
            col.names = T, 
            row.names = F)



# SCENIC
MCF7_Scenic <- read.csv("MCF7_SCENIC_ESR1_scores.tsv", sep = "\t") #read in scores from SCENIC
MCF7_ESR1_SCENIC_Norm <- scale(MCF7_Scenic$MCF7_ESR1_SCENIC, scale=T, center=T)
names(MCF7_ESR1_SCENIC_Norm) <- rownames(MCF7_Scenic)
MCF7_E2$Scenic <- MCF7_ESR1_SCENIC_Norm
library(TITAN)
library(tidyverse)
MCF7cols <- brewer.pal(9,"BuGn")
MCF7cols <- MCF7cols[3:8]
MCF7cols
MCF7_All <- MCF7_E2@meta.data %>%
  mutate(TITAN_ESR1 = Topic_8 - median(Topic_8[mapLabel == "Control"]),
         counClust_ESR1 = countClustCluster18 - median(countClustCluster18[mapLabel == "Control"]),
         Cluster1Top50DE = Cluster1 - median(Cluster1[mapLabel == "Control"]),
         SCENIC_ESR1 = Scenic - median(Scenic[mapLabel == "Control"], na.rm = T)) %>%
  select(hash.ID,TITAN_ESR1, counClust_ESR1, Cluster1Top50DE, SCENIC_ESR1) %>%
  gather("Tool", "ZScore", 2:5)
ggplot(MCF7_All, aes(x = Tool, y = ZScore, fill = hash.ID)) +
  geom_boxplot(position = "dodge") +
  theme_minimal() +
  ggtitle("MCF7") +
  scale_fill_manual(values = MCF7cols)


Test <- MCF7_E2@meta.data %>%
  mutate(TITAN_ESR1 = Topic_8 - median(Topic_8[mapLabel == "Control"]),
         counClust_ESR1 = countClustCluster14 - median(countClustCluster14[mapLabel == "Control"]),
         Cluster1Top50DE = Cluster1 - median(Cluster1[mapLabel == "Control"]),
         SCENIC_ESR1 = Scenic - median(Scenic[mapLabel == "Control"], na.rm = T)) %>%
  select(hash.ID,TITAN_ESR1, counClust_ESR1, Cluster1Top50DE, SCENIC_ESR1) %>%
  gather("Tool", "ZScore", 2:5)
ggplot(Test, aes(x = Tool, y = ZScore, fill = hash.ID)) +
  geom_boxplot(position = "dodge") +
  theme_minimal() +
  ggtitle("MCF7")


MCF7_E2.countClust <- MCF7_E2@meta.data %>% select(starts_with("countClust"))

MCF7_E2.TITAN <- MCF7_E2@meta.data %>% select(starts_with("Topic_"))
library(pheatmap)
pheatmap(cor(MCF7_E2.countClust, 
             MCF7_E2.TITAN, 
             method = "spearman"),
         cluster_cols = T, 
         cluster_rows = T, 
         display_numbers = T, fontsize_number = 10,
         main = "MCF7")


### MCF7 Viper
MCF7_VIPER <- read.delim("MCF7_VIPER_ESR1_scores.tsv", #load in Viper scores
                         header=F, 
                         sep = "\t", row.names = 1)
MCF7_VIPER_cellnames <- rownames(MCF7_VIPER)
MCF7cols <- brewer.pal(9,"BuGn")
MCF7cols <- MCF7cols[4:8]
MCF7_E2$VIPER_UnScaled <- MCF7_VIPER[,1,drop = F]
MCF7_VIPER <- scale(MCF7_VIPER[,1,drop = F], center=T, scale = T)
rownames(MCF7_VIPER) <- MCF7_VIPER_cellnames
MCF7_E2$VIPER <- MCF7_VIPER[,1]
ggplot(MCF7_E2@meta.data %>% filter(hash.ID != "M7-CTL"), aes(x = hash.ID, y = VIPER, fill = hash.ID)) +
  geom_boxplot(position = "dodge") +
  theme_minimal() +
  ggtitle("MCF7") +
  ylim(-3, 5) +
  scale_fill_manual(values = MCF7cols)
MCF7_All <- MCF7_E2@meta.data %>%
  mutate(TITAN_ESR1 = Topic_8 - median(Topic_8[mapLabel == "Control"]),
         counClust_ESR1 = countClustCluster18 - median(countClustCluster18[mapLabel == "Control"]),
         Cluster1Top50DE = Cluster1 - median(Cluster1[mapLabel == "Control"]),
         SCENIC_ESR1 = Scenic - median(Scenic[mapLabel == "Control"], na.rm = T),
         VIPER_ESR1 = VIPER) %>%
  select(hash.ID,TITAN_ESR1, counClust_ESR1, Cluster1Top50DE, SCENIC_ESR1, VIPER_ESR1) %>%
  gather("Tool", "ZScore", 2:6)

MCF7cols <- brewer.pal(9,"BuGn")
MCF7cols <- MCF7cols[3:8]
ggplot(MCF7_All, aes(x = Tool, y = ZScore, fill = hash.ID)) +
  geom_boxplot(position = "dodge") +
  theme_minimal() +
  ggtitle("MCF7") +
  scale_fill_manual(values = MCF7cols)
ggplot(MCF7_All, aes(x = Tool, y = ZScore, fill = hash.ID)) +
  geom_boxplot(position = "dodge") +
  theme_minimal() +
  ggtitle("MCF7") 


MCF7_All <- MCF7_E2@meta.data %>%
  mutate(TITAN_ESR1 = Topic_8 - median(Topic_8[mapLabel == "Control"]),
         counClust_ESR1 = countClustCluster18 - median(countClustCluster18[mapLabel == "Control"]),
         Cluster1Top50DE = Cluster1 - median(Cluster1[mapLabel == "Control"]),
         SCENIC_ESR1 = Scenic - median(Scenic[mapLabel == "Control"], na.rm = T),
         VIPER_ESR1 = VIPER_UnScaled) %>%
  select(hash.ID,TITAN_ESR1, counClust_ESR1, Cluster1Top50DE, SCENIC_ESR1, VIPER_ESR1) %>%
  gather("Tool", "ZScore", 2:6)
ggplot(MCF7_All, aes(x = Tool, y = ZScore, fill = hash.ID)) +
  geom_boxplot(position = "dodge") +
  theme_minimal() +
  ggtitle("MCF7") +
  scale_fill_manual(values = MCF7cols)
ggplot(MCF7_All, aes(x = Tool, y = ZScore, fill = hash.ID)) +
  geom_boxplot(position = "dodge") +
  theme_minimal() +
  ggtitle("MCF7") 
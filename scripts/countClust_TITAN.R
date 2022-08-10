#Author: Alex Chitsazan
## COUNT CLUST
library(devtools)
install_github("TaddyLab/maptpx")
install_github('kkdey/CountClust')
library(CountClust)
library(Seurat)
#### MCF7 is used in this example, but the process is the same for any Seurat Object
#### MCF7_E2 is the Seurat object containing the MCF7 E2 treatment series data
MCF7_E2.counts <- GetAssayData(MCF7_E2, assay = "RNA", slot = "counts")
# MCF7_E2.counts <- as.matrix(MCF7_E2.counts)
MCF7_E2.metaData <- MCF7_E2@meta.data
MCF7_E2.gene_names <- rownames(MCF7_E2)

FitGoM(t(as.matrix(MCF7_E2.counts)),
       K=20, tol=0.1,
       path_rda="MCF_E2.rda")
load("MCF_E2.rda")
MCF7_E2.topicClus <- Topic_clus
MCF7_E2.topicClus$D
omega <- MCF7_E2.topicClus$omega
colnames(omega) <- c(1:NCOL(omega))
clusterlabels  <- MCF7_E2.metaData$mapLabel
library(CountClust)
annotation <- data.frame(
  sample_id = paste0("X", 1:length(clusterlabels)),
  tissue_label =  clusterlabels)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
StructureGGplot(omega = omega, ### CHANGE COLORS SO IT"S BRIGHT AND GREY
                annotation= annotation, 
                palette = gg_color_hue(20),
                yaxis_label = "",
                order_sample = TRUE,
                split_line = list(split_lwd = .4,
                                  split_col = "white"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))
clusterlabels  <- MCF7_E2.metaData$integrated_snn_res.0.3
library(CountClust)
annotation <- data.frame(
  sample_id = paste0("X", 1:length(clusterlabels)),
  tissue_label =  clusterlabels)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
StructureGGplot(omega = omega,
                annotation= annotation, 
                palette = gg_color_hue(20),
                yaxis_label = "",
                order_sample = TRUE,
                split_line = list(split_lwd = .4,
                                  split_col = "white"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))



colnames(omega) <- paste0("countClustCluster", 1:20)
omegaToTitan <- as.data.frame(omega)
MCF7_E2@meta.data <- cbind(MCF7_E2@meta.data, omegaToTitan)
test <- t(scale(t(omegaToTitan), center=TRUE, scale = T))
colnames(test) <-paste0("scaled", colnames(test))
MCF7_E2@meta.data <- cbind(MCF7_E2@meta.data, test)
dev.off()
addTopicsToSeuratObject()
VlnPlot(MCF7_E2, "countClustCluster18", group.by = "mapLabel", pt.size = 0.1) +
  geom_boxplot(size=0.1) +
  ggtitle("MCF7 E2 New ESR1 countClust Cluster18")
VlnPlot(MCF7_E2, "scaledcountClustCluster18", group.by = "mapLabel", pt.size = 0.1) +
  geom_boxplot(size=0.1) +
  ggtitle("MCF7 E2 New ESR1 countCkyst Cluster18")
MCF7_E2@meta.data %>%
  group_by(mapLabel) %>%
  summarize(NormTopicScore = mean(countClustCluster18)) %>%
  mutate(TITANFC = NormTopicScore / NormTopicScore[1]) 
HeatmapTopic
MCF7_E2.topicClus$X
simple_triplet_matrix_sparse <-  sparseMatrix(i=MCF7_E2.topicClus$X$i, 
                                              j=MCF7_E2.topicClus$X$j, 
                                              x=MCF7_E2.topicClus$X$v,
                                              dims=c(MCF7_E2.topicClus$X$nrow, 
                                                     MCF7_E2.topicClus$X$ncol))
simple_triplet_matrix_sparse[1:6,1:50]

png(paste("MCF7_E2", "Heatmap_TimePoints_UnNormalized.png", sep = "_"), height = 8.5, width = 11, units = "in", res = 500)
TITAN::HeatmapTopic(Object = MCF7_E2,
                    topics =  omegaToTitan,
                    AnnoVector = MCF7_E2@meta.data$mapLabel,
                    AnnoName = "Time Point", clusterTopics = T)
dev.off()

library(Seurat)
library(TITAN)
library(tidyverse)

swarbrick <- ReadMtx("swarbrick/matrix.mtx.gz", "swarbrick/barcodes.tsv.gz", "swarbrick/features.tsv.gz")
SB_SO <- CreateSeuratObject(counts=swarbrick, min.cells = 5, min.features = 200)

SB_SO[["percent.mt"]] <- PercentageFeatureSet(object = SB_SO, pattern = "^MT-")

epi_cells <- read.table("swarbrick/subset_object_Cancer_epithelial_subcluster.tsv", header = T, sep = "\t")

epi_cells <- epi_cells[-1,]
rownames(epi_cells) <- epi_cells$NAME
epi_cells <- epi_cells[,2:3]
colnames(epi_cells) <- c("paperDim_1", "paperDim_2")
epi_cells$paperDim_1 <- as.numeric(epi_cells$paperDim_1)
epi_cells$paperDim_2 <- as.numeric(epi_cells$paperDim_2)

epi_SB <- subset(SB_SO, cells = rownames(epi_cells))

epi_SB[["paperDim"]] <- CreateDimReducObject(
  embeddings = as.matrix(epi_cells),
  key = "paperDim_",
  assay = "RNA",
  global = TRUE
)

pdf("swarbrick/swarbrick_origin_UMAP.pdf")
DimPlot(epi_SB, reduction = "paperDim", group.by = "orig.ident")
dev.off()

T47D_model <- readRDS("/home/groups/CEDAR/doe/projects/LDA/Model_PEPE_T47D_20T_CLR_5000Variable_M10.rds")

SO <- ImputeAndAddTopics(epi_SB, T47D_model)

SO$origin <- SO$orig.ident

pdf("swarbrick/swarbrick_T47DimputedTopics_heatmap.pdf")
HeatmapSortByTopic(SO, SO@meta.data[,5:24], sortByTopic="ImputedTopic_11", SO$origin, "origin")
dev.off()

pdf("swarbrick/swarbrick_T47DimputedTopics_heatmap_bySample.pdf")
HeatmapTopic(SO, SO@meta.data[,5:24], SO$origin, "origin", clusterTopics = T)
dev.off()

Idents(SO) <- "origin"

pdf("swarbrick/swarbrick_T47DimputedTopic11_boxplot.pdf")
VlnPlot(SO, features = "ImputedTopic_11") + geom_boxplot()
dev.off()

pdf("swarbrick/swarbrick_T47DimputedTopic17_boxplot.pdf")
VlnPlot(SO, features = "ImputedTopic_17") + geom_boxplot()
dev.off()

pdf("swarbrick/swarbrick_T47DimputedTopic_RidgePlot.pdf")
RidgePlot(SO, features = c("ImputedTopic_11", "ImputedTopic_17"))
dev.off()

pdf("swarbrick/swarbrick_T47DimputedTopic11_FeaturePlot.pdf")
FeaturePlot(SO, features = "ImputedTopic_11", reduction = "paperDim")
dev.off()

pdf("swarbrick/swarbrick_T47DimputedTopic17_FeaturePlot.pdf")
FeaturePlot(SO, features = "ImputedTopic_17", reduction = "paperDim")
dev.off()

saveRDS(SO, "swarbrick/swarbrick_epi_wImputedTopics_SO.rds")
SO <- readRDS("swarbrick/swarbrick_epi_wImputedTopics_SO.rds")

pdf("swarbrick/swarbrick_T47DimputedTopic11vs17_DimPlot.pdf")
DimPlot(SO, reduction = "imputedLDA", dims = c(11,17))
dev.off()



############ swarbrick by subtype #############################

SO$subType <- "ER+"
SO$subType[which(SO$origin %in% c("CID44041", "CID4465", "CID4495", "CID44971", "CID44991", "CID4513", "CID4515", "CID4523", "CID3946", "CID3963"))] <- "TNBC"
SO$subType[which(SO$origin %in% c("CID3586", "CID3921", "CID45171", "CID3838", "CID4066"))] <- "HER2+"

pdf("swarbrick/swarbrick_T47DimputedTopics_heatmap_bySubType_sorted.pdf")
HeatmapSortByTopic(SO, SO@meta.data[,5:24], sortByTopic="ImputedTopic_11", SO$subType, "subType")
dev.off()

pdf("swarbrick/swarbrick_T47DimputedTopics_heatmap_bySubType.pdf")
HeatmapTopic(SO, SO@meta.data[,5:24], SO$subType, "subType", clusterTopics = T)
dev.off()


saveRDS(SO, "swarbrick/swarbrick_epi_wImputedTopics_SO.rds")
SO <- readRDS("swarbrick/swarbrick_epi_wImputedTopics_SO.rds")

pdf("swarbrick/swarbrick_T47DimputedTopic11vs17_DimPlot_bySubType.pdf")
DimPlot(SO, reduction = "imputedLDA", dims = c(11,17))
dev.off()


SO <- readRDS("swarbrick_epi_wImputedTopics_SO.rds")

runLDA(SO, ntopics = seq(10,100, by=10), parallel = T, outDir = "swarbrick_Models", cores = 10)

pdf("swarbrick_Models/TITAN_elbowPlot.pdf")
LDAelbowPlot("swarbrick_Models/", SO)
dev.off()

swarbrick_model <- readRDS("swarbrick_Models/Model_40topics.rds")

SO <- addTopicsToSeuratObject(swarbrick_model, SO)

pdf("swarbrick_TopicHeatmap.pdf")
HeatmapTopic(SO, topics = SO@meta.data[,26:65], SO$origin, "origin", clusterTopics = T)
dev.off()

pdf("swarbrick_both36ESR1and39FOXM1_FeaturePlot.pdf", width=15)
FeaturePlot(SO, features = c("Topic_36", "Topic_39"), reduction = "paperDim", blend = T)
dev.off()


pdf("swarbrick_ownTopics_RidgePlot36v39.pdf")
RidgePlot(SO, features = c("Topic_36", "Topic_39"))
dev.off()

SO <- RunUMAP(SO, reduction = "lda", dims = 1:40)

pdf("swarbrick_TITANumap.pdf")
DimPlot(SO, reduction = "umap", group.by = "origin")
dev.off()

saveRDS(SO, "swarbrick_wAllTopics.rds")

########################## by subtype ################################

SO <- readRDS("swarbrick_wAllTopics.rds")

SO$subType <- "ER+"
SO$subType[which(SO$origin %in% c("CID44041", "CID4465", "CID4495", "CID44971", "CID44991", "CID4513", "CID4515", "CID4523", "CID3946", "CID3963"))] <- "TNBC"
SO$subType[which(SO$origin %in% c("CID3586", "CID3921", "CID45171", "CID3838", "CID4066"))] <- "HER2+"

pdf("swarbrick_TopicHeatmap_bySubType.pdf")
HeatmapTopic(SO, topics = SO@meta.data[,26:65], SO$subType, "subType", clusterTopics = T)
dev.off()


p1 <- pheatmap(SO@meta.data[,26:65][order(SO@meta.data[,c("subType")], SO@meta.data[,"origin"]),],
               hclustfun = function(x) hclust(x, method="ward.D2"),
               scale = "row",
               cluster_cols = T,
               cluster_rows = F,show_rownames = F,
               col=colorRampPalette(rev(brewer.pal(11, "RdBu"))[c(1:4,8:11)])(256),
               annotation_row = SO@meta.data[,c("subType", "origin")],
               annotation_names_row = T,
               #annotation_colors = anno_colors,
               cex=1)

pdf("swarbrick_TopicHeatmap_bySubTypeANDorigin.pdf")
p1
dev.off()

Idents(SO) <- "subType"

pdf("swarbrick_FOXM1topic39_FeaturePlot.pdf")
FeaturePlot(SO, features = "Topic_39", reduction = "paperDim")
dev.off()

pdf("swarbrick_ESR1topic36_FeaturePlot.pdf")
FeaturePlot(SO, features = "Topic_36", reduction = "paperDim")
dev.off()


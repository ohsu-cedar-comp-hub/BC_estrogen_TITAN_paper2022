library(Seurat)
library(tidyverse)
library(TITAN)

plot2Topics <- function(SeuratObj, model, topics, GOI, plot_split, IntIdents, cols = NULL) {
  SO <- addTopicsToSeuratObject(Object = SeuratObj, model = model)
  if (length(IntIdents) > 0) {
    keeptime <- (unname(SO$ADT_maxID) %in% IntIdents)
    keepclus <- (unname(SO$seurat_clusters) %in% IntIdents)
    keepall <- (keeptime | keepclus)
    SO <- SO[,keepall]
  }
  
  SO$seurat_clusters <- factor(SO$seurat_clusters)
  SO$ADT_maxID <- factor(SO$ADT_maxID)
  if (plot_split == "Timepoint") {
    if (GOI[1] %in% rownames(SO)){
      if (length(GOI) == 1) {
        fp <- FeaturePlot(SO, features = GOI[1], dims = topics, reduction = "lda", split.by = "ADT_maxID", min.cutoff = 'q1', combine = F, cols = brewer.pal(n = 2, name = cols))
        print(CombinePlots(fp, ncol = 2))
      }
      if (length(GOI) == 2) {
        fp <- FeaturePlot(SO, features = GOI, dims = topics, reduction = "lda", split.by = "ADT_maxID", min.cutoff = 'q1', combine = F, blend = T, cols = brewer.pal(n = 3, name = cols))
        print(CombinePlots(fp, ncol = 2))
      }
    }
    else {
      if (GOI[1] == "CC"){
        fp <- DimPlot(SO, reduction = "lda", dims = topics, group.by = "Phase", split.by = "ADT_maxID", ncol = 2, cols = cols)
        print(fp)
      } 
      if (GOI[1] == "Clone") {
        fp <- DimPlot(SO, reduction = "lda", dims = topics, group.by = "Clone", split.by = "ADT_maxID", ncol = 2, cols = cols)
        print(fp)
      }
      if (GOI[1] == "Time") {
        fp <- DimPlot(SO, reduction = "lda", dims = topics, group.by = "ADT_maxID", split.by = "ADT_maxID", ncol = 6, cols = cols) + theme(text = element_text(size = 5), axis.text = element_text(size = 5), axis.line = element_line(colour = "black", size = 0.1), axis.ticks = element_line(colour = "black", size = 0.1)) +
          guides(shape = guide_legend(override.aes = list(size = 1)),
                 color = guide_legend(override.aes = list(size = 1))) +
          theme(legend.title = element_text(size = 5), 
                legend.text  = element_text(size = 5),
                legend.key.size = unit(0.1, "lines"))
        print(fp)
      }
      if (GOI[1] != "CC" & GOI[1] != "Clone" & GOI[1] != "Time") {
        fp <- DimPlot(SO, reduction = "lda", dims = topics, split.by = "ADT_maxID", ncol = 2, cols = cols)
        print(fp)
      }
    }
  } 
  if (plot_split == "Cluster") {
    if (GOI[1] %in% rownames(SO)){
      if (length(GOI) == 1) {
        fp <- FeaturePlot(SO, features = GOI[1], dims = topics, reduction = "lda", split.by = "seurat_clusters", min.cutoff = 'q1', combine = F, cols = brewer.pal(n = 2, name = cols))
        print(CombinePlots(fp, ncol = 2))
      }
      if (length(GOI) == 2) {
        fp <- FeaturePlot(SO, features = GOI, dims = topics, reduction = "lda", split.by = "seurat_clusters", min.cutoff = 'q1', combine = F, blend = T, cols = brewer.pal(n = 3, name = cols))
        print(CombinePlots(fp, ncol = 2))
      }
    }
    else {
      if (GOI[1] == "CC"){
        fp <- DimPlot(SO, reduction = "lda", dims = topics, group.by = "Phase", split.by = "seurat_clusters", ncol = 2, cols = cols)
        print(fp)
      }
      if (GOI[1] == "Clone") {
        fp <- DimPlot(SO, reduction = "lda", dims = topics, group.by = "Clone", split.by = "seurat_clusters", ncol = 2, cols = cols)
        print(fp)
      }
      if (GOI[1] == "Time") {
        fp <- DimPlot(SO, reduction = "lda", dims = topics, group.by = "ADT_maxID", split.by = "seurat_clusters", ncol = 2, cols = cols)
        print(fp)
      }
      if (GOI[1] != "CC" & GOI[1] != "Clone" & GOI[1] != "Time") {
        fp <- DimPlot(SO, reduction = "lda", dims = topics, split.by = "seurat_clusters", ncol = 2, cols = cols)
        print(fp)
      }
    }
  } 
  if (plot_split != "Timepoint" & plot_split != "Cluster") {
    if (GOI[1] %in% rownames(SO)){
      if (length(GOI) == 1) {
        fp <- FeaturePlot(SO, features = GOI[1], dims = topics, reduction = "lda", min.cutoff = 'q1', cols = brewer.pal(n = 2, name = cols))
        print(fp)
      }
      if (length(GOI) == 2) {
        fp <- FeaturePlot(SO, features = GOI, dims = topics, reduction = "lda", min.cutoff = 'q1', blend = T, combine = F, cols = brewer.pal(n = 3, name = cols))
        print(CombinePlots(fp, ncol = 2))
      }
    } else {
      if (GOI[1] == "CC"){
        fp <- DimPlot(SO, reduction = "lda", dims = topics, group.by = "Phase", ncol = 2, cols = cols)
        print(fp)
      }
      if (GOI[1] == "Clone") {
        fp <- DimPlot(SO, reduction = "lda", dims = topics, group.by = "Clone", ncol = 2, cols = cols)
        print(fp)
      }
      if (GOI[1] == "Time") {
        fp <- DimPlot(SO, reduction = "lda", dims = topics, group.by = "ADT_maxID", ncol = 2, cols = cols)
        print(fp)
      }
      if (GOI[1] != "CC" & GOI[1] != "Clone" & GOI[1] != "Time") {
        fp <- DimPlot(SO, reduction = "lda", dims = topics, ncol = 2, cols = cols)
        print(fp)
      }
    }
  }
}

Top_VlnPlot <- function(Object, model, topic, cols) {
  fp <- VlnPlot(addTopicsToSeuratObject(Object = Object, model = model), 
                features = paste("Topic", topic, sep = "_"),
                cols = cols, 
                group.by = "ADT_maxID", pt.size = 0) + geom_boxplot(lwd = 0.1, outlier.size = 0) + theme(text = element_text(size = 5), axis.text = element_text(size = 5), axis.line = element_line(colour = "black", size = 0.1), axis.ticks = element_line(colour = "black", size = 0.1)) +
    guides(shape = guide_legend(override.aes = list(size = 1)),
           color = guide_legend(override.aes = list(size = 1))) +
    theme(legend.title = element_text(size = 5), 
          legend.text  = element_text(size = 5),
          legend.key.size = unit(0.1, "lines"),
          axis.title.x = element_blank())
  fp$layers[[1]]$aes_params$size = 0.1
  print(fp)
}

Seurat_object <- readRDS("ZR_ER_SO.rds")
LDA_model <- readRDS("Model_ZR_20T_CLR_5000Variable_M10.rds")

Seurat_object <- addTopicsToSeuratObject(LDA_model, Seurat_object)

ZRcols <- brewer.pal(9,"PuRd")
ZRcols <- ZRcols[3:8]

pdf("ZR_FeatureTopic720.pdf", height = 1.5, width = 6)
plot2Topics(Seurat_object, LDA_model, c(7,20), "Time", "Timepoint", IntIdents = NULL, cols = ZRcols)
dev.off()

pdf("ZR_Top7.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 7, ZRcols)
dev.off()

pdf("ZR_Top20.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 20, ZRcols)
dev.off()

pdf("ZR_FeatureTopic45.pdf", height = 1.5, width = 6)
plot2Topics(Seurat_object, LDA_model, c(4,5), "Time", "Timepoint", IntIdents = NULL, cols = ZRcols)
dev.off()

pdf("ZR_Top4.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 4, ZRcols)
dev.off()

pdf("ZR_Top5.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 5, ZRcols)
dev.off()

pdf("ZR_FeatureTopic1710.pdf", height = 1.5, width = 6)
plot2Topics(Seurat_object, LDA_model, c(17,10), "Time", "Timepoint", IntIdents = NULL, cols = ZRcols)
dev.off()

pdf("ZR_Top17.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 17, ZRcols)
dev.off()

pdf("ZR_Top10.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 10, ZRcols)
dev.off()

pdf("ZR_FeatureTopic13.pdf", height = 1.5, width = 6)
plot2Topics(Seurat_object, LDA_model, c(1,3), "Time", "Timepoint", IntIdents = NULL, cols = ZRcols)
dev.off()

pdf("ZR_Top1.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 1, ZRcols)
dev.off()

pdf("ZR_Top3.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 3, ZRcols)
dev.off()

pdf("ZR_Topic_8_13_FeaturePlot.pdf")
FeaturePlot(Seurat_object, reduction = "umap", features = c("Topic_8", "Topic_13"), blend = T, cols = c("orange","purple"))[[3]]
dev.off()

# Figure 2G RidgePlot

subSeurat_object <- subset(SO, hash.ID %in% c("ZR-CTL", "ZR-E-48"))

ggplot(data=subSeurat_object@meta.data, aes(x=Topic_of_Choice, fill = hash.ID)) + geom_density(alpha = 0.5) + 
  scale_fill_manual(values=brewer.pal(n=2, name = 'Pastel1')[c(1,2)]) + theme_bw()

Seurat_object <- readRDS("newE2_MCF7.rds")
LDA_model <- readRDS("Model_newMCF7_E2_20T_CLR_5000Variable_M10.rds")

Seurat_object <- addTopicsToSeuratObject(LDA_model, Seurat_object)

MCF7cols <- brewer.pal(9,"BuGn")
MCF7cols <- MCF7cols[3:8]

pdf("MCF7_FeatureTopic23.pdf", height = 1.5, width = 6)
plot2Topics(Seurat_object, LDA_model, c(2,3), "Time", "Timepoint", IntIdents = NULL, cols = MCF7cols)
dev.off()

pdf("MCF7_FeatureTopic28.pdf", height = 1.5, width = 6)
plot2Topics(Seurat_object, LDA_model, c(2,8), "Time", "Timepoint", IntIdents = NULL, cols = MCF7cols)
dev.off()

pdf("MCF7_HMMR_TFF3_FeaturePlot.pdf")
FeaturePlot(Seurat_object, features = c("HMMR", "TFF3"), reduction = "lda", dims = c(2,8), blend = T, cols = c("lightgrey", "blue", "red"), combine = F, pt.size = .9)[[3]]
dev.off()

pdf("MCF7_HMMR_PGR_FeaturePlot.pdf")
FeaturePlot(Seurat_object, features = c("HMMR", "PGR"), reduction = "lda", dims = c(2,8), blend = T, cols = c("lightgrey", "blue", "red"), combine = F, pt.size = .9)[[3]]
dev.off()

pdf("MCF7_Top2.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 2, MCF7cols)
dev.off()

pdf("MCF7_Top3.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 3, MCF7cols)
dev.off()

pdf("MCF7_FeatureTopic98.pdf", height = 1.5, width = 6)
plot2Topics(Seurat_object, LDA_model, c(9,8), "Time", "Timepoint", IntIdents = NULL, cols = MCF7cols)
dev.off()

pdf("MCF7_Top9.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 9, MCF7cols)
dev.off()

pdf("MCF7_Top8_pubColors.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 8, MCF7cols)
dev.off()

pdf("MCF7_Top2_pubColors.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 2, MCF7cols)
dev.off()

pdf("MCF7_FeatureTopic67.pdf", height = 1.5, width = 6)
plot2Topics(Seurat_object, LDA_model, c(6,7), "Time", "Timepoint", IntIdents = NULL, cols = MCF7cols)
dev.off()

pdf("MCF7_FeatureTopic28.pdf", height = 1.5, width = 6)
plot2Topics(Seurat_object, LDA_model, c(2,8), "Time", "Timepoint", IntIdents = NULL, cols = MCF7cols)
dev.off()

pdf("MCF7_Top6.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 6, MCF7cols)
dev.off()

pdf("MCF7_Top7.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 7, MCF7cols)
dev.off()

pdf("MCF7_FeatureTopic1112.pdf", height = 1.5, width = 6)
plot2Topics(Seurat_object, LDA_model, c(11,12), "Time", "Timepoint", IntIdents = NULL, cols = MCF7cols)
dev.off()

pdf("MCF7_Top11.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 11, MCF7cols)
dev.off()

pdf("MCF7_Top12.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 12, MCF7cols)
dev.off()

pdf("MCF7_Topic_8_2_FeaturePlot.pdf")
FeaturePlot(Seurat_object, reduction = "umap", features = c("Topic_8", "Topic_2"), blend = T, cols = c("orange","purple"))[[3]]
dev.off()

# Figure 2G RidgePlot

subSeurat_object <- subset(SO, hash.ID %in% c("M7-CTL", "M7-E-48"))

ggplot(data=subSeurat_object@meta.data, aes(x=Topic_of_Choice, fill = hash.ID)) + geom_density(alpha = 0.5) + 
  scale_fill_manual(values=brewer.pal(n=2, name = 'Pastel1')[c(1,2)]) + theme_bw()


Seurat_object <- readRDS("E2_T47D.rds")
LDA_model <- readRDS("Model_T47D_20T_CLR_5000Variable_M10.rds")

Seurat_object <- addTopicsToSeuratObject(LDA_model, Seurat_object)

T47Dcols <- brewer.pal(9,"Greens")
T47Dcols <- T47Dcols[3:8]

pdf("T47D_FeatureTopic420_greens.pdf", height = 1.5, width = 6)
plot2Topics(Seurat_object, LDA_model, c(4,20), "Time", "Timepoint", IntIdents = NULL, cols = T47Dcols)
dev.off()

pdf("T47D_Top4.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 4, T47Dcols)
dev.off()

pdf("T47D_Top20.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 20, T47Dcols)
dev.off()

pdf("T47D_FeatureTopic1614.pdf", height = 1.5, width = 6)
plot2Topics(Seurat_object, LDA_model, c(16,14), "Time", "Timepoint", IntIdents = NULL, cols = T47Dcols)
dev.off()

pdf("T47D_Top16.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 16, T47Dcols)
dev.off()

pdf("T47D_Top14.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 14, T47Dcols)
dev.off()

pdf("T47D_FeatureTopic817.pdf", height = 1.5, width = 6)
plot2Topics(Seurat_object, LDA_model, c(8,17), "Time", "Timepoint", IntIdents = NULL, cols = T47Dcols)
dev.off()

pdf("T47D_Top8.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 8, T47Dcols)
dev.off()

pdf("T47D_Top17.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 17, T47Dcols)
dev.off()

pdf("T47D_FeatureTopic157.pdf", height = 1.5, width = 6)
plot2Topics(Seurat_object, LDA_model, c(15,7), "Time", "Timepoint", IntIdents = NULL, cols = T47Dcols)
dev.off()

pdf("T47D_Top15.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 15, T47Dcols)
dev.off()

pdf("T47D_Top7.pdf", height = 1.5, width = 2)
Top_VlnPlot(Seurat_object, LDA_model, 7, T47Dcols)
dev.off()

pdf("T47D_Topic_20_4_FeaturePlot.pdf")
FeaturePlot(Seurat_object, reduction = "umap", features = c("Topic_20", "Topic_4"), blend = T, cols = c("orange","purple"))[[3]]
dev.off()

# Figure 2G RidgePlot

subSeurat_object <- subset(SO, hash.ID %in% c("TD-CTL", "TD-E-48"))

ggplot(data=subSeurat_object@meta.data, aes(x=Topic_of_Choice, fill = hash.ID)) + geom_density(alpha = 0.5) + 
  scale_fill_manual(values=brewer.pal(n=2, name = 'Pastel1')[c(1,2)]) + theme_bw()


clone_SO <- readRDS("clones_MCF7.rds")
MCF7_E2model <- readRDS("Model_newMCF7_E2_20T_CLR_5000Variable_M10.rds")

clone_SO <- ImputeAndAddTopics(clone_SO, MCF7_E2model)

clone_SO$hash.ID <- factor(x = clone_SO$hash.ID, levels = c('M7-C2-C', 'M7-C2-E24', 'M7-C2-E48', 'M7-C4-C', 'M7-C4-E24', 'M7-C4-E48', 'M7-CTL', 'M7-E-24', 'M7-E-48'))

Idents(clone_SO) <- "hash.ID"

pdf("Clone_Impute_Top2.pdf")
VlnPlot(clone_SO, features = "Topic_2", ncol = 2, pt.size = 0, cols = c('#66B2FF','#0080FF','#004C99','#FFB266','#FF8000','#994C00','#66FF66','#00CC00','#006600')) + geom_boxplot()
dev.off()

pdf("Clone_Impute_Top8.pdf")
VlnPlot(clone_SO, features = "Topic_8", ncol = 2, pt.size = 0, cols = c('#66B2FF','#0080FF','#004C99','#FFB266','#FF8000','#994C00','#66FF66','#00CC00','#006600')) + geom_boxplot()
dev.off()


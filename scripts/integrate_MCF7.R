library(Seurat)
library(tidyverse)

SO1 <- readRDS("MCF7SO_A.rds")
SO2 <- readRDS("MCF7SO_B.rds")
SO3 <- readRDS("MCF7SO_C.rds")

SO.anchors <- FindIntegrationAnchors(object.list = c(SO1,SO2,SO3), dims = 1:40)
SO.integrated <- IntegrateData(anchorset = SO.anchors, dims = 1:40)

cc_genes  <- read.table("<location_of_CCgenes>/regev_lab_cell_cycle_genes.txt")
s.genes   <- cc_genes$V1[1:43]
g2m.genes <- cc_genes$V1[44:97]

SO.integrated <- CellCycleScoring(SO.integrated, s.features = s.genes, g2m.features = g2m.genes)

SO.integrated <- ScaleData(SO.integrated, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(SO.integrated))

SO.integrated <- RunPCA(SO.integrated, npcs = 40, verbose = FALSE)
SO.integrated <- FindNeighbors(object=SO.integrated, dims=1:40)
SO.integrated <- FindClusters(object=SO.integrated, resolution = 0.3)
SO.integrated <- RunUMAP(SO.integrated, reduction = "pca", dims = 1:40)

### BELOW IS WHERE THE PROPER LABELS ARE ASSIGNED (there were barcoding issues that led to labels being swapped)

SO.integrated@meta.data$mapLabel <- "HERE"
for (i in seq(1,length(SO.integrated$ADT_classification))) {
  if (SO.integrated$ADT_classification[i] == "M7-CTL") {
    SO.integrated@meta.data$mapLabel[i] <- "Control"
  }
  if (SO.integrated$ADT_classification[i] == "M7-E-03") {
    SO.integrated@meta.data$mapLabel[i] <- "E-3"
  }
  if (SO.integrated$ADT_classification[i] == "M7-E-06") {
    SO.integrated@meta.data$mapLabel[i] <- "E-6"
  }
  if (SO.integrated$ADT_classification[i] == "M7-E-24") {
    SO.integrated@meta.data$mapLabel[i] <- "E-24"
  }
  if (SO.integrated$ADT_classification[i] == "M7-E-48") {
    SO.integrated@meta.data$mapLabel[i] <- "E-48"
  }
  if (SO.integrated$ADT_classification[i] == "M7-E-72") {
    SO.integrated@meta.data$mapLabel[i] <- "E-72"
  }
  if (SO.integrated$ADT_classification[i] == "M7-ED-03") {
    SO.integrated@meta.data$mapLabel[i] <- "TE-3"
  }
  if (SO.integrated$ADT_classification[i] == "M7-ED-48") {
    SO.integrated@meta.data$mapLabel[i] <- "TE-48"
  }
  if (SO.integrated$ADT_classification[i] == "M7-EP-03") {
    SO.integrated@meta.data$mapLabel[i] <- "EP-3"
  }
  if (SO.integrated$ADT_classification[i] == "M7-EP-48") {
    SO.integrated@meta.data$mapLabel[i] <- "EP-48"
  }
  if (SO.integrated$ADT_classification[i] == "M7-D-03") {
    SO.integrated@meta.data$mapLabel[i] <- "T-3"
  }
  if (SO.integrated$ADT_classification[i] == "M7-D-48") {
    SO.integrated@meta.data$mapLabel[i] <- "T-48"
  }
  if (SO.integrated$ADT_classification[i] == "M7-C2-C") {
    SO.integrated@meta.data$mapLabel[i] <- "Clone2-Control"
  }
  if (SO.integrated$ADT_classification[i] == "M7-C2-E24") {
    SO.integrated@meta.data$mapLabel[i] <- "Clone2-E24"
  }
  if (SO.integrated$ADT_classification[i] == "M7-C2-E48") {
    SO.integrated@meta.data$mapLabel[i] <- "Clone2-E48"
  }
  if (SO.integrated$ADT_classification[i] == "M7-C4-C") {
    SO.integrated@meta.data$mapLabel[i] <- "Clone4-Control"
  }
  if (SO.integrated$ADT_classification[i] == "M7-C4-E24") {
    SO.integrated@meta.data$mapLabel[i] <- "Clone4-E24"
  }
  if (SO.integrated$ADT_classification[i] == "M7-C4-E48") {
    SO.integrated@meta.data$mapLabel[i] <- "Clone4-E48"
  }
  if (SO.integrated$ADT_classification[i] == "M7-P-03") {
    SO.integrated@meta.data$mapLabel[i] <- "P-3"
  }
  if (SO.integrated$ADT_classification[i] == "M7-P-48") {
    SO.integrated@meta.data$mapLabel[i] <- "P-48"
  }
}

DimPlot(SO.integrated, reduction = "umap", group.by = "seurat_clusters")
DimPlot(SO.integrated, reduction = "umap", group.by = "mapLabel")
DimPlot(SO.integrated, reduction = "umap", group.by = "Phase")


PEPE_MCF7 <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "EP-3" | mapLabel == "EP-48" | mapLabel == "E-3" | mapLabel == "E-48" | mapLabel == "P-3" | mapLabel == "P-48")
E2_MCF7 <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "E-3" | mapLabel == "E-6" | mapLabel == "E-24" | mapLabel == "E-48" | mapLabel == "E-72")
E2PG_MCF7 <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "E-3" | mapLabel == "E-48" | mapLabel == "P-3" | mapLabel == "P-48")
clones_MCF7 <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "E-24" | mapLabel == "E-48" | mapLabel == "Clone2-Control" | mapLabel == "Clone2-E24" | mapLabel == "Clone2-E48" | mapLabel == "Clone4-Control" | mapLabel == "Clone4-E24" | mapLabel == "Clone4-E48")


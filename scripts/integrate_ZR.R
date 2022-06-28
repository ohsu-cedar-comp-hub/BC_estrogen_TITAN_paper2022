library(Seurat)
library(tidyverse)

SOZR1 <- readRDS("ZR_ER_SO.rds")
SOZR2 <- readRDS("ZR_2.rds")
SOZR3 <- readRDS("ZR_3.rds")

ZR.anchors <- FindIntegrationAnchors(object.list = c(SOZR1, SOZR2, SOZR3), dims = 1:40)
SO.integrated <- IntegrateData(anchorset = ZR.anchors, dims = 1:40)

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
  if (SO.integrated$ADT_classification[i] == "ZR-CTL") {
    SO.integrated@meta.data$mapLabel[i] <- "Control"
  }
  if (SO.integrated$ADT_classification[i] == "ZR-E-03") {
    SO.integrated@meta.data$mapLabel[i] <- "E-3"
  }
  if (SO.integrated$ADT_classification[i] == "ZR-E-06") {
    SO.integrated@meta.data$mapLabel[i] <- "E-6"
  }
  if (SO.integrated$ADT_classification[i] == "ZR-E-24") {
    SO.integrated@meta.data$mapLabel[i] <- "E-24"
  }
  if (SO.integrated$ADT_classification[i] == "ZR-E-48") {
    SO.integrated@meta.data$mapLabel[i] <- "E-48"
  }
  if (SO.integrated$ADT_classification[i] == "ZR-E-72") {
    SO.integrated@meta.data$mapLabel[i] <- "E-72"
  }
  if (SO.integrated$ADT_classification[i] == "ZR-ED-04") {
    SO.integrated@meta.data$mapLabel[i] <- "T-48"
  }
  if (SO.integrated$ADT_classification[i] == "ZR-ED-49") {
    SO.integrated@meta.data$mapLabel[i] <- "TE-48"
  }
  if (SO.integrated$ADT_classification[i] == "ZR-EP-04") {
    SO.integrated@meta.data$mapLabel[i] <- "P-48"
  }
  if (SO.integrated$ADT_classification[i] == "ZR-EP-49") {
    SO.integrated@meta.data$mapLabel[i] <- "EP-48"
  }
  if (SO.integrated$ADT_classification[i] == "ZR-D-04") {
    SO.integrated@meta.data$mapLabel[i] <- "T-3"
  }
  if (SO.integrated$ADT_classification[i] == "ZR-D-49") {
    SO.integrated@meta.data$mapLabel[i] <- "TE-3"
  }
  if (SO.integrated$ADT_classification[i] == "ZR-IC-03") {
    SO.integrated@meta.data$mapLabel[i] <- "F-3"
  }
  if (SO.integrated$ADT_classification[i] == "ZR-IC-48") {
    SO.integrated@meta.data$mapLabel[i] <- "F-48"
  }
  if (SO.integrated$ADT_classification[i] == "ZR-T-03") {
    SO.integrated@meta.data$mapLabel[i] <- "TAM-3"
  }
  if (SO.integrated$ADT_classification[i] == "ZR-T-48") {
    SO.integrated@meta.data$mapLabel[i] <- "TAM-48"
  }
  if (SO.integrated$ADT_classification[i] == "ZR-P-04") {
    SO.integrated@meta.data$mapLabel[i] <- "P-3"
  }
  if (SO.integrated$ADT_classification[i] == "ZR-P-49") {
    SO.integrated@meta.data$mapLabel[i] <- "EP-3"
  }
}

DimPlot(SO.integrated, reduction = "umap", group.by = "seurat_clusters")
DimPlot(SO.integrated, reduction = "umap", group.by = "mapLabel")
DimPlot(SO.integrated, reduction = "umap", group.by = "Phase")


EP_ZR <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "EP-3" | mapLabel == "EP-48" | mapLabel == "E-3" | mapLabel == "E-48")
TE_ZR <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "TE-3" | mapLabel == "TE-48" | mapLabel == "E-3" | mapLabel == "E-48")
EPTE_ZR <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "EP-3" | mapLabel == "EP-48" | mapLabel == "E-3" | mapLabel == "E-48" | mapLabel == "TE-3" | mapLabel == "TE-48")
TAMF_ZR <- subset(SO.integrated, mapLabel == "E-48" | mapLabel == "TAM-3" | mapLabel == "TAM-48" | mapLabel == "F-3" | mapLabel == "F-48")
PEPE_ZR <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "EP-3" | mapLabel == "EP-48" | mapLabel == "E-3" | mapLabel == "E-48" | mapLabel == "P-3" | mapLabel == "P-48")
TETE_ZR <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "TE-3" | mapLabel == "TE-48" | mapLabel == "E-3" | mapLabel == "E-48" | mapLabel == "T-3" | mapLabel == "T-48")
TAME2_ZR <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "E-3" | mapLabel == "E-48" | mapLabel == "TAM-3" | mapLabel == "TAM-48")
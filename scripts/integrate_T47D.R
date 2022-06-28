library(Seurat)
library(tidyverse)

SO1 <- readRDS("E2_T47D/T47D_SO.rds")
SO2 <- readRDS("PG_TE_T47D/T47D_SO.rds")
SO3 <- readRDS("TAM_PG_IC_T47D/T47D_SO.rds")

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
  if (SO.integrated$ADT_classification[i] == "TD-CTL") {
    SO.integrated@meta.data$mapLabel[i] <- "Control"
  }
  if (SO.integrated$ADT_classification[i] == "TD-E-03") {
    SO.integrated@meta.data$mapLabel[i] <- "E-3"
  }
  if (SO.integrated$ADT_classification[i] == "TD-E-06") {
    SO.integrated@meta.data$mapLabel[i] <- "E-6"
  }
  if (SO.integrated$ADT_classification[i] == "TD-E-24") {
    SO.integrated@meta.data$mapLabel[i] <- "E-24"
  }
  if (SO.integrated$ADT_classification[i] == "TD-E-48") {
    SO.integrated@meta.data$mapLabel[i] <- "E-48"
  }
  if (SO.integrated$ADT_classification[i] == "TD-E-72") {
    SO.integrated@meta.data$mapLabel[i] <- "E-72"
  }
  if (SO.integrated$ADT_classification[i] == "TD-ED-03") {
    SO.integrated@meta.data$mapLabel[i] <- "T-48"
  }
  if (SO.integrated$ADT_classification[i] == "TD-ED-48") {
    SO.integrated@meta.data$mapLabel[i] <- "TE-48"
  }
  if (SO.integrated$ADT_classification[i] == "TD-EP-03") {
    SO.integrated@meta.data$mapLabel[i] <- "P-48"
  }
  if (SO.integrated$ADT_classification[i] == "TD-EP-48") {
    SO.integrated@meta.data$mapLabel[i] <- "EP-48"
  }
  if (SO.integrated$ADT_classification[i] == "TD-D-03") {
    SO.integrated@meta.data$mapLabel[i] <- "T-3"
  }
  if (SO.integrated$ADT_classification[i] == "TD-D-48") {
    SO.integrated@meta.data$mapLabel[i] <- "TE-3"
  }
  if (SO.integrated$ADT_classification[i] == "TD-IC-03") {
    SO.integrated@meta.data$mapLabel[i] <- "F-3"
  }
  if (SO.integrated$ADT_classification[i] == "TD-IC-48") {
    SO.integrated@meta.data$mapLabel[i] <- "F-48"
  }
  if (SO.integrated$ADT_classification[i] == "TD-T-03") {
    SO.integrated@meta.data$mapLabel[i] <- "TAM-3"
  }
  if (SO.integrated$ADT_classification[i] == "TD-T-48") {
    SO.integrated@meta.data$mapLabel[i] <- "TAM-48"
  }
  if (SO.integrated$ADT_classification[i] == "TD-P-03") {
    SO.integrated@meta.data$mapLabel[i] <- "P-3"
  }
  if (SO.integrated$ADT_classification[i] == "TD-P-48") {
    SO.integrated@meta.data$mapLabel[i] <- "EP-3"
  }
}

DimPlot(SO.integrated, reduction = "umap", group.by = "seurat_clusters")
DimPlot(SO.integrated, reduction = "umap", group.by = "mapLabel")
DimPlot(SO.integrated, reduction = "umap", group.by = "Phase")

TAM_T47D <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "TAM-3" | mapLabel == "TAM-48")
Fulv_T47D <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "F-3" | mapLabel == "F-48")
PG_T47D <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "P-3" | mapLabel == "P-48")
EP_T47D <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "EP-3" | mapLabel == "EP-48" | mapLabel == "E-3" | mapLabel == "E-48")
DHT_T47D <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "T-3" | mapLabel == "T-48")
TE_T47D <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "TE-3" | mapLabel == "TE-48" | mapLabel == "E-3" | mapLabel == "E-48")
EPTE_T47D <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "EP-3" | mapLabel == "EP-48" | mapLabel == "E-3" | mapLabel == "E-48" | mapLabel == "TE-3" | mapLabel == "TE-48")
TAMF_T47D <- subset(SO.integrated, mapLabel == "E-48" | mapLabel == "TAM-3" | mapLabel == "TAM-48" | mapLabel == "F-3" | mapLabel == "F-48")
PEPE_T47D <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "EP-3" | mapLabel == "EP-48" | mapLabel == "E-3" | mapLabel == "E-48" | mapLabel == "P-3" | mapLabel == "P-48")
TETE_T47D <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "TE-3" | mapLabel == "TE-48" | mapLabel == "E-3" | mapLabel == "E-48" | mapLabel == "T-3" | mapLabel == "T-48")
TAME2_T47D <- subset(SO.integrated, mapLabel == "Control" | mapLabel == "E-3" | mapLabel == "E-48" | mapLabel == "TAM-3" | mapLabel == "TAM-48")


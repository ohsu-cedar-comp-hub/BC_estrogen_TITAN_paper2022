# This script shows the general pipeline that was used for the pre-processing and basic Seurat analysis of all 
# of the single cell RNA sequencing data

# This script should be run after the raw data has been run through CellRanger

# As with cellRanger, this should be run separately for each cell line or organoid
 
nFeature_min = #800 for T47D E2   #700 for T47D PG_TE_T47D   #no min for T47D TAM_PG_IC_T47D  #700 for T47D PEPE    #1000 for ZR_E2(ZR1)   #0 for ZR_23(ZR2 and ZR3)      #500 for MCF7_E2   #500 for MCF7_PD   #500 for MCF7_CL   #1000 for TAMR  #no min for HCI003  #no min for HCI011 (C,E,and E2PG)  
nFeature_max = #3500 for T47D E2  #3000 for T47D PG_TE_T47D  #3000 for T47D TAM_PG_IC_T47D    #3000 for T47D PEPE   #4000 for ZR_E2(ZR1)   #3600 for ZR_23(ZR2 and ZR3)   #2500 for MCF7_E2  #2500 for MCF7_PD  #2500 for MCF7_CL  #4000 for TAMR  #2000 for HCI003    #no min for HCI011 (C,E,and E2PG)  
max_mt       = #30 for T47D E2    #30 for T47D PG_TE_T47D    #28 for T47D TAM_PG_IC_T47D      #30 for T47D PEPE     #13 for ZR_E2(ZR1)     #20 for ZR_23(ZR2 and ZR3)     #15 for MCF7_E2    #15 for MCF7_PD    #15 for MCF7_CL    #20 for TAMR    #10 for HCI003      #20 for HCI011 (C,E,and E2PG)
PCA_dims     = #1:40 for T47D E2  #1:40 for T47D PG_TE_T47D  #1:40 for T47D TAM_PG_IC_T47D    #1:40 for T47D PEPE   #1:40 for ZR_E2(ZR1)   #1:40 for ZR_23(ZR2 and ZR3)   #1:40 for MCF7_E2  #1:40 for MCF7_PD  #1:40 for MCF7_CL  #1:40 for TAMR  #1:40 for HCI003    #1:20 for HCI011 (C,E,and E2PG)
clus_res     = #0.3 for T47D E2   #0.3 for T47D PG_TE_T47D   #0.4 for T47D TAM_PG_IC_T47D     #0.3 for T47D PEPE    #0.3 for ZR_E2(ZR1)    #0.3 for ZR_23(ZR2 and ZR3)    #0.3 for MCF7_E2   #0.3 for MCF7_PD   #0.3 for MCF7_CL   #0.3 for TAMR   #0.5 for HCI003     #0.2 for HCI011 (C,E,and E2PG)
n_neighbors  = #20 for T47D E2    #20 for T47D PG_TE_T47D    #20 for T47D TAM_PG_IC_T47D      #20 for T47D PEPE     #40 for ZR_E2(ZR1)     #40 for ZR_23(ZR2 and ZR3)     #20 for MCF7_E2    #20 for MCF7_PD    #20 for MCF7_CL    #20 for TAMR    #50 for HCI003      #30 for HCI011 (C,E,and E2PG)
min_dist     = #0.2 for T47D E2   #0.2 for T47D PG_TE_T47D   #0.2 for T47D TAM_PG_IC_T47D     #0.2 for T47D PEPE    #0.3 for ZR_E2(ZR1)    #0.3 for ZR_23(ZR2 and ZR3)    #0.2 for MCF7_E2   #0.2 for MCF7_PD   #0.2 for MCF7_CL   #0.2 for TAMR   #0.4 for HCI003     #0.3 for HCI011 (C,E,and E2PG)
nCount_max   = #60000 for HCI003  #50000 for HCI011 (C,E,and E2PG)
  
  
library(Seurat)
library(tidyverse)

raw.counts <- Read10X(data.dir = "<path_to_cellRanger_output>/outs/filtered_feature_bc_matrix")

### If the data has been barcoded to separate out the different treatments/timepoints
### then this next section should be run to add the barcode counts to the Seurat object

barcode_ID <- raw.counts[["Antibody Capture"]]

SO <- CreateSeuratObject(counts=raw.counts$'Gene Expression', min.cells = 5, min.features = 200)
ADTassay <- (counts = barcode_ID)
missing_link <- setdiff(colnames(x = CreateAssayObject(counts = barcode_ID)), colnames(x = SO))
ADTassay <- ADTassay[ , -which(colnames(ADTassay) %in% missing_link)]
SO[["ADT"]] <- CreateAssayObject(ADTassay)

### This next section should be run for all datasets

SO[["percent.mt"]] <- PercentageFeatureSet(object = SO, pattern = "^MT-")

VlnPlot(object = SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

SO <- subset(x = SO, subset = nFeature_RNA > nFeature_min & nFeature_RNA < nFeature_max & percent.mt < max_mt) # & nCount_RNA < nCount_max (only used on organoids)

SO <- NormalizeData(SO)
SO <- FindVariableFeatures(object = SO, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(SO)
SO <- ScaleData(object = SO, features = all.genes)

SO <- RunPCA(SO, features = VariableFeatures(object = SO))
ElbowPlot(SO, ndims = 50)

### This next section should only be run if an ADT assay is present and will identify the cells

SO <- NormalizeData(SO, assay = "ADT", normalization.method = "CLR")
SO <- ScaleData(SO, assay = "ADT")

SO <- HTODemux(SO, assay = "ADT", positive.quantile = 0.99)

table(SO$ADT_classification.global)
HTOHeatmap(SO, assay = "ADT")

SO <- subset(SO, ADT_classification.global == "Singlet")

### The rest of the script should be run on all datasets

##########################################################################################################
##### PLEASE READ. VERY IMPORTANT TO INTEGRATE FIRST TO GET CORRECT LABELS
##### For MCF7, T47D, and ZR, integration scripts were run and different treatmant combo objects were 
##### created prior to this step, then those objects were used for the rest of the steps
##########################################################################################################

cc.genes <- readLines(con = "<location_of_CCgenes>/regev_lab_cell_cycle_genes.txt") #This text file can be downloaded from the Seurat website

s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

SO <- CellCycleScoring(SO, s.features = s.genes, g2m.features = g2m.genes)
SO <- ScaleData(SO, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(SO))

SO <- RunPCA(object=SO, features=VariableFeatures(object = SO))
SO <- FindNeighbors(object=SO, dims=PCA_dims)
SO <- FindClusters(object=SO, resolution = clus_res)
SO <- RunUMAP(object=SO, dims=PCA_dims, n.neighbors = n_neighbors, min.dist = min_dist)

DimPlot(SO, reduction = "umap", group.by = "seurat_clusters")
DimPlot(SO, reduction = "umap", group.by = "ADT_classification")
DimPlot(SO, reduction = "umap", group.by = "Phase")

cluster.markers <- FindAllMarkers(object=SO, only.pos=T, logfc.threshold=0.2)
top.markers <- cluster.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)

meta.data <- SO@meta.data
condition <- as.character(meta.data$seurat_clusters)

counts <- group_by(meta.data, seurat_clusters, ADT_classification) %>% summarise(count=n())

ggplot(counts, aes(y=count,x=ADT_classification,fill=seurat_clusters)) + geom_bar(position='fill',stat='identity')

top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DoHeatmap(SO, features= top10$gene, group.by = "seurat_clusters") + NoLegend()


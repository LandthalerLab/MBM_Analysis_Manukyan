library(VoltRon)
library(Seurat)
library(harmony)
library(dplyr)
library(ComplexHeatmap)
library(xlsx)

####
# data ####
####

# import data
brain_vr <- readRDS("../data/VoltRonData/brain_all_merged.rds")

# score markers (Interferon, Immune etc.)
score_markers <- read.xlsx("../data/Supplementary Information/score_markers.xlsx", sheetName = "Sheet1")
score_markers <- as.list(score_markers)
score_markers <- lapply(score_markers, function(x) x[!is.na(x)])

####
# overlay images ####
####

# get clustered data
brain_vr_Pat_90_P <- subset(brain_vr, assays = c("Assay14"))

# get image data and register
brain_image_Pat_90_P  <- importImageData(image = "../data/images/H&E/TMA.TIFF/Xenium-5.tif", tile.size = 10)

# register and save images, use 2016
xen_reg1 <- registerSpatialData(object_list = list(brain_vr_Pat_90_P, brain_image_Pat_90_P))
vrImages(brain_vr_Pat_90_P[["Assay14"]], name = "DAPI", channel = "H&E") <- vrImages(xen_reg1$registered_spat[[2]], name = "main_reg")

####
# clustering ####
####

# clustering with Seurat
brain_seu_Pat_90_P <- VoltRon::as.Seurat(brain_vr_Pat_90_P, cell.assay = "Xenium", type = "image")
brain_seu_Pat_90_P <- NormalizeData(brain_seu_Pat_90_P, scale.factor = 100)

# dimensionality reduction
brain_seu_Pat_90_P <- ScaleData(brain_seu_Pat_90_P)
brain_seu_Pat_90_P <- RunPCA(brain_seu_Pat_90_P, features = Features(brain_seu_Pat_90_P), npcs = 20)
brain_seu_Pat_90_P <- RunUMAP(brain_seu_Pat_90_P, dims = 1:20)
DimPlot(brain_seu_Pat_90_P, reduction = "umap", group.by = "Sample", split.by = "Sample")

# clustering
brain_seu_Pat_90_P <- FindNeighbors(brain_seu_Pat_90_P, dims = 1:20)
brain_seu_Pat_90_P <- FindClusters(brain_seu_Pat_90_P, resolution = c(0.6,0.8,1.0,1.1,1.2))

# visualize
g_list <- list()
g_list[["Xenium_snn_res.0.6"]] <- DimPlot(brain_seu_Pat_90_P, reduction = "umap", group.by = "Xenium_snn_res.0.6", label = T)
g_list[["Xenium_snn_res.0.8"]] <- DimPlot(brain_seu_Pat_90_P, reduction = "umap", group.by = "Xenium_snn_res.0.8", label = T)
g_list[["Xenium_snn_res.1"]] <- DimPlot(brain_seu_Pat_90_P, reduction = "umap", group.by = "Xenium_snn_res.1", label = T)
g_list[["Xenium_snn_res.1.1"]] <- DimPlot(brain_seu_Pat_90_P, reduction = "umap", group.by = "Xenium_snn_res.1.1", label = T)
g_list[["Xenium_snn_res.1.2"]] <- DimPlot(brain_seu_Pat_90_P, reduction = "umap", group.by = "Xenium_snn_res.1.2", label = T)
ggarrange(plotlist = g_list, ncol = 3, nrow = 2)

# res 1.1 looks good
brain_vr_Pat_90_P$cluster <- as.character(brain_seu_Pat_90_P$Xenium_snn_res.1.1)

####
# marker analysis ####
####

####
## check resolution 1.1 ####
####

# marker analysis
brain_seu_Pat_90_P <- as.Seurat(brain_vr_Pat_90_P, type = "image", cell.assay = "Xenium")
Idents(brain_seu_Pat_90_P) <- "Xenium_snn_res.1.1"
markers <- FindAllMarkers(brain_seu_Pat_90_P)
top_markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.8, p_val_adj < 0.05, pct.1 > 0.5)

####
# annotation ####
####

# select a top_markers from marker analysis section
top_markers <- top_markers %>%
  group_by(cluster) %>%
  mutate(temp = gsub(" ", "", paste(Cell.type.label, collapse = ","))) %>%
  mutate(FinalAnnotation = names(table(strsplit(temp, split = ",")[[1]]))[which.max(table(strsplit(temp, split = ",")[[1]]))]) %>%
  select(-temp)

####
## Manual annotation, ignore Final ####
####

# cluster 0 Tumorcells (CHL1/MET+)
# cluster 1 Tumorcells? (FBLN1+)
# cluster 2 Tumorcells? (TGFBI/SOX9+)
# cluster 3 Tumorcells (GAS2L3/) (subcluster ? SNTB2/BZW2)
# cluster 4 Tumorcells (TOP2A/MKI67/CDK1/CCNB2/CENPF)
# cluster 5 TAMs_1
# cluster 6 TAMs_2
# cluster 7 VLMC_1 (DCN+)
# cluster 8 TAP1+ cells
# cluster 9 VLMC_2 (DCN+)
# cluster 10 undefined?
# cluster 11 Tumorcells (MET/POSTN+) 
# cluster 12 endothelialcells_1
# cluster 13 endothelialcells_2 
# cluster 14 CD8 t cells
# cluster 15 t cells
# cluster 16 CYTIP+ cells

# get final annotations
annotations <- c(
  "Tumorcells_1 (CHL1/MET+)",
  "Tumorcells_2? (FBLN1+)",
  "Tumorcells_3? (TGFBI/SOX9+)",
  "Tumorcells_4 (GAS2L3/) (subcluster ? SNTB2/BZW2)",
  "Tumorcells_5 (TOP2A/MKI67/CDK1/CCNB2/CENPF)",
  "TAMs_1",
  "TAMs_2",
  "VLMC_1 (DCN+)",
  "TAP1+ cells",
  "VLMC_2 (DCN+)",
  "undefined?",
  "Tumorcells_6 (MET/POSTN+)",
  "endothelialcells_1",
  "endothelialcells_2",
  "CD8 t cells",
  "t cells",
  "CYTIP+ cells"
)

# mark cell types in marker table
top_markers$FinalAnnotation <- annotations[match(as.numeric(as.character(top_markers$cluster))+1, 1:length(annotations))]

# mark cell types in VoltRon and Seurat object
brain_vr_Pat_90_P$CellType <- annotations[match(as.numeric(as.character(brain_vr_Pat_90_P$cluster))+1, 1:length(annotations))]
brain_seu_Pat_90_P$CellType <- annotations[match(as.numeric(as.character(brain_vr_Pat_90_P$cluster))+1, 1:length(annotations))]

# transfer umap
vrEmbeddings(brain_vr_Pat_90_P, type = "umap", overwrite = TRUE) <- Embeddings(brain_seu_Pat_90_P, reduction = "umap")

####
## Update annotation (Major Cell Type) ####
####

temp <- brain_vr_Pat_90_P$CellType
temp[temp == "CD8 t cells"] <- "Immune Cells"
temp[temp == "t cells"] <- "Immune Cells"
temp[temp == "CYTIP+ cells"] <- "Immune Cells"
temp[grep("TAMs", temp)] <- "TAMs"
temp[grep("endothelialcells", temp)] <- "Endothelial Cells"
temp[grep("VLMC", temp)] <- "VLMC/CAF"
temp[grep("Tumorcells_5", temp)] <- "Proliferating Tumor Cells 1"
temp[grep("Tumorcells_1", temp)] <- "Tumor Cells 1"
temp[grep("Tumorcells_2", temp)] <- "Tumor Cells 2"
temp[grep("Tumorcells_3", temp)] <- "Tumor Cells 3"
temp[grep("Tumorcells_4", temp)] <- "Tumor Cells 4"
temp[grep("Tumorcells_6", temp)] <- "Tumor Cells 5"
temp[temp == "undefined?"] <- "undefined"
temp[temp == "TAP1+ cells"] <- "Tumor Cells 6"
brain_vr_Pat_90_P$MajorCellType <- temp
# saveRDS(brain_vr_Pat_90_P, file = "../data/VoltRonData/brain_vr_Pat_90_P_annotated.rds")

####
# annotate Hot tumors ####
####

####
## get immunescores ####
####

# get score per modules using Seurat
brain_vr_Pat_90_P_assay1 <- brain_vr_Pat_90_P
brain_vr_Pat_90_P_assay1_seu <- VoltRon::as.Seurat(brain_vr_Pat_90_P_assay1, cell.assay = "Xenium", type = "image")
brain_vr_Pat_90_P_assay1_seu <- NormalizeData(brain_vr_Pat_90_P_assay1_seu, scale.factor = 1000)
for(i in 1:length(score_markers)){
  cur_markers <- intersect(score_markers[[i]], SeuratObject::Features(brain_vr_Pat_90_P_assay1_seu))
  brain_vr_Pat_90_P_assay1_seu <- AddModuleScore(brain_vr_Pat_90_P_assay1_seu, features = list(cur_markers), name = names(score_markers)[i], ctrl = 5)  
}

# transfer module scores to voltron objects
temp <- rep(NA, length(vrSpatialPoints(brain_vr_Pat_90_P)))
names(temp) <- vrSpatialPoints(brain_vr_Pat_90_P)
temp[names(brain_vr_Pat_90_P_assay1_seu$ImmuneScores1)] <- brain_vr_Pat_90_P_assay1_seu$ImmuneScores1
brain_vr_Pat_90_P$ImmuneScores <- temp

# visualize score
vrSpatialFeaturePlot(brain_vr_Pat_90_P, features = c("ImmuneScores"),          
                     plot.segments = FALSE, pt.size = 1, n.tile = 300, alpha = 1)
vrSpatialPlot(brain_vr_Pat_90_P, group.by = "annotation")

####
## annotate hot tumor ####
####

brain_vr_Pat_90_P <- annotateSpatialData(brain_vr_Pat_90_P, use.image = TRUE, channel = "H&E")
# saveRDS(brain_vr_Pat_90_P, file = "../data/VoltRonData/brain_vr_Pat_90_P_annotated.rds")

####
# save as AnnData ####
####

# interactive TissUUmaps tool
as.AnnData(brain_vr_Pat_90_P, file = "../data/AnnDataData/brain_anndata_Pat_90_P.h5ad", flip_coordinates = TRUE, 
           assay = "Xenium", method = "anndata", name = "DAPI", channel = "H&E")

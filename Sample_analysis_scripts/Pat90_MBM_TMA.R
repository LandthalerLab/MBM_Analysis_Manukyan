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
brain_vr_Pat_90_MBM_TMA <- subset(brain_vr, assays = c("Assay15", "Assay16"))

# get image data and register
brain_image_Pat_90_MBM_TMA_A <- importImageData(image = "../data/images/H&E/TMA.TIFF/Xenium-13.tif", tile.size = 10)
brain_image_Pat_90_MBM_TMA_B <- importImageData(image = "../data/images/H&E/TMA.TIFF/Xenium-9.tif", tile.size = 10)

# register and save images
xen_reg1 <- registerSpatialData(object_list = list(subset(brain_vr_Pat_90_MBM_TMA, samples = "Pat_90_MBM_TMA_23_D_A"),
                                                   brain_image_Pat_90_MBM_TMA_A))
vrImages(brain_vr_Pat_90_MBM_TMA[["Assay16"]], name = "DAPI", channel = "H&E") <- vrImages(xen_reg1$registered_spat[[2]], name = "main_reg")

xen_reg2 <- registerSpatialData(object_list = list(subset(brain_vr_Pat_90_MBM_TMA, samples = "Pat_90_MBM_TMA_23_D_B"),
                                                   brain_image_Pat_90_MBM_TMA_B))
vrImages(brain_vr_Pat_90_MBM_TMA[["Assay15"]], name = "DAPI", channel = "H&E") <- vrImages(xen_reg2$registered_spat[[2]], name = "main_reg")


####
# clustering ####
####

# clustering with Seurat
brain_seu_Pat_90_MBM_TMA <- VoltRon::as.Seurat(brain_vr_Pat_90_MBM_TMA, cell.assay = "Xenium", type = "image")
brain_seu_Pat_90_MBM_TMA <- NormalizeData(brain_seu_Pat_90_MBM_TMA, scale.factor = 100)

# dimensionality reduction
brain_seu_Pat_90_MBM_TMA <- ScaleData(brain_seu_Pat_90_MBM_TMA)
brain_seu_Pat_90_MBM_TMA <- RunPCA(brain_seu_Pat_90_MBM_TMA, features = Features(brain_seu_Pat_90_MBM_TMA), npcs = 20)
brain_seu_Pat_90_MBM_TMA <- RunUMAP(brain_seu_Pat_90_MBM_TMA, dims = 1:20)
DimPlot(brain_seu_Pat_90_MBM_TMA, reduction = "umap", group.by = "Sample", split.by = "Sample")

# clustering
brain_seu_Pat_90_MBM_TMA <- FindNeighbors(brain_seu_Pat_90_MBM_TMA, dims = 1:20)
brain_seu_Pat_90_MBM_TMA <- FindClusters(brain_seu_Pat_90_MBM_TMA, resolution = c(0.6,0.8,1.0,1.1,1.2))

# visualize
g_list <- list()
g_list[["Xenium_snn_res.0.6"]] <- DimPlot(brain_seu_Pat_90_MBM_TMA, reduction = "umap", group.by = "Xenium_snn_res.0.6", label = T)
g_list[["Xenium_snn_res.0.8"]] <- DimPlot(brain_seu_Pat_90_MBM_TMA, reduction = "umap", group.by = "Xenium_snn_res.0.8", label = T)
g_list[["Xenium_snn_res.1"]] <- DimPlot(brain_seu_Pat_90_MBM_TMA, reduction = "umap", group.by = "Xenium_snn_res.1", label = T)
g_list[["Xenium_snn_res.1.1"]] <- DimPlot(brain_seu_Pat_90_MBM_TMA, reduction = "umap", group.by = "Xenium_snn_res.1.1", label = T)
g_list[["Xenium_snn_res.1.2"]] <- DimPlot(brain_seu_Pat_90_MBM_TMA, reduction = "umap", group.by = "Xenium_snn_res.1.2", label = T)
ggarrange(plotlist = g_list, ncol = 3, nrow = 2)

# res 1.1 looks good
brain_vr_Pat_90_MBM_TMA$cluster <- as.character(brain_seu_Pat_90_MBM_TMA$Xenium_snn_res.1.1)

####
# marker analysis ####
####

####
## check resolution 1.1 ####
####

# marker analysis
# brain_seu_Pat_90_MBM_TMA <- as.Seurat(brain_vr_Pat_90_MBM_TMA, type = "image", cell.assay = "Xenium")
brain_seu_Pat_90_MBM_TMA <- readRDS("../data/SeuratData/brain_seu_Pat_90_MBM_TMA_clustered.rds")
Idents(brain_seu_Pat_90_MBM_TMA) <- "Xenium_snn_res.1.1"
markers <- FindAllMarkers(brain_seu_Pat_90_MBM_TMA)
top_markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.8, p_val_adj < 0.05, pct.1 > 0.5)

# check markers
g1 <- FeaturePlot(brain_seu_Pat_90_MBM_TMA, features = c("MET"), reduction = "umap")
g2 <- DimPlot(brain_seu_Pat_90_MBM_TMA, group.by = "Xenium_snn_res.1.1", reduction = "umap", label = TRUE)
g1 | g2


####
# annotation ####
####

# select a top_markers from marker analysis section
top_markers <- top_markers %>%
  group_by(cluster) %>%
  mutate(temp = gsub(" ", "", paste(Cell.type.label, collapse = ","))) %>%
  mutate(FinalAnnotation = names(table(strsplit(temp, split = ",")[[1]]))[which.max(table(strsplit(temp, split = ",")[[1]]))]) %>%
  select(-temp)
top_markers_celltype <- top_markers %>%
  select(cluster, FinalAnnotation) %>%
  distinct() %>%
  group_by(FinalAnnotation) %>%
  mutate(FinalAnnotation = paste(FinalAnnotation, 1:n(), sep = "_"))

####
## Manual annotation, ignore Final ####
####

# cluster 0 Tumorcells_1 (CENPF+, BZW2+, TOP2A+)
# cluster 1 Neurons
# cluster 2 NES+ cells
# cluster 3 Neurons
# cluster 4 Tumorcells_2 (MET+)
# cluster 5 CD8 T cells
# cluster 6 Tumorcells_3 (MET+, CENPF+, BZW2+, TOP2A+)
# cluster 7 Tumorcells_4 (MET+)
# cluster 8 TAMs_1
# cluster 9 endothelial cells
# cluster 10 TAMs_2
# cluster 11 TAMs_3
# cluster 12 CD8 T cells
# cluster 13 Neurons
# cluster 14 ITGAX+ cells
# cluster 15 Tumorcells_5 (MET+, ABCB5+)
# cluster 16 TAMs_4
# cluster 17 Tumorcells_6 (ABCB5+)

# get final annotations
annotations <- c(
  "Tumorcells_1 (CENPF+, BZW2+, TOP2A+)",
  "Neurons_1",
  "NES+ cells",
  "Neurons_2",
  "Tumorcells_2 (MET+)",
  "CD8 T cells",
  "Tumorcells_3 (MET+, CENPF+, BZW2+, TOP2A+)",
  "Tumorcells_4 (MET+)",
  "TAMs_1",
  "endothelial cells",
  "TAMs_2",
  "TAMs_3",
  "CD8 T cells",
  "Neurons_3 ??",
  "ITGAX+ cells",
  "Tumorcells_5 (MET+, ABCB5+)",
  "TAMs_4",
  "Tumorcells_6 (ABCB5+)"
)

# mark cell types in marker table
top_markers$FinalAnnotation <- annotations[match(as.numeric(as.character(top_markers$cluster))+1, 1:length(annotations))]

# mark cell types in VoltRon and Seurat object
brain_vr_Pat_90_MBM_TMA$CellType <- annotations[match(as.numeric(as.character(brain_vr_Pat_90_MBM_TMA$cluster))+1, 1:length(annotations))]
brain_seu_Pat_90_MBM_TMA$CellType <- annotations[match(as.numeric(as.character(brain_vr_Pat_90_MBM_TMA$cluster))+1, 1:length(annotations))]

# update umap
vrEmbeddings(brain_vr_Pat_90_MBM_TMA, type = "umap", overwrite = TRUE) <- Embeddings(brain_seu_Pat_90_MBM_TMA, reduction = "umap")
vrEmbeddings(brain_vr_Pat_90_MBM_TMA, type = "pca_seurat", overwrite = TRUE) <- Embeddings(brain_seu_Pat_90_MBM_TMA, reduction = "pca")

####
## Update annotation (Major Cell Type) ####
####

temp <- brain_vr_Pat_90_MBM_TMA$CellType
temp[temp == "CD8 T cells"] <- "Immune Cells"
temp[grep("Tumorcells_4", temp)] <- "Tumor Cells 1"
temp[grep("Tumorcells_5", temp)] <- "Tumor Cells 2"
temp[grep("Tumorcells_6", temp)] <- "Tumor Cells 3"
temp[grep("endothelial cells", temp)] <- "Endothelial Cells"
temp[grep("Neurons", temp)] <- "undefined"
temp[grep("ITGAX", temp)] <- "TAMs"
temp[grep("TAMs", temp)] <- "TAMs"
temp[grep("Tumorcells_3", temp)] <- "Tumor Cells 4 (NES+)"
temp[grep("Tumorcells_2", temp)] <- "Tumor Cells 5 (NES+)"
temp[grep("Tumorcells_1", temp)] <- "Proliferating Tumor Cells 1 (NES+)"
brain_vr_Pat_90_MBM_TMA$MajorCellType <- temp
# saveRDS(brain_vr_Pat_90_MBM_TMA, file = "../data/VoltRonData/brain_vr_Pat_90_MBM_TMA_annotated.rds")

####
# Sample Embeddings ####
####

sample_names <- unique(brain_vr_Pat_90_MBM_TMA$Sample)
brain_vr_Pat_90_MBM_TMA_xenium <- subset(brain_vr_Pat_90_MBM_TMA, assay = "Xenium")
`%notin%` <- Negate("%in%")
brain_vr_Pat_90_MBM_TMA_xenium <- subset(brain_vr_Pat_90_MBM_TMA_xenium, MajorCellType %notin% c("Low Count Cells", "undefined"))
for(samp in sample_names){
  print(samp)
  brain_vr_Pat_90_MBM_TMA_sub <- subset(brain_vr_Pat_90_MBM_TMA_xenium, samples = samp)
  brain_vr_Pat_90_MBM_TMA_subseu <- as.Seurat(brain_vr_Pat_90_MBM_TMA_sub, cell.assay = "Xenium", type = "image")
  brain_vr_Pat_90_MBM_TMA_subseu <- NormalizeData(brain_vr_Pat_90_MBM_TMA_subseu, scale.factor = 100)
  brain_vr_Pat_90_MBM_TMA_subseu <- ScaleData(brain_vr_Pat_90_MBM_TMA_subseu)
  selected_features <- Features(brain_vr_Pat_90_MBM_TMA_subseu)[!grepl("^S2-|^HSV1", Features(brain_vr_Pat_90_MBM_TMA_subseu))]
  brain_vr_Pat_90_MBM_TMA_subseu <- RunPCA(brain_vr_Pat_90_MBM_TMA_subseu, features = selected_features, npcs = 20)
  brain_vr_Pat_90_MBM_TMA_subseu <- RunUMAP(brain_vr_Pat_90_MBM_TMA_subseu, dims = 1:20)
  vrEmbeddings(brain_vr_Pat_90_MBM_TMA, type = paste0("umap_", samp), overwrite = TRUE) <- Embeddings(brain_vr_Pat_90_MBM_TMA_subseu, reduction = "umap")
}

vrEmbeddingPlot(brain_vr_Pat_90_MBM_TMA, embedding = "umap_Pat_90_MBM_TMA_23_D_B", group.by = "MajorCellType", label = T)

####
# annotate Hot tumors ####
####

####
## get immunescores ####
####

# get score per modules using Seurat
brain_vr_Pat_90_MBM_TMA_assay1 <- subset(brain_vr_Pat_90_MBM_TMA, assays = "Assay15")
brain_vr_Pat_90_MBM_TMA_assay1_seu <- VoltRon::as.Seurat(brain_vr_Pat_90_MBM_TMA_assay1, cell.assay = "Xenium", type = "image")
brain_vr_Pat_90_MBM_TMA_assay1_seu <- NormalizeData(brain_vr_Pat_90_MBM_TMA_assay1_seu, scale.factor = 1000)
for(i in 1:length(score_markers)){
  cur_markers <- intersect(score_markers[[i]], SeuratObject::Features(brain_vr_Pat_90_MBM_TMA_assay1_seu))
  brain_vr_Pat_90_MBM_TMA_assay1_seu <- AddModuleScore(brain_vr_Pat_90_MBM_TMA_assay1_seu, features = list(cur_markers), name = names(score_markers)[i], ctrl = 5)  
}

# transfer module scores to voltron objects
temp <- rep(NA, length(vrSpatialPoints(brain_vr_Pat_90_MBM_TMA)))
names(temp) <- vrSpatialPoints(brain_vr_Pat_90_MBM_TMA)
temp[names(brain_vr_Pat_90_MBM_TMA_assay1_seu$ImmuneScores1)] <- brain_vr_Pat_90_MBM_TMA_assay1_seu$ImmuneScores1
brain_vr_Pat_90_MBM_TMA$ImmuneScores <- temp

# visualize score and annotations
g1 <- vrSpatialFeaturePlot(brain_vr_Pat_90_MBM_TMA, assay = "Assay15", features = c("ImmuneScores"),          
                           plot.segments = FALSE, pt.size = 1, n.tile = 200, alpha = 1, ncol = 3)
vrMainAssay(brain_vr_Pat_90_MBM_TMA) <- "ROIAnnotation"
brain_vr_Pat_90_MBM_TMA$Annotation <- "hot"
vrImages(brain_vr_Pat_90_MBM_TMA[["Assay17"]], channel = "H&E") <- vrImages(brain_vr_Pat_90_MBM_TMA, assay = "Assay15", channel = "H&E")
g2 <- vrSpatialPlot(brain_vr_Pat_90_MBM_TMA, assay = "ROIAnnotation", group.by = "Annotation",
                    plot.segments = TRUE, pt.size = 1, alpha = 1, ncol = 3, channel = "H&E")
g1 | g2

####
## annotate hot tumor ####
####

brain_vr_Pat_90_MBM_TMA <- annotateSpatialData(brain_vr_Pat_90_MBM_TMA, assay = "Assay15", use.image = TRUE, channel = "H&E")
# saveRDS(brain_vr_Pat_90_MBM_TMA, file = "../data/VoltRonData/brain_vr_Pat_90_MBM_TMA_annotated.rds")

####
# save as AnnData ####
####

# interactive TissUUmaps tool
as.AnnData(brain_vr_Pat_90_MBM_TMA, file = "../data/AnnDataData/brain_anndata_Pat_90_MBM_TMA.h5ad", flip_coordinates = TRUE, 
           assay = "Xenium", method = "anndata", name = "DAPI", channel = "H&E")

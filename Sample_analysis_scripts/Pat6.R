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

####
# overlay images ####
####

# get clustered data
brain_vr_Pat6 <- subset(brain_vr, assays = c("Assay1"))
brain_vr_Pat6 <- modulateImage(brain_vr_Pat6, brightness = 200)

# get image data and register
brain_image_Pat6  <- importImageData(image = "../data/images/H&E/TMA.TIFF/Xenium-4.tif", tile.size = 10)

# register and save images, use default size
# xen_reg1 <- registerSpatialData(object_list = list(brain_vr_Pat6, brain_image_Pat6))
xen_reg1 <- registerSpatialData(object_list = list(brain_vr_Pat6, brain_image_Pat6), keypoints = readRDS("../data/AuxiliaryData/Pat6_keypoints.rds"))
vrImages(brain_vr_Pat6[["Assay1"]], name = "DAPI", channel = "H&E") <- vrImages(xen_reg1$registered_spat[[2]], name = "main_reg")

####
# clustering ####
####

# clustering with Seurat
brain_seu_Pat6 <- VoltRon::as.Seurat(brain_vr_Pat6, cell.assay = "Xenium", type = "image")
brain_seu_Pat6 <- NormalizeData(brain_seu_Pat6, scale.factor = 100)

# dimensionality reduction
brain_seu_Pat6 <- ScaleData(brain_seu_Pat6)
brain_seu_Pat6 <- RunPCA(brain_seu_Pat6, features = Features(brain_seu_Pat6), npcs = 20)
brain_seu_Pat6 <- RunUMAP(brain_seu_Pat6, dims = 1:20)
DimPlot(brain_seu_Pat6, reduction = "umap", group.by = "Sample", split.by = "Sample")

# clustering
brain_seu_Pat6 <- FindNeighbors(brain_seu_Pat6, dims = 1:20)
brain_seu_Pat6 <- FindClusters(brain_seu_Pat6, resolution = c(0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8, 0.9, 1.0,1.1))

# visualize
g_list <- list()
g_list[["Xenium_snn_res.0.1"]] <- DimPlot(brain_seu_Pat6, reduction = "umap", group.by = "Xenium_snn_res.0.1", label = T)
g_list[["Xenium_snn_res.0.2"]] <- DimPlot(brain_seu_Pat6, reduction = "umap", group.by = "Xenium_snn_res.0.2", label = T)
g_list[["Xenium_snn_res.0.3"]] <- DimPlot(brain_seu_Pat6, reduction = "umap", group.by = "Xenium_snn_res.0.3", label = T)
g_list[["Xenium_snn_res.0.4"]] <- DimPlot(brain_seu_Pat6, reduction = "umap", group.by = "Xenium_snn_res.0.4", label = T)
g_list[["Xenium_snn_res.0.5"]] <- DimPlot(brain_seu_Pat6, reduction = "umap", group.by = "Xenium_snn_res.0.5", label = T)
g_list[["Xenium_snn_res.0.6"]] <- DimPlot(brain_seu_Pat6, reduction = "umap", group.by = "Xenium_snn_res.0.6", label = T)
g_list[["Xenium_snn_res.0.7"]] <- DimPlot(brain_seu_Pat6, reduction = "umap", group.by = "Xenium_snn_res.0.7", label = T)
g_list[["Xenium_snn_res.0.8"]] <- DimPlot(brain_seu_Pat6, reduction = "umap", group.by = "Xenium_snn_res.0.8", label = T)
g_list[["Xenium_snn_res.0.9"]] <- DimPlot(brain_seu_Pat6, reduction = "umap", group.by = "Xenium_snn_res.0.9", label = T)
g_list[["Xenium_snn_res.1"]] <- DimPlot(brain_seu_Pat6, reduction = "umap", group.by = "Xenium_snn_res.1", label = T)
g_list[["Xenium_snn_res.1.1"]] <- DimPlot(brain_seu_Pat6, reduction = "umap", group.by = "Xenium_snn_res.1.1", label = T)
ggarrange(plotlist = g_list, ncol = 3, nrow = 4)

# res 0.5 looks good
brain_vr_Pat6$cluster <- as.character(brain_seu_Pat6$Xenium_snn_res.0.5)

####
# marker analysis ####
####

####
## check resolution 1.1 ####
####

# marker analysis
brain_seu_Pat6 <- as.Seurat(brain_vr_Pat6, type = "image", cell.assay = "Xenium")
Idents(brain_seu_Pat6) <- "Xenium_snn_res.0.5"
markers <- FindAllMarkers(brain_seu_Pat6)
top_markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.8, p_val_adj < 0.05, pct.1 > 0.5)

####
# Annotation ####
####

####
## Manual annotation, ignore Final ####
####

# cluster 0 0 ? 
# cluster 1 Oligodendrocytes
# cluster 2 microglia/macrophages
# cluster 3 T cells
# cluster 4 Astrocytes
# cluster 5 5?
# cluster 6 OPC
# cluster 7 VLMC (DCN+)
# cluster 8 8?
# cluster 9 IGFBP3+ cells

# get final annotations
annotations <- c(
  "0?",  
  "Oligodendrocytes",
  "microglia/macrophages",
  "T cells",
  "Astrocytes",
  "5?",
  "OPC",
  "VLMC (DCN+)",
  "8?",
  "IGFBP3+ cells"
)

# mark cell types in marker table
top_markers$FinalAnnotation <- annotations[match(as.numeric(as.character(top_markers$cluster))+1, 1:length(annotations))]

# mark cell types in VoltRon and Seurat object
brain_vr_Pat6$CellType <- annotations[match(as.numeric(as.character(brain_vr_Pat6$cluster))+1, 1:length(annotations))]
brain_seu_Pat6$CellType <- annotations[match(as.numeric(as.character(brain_vr_Pat6$cluster))+1, 1:length(annotations))]

# update umap
vrEmbeddings(brain_vr_Pat6, type = "umap", overwrite = TRUE) <- Embeddings(brain_seu_Pat6, reduction = "umap")

####
## Update annotation (Major Cell Type) ####
####

temp <- brain_vr_Pat6$CellType
temp[grep("endothelialcells", temp)] <- "Endothelial Cells"
temp[grep("VLMC", temp)] <- "VLMC/CAF"
temp[grep("microglia", temp)] <- "TAMs"
temp[grep("T cells", temp)] <- "Immune Cells"
temp[grep("OPC", temp)] <- "Oligodendrocytes"
temp[temp=="5?"] <- "Low Count Cells"
temp[temp=="0?"] <- "undefined"
brain_vr_Pat6$MajorCellType <- temp
# saveRDS(brain_vr_Pat6, file = "../data/VoltRonData/brain_vr_Pat6_annotated.rds")

####
# save as AnnData ####
####

# interactive TissUUmaps tool
as.AnnData(brain_vr_Pat6, file = "../data/AnnDataData/brain_anndata_Pat6.h5ad", flip_coordinates = TRUE, 
           assay = "Xenium", method = "anndata", name = "DAPI", channel = "H&E")
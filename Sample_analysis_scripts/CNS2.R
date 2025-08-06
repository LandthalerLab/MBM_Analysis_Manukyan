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
brain_vr_CNS2 <- subset(brain_vr, assays = c("Assay9"))

# get image data and register
brain_image_CNS2  <- importImageData(image = "../data/images/H&E/TMA.TIFF/Xenium-2.tif", tile.size = 10)

# register and save images, use 2022 vs 1992
brain_vr_CNS2 <- modulateImage(brain_vr_CNS2, brightness = 300)
# xen_reg1 <- registerSpatialData(object_list = list(brain_vr_CNS2, brain_image_CNS2))
xen_reg1 <- registerSpatialData(object_list = list(brain_vr_CNS2, brain_image_CNS2), 
                                keypoints = readRDS(file = "../data/AuxiliaryData/CNS2_keypoints.rds")[[1]])
vrImages(brain_vr_CNS2[["Assay9"]], name = "DAPI", channel = "H&E") <- vrImages(xen_reg1$registered_spat[[2]], name = "main_reg")

####
# clustering ####
####

# clustering with Seurat
brain_seu_CNS2 <- VoltRon::as.Seurat(brain_vr_CNS2, cell.assay = "Xenium", type = "image")
brain_seu_CNS2 <- NormalizeData(brain_seu_CNS2, scale.factor = 100)

# dimensionality reduction
brain_seu_CNS2 <- ScaleData(brain_seu_CNS2)
brain_seu_CNS2 <- RunPCA(brain_seu_CNS2, features = Features(brain_seu_CNS2), npcs = 20)
brain_seu_CNS2 <- RunUMAP(brain_seu_CNS2, dims = 1:20)
DimPlot(brain_seu_CNS2, reduction = "umap", group.by = "Sample", split.by = "Sample")

# clustering
brain_seu_CNS2 <- FindNeighbors(brain_seu_CNS2, dims = 1:20)
brain_seu_CNS2 <- FindClusters(brain_seu_CNS2, resolution = c(0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,1.0,1.1))

# visualize
g_list <- list()
g_list[["Xenium_snn_res.0.1"]] <- DimPlot(brain_seu_CNS2, reduction = "umap", group.by = "Xenium_snn_res.0.1", label = T)
g_list[["Xenium_snn_res.0.2"]] <- DimPlot(brain_seu_CNS2, reduction = "umap", group.by = "Xenium_snn_res.0.2", label = T)
g_list[["Xenium_snn_res.0.3"]] <- DimPlot(brain_seu_CNS2, reduction = "umap", group.by = "Xenium_snn_res.0.3", label = T)
g_list[["Xenium_snn_res.0.4"]] <- DimPlot(brain_seu_CNS2, reduction = "umap", group.by = "Xenium_snn_res.0.4", label = T)
g_list[["Xenium_snn_res.0.5"]] <- DimPlot(brain_seu_CNS2, reduction = "umap", group.by = "Xenium_snn_res.0.5", label = T)
g_list[["Xenium_snn_res.0.6"]] <- DimPlot(brain_seu_CNS2, reduction = "umap", group.by = "Xenium_snn_res.0.6", label = T)
g_list[["Xenium_snn_res.0.7"]] <- DimPlot(brain_seu_CNS2, reduction = "umap", group.by = "Xenium_snn_res.0.7", label = T)
g_list[["Xenium_snn_res.1"]] <- DimPlot(brain_seu_CNS2, reduction = "umap", group.by = "Xenium_snn_res.1", label = T)
g_list[["Xenium_snn_res.1.1"]] <- DimPlot(brain_seu_CNS2, reduction = "umap", group.by = "Xenium_snn_res.1.1", label = T)
ggarrange(plotlist = g_list, ncol = 3, nrow = 3)

# res 0.5 looks good
brain_vr_CNS2$cluster <- as.character(brain_seu_CNS2$Xenium_snn_res.0.5)

####
# marker analysis ####
####

####
## check resolution 0.5 ####
####

# marker analysis
brain_seu_CNS2 <- as.Seurat(brain_vr_CNS2, type = "image", cell.assay = "Xenium")
Idents(brain_seu_CNS2) <- "Xenium_snn_res.0.5"
markers <- FindAllMarkers(brain_seu_CNS2)
top_markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.8, p_val_adj < 0.05, pct.1 > 0.5)

####
## Manual annotation, ignore Final ####
####

# cluster 0 Astocytes (GFAP+)
# cluster 1 Oligodendrocytes
# cluster 2 Microglia/Macrophages_1, Macrophages new (Josefine)
# cluster 3 Astrocytes (CHI3L1+))
# cluster 4 Neurons (STMN2+)
# cluster 5 Microglia/Macrophages_2, Microglia new (Josefine)
# cluster 6 Astorcytes (GJA1, AQP4+)
# cluster 7 endothelialcells

# get final annotations
annotations <- c(
  "Astocytes (GFAP+)",
  "Oligodendrocytes",
  "Macrophages",
  "Astrocytes (CHI3L1+))",
  "Neurons (STMN2+)",
  "Microglia",
  "Astorcytes (GJA1, AQP4+)",
  "endothelialcells"
)

# mark cell types in marker table
top_markers$FinalAnnotation <- annotations[match(as.numeric(as.character(top_markers$cluster))+1, 1:length(annotations))]

# mark cell types in VoltRon and Seurat object
brain_vr_CNS2$CellType <- annotations[match(as.numeric(as.character(brain_vr_CNS2$cluster))+1, 1:length(annotations))]
brain_seu_CNS2$CellType <- annotations[match(as.numeric(as.character(brain_vr_CNS2$cluster))+1, 1:length(annotations))]

# update umap
vrEmbeddings(brain_vr_CNS2, type = "umap", overwrite = TRUE) <- Embeddings(brain_seu_CNS2, reduction = "umap")

####
## Update annotation (Major Cell Type) ####
####

temp <- brain_vr_CNS2$CellType
temp[grep("endothelialcells", temp)] <- "Endothelial Cells"
temp[grep("VLMC", temp)] <- "VLMC/CAF"
temp[grep("Macrophages|Microglia", temp)] <- "TAMs"
temp[temp == "Astocytes (GFAP+)"] <- "Astrocytes"
temp[temp == "Astorcytes (GJA1, AQP4+)"] <- "Astrocytes"
temp[temp == "Astrocytes (CHI3L1+))"] <- "Astrocytes"
temp[temp == "Neurons (STMN2+)"] <- "Neurons"
brain_vr_CNS2$MajorCellType <- temp
# saveRDS(brain_vr_CNS2, file = "../data/VoltRonData/brain_vr_CNS2_annotated.rds")

####
# save as AnnData ####
####

# interactive TissUUmaps tool
as.AnnData(brain_vr_CNS2, file = "../data/AnnDataData/brain_anndata_CNS2.h5ad", flip_coordinates = TRUE, 
           assay = "Xenium", method = "anndata", name = "DAPI", channel = "H&E")

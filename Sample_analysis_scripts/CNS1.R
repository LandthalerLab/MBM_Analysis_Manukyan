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
brain_vr_CNS1 <- subset(brain_vr, assays = c("Assay12"))

# get image data and register
brain_image_CNS1  <- importImageData(image = "../data/images/H&E/TMA.TIFF/Xenium-14.tif", tile.size = 10)

# register and save images, use 1092 X 1068
xen_reg1 <- registerSpatialData(object_list = list(brain_vr_CNS1, brain_image_CNS1))
vrImages(brain_vr_CNS1[["Assay12"]], name = "DAPI", channel = "H&E") <- vrImages(xen_reg1$registered_spat[[2]], name = "main_reg")

####
# clustering ####
####

# clustering with Seurat
brain_seu_CNS1 <- VoltRon::as.Seurat(brain_vr_CNS1, cell.assay = "Xenium", type = "image")
brain_seu_CNS1 <- NormalizeData(brain_seu_CNS1, scale.factor = 100)

# dimensionality reduction
brain_seu_CNS1 <- ScaleData(brain_seu_CNS1)
brain_seu_CNS1 <- RunPCA(brain_seu_CNS1, features = Features(brain_seu_CNS1), npcs = 20)
brain_seu_CNS1 <- RunUMAP(brain_seu_CNS1, dims = 1:20)
DimPlot(brain_seu_CNS1, reduction = "umap", group.by = "Sample", split.by = "Sample")

# clustering
brain_seu_CNS1 <- FindNeighbors(brain_seu_CNS1, dims = 1:20)
brain_seu_CNS1 <- FindClusters(brain_seu_CNS1, resolution = c(0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8, 0.9, 1.0,1.1))

# visualize
g_list <- list()
g_list[["Xenium_snn_res.0.1"]] <- DimPlot(brain_seu_CNS1, reduction = "umap", group.by = "Xenium_snn_res.0.1", label = T)
g_list[["Xenium_snn_res.0.2"]] <- DimPlot(brain_seu_CNS1, reduction = "umap", group.by = "Xenium_snn_res.0.2", label = T)
g_list[["Xenium_snn_res.0.3"]] <- DimPlot(brain_seu_CNS1, reduction = "umap", group.by = "Xenium_snn_res.0.3", label = T)
g_list[["Xenium_snn_res.0.4"]] <- DimPlot(brain_seu_CNS1, reduction = "umap", group.by = "Xenium_snn_res.0.4", label = T)
g_list[["Xenium_snn_res.0.5"]] <- DimPlot(brain_seu_CNS1, reduction = "umap", group.by = "Xenium_snn_res.0.5", label = T)
g_list[["Xenium_snn_res.0.6"]] <- DimPlot(brain_seu_CNS1, reduction = "umap", group.by = "Xenium_snn_res.0.6", label = T)
g_list[["Xenium_snn_res.0.7"]] <- DimPlot(brain_seu_CNS1, reduction = "umap", group.by = "Xenium_snn_res.0.7", label = T)
g_list[["Xenium_snn_res.0.8"]] <- DimPlot(brain_seu_CNS1, reduction = "umap", group.by = "Xenium_snn_res.0.8", label = T)
g_list[["Xenium_snn_res.0.9"]] <- DimPlot(brain_seu_CNS1, reduction = "umap", group.by = "Xenium_snn_res.0.9", label = T)
g_list[["Xenium_snn_res.1"]] <- DimPlot(brain_seu_CNS1, reduction = "umap", group.by = "Xenium_snn_res.1", label = T)
g_list[["Xenium_snn_res.1.1"]] <- DimPlot(brain_seu_CNS1, reduction = "umap", group.by = "Xenium_snn_res.1.1", label = T)
ggarrange(plotlist = g_list, ncol = 3, nrow = 4)

# res 0.5 looks good
brain_vr_CNS1$cluster <- as.character(brain_seu_CNS1$Xenium_snn_res.0.8)

####
# marker analysis ####
####

####
## check resolution 0.8 ####
####

# marker analysis
brain_seu_CNS1 <- as.Seurat(brain_vr_CNS1, type = "image", cell.assay = "Xenium")
Idents(brain_seu_CNS1) <- "Xenium_snn_res.0.8"
markers <- FindAllMarkers(brain_seu_CNS1)
top_markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.8, p_val_adj < 0.05, pct.1 > 0.5)

####
## Manual annotation, ignore Final ####
####

# cluster 0 VLMC (DCN+)
# cluster 1 Astocytes (GFAP+)
# cluster 2 Oligodendrocytes
# cluster 3 microglia/macrophages_1
# cluster 4 microglia/macrophages_2
# cluster 5 endothelialcells_1
# cluster 6 Astorcytes (GJA1, AQP4, CHI3L1+)
# cluster 7 microglia/macrophages_3 
# cluster 8 T cells
# cluster 9 Astrocytes (? CHL1+)
# cluster 10 endothelialcells_2

# get final annotations
annotations <- c(
  "VLMC (DCN+)",
  "Astocytes (GFAP+)",
  "Oligodendrocytes",
  "microglia/macrophages_1",
  "microglia/macrophages_2",
  "endothelialcells_1",
  "Astorcytes (GJA1, AQP4, CHI3L1+)",
  "microglia/macrophages_3",
  "T cells",
  "Astrocytes (? CHL1+)",
  "endothelialcells_2"
)

# mark cell types in marker table
top_markers$FinalAnnotation <- annotations[match(as.numeric(as.character(top_markers$cluster))+1, 1:length(annotations))]

# mark cell types in VoltRon and Seurat object
brain_vr_CNS1$CellType <- annotations[match(as.numeric(as.character(brain_vr_CNS1$cluster))+1, 1:length(annotations))]
brain_seu_CNS1$CellType <- annotations[match(as.numeric(as.character(brain_vr_CNS1$cluster))+1, 1:length(annotations))]

# update umap
vrEmbeddings(brain_vr_CNS1, type = "umap", overwrite = TRUE) <- Embeddings(brain_seu_CNS1, reduction = "umap")

####
## Update annotation (Major Cell Type) ####
####

temp <- brain_vr_CNS1$CellType
temp[grep("endothelialcells", temp)] <- "Endothelial Cells"
temp[grepl("T cells", temp)] <- "Immune Cells"
temp[grep("VLMC", temp)] <- "VLMC/CAF"
temp[grep("macrophages", temp)] <- "TAMs"
temp[temp == "Astocytes (GFAP+)"] <- "Astrocytes"
temp[temp == "Astorcytes (GJA1, AQP4, CHI3L1+)"] <- "Astrocytes"
temp[temp == "Astrocytes (? CHL1+)"] <- "Astrocytes"
brain_vr_CNS1$MajorCellType <- temp
# saveRDS(brain_vr_CNS1, file = "../data/VoltRonData/brain_vr_CNS1_annotated.rds")

####
# save as AnnData ####
####

# interactive TissUUmaps tool
as.AnnData(brain_vr_CNS1, file = "../data/AnnDataData/brain_anndata_CNS1.h5ad", flip_coordinates = TRUE, 
           assay = "Xenium", method = "anndata", name = "DAPI", channel = "H&E")
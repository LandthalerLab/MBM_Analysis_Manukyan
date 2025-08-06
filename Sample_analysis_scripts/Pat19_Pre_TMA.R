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
brain_vr_Pat19_Pre_TMA <- subset(brain_vr, assays = c("Assay13"))

# get image data and register
brain_image_Pat19_Pre_TMA  <- importImageData(image = "../data/images/H&E/TMA.TIFF/Xenium-1.tif", tile.size = 10)

# register and save images, use 2022 vs 1992
xen_reg1 <- registerSpatialData(object_list = list(brain_vr_Pat19_Pre_TMA, brain_image_Pat19_Pre_TMA))
vrImages(brain_vr_Pat19_Pre_TMA[["Assay13"]], name = "DAPI", channel = "H&E") <- vrImages(xen_reg1$registered_spat[[2]], name = "main_reg")

####
# clustering ####
####

# clustering with Seurat
brain_seu_Pat19_Pre_TMA <- VoltRon::as.Seurat(brain_vr_Pat19_Pre_TMA, cell.assay = "Xenium", type = "image")
brain_seu_Pat19_Pre_TMA <- NormalizeData(brain_seu_Pat19_Pre_TMA, scale.factor = 100)

# dimensionality reduction
brain_seu_Pat19_Pre_TMA <- ScaleData(brain_seu_Pat19_Pre_TMA)
brain_seu_Pat19_Pre_TMA <- RunPCA(brain_seu_Pat19_Pre_TMA, features = Features(brain_seu_Pat19_Pre_TMA), npcs = 30)
brain_seu_Pat19_Pre_TMA <- RunUMAP(brain_seu_Pat19_Pre_TMA, dims = 1:10)
DimPlot(brain_seu_Pat19_Pre_TMA, reduction = "umap", group.by = "Sample", split.by = "Sample")

# clustering
brain_seu_Pat19_Pre_TMA <- FindNeighbors(brain_seu_Pat19_Pre_TMA, dims = 1:10)
brain_seu_Pat19_Pre_TMA <- FindClusters(brain_seu_Pat19_Pre_TMA, resolution = c(0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,1.0,1.1))

# visualize
g_list <- list()
g_list[["Xenium_snn_res.0.1"]] <- DimPlot(brain_seu_Pat19_Pre_TMA, reduction = "umap", group.by = "Xenium_snn_res.0.1", label = T)
g_list[["Xenium_snn_res.0.2"]] <- DimPlot(brain_seu_Pat19_Pre_TMA, reduction = "umap", group.by = "Xenium_snn_res.0.2", label = T)
g_list[["Xenium_snn_res.0.3"]] <- DimPlot(brain_seu_Pat19_Pre_TMA, reduction = "umap", group.by = "Xenium_snn_res.0.3", label = T)
g_list[["Xenium_snn_res.0.4"]] <- DimPlot(brain_seu_Pat19_Pre_TMA, reduction = "umap", group.by = "Xenium_snn_res.0.4", label = T)
g_list[["Xenium_snn_res.0.5"]] <- DimPlot(brain_seu_Pat19_Pre_TMA, reduction = "umap", group.by = "Xenium_snn_res.0.5", label = T)
g_list[["Xenium_snn_res.0.6"]] <- DimPlot(brain_seu_Pat19_Pre_TMA, reduction = "umap", group.by = "Xenium_snn_res.0.6", label = T)
g_list[["Xenium_snn_res.0.7"]] <- DimPlot(brain_seu_Pat19_Pre_TMA, reduction = "umap", group.by = "Xenium_snn_res.0.7", label = T)
g_list[["Xenium_snn_res.1"]] <- DimPlot(brain_seu_Pat19_Pre_TMA, reduction = "umap", group.by = "Xenium_snn_res.1", label = T)
g_list[["Xenium_snn_res.1.1"]] <- DimPlot(brain_seu_Pat19_Pre_TMA, reduction = "umap", group.by = "Xenium_snn_res.1.1", label = T)
ggarrange(plotlist = g_list, ncol = 3, nrow = 3)

# res 0.5 looks good
brain_vr_Pat19_Pre_TMA$cluster <- as.character(brain_seu_Pat19_Pre_TMA$Xenium_snn_res.0.5)

####
# marker analysis ####
####

####
## check resolution 0.5 ####
####

# marker analysis
brain_seu_Pat19_Pre_TMA <- as.Seurat(brain_vr_Pat19_Pre_TMA, type = "image", cell.assay = "Xenium")
Idents(brain_seu_Pat19_Pre_TMA) <- "Xenium_snn_res.0.5"
markers <- FindAllMarkers(brain_seu_Pat19_Pre_TMA)
top_markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.8, p_val_adj < 0.05, pct.1 > 0.3)

####
# subclustering of res 0.5 clus 0,1,2,7 ####
####

brain_seu_Pat19_Pre_TMA_subset <- subset(brain_seu_Pat19_Pre_TMA, subset = Xenium_snn_res.0.5 %in% c("0", "1", "2", "7"))
brain_seu_Pat19_Pre_TMA_subset <- NormalizeData(brain_seu_Pat19_Pre_TMA_subset, scale.factor = 100)
brain_seu_Pat19_Pre_TMA_subset <- ScaleData(brain_seu_Pat19_Pre_TMA_subset)
brain_seu_Pat19_Pre_TMA_subset <- RunPCA(brain_seu_Pat19_Pre_TMA_subset, features = Features(brain_seu_Pat19_Pre_TMA_subset), npcs = 30)
ElbowPlot(brain_seu_Pat19_Pre_TMA_subset)
brain_seu_Pat19_Pre_TMA_subset <- RunUMAP(brain_seu_Pat19_Pre_TMA_subset, dims = 1:8)
DimPlot(brain_seu_Pat19_Pre_TMA_subset, reduction = "umap", group.by = "Sample", split.by = "Sample")
brain_seu_Pat19_Pre_TMA_subset <- FindNeighbors(brain_seu_Pat19_Pre_TMA_subset, dims = 1:8)
brain_seu_Pat19_Pre_TMA_subset <- FindClusters(brain_seu_Pat19_Pre_TMA_subset, resolution = c(0.1, 0.2, 0.3, 0.4,0.5,0.6))

# visualize
g_list <- list()
g_list[["Xenium_snn_res.0.1"]] <- DimPlot(brain_seu_Pat19_Pre_TMA_subset, reduction = "umap", group.by = "Xenium_snn_res.0.1", label = T)
g_list[["Xenium_snn_res.0.2"]] <- DimPlot(brain_seu_Pat19_Pre_TMA_subset, reduction = "umap", group.by = "Xenium_snn_res.0.2", label = T)
g_list[["Xenium_snn_res.0.3"]] <- DimPlot(brain_seu_Pat19_Pre_TMA_subset, reduction = "umap", group.by = "Xenium_snn_res.0.3", label = T)
g_list[["Xenium_snn_res.0.4"]] <- DimPlot(brain_seu_Pat19_Pre_TMA_subset, reduction = "umap", group.by = "Xenium_snn_res.0.4", label = T)
g_list[["Xenium_snn_res.0.5"]] <- DimPlot(brain_seu_Pat19_Pre_TMA_subset, reduction = "umap", group.by = "Xenium_snn_res.0.5", label = T)
g_list[["Xenium_snn_res.0.6"]] <- DimPlot(brain_seu_Pat19_Pre_TMA_subset, reduction = "umap", group.by = "Xenium_snn_res.0.6", label = T)
ggarrange(plotlist = g_list, ncol = 3, nrow = 2)

####
## check resolution 0.5 ####
####

# marker analysis
Idents(brain_seu_Pat19_Pre_TMA_subset) <- "Xenium_snn_res.0.5"
markers <- FindAllMarkers(brain_seu_Pat19_Pre_TMA_subset)
top_markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.8, p_val_adj < 0.05, pct.1 > 0.3)

# investigate clusters
g1 <- FeaturePlot(brain_seu_Pat19_Pre_TMA_subset, features = c("ABCB5", "TOP2A", "ERBB3"), reduction = "umap")
g2 <- DimPlot(brain_seu_Pat19_Pre_TMA_subset, reduction = "umap", group.by = "Xenium_snn_res.0.5", label = TRUE)
g1 | g2

g1 <- FeaturePlot(brain_seu_Pat19_Pre_TMA, features = c("ABCB5", "TOP2A", "ERBB3"), reduction = "umap")
g2 <- DimPlot(brain_seu_Pat19_Pre_TMA, reduction = "umap", group.by = "Xenium_snn_res.0.5", label = TRUE)
g1 | g2

# appoint clusters
temp <- as.character(brain_seu_Pat19_Pre_TMA$Xenium_snn_res.0.5)
names(temp) <- names(brain_seu_Pat19_Pre_TMA$Xenium_snn_res.0.5)
temp[names(brain_seu_Pat19_Pre_TMA_subset$Xenium_snn_res.0.5)] <- paste0("subcluster_", as.character(brain_seu_Pat19_Pre_TMA_subset$Xenium_snn_res.0.5))
brain_seu_Pat19_Pre_TMA$combined_clusters <- temp
DimPlot(brain_seu_Pat19_Pre_TMA, reduction = "umap", group.by = "combined_clusters", label = TRUE)
VlnPlot(brain_seu_Pat19_Pre_TMA, features = "nCount_Xenium", group.by = "combined_clusters")

# save files
brain_vr_Pat19_Pre_TMA$cluster <- brain_seu_Pat19_Pre_TMA$combined_clusters

####
# Manual annotation ####
####

####
## Update annotation ####
####

# cluster subcluster_0 undefined
# cluster subcluster_1 Tumorcells_1 (ABCB5, ERRB3)
# cluster subcluster_2 Tumorcells_2 (ABCB5)
# cluster subcluster_3 Tumorcells_1 (ABCB5, ERRB3)
# cluster subcluster_4 Tumorcells_3 (TOP2A, MKI67, CENPF)
# cluster 3 TAMs_1
# cluster 4 endothelial cells
# cluster 5 MX1+ IFI6+ IFI44L+ cells
# cluster 6 TAMs_2

# get final annotations
annotations <- list(`subcluster_0` = "undefined", 
                    `subcluster_1` = "Tumorcells_1 (ABCB5, ERRB3)",
                    `subcluster_2` = "Tumorcells_2 (ABCB5)",
                    `subcluster_3` = "Tumorcells_1 (ABCB5, ERRB3)",
                    `subcluster_4` = "Tumorcells_3 (TOP2A, MKI67, CENPF)",
                    `3` = "TAMs_1",
                    `4` = "endothelial cells",
                    `5` = "MX1+ IFI6+ IFI44L+ cells",
                    `6` = "TAMs_2")

# mark cell types in VoltRon and Seurat object
brain_vr_Pat19_Pre_TMA$CellType <- unlist(annotations[brain_vr_Pat19_Pre_TMA$cluster])
temp <- unlist(annotations[brain_seu_Pat19_Pre_TMA$combined_clusters])
names(temp) <- NULL
brain_seu_Pat19_Pre_TMA$CellType <- temp

# update umap
vrEmbeddings(brain_vr_Pat19_Pre_TMA, type = "umap", overwrite = TRUE) <- Embeddings(brain_seu_Pat19_Pre_TMA, reduction = "umap")

####
## Update annotation (Major Cell Type) ####
####

temp <- brain_vr_Pat19_Pre_TMA$CellType
temp[grep("endothelial", temp)] <- "Endothelial Cells"
temp[grep("TAMs", temp)] <- "TAMs"
temp[temp == "MX1+ IFI6+ IFI44L+ cells"] <- "IFN Cells"
temp[temp == "Tumorcells_3 (TOP2A, MKI67, CENPF)"] <- "Proliferating Tumor Cells 1"
temp[grepl("Tumorcells_1", temp)] <- "Tumor Cells 1"
temp[grepl("Tumorcells_2", temp)] <- "Tumor Cells 2"
brain_vr_Pat19_Pre_TMA$MajorCellType <- temp
# saveRDS(brain_vr_Pat19_Pre_TMA, file = "../data/VoltRonData/brain_vr_Pat19_Pre_TMA_annotated.rds")

####
# annotate Hot tumors ####
####

####
## get immunescores ####
####

# get score per modules using Seurat
brain_vr_Pat19_Pre_TMA_assay1_seu <- VoltRon::as.Seurat(brain_vr_Pat19_Pre_TMA, cell.assay = "Xenium", type = "image")
brain_vr_Pat19_Pre_TMA_assay1_seu <- NormalizeData(brain_vr_Pat19_Pre_TMA_assay1_seu, scale.factor = 1000)
for(i in 1:length(score_markers)){
  cur_markers <- intersect(score_markers[[i]], SeuratObject::Features(brain_vr_Pat19_Pre_TMA_assay1_seu))
  brain_vr_Pat19_Pre_TMA_assay1_seu <- AddModuleScore(brain_vr_Pat19_Pre_TMA_assay1_seu, features = list(cur_markers), name = names(score_markers)[i], ctrl = 5)  
}

# transfer module scores to voltron objects
temp <- rep(NA, length(vrSpatialPoints(brain_vr_Pat19_Pre_TMA)))
names(temp) <- vrSpatialPoints(brain_vr_Pat19_Pre_TMA)
temp[names(brain_vr_Pat19_Pre_TMA_assay1_seu$ImmuneScores1)] <- brain_vr_Pat19_Pre_TMA_assay1_seu$ImmuneScores1
brain_vr_Pat19_Pre_TMA$ImmuneScores <- temp

# visualize score
vrSpatialFeaturePlot(brain_vr_Pat19_Pre_TMA, assay = "Xenium", features = c("ImmuneScores"),          
                     plot.segments = FALSE, pt.size = 1, n.tile = 300, alpha = 1)

####
## annotate hot tumor ####  
####

brain_vr_Pat19_Pre_TMA <- annotateSpatialData(brain_vr_Pat19_Pre_TMA, use.image = TRUE, channel = "H&E")
vrSpatialPlot(brain_vr_Pat19_Pre_TMA, group.by = "annotation", n.tile = 300)

####
# save as AnnData ####
####

# interactive TissUUmaps tool
as.AnnData(brain_vr_Pat19_Pre_TMA, file = "../data/AnnDataData/brain_anndata_Pat19_Pre_TMA.h5ad", flip_coordinates = TRUE, 
           name = "DAPI", channel = "H&E")
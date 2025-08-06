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

# xenium annotation
xenium_markers <- read.xlsx("../../RadkeRedmerSpatial/data/Supplementary Information/xenium_brain_annotation.xlsx", sheetIndex = 1)

####
# overlay images ####
####

# get clustered data
brain_vr_Pat3 <- subset(brain_vr, assays = c("Assay2"))

# get image data and register
brain_image_Pat3  <- importImageData(image = "../data/images/H&E/TMA.TIFF/Xenium-8.tif", tile.size = 10)

# register and save images, use default size
xen_reg1 <- registerSpatialData(object_list = list(brain_vr_Pat3, brain_image_Pat3))
vrImages(brain_vr_Pat3[["Assay2"]], name = "DAPI", channel = "H&E") <- vrImages(xen_reg1$registered_spat[[2]], name = "main_reg")

# save data
saveRDS(brain_vr_Pat3, file = "../data/VoltRonData/brain_vr_Pat3_withHE.rds")

####
# clustering ####
####

# read registered VoltRon objects
# brain_vr_Pat3 <- readRDS(file = "../data/VoltRonData/brain_vr_Pat3_withHE.rds")

# clustering with Seurat
brain_seu_Pat3 <- VoltRon::as.Seurat(brain_vr_Pat3, cell.assay = "Xenium", type = "image")
brain_seu_Pat3 <- NormalizeData(brain_seu_Pat3, scale.factor = 100)

# dimensionality reduction
brain_seu_Pat3 <- ScaleData(brain_seu_Pat3)
brain_seu_Pat3 <- RunPCA(brain_seu_Pat3, features = Features(brain_seu_Pat3), npcs = 20)
brain_seu_Pat3 <- RunUMAP(brain_seu_Pat3, dims = 1:20)
DimPlot(brain_seu_Pat3, reduction = "umap", group.by = "Sample", split.by = "Sample")

# clustering
brain_seu_Pat3 <- FindNeighbors(brain_seu_Pat3, dims = 1:20)
brain_seu_Pat3 <- FindClusters(brain_seu_Pat3, resolution = c(0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8, 0.9, 1.0,1.1))

# visualize
g_list <- list()
g_list[["Xenium_snn_res.0.1"]] <- DimPlot(brain_seu_Pat3, reduction = "umap", group.by = "Xenium_snn_res.0.1", label = T)
g_list[["Xenium_snn_res.0.2"]] <- DimPlot(brain_seu_Pat3, reduction = "umap", group.by = "Xenium_snn_res.0.2", label = T)
g_list[["Xenium_snn_res.0.3"]] <- DimPlot(brain_seu_Pat3, reduction = "umap", group.by = "Xenium_snn_res.0.3", label = T)
g_list[["Xenium_snn_res.0.4"]] <- DimPlot(brain_seu_Pat3, reduction = "umap", group.by = "Xenium_snn_res.0.4", label = T)
g_list[["Xenium_snn_res.0.5"]] <- DimPlot(brain_seu_Pat3, reduction = "umap", group.by = "Xenium_snn_res.0.5", label = T)
g_list[["Xenium_snn_res.0.6"]] <- DimPlot(brain_seu_Pat3, reduction = "umap", group.by = "Xenium_snn_res.0.6", label = T)
g_list[["Xenium_snn_res.0.7"]] <- DimPlot(brain_seu_Pat3, reduction = "umap", group.by = "Xenium_snn_res.0.7", label = T)
g_list[["Xenium_snn_res.0.8"]] <- DimPlot(brain_seu_Pat3, reduction = "umap", group.by = "Xenium_snn_res.0.8", label = T)
g_list[["Xenium_snn_res.0.9"]] <- DimPlot(brain_seu_Pat3, reduction = "umap", group.by = "Xenium_snn_res.0.9", label = T)
g_list[["Xenium_snn_res.1"]] <- DimPlot(brain_seu_Pat3, reduction = "umap", group.by = "Xenium_snn_res.1", label = T)
g_list[["Xenium_snn_res.1.1"]] <- DimPlot(brain_seu_Pat3, reduction = "umap", group.by = "Xenium_snn_res.1.1", label = T)
ggarrange(plotlist = g_list, ncol = 3, nrow = 4)

# res 1.1 looks good
brain_vr_Pat3$cluster <- as.character(brain_seu_Pat3$Xenium_snn_res.1.1)

####
# marker analysis ####
####

####
## check resolution 1.1 ####
####

# marker analysis
brain_seu_Pat3 <- as.Seurat(brain_vr_Pat3, type = "image", cell.assay = "Xenium")
Idents(brain_seu_Pat3) <- "Xenium_snn_res.1.1"
markers <- FindAllMarkers(brain_seu_Pat3)
markers <- markers %>% left_join(xenium_markers, by = c("gene" = "Genes"))
top_markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.8, p_val_adj < 0.05, pct.1 > 0.5)

####
# annotation ####
####

####
## Manual annotation, ignore Final ####
####

# cluster 0 undefined, previously 0_mix (subcluster?)
# cluster 1 undefined (TOP2A/CENPF+), previously Tumorcells_1 (TOP2A/CENPF+)
# cluster 2 TAMs_1
# cluster 3 undefined (TOP2A/CENPF+), previously Tumorcells_1 (TOP2A/CENPF+)
# cluster 4 Astrocytes_1
# cluster 5 Astrocytes_2
# cluster 6 Oligodendrocytes
# cluster 7 Neurons_1
# cluster 8 endothelialcells
# cluster 9 Neurons_2
# cluster 10 Neurons_3
# cluster 11 Astrocytes_3
# cluster 12 immune cells
# cluster 13 Astrocytes_4
# cluster 14 VLMC (DCN+)
# cluster 15 OPC

# get final annotations
annotations <- c(
  "undefined", # 0_mix (subcluster?)", # Josi: highlights Tumor cells (malignant B-cells)
  "undefined", # "Tumorcells_1 (TOP2A/CENPF+)"
  "TAMs_1",
  "undefined", # "Tumorcells_1 (TOP2A/CENPF+)"
  "Astrocytes_1",
  "Astrocytes_2",
  "Oligodendrocytes",
  "Neurons_1",
  "endothelialcells",
  "Neurons_2",
  "Neurons_3",
  "Astrocytes_3",
  "immune cells",
  "Astrocytes_4",
  "VLMC (DCN+)",
  "OPC"
)

# mark cell types in marker table
top_markers$FinalAnnotation <- annotations[match(as.numeric(as.character(top_markers$cluster))+1, 1:length(annotations))]

# mark cell types in VoltRon and Seurat object
brain_vr_Pat3$CellType <- annotations[match(as.numeric(as.character(brain_vr_Pat3$cluster))+1, 1:length(annotations))]
brain_seu_Pat3$CellType <- annotations[match(as.numeric(as.character(brain_vr_Pat3$cluster))+1, 1:length(annotations))]

# update umap
vrEmbeddings(brain_vr_Pat3, type = "umap", overwrite = TRUE) <- Embeddings(brain_seu_Pat3, reduction = "umap")

####
## Update annotation (Major Cell Type) ####
####

temp <- brain_vr_Pat3$CellType
temp[temp %in% c("Astrocytes_2","Astrocytes_3","Astrocytes_4")] <- "Astrocytes"
temp[temp == "immune cells"] <- "Immune Cells"
temp[grep("endothelialcells", temp)] <- "Endothelial Cells"
temp[temp == "Tumorcells_1 (TOP2A/CENPF+)"] <- "Proliferating Tumor Cells 1"
# temp[grep("TAMs_1|0_mix", temp)] <- "TAMs"
temp[grep("TAMs_1", temp)] <- "TAMs"
temp[grep("0_mix", temp)] <- "undefined"
temp[grep("VLMC", temp)] <- "VLMC/CAF"
temp[grep("Neurons", temp)] <- "Neurons"
temp[temp == "OPC"] <- "Oligodendrocytes"
brain_vr_Pat3$MajorCellType <- temp
# saveRDS(brain_vr_Pat3, file = "../data/VoltRonData/brain_vr_Pat3_annotated.rds")

####
# annotate Hot tumors ####
####

####
## get immunescores ####
####

# get score per modules using Seurat
brain_vr_Pat3_assay1 <- brain_vr_Pat3
brain_vr_Pat3_assay1_seu <- VoltRon::as.Seurat(brain_vr_Pat3_assay1, cell.assay = "Xenium", type = "image")
brain_vr_Pat3_assay1_seu <- NormalizeData(brain_vr_Pat3_assay1_seu, scale.factor = 1000)
for(i in 1:length(score_markers)){
  cur_markers <- intersect(score_markers[[i]], SeuratObject::Features(brain_vr_Pat3_assay1_seu))
  brain_vr_Pat3_assay1_seu <- AddModuleScore(brain_vr_Pat3_assay1_seu, features = list(cur_markers), name = names(score_markers)[i], ctrl = 5)  
}

# transfer module scores to voltron objects
temp <- rep(NA, length(vrSpatialPoints(brain_vr_Pat3)))
names(temp) <- vrSpatialPoints(brain_vr_Pat3)
temp[names(brain_vr_Pat3_assay1_seu$ImmuneScores1)] <- brain_vr_Pat3_assay1_seu$ImmuneScores1
brain_vr_Pat3$ImmuneScores <- temp

# visualize score
vrSpatialFeaturePlot(brain_vr_Pat3, features = c("ImmuneScores"),          
                     plot.segments = FALSE, pt.size = 1, n.tile = 300, alpha = 1, spatial = "DAPI", channel = "H&E")

####
## annotate hot tumor ####
####

brain_vr_Pat3 <- annotateSpatialData(brain_vr_Pat3, use.image = TRUE, channel = "H&E")
# saveRDS(brain_vr_Pat3, file = "../data/VoltRonData/brain_vr_Pat3_annotated.rds")

####
# save as AnnData ####
####

# interactive TissUUmaps tool
as.AnnData(brain_vr_Pat3, file = "../data/AnnDataData/brain_anndata_Pat3.h5ad", flip_coordinates = TRUE, 
           assay = "Xenium", method = "anndata", name = "DAPI", channel = "H&E")

# interactive TissUUmaps tool
brain_vr_Pat3$subcluster <- brain_seu_Pat3$subcluster
as.AnnData(brain_vr_Pat3, file = "../data/AnnDataData/brain_anndata_Pat3_temp.h5ad", flip_coordinates = TRUE, 
           assay = "Xenium", method = "anndata", name = "DAPI", channel = "H&E")


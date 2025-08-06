library(VoltRon)
library(Seurat)
library(harmony)
library(dplyr)
library(ComplexHeatmap)
library(xlsx)
library(patchwork)
library(uwot)

####
# Import data ####
####

brain_vr <- readRDS("../data/VoltRonData/brain_all_merged.rds")

# get subset
brain_vr_Pat23 <- subset(brain_vr, assays = c("Assay6", "Assay7", "Assay8"))

# score markers (Interferon, Immune etc.)
score_markers <- read.xlsx("../data/Supplementary Information/score_markers.xlsx", sheetName = "Sheet1")
score_markers <- as.list(score_markers)
score_markers <- lapply(score_markers, function(x) x[!is.na(x)]) 

####
# overlay images ####
####

# get image data and register
brain_image_Pat23_A <- importImageData(image = "../data/images/H&E/TMA.TIFF/Xenium-15.tif", tile.size = 10)
xen_reg <- registerSpatialData(object_list = list(subset(brain_vr_Pat23, samples = "Pat23_19_1A_A"),
                                                  brain_image_Pat23_A), keypoints = readRDS("../data/AuxiliaryData/"))
vrImages(brain_vr_Pat23[["Assay8"]], name = "DAPI", channel = "H&E") <- vrImages(xen_reg$registered_spat[[2]], name = "image_1_reg")
vrSpatialPlot(brain_vr_Pat23, assay = "Assay8", plot.segments = TRUE, group.by = "Sample", ncol = 2, background = c("DAPI", "H&E"))
brain_vr_Pat23_A <- subset(brain_vr_Pat23, assays = "Assay8")
saveRDS(brain_vr_Pat23_A, file = "../data/VoltRonData/brain_vr_Pat23_A_withHE.rds")

brain_image_Pat23_B <- importImageData(image = "../data/images/H&E/TMA.TIFF/Xenium-11.tif", tile.size = 10)
xen_reg <- registerSpatialData(object_list = list(subset(brain_vr_Pat23, samples = "Pat23_19_1A_B"),
                                                  brain_image_Pat23_B))
vrImages(brain_vr_Pat23[["Assay7"]], name = "DAPI", channel = "H&E") <- vrImages(xen_reg$registered_spat[[2]], name = "image_1_reg")
vrSpatialPlot(brain_vr_Pat23, assay = "Assay7", plot.segments = TRUE, group.by = "Sample", ncol = 2, background = c("DAPI", "H&E"))

brain_image_Pat23_C <- importImageData(image = "../data/images/H&E/TMA.TIFF/Xenium-7.tif", tile.size = 10)
xen_reg <- registerSpatialData(object_list = list(subset(brain_vr_Pat23, samples = "Pat23_19_1A_C"),
                                                  brain_image_Pat23_C))
vrImages(brain_vr_Pat23[["Assay6"]], name = "DAPI", channel = "H&E") <- vrImages(xen_reg$registered_spat[[2]], name = "image_1_reg")
vrSpatialPlot(brain_vr_Pat23, assay = "Assay6", plot.segments = TRUE, group.by = "Sample", ncol = 2, background = c("DAPI", "H&E"))

saveRDS(brain_vr_Pat23, file = "../data/VoltRonData/brain_vr_Pat23_withHE.rds")

####
# clustering ####
####

# clustering with Seurat
brain_seu_Pat23 <- as.Seurat(brain_vr_Pat23, cell.assay = "Xenium", type = "image")
brain_seu_Pat23 <- NormalizeData(brain_seu_Pat23, scale.factor = 100)

# dimensionality reduction
brain_seu_Pat23 <- ScaleData(brain_seu_Pat23)
selected_features <- Features(brain_seu_Pat23)[!grepl("^S2-|^HSV1", Features(brain_seu_Pat23))]
brain_seu_Pat23 <- RunPCA(brain_seu_Pat23, features = selected_features, npcs = 20)
brain_seu_Pat23 <- RunUMAP(brain_seu_Pat23, dims = 1:20)
DimPlot(brain_seu_Pat23, reduction = "umap", group.by = "Sample")

# clustering
brain_seu_Pat23 <- FindNeighbors(brain_seu_Pat23, dims = 1:20)
brain_seu_Pat23 <- FindClusters(brain_seu_Pat23, resolution = c(0.6,0.8,1.0,1.1,1.2))

# visualize, looks like 1.1 is nice
g_list <- list()
g_list[["Xenium_snn_res.0.6"]] <- DimPlot(brain_seu_Pat23, reduction = "umap", group.by = "Xenium_snn_res.0.6", label = T)
g_list[["Xenium_snn_res.0.8"]] <- DimPlot(brain_seu_Pat23, reduction = "umap", group.by = "Xenium_snn_res.0.8", label = T)
g_list[["Xenium_snn_res.1"]] <- DimPlot(brain_seu_Pat23, reduction = "umap", group.by = "Xenium_snn_res.1", label = T)
g_list[["Xenium_snn_res.1.1"]] <- DimPlot(brain_seu_Pat23, reduction = "umap", group.by = "Xenium_snn_res.1.1", label = T)
g_list[["Xenium_snn_res.1.2"]] <- DimPlot(brain_seu_Pat23, reduction = "umap", group.by = "Xenium_snn_res.1.2", label = T)
ggarrange(plotlist = g_list, ncol = 3, nrow = 2)

# looks like res 1.1 is nice
DimPlot(brain_seu_Pat23, reduction = "umap", group.by = "Xenium_snn_res.1.1", label = T)
brain_vr_Pat23$cluster <- as.character(brain_seu_Pat23$Xenium_snn_res.1.1)

####
# marker analysis ####
####

####
## check resolution 1.1 ####
####

# marker analysis
Idents(brain_seu_Pat23) <- "Xenium_snn_res.1.1"
markers <- FindAllMarkers(brain_seu_Pat23, )
top_markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.8, p_val_adj < 0.05, pct.1 > 0.5)
unique_markers <- top_markers %>%
  dplyr::filter(avg_log2FC > 1) %>%
  group_by(gene) %>% mutate(n_gene = n()) %>%
  filter(n_gene < 2) %>% select(-n_gene)
top_markers$is.unique <- ifelse(top_markers$gene %in% unique_markers$gene, "Yes", "No")

####
# annotation ####
####

####
## Manual annotation, ignore Final ####
####

# cluster 0 Tumorcells_1 (SOX11+)
# cluster 1 Tumorcells_2 (SOX9+)
# cluster 2 Astrocytes (ELOVL2,TRIL,CHL1+)
# cluster 3 Tumorcells_3 (IGFBP3,IGFBP5+)
# cluster 4 Neurons (RELN+)
# cluster 5 5_Mix? (cluster again)
# cluster 6 6_Mix? (cluster again)
# cluster 7 Reactive glia (CD163+)
# cluster 8 Tumorcells_4 (SERPINA3+)
# cluster 9 Tumorcells_5 (CCNB2, TOP2A, MKI67+)
# cluster 10 Astrocytes (GLIS3+?)
# cluster 11 TAMs_1 (CD68+)
# cluster 12 12? (could be trash)
# cluster 13 endothelialcells_1
# cluster 14 immunecells
# cluster 15 Tumorcells_6 (SERPINA3+)
# cluster 16 endothelialcells_2
# cluster 17 ISG15+ OAS3+ cells
# cluster 18 VLMC (DCN+)
# cluster 19 TAMs_2 (ITGAX+)
# cluster 20 Oligodendroytes (ERMN, MAG)
# cluster 21 Astrocytes (AQP4+)

# get final annotations
annotations <- c(
   "Tumorcells_1 (SOX11+)", "Tumorcells_2 (SOX9+)", "Astrocytes (ELOVL2,TRIL,CHL1+)",
   "Tumorcells_3 (IGFBP3,IGFBP5+)", "Neurons (RELN+)", "5_Mix? (cluster again)", "6_Mix? (cluster again)",
   "Reactive glia (CD163+)", "Tumorcells_4 (SERPINA3+)", "Tumorcells_5 (CCNB2, TOP2A, MKI67+)",
   "Astrocytes (GLIS3+?)", "TAMs_1 (CD68+)", "12? (could be trash)",
   "endothelialcells_1", "immunecells", "Tumorcells_6 (SERPINA3+)", "endothelialcells_2",
   "ISG15+ OAS3+ cells", "VLMC (DCN+)", "TAMs_2 (ITGAX+)", "Oligodendroytes (ERMN, MAG)", "Astrocytes (AQP4+)"
)

# mark cell types in marker table
top_markers$FinalAnnotation <- annotations[match(as.numeric(as.character(top_markers$cluster))+1, 1:length(annotations))]

# mark cell types in VoltRon and Seurat object
brain_vr_Pat23$CellType <- annotations[match(as.numeric(brain_vr_Pat23$cluster)+1, 1:length(annotations))]
brain_seu_Pat23$CellType <- annotations[match(as.numeric(brain_vr_Pat23$cluster)+1, 1:length(annotations))]

####
## Update Annotation ####
####

# update tumor annotation
top_markers <- read.xlsx(file = "../data/MarkerData/brain_seu_Pat23_marker_analysis.xlsx", sheetName = "results")

oldnames <- annotations[grepl("Tumorcells", annotations)]
newnames <- c("Tumorcells_1 (MKI67/TOP2A/CENPF+)", "SOX9+ cells (tumor?)", "CAV1/AXL+ cells (tumor?)", "SERPINA3+ cells",
              "Tumorcells_5 (MKI67/TOP2A/CENPF/CCNB2/GAS2L3+)", "Tumorcells_6 (GAS2L3+)")
brain_vr_Pat23$CellType_TumorAnnotation <- brain_vr_Pat23$CellType
brain_seu_Pat23$CellType_TumorAnnotation <- brain_seu_Pat23$CellType
for(i in 1:length(oldnames)){
  brain_vr_Pat23$CellType_TumorAnnotation[brain_vr_Pat23$CellType_TumorAnnotation == oldnames[i]] <- newnames[i]
  brain_seu_Pat23$CellType_TumorAnnotation[brain_vr_Pat23$CellType_TumorAnnotation == oldnames[i]] <- newnames[i]
}

# update umap
vrEmbeddings(brain_vr_Pat23, type = "umap", overwrite = TRUE) <- Embeddings(brain_seu_Pat23, reduction = "umap")
vrEmbeddings(brain_vr_Pat23, type = "pca_seurat", overwrite = TRUE) <- Embeddings(brain_seu_Pat23, reduction = "pca")

####
## Update annotation (Major Cell Type) ####
####

temp <- brain_vr_Pat23$CellType_TumorAnnotation
temp[temp == "ISG15+ OAS3+ cells"] <- "IFN Cells"
temp[temp == "immunecells"] <- "Immune Cells"
temp[temp == "Astrocytes (AQP4+)"] <- "Astrocytes"
temp[temp == "Astrocytes (GLIS3+?)"] <- "Astrocytes" # 
temp[temp == "Oligodendroytes (ERMN, MAG)"] <- "Oligodendrocytes"
temp[grep("TAMs", temp)] <- "TAMs"
temp[temp == "Reactive glia (CD163+)"] <- "TAMs"
temp[grep("endothelialcells", temp)] <- "Endothelial Cells"
temp[grep("VLMC", temp)] <- "VLMC/CAF"
temp[temp == "Tumorcells_5 (MKI67/TOP2A/CENPF/CCNB2/GAS2L3+)"] <- "Proliferating Tumor Cells 1"
temp[temp == "Tumorcells_1 (MKI67/TOP2A/CENPF+)"] <- "Proliferating Tumor Cells 2"
temp[temp == "Tumorcells_6 (GAS2L3+)"] <- "Tumor Cells 1"
temp[temp == "SOX9+ cells (tumor?)"] <- "Tumor Cells 2"
temp[temp == "SERPINA3+ cells"] <- "Tumor Cells 3"
temp[temp == "CAV1/AXL+ cells (tumor?)"] <- "Tumor Cells 4"
temp[temp == "Neurons (RELN+)"] <- "Tumor Cells 5"
temp[temp == "Astrocytes (ELOVL2,TRIL,CHL1+)"] <- "Tumor Cells 6"
temp[temp == "12? (could be trash)"] <- "Tumor Cells 7" # These were originally TAMs, see if you need change them later
temp[temp == "6_Mix? (cluster again)"] <- "Tumor Cells 8" # These were originally TAMs, see if you need change them later
temp[temp == "5_Mix? (cluster again)"] <- "Low Count Cells"
brain_vr_Pat23$MajorCellType <- temp
# saveRDS(brain_vr_Pat23, file = "../data/VoltRonData/brain_vr_Pat23_annotated.rds")

####
# Sample Embeddings ####
####

sample_names <- unique(brain_vr_Pat23$Sample)
brain_vr_Pat23_xenium <- subset(brain_vr_Pat23, assay = "Xenium")
`%notin%` <- Negate("%in%")
brain_vr_Pat23_xenium <- subset(brain_vr_Pat23_xenium, MajorCellType %notin% "Low Count Cells")
for(samp in sample_names){
  print(samp)
  brain_vr_Pat23_sub <- subset(brain_vr_Pat23_xenium, samples = samp)
  brain_vr_Pat23_subseu <- as.Seurat(brain_vr_Pat23_sub, cell.assay = "Xenium", type = "image")
  brain_vr_Pat23_subseu <- NormalizeData(brain_vr_Pat23_subseu, scale.factor = 100)
  brain_vr_Pat23_subseu <- ScaleData(brain_vr_Pat23_subseu)
  selected_features <- Features(brain_vr_Pat23_subseu)[!grepl("^S2-|^HSV1", Features(brain_vr_Pat23_subseu))]
  brain_vr_Pat23_subseu <- RunPCA(brain_vr_Pat23_subseu, features = selected_features, npcs = 20)
  brain_vr_Pat23_subseu <- RunUMAP(brain_vr_Pat23_subseu, dims = 1:20)
  vrEmbeddings(brain_vr_Pat23, type = paste0("umap_", samp), overwrite = TRUE) <- Embeddings(brain_vr_Pat23_subseu, reduction = "umap")
}

vrEmbeddingPlot(brain_vr_Pat23, embedding = "umap_Pat23_19_1A_C", group.by = "MajorCellType")

####
# annotate Hot tumors ####
####

####
## get immunescores ####
####

# get score per modules using Seurat
brain_vr_Pat23_assay1_seu <- VoltRon::as.Seurat(brain_vr_Pat23, cell.assay = "Xenium", type = "image")
brain_vr_Pat23_assay1_seu <- NormalizeData(brain_vr_Pat23_assay1_seu, scale.factor = 1000)
for(i in 1:length(score_markers)){
  cur_markers <- intersect(score_markers[[i]], SeuratObject::Features(brain_vr_Pat23_assay1_seu))
  brain_vr_Pat23_assay1_seu <- AddModuleScore(brain_vr_Pat23_assay1_seu, features = list(cur_markers), name = names(score_markers)[i], ctrl = 5)  
}

# transfer module scores to voltron objects
temp <- rep(NA, length(vrSpatialPoints(brain_vr_Pat23)))
names(temp) <- vrSpatialPoints(brain_vr_Pat23)
temp[names(brain_vr_Pat23_assay1_seu$ImmuneScores1)] <- brain_vr_Pat23_assay1_seu$ImmuneScores1
brain_vr_Pat23$ImmuneScores <- temp

# visualize score and annotations
g1 <- vrSpatialFeaturePlot(brain_vr_Pat23, assay = "Xenium", features = c("ImmuneScores"),          
                     plot.segments = FALSE, pt.size = 1, n.tile = 200, alpha = 1, ncol = 3)
vrMainAssay(brain_vr_Pat23) <- "ROIAnnotation"
brain_vr_Pat23$Annotation <- "hot"
vrImages(brain_vr_Pat23[["Assay9"]], channel = "H&E") <- vrImages(brain_vr_Pat23, assay = "Assay6", channel = "H&E")
vrImages(brain_vr_Pat23[["Assay10"]], channel = "H&E") <- vrImages(brain_vr_Pat23, assay = "Assay7", channel = "H&E")
vrImages(brain_vr_Pat23[["Assay11"]], channel = "H&E") <- vrImages(brain_vr_Pat23, assay = "Assay8", channel = "H&E")
g2 <- vrSpatialPlot(brain_vr_Pat23, assay = "ROIAnnotation", group.by = "Annotation",
                    plot.segments = TRUE, pt.size = 1, alpha = 1, ncol = 3, channel = "H&E")
g1 / g2

####
## annotate hot tumor ####  
####

brain_vr_Pat23 <- annotateSpatialData(brain_vr_Pat23, assay = "Assay6", use.image = TRUE, channel = "H&E")
vrSpatialPlot(brain_vr_Pat23, group.by = "annotation", n.tile = 300)

brain_vr_Pat23 <- annotateSpatialData(brain_vr_Pat23, assay = "Assay7", use.image = TRUE, channel = "H&E", annotation_assay = "ROIAnnotation")
brain_vr_Pat23 <- annotateSpatialData(brain_vr_Pat23, assay = "Assay8", use.image = TRUE, channel = "H&E", annotation_assay = "ROIAnnotation")
brain_vr_Pat23$annotation <- ifelse(grepl("hot", brain_vr_Pat23$annotation.1), brain_vr_Pat23$annotation.1, brain_vr_Pat23$annotation)
brain_vr_Pat23$annotation <- ifelse(grepl("hot", brain_vr_Pat23$annotation.2), brain_vr_Pat23$annotation.2, brain_vr_Pat23$annotation)

# saveRDS(brain_vr_Pat23, file = "../data/VoltRonData/brain_vr_Pat23_annotated.rds")

####
# save as AnnData ####
####

# interactive TissUUmaps tool
as.AnnData(brain_vr_Pat23, file = "../data/AnnDataData/brain_anndata_Pat23.h5ad", flip_coordinates = TRUE, 
           assay = "Xenium", method = "anndata", name = "DAPI", channel = "H&E")

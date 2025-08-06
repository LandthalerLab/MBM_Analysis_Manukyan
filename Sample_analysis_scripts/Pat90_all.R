library(VoltRon)
library(Seurat)
library(harmony)
library(dplyr)
library(ComplexHeatmap)
library(xlsx)

####
# Import data ####
####

# import (change this path to where the xenium output is)
Pat90_section <- importXenium("../../../data/RadkeRedmerSpatial_Brain_05112023/20240124__160424__Landthaler_human_brain_24012024/output-XETG00046__0020705__Region_1__20240124__161116/",
                              resolution_level = 2, overwrite_resolution = TRUE, sample_name = "Pat90")
Pat90_MBM_section <- importXenium("../../../data/RadkeRedmerSpatial_Brain_05112023/20240124__160424__Landthaler_human_brain_24012024/output-XETG00046__0020705__Region_2__20240124__161116/",
                              resolution_level = 2, overwrite_resolution = TRUE, sample_name = "Pat90_MBM")

# gene set score markers
score_markers <- read.xlsx("../data/Supplementary Information/score_markers.xlsx", sheetName = "Sheet1")
score_markers <- as.list(score_markers)
score_markers <- lapply(score_markers, function(x) x[!is.na(x)])

####
# overlay images ####
####

# width 
imgdata <- importImageData("../data/Images/H&E/sections/MBM_HE_1.2.24/Xenium ID0020705_sec2.tif")
xen_reg <- registerSpatialData(object_list = list(Pat90_section, imgdata))
vrImages(Pat90_section[["Assay1"]], name = "main", channel = "H&E") <- vrImages(xen_reg$registered_spat[[2]], name = "main_reg")

# width 3145 vs 3291
imgdata2 <- importImageData("../data/Images/H&E/sections/MBM_HE_1.2.24/Xenium ID0020705_sec1.tif")
xen_reg <- registerSpatialData(object_list = list(Pat90_MBM_section, imgdata2))
vrImages(Pat90_MBM_section[["Assay1"]], name = "main", channel = "H&E") <- vrImages(xen_reg$registered_spat[[2]], name = "main_reg")

####
# clustering ####
####

# merge samples 
brain_vr_Pat90_all <- merge(Pat90_section, Pat90_MBM_section)

# clustering with Seurat
brain_seu_Pat90_all <- VoltRon::as.Seurat(brain_vr_Pat90_all, cell.assay = "Xenium", type = "image")
brain_seu_Pat90_all <- subset(brain_seu_Pat90_all, subset = nCount_Xenium > 5)
brain_seu_Pat90_all <- NormalizeData(brain_seu_Pat90_all, scale.factor = 100)

# dimensionality reduction
brain_seu_Pat90_all <- ScaleData(brain_seu_Pat90_all)
brain_seu_Pat90_all <- RunPCA(brain_seu_Pat90_all, features = Features(brain_seu_Pat90_all), npcs = 30)
brain_seu_Pat90_all <- RunUMAP(brain_seu_Pat90_all, dims = 1:30)
DimPlot(brain_seu_Pat90_all, reduction = "umap", group.by = "Sample")
DimPlot(brain_seu_Pat90_all, reduction = "umap", group.by = "Sample", split.by = "Sample")

# clustering
brain_seu_Pat90_all <- FindNeighbors(brain_seu_Pat90_all, dims = 1:30)
brain_seu_Pat90_all <- FindClusters(brain_seu_Pat90_all, resolution = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3))

# visualize
g_list <- list()
g_list[["Xenium_snn_res.0.5"]] <- DimPlot(brain_seu_Pat90_all, reduction = "umap", group.by = "Xenium_snn_res.0.5", label = T) + NoLegend()
g_list[["Xenium_snn_res.0.6"]] <- DimPlot(brain_seu_Pat90_all, reduction = "umap", group.by = "Xenium_snn_res.0.6", label = T) + NoLegend()
g_list[["Xenium_snn_res.0.7"]] <- DimPlot(brain_seu_Pat90_all, reduction = "umap", group.by = "Xenium_snn_res.0.7", label = T) + NoLegend()
g_list[["Xenium_snn_res.0.8"]] <- DimPlot(brain_seu_Pat90_all, reduction = "umap", group.by = "Xenium_snn_res.0.8", label = T) + NoLegend()
g_list[["Xenium_snn_res.0.9"]] <- DimPlot(brain_seu_Pat90_all, reduction = "umap", group.by = "Xenium_snn_res.0.9", label = T) + NoLegend()
g_list[["Xenium_snn_res.1"]] <- DimPlot(brain_seu_Pat90_all, reduction = "umap", group.by = "Xenium_snn_res.1", label = T) + NoLegend()
g_list[["Xenium_snn_res.1.1"]] <- DimPlot(brain_seu_Pat90_all, reduction = "umap", group.by = "Xenium_snn_res.1.1", label = T) + NoLegend()
g_list[["Xenium_snn_res.1.2"]] <- DimPlot(brain_seu_Pat90_all, reduction = "umap", group.by = "Xenium_snn_res.1.2", label = T) + NoLegend()
g_list[["Xenium_snn_res.1.3"]] <- DimPlot(brain_seu_Pat90_all, reduction = "umap", group.by = "Xenium_snn_res.1.3", label = T) + NoLegend()
ggarrange(plotlist = g_list, ncol = 4, nrow = 3)

# clusters, res = 1
# get clusters to VoltRon
Clusters <- brain_seu_Pat90_all$Xenium_snn_res.1
Clusters <- as.character(Clusters[vrSpatialPoints(brain_vr_Pat90_all)])
Clusters[is.na(Clusters)] <- "undefined"
brain_vr_Pat90_all$Clusters <- Clusters

# save embeddings
vrEmbeddings(brain_vr_Pat90_all, type = "umap") <- Embeddings(brain_seu_Pat90_all, reduction = "umap")

####
# Marker Analysis ####
####

# Xenium_snn_res.1 looks good
Idents(brain_seu_Pat90_all) <- "Xenium_snn_res.1"
markers <- FindAllMarkers(brain_seu_Pat90_all)
top_markers <- markers %>% filter(pct.1 > 0.5, avg_log2FC > 0.5, p_val_adj < 0.05) %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 20)
unique_markers <- top_markers %>%
  dplyr::filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
  group_by(gene) %>% mutate(n_gene = n()) %>%
  filter(n_gene < 2) %>% select(-n_gene)

####
## Manual annotation ####
####

# cluster 0 Tumorcells_1 
# cluster 1 1? (TAP1+)
# cluster 2 SLC17A6+ cells
# cluster 3 3?
# cluster 4 T cells 1
# cluster 5 NES+ cells (necrotic?)
# cluster 6 Tumorcells_9 (BZWZ+)
# cluster 7 Neurons_1
# cluster 8 T cells 2 
# cluster 9 TAMs_1
# cluster 10 Tumorcells_10 (BZWZ+)
# cluster 11 Tumorcells_2
# cluster 12 Tumorcells_3
# cluster 13 13??
# cluster 14 endothelialcells
# cluster 15 endothelialcells (L_Endo?) 
# cluster 16 TAMs_2
# cluster 17 T cells 3
# cluster 18 Tumorcells_4 (MET+)
# cluster 19 Tumorcells_5 (LOX+ (%50), subcluster later)
# cluster 20 Astrocytes_1 (GFAP+)
# cluster 21 Tumorcells_6
# cluster 22 VLMC_1 (DCN+)
# cluster 23 23?
# cluster 24 IFI6+ IFI44L+ cells
# cluster 25 VLMC_2 (DCN+?)
# cluster 26 26?
# cluster 27 TAMs_3 (ITGAX+)
# cluster 28 Tumorcells_7 (MET+)
# cluster 29 endothelialcells + T cells (subcluster!)
# cluster 30 Tumorcells_8 (CDH1+)
# cluster 31 31? 
# cluster 32 32?
# cluster 33 PSENEN+ (Neuron?) cells
# cluster 34 VAMP1+ (Neuron?) cells
# cluster 35 PSEN1+ (Neuron?) cells
# cluster 36 SYNPR+ (Neuron?) cells
# cluster 37 PTPRZ1+ (Neuron?) cells
# cluster 38 undefined
# cluster 39 undefined

annotations <- c(
  "Tumorcells_1",
  "1? (TAP1+)",
  "SLC17A6+ cells",
  "3?",
  "T cells 1",
  "NES+ cells (necrotic?)",
  "Tumorcells_9 (BZWZ+)",
  "Neurons_1",
  "T cells 2",
  "TAMs_1",
  "Tumorcells_10 (BZWZ+)",
  "Tumorcells_2",
  "Tumorcells_3",
  "13??",
  "endothelialcells",
  "endothelialcells (L_Endo?)",
  "TAMs_2",
  "T cells 3",
  "Tumorcells_4 (MET+)",
  "Tumorcells_5 (LOX+ (%50), subcluster later)",
  "Astrocytes_1 (GFAP+)",
  "Tumorcells_6",
  "VLMC_1 (DCN+)",
  "23?",
  "IFI6+ IFI44L+ cells",
  "VLMC_2 (DCN+?)",
  "26?",
  "TAMs_3 (ITGAX+)",
  "Tumorcells_7 (MET+)",
  "endothelialcells + T cells (subcluster!)",
  "Tumorcells_8 (CDH1+)",
  "31?",
  "32?",
  "PSENEN+ (Neuron?) cells",
  "VAMP1+ (Neuron?) cells",
  "PSEN1+ (Neuron?) cells",
  "SYNPR+ (Neuron?) cells",
  "PTPRZ1+ (Neuron?) cells",
  "undefined",
  "undefined"
)

# mark cell types in marker table
top_markers$FinalAnnotation <- annotations[match(as.numeric(top_markers$cluster)+1, 1:length(annotations))]

# mark cell types in VoltRon object
brain_vr_Pat90_all$CellType <- annotations[match(as.numeric(brain_vr_Pat90_all$Clusters)+1, 1:length(annotations))]

# annotate Seurat object
brain_seu_Pat90_all$CellType <- annotations[match(as.numeric(as.character(brain_seu_Pat90_all$Xenium_snn_res.1))+1, 1:length(annotations))]

####
## Update T cell annotation ####
####

brain_vr_Pat90_all$CellType_withtcells <- brain_vr_Pat90_all$CellType
brain_vr_Pat90_all$CellType_withtcells[match(names(brain_seu_Pat90_all$CellType_withTcells), rownames(Metadata(brain_vr_Pat90_all)))] <- 
  brain_seu_Pat90_all$CellType_withTcells

####
## Update annotation (Major Cell Type) ####
####

temp <- brain_vr_Pat90_all$CellType_withtcells
temp[temp == "subcluster VLMC_1 (DCN+)"] <- "Immune Cells"
temp[temp %in% c("CD4 T cells","CD8 T cells", "NK Cells?", "NK Cells? (2)", "Dendritic cells? (CTSH, CCR7, CD83)", "subcluster CTSH+ Cells", "subcluster NES+ Cells")] <- "Immune Cells"
temp[grepl("VLMC", temp)] <- "VLMC/CAF"
temp[temp == "Tumorcells_5 (LOX+ (%50), subcluster later)"] <- "VLMC/CAF"
temp[grepl("endothelial", temp)] <- "Endothelial Cells"
temp[temp == "Tumorcells_8 (CDH1+)"] <- "Epithelial Cells"
temp[temp == "IFI6+ IFI44L+ cells"] <- "IFN Cells"
temp[temp == "23?"] <- "IFN Cells"
temp[grepl("TAM", temp)] <- "TAMs"
temp[temp == "Astrocytes_1 (GFAP+)"] <- "Astrocytes (GFAP+)"
temp[temp == "Tumorcells_7 (MET+)"] <- "Proliferating Tumor Cells 1"
temp[temp == "Tumorcells_2"] <- "Proliferating Tumor Cells 2"
temp[temp == "Tumorcells_3"] <- "Proliferating Tumor Cells 3"
temp[temp == "subcluster Tumorcells (TOP2A, MKI67)"] <- "Immune Cells"
temp[temp == "32?"] <- "Proliferating Tumor Cells 4"
temp[grepl("(Neuron\\?)", temp)] <- "Tumor Cells 1"
temp[grepl("MET+", temp)] <- "Tumor Cells 2"
temp[temp %in% c("SLC17A6+ cells","Neurons_1")] <- "Tumor Cells 3"
temp[temp == "1? (TAP1+)"] <- "Tumor Cells 4"
temp[temp == "Tumorcells_9 (BZWZ+)"] <- "Tumor Cells 5"
temp[temp == "Tumorcells_10 (BZWZ+)"] <- "Tumor Cells 6"
temp[temp == "Tumorcells_1"] <- "Tumor Cells 7"
temp[temp == "Tumorcells_6"] <- "Tumor Cells 8"
temp[temp == "NES+ cells (necrotic?)"] <- "NES+ Cells"
temp[temp %in% c("31?", "26?", "13??", "3?")] <- "Low Count Cells"
brain_vr_Pat90_all$MajorCellType <- temp
# saveRDS(brain_vr_Pat90_all, file = "../data/VoltRonData/Pat90_all_section_annotated.rds")

####
# Make new filtered UMAP ####
####

# transfer labels
tmp <- setNames(brain_vr_Pat90_all$MajorCellType, vrSpatialPoints(brain_vr_Pat90_all))
brain_seu_Pat90_all$MajorCellType <- tmp[colnames(brain_seu_Pat90_all)]

# subset for umap
brain_seu_Pat90_all_sub <- brain_seu_Pat90_all[,!brain_seu_Pat90_all$MajorCellType %in% c("undefined", "Low Count Cells", "NES+ Cells")]
brain_seu_Pat90_all_sub <- NormalizeData(brain_seu_Pat90_all_sub, scale.factor = 1000)
brain_seu_Pat90_all_sub <- ScaleData(brain_seu_Pat90_all_sub)
brain_seu_Pat90_all_sub <- RunPCA(brain_seu_Pat90_all_sub, features = Features(brain_seu_Pat90_all_sub), npcs = 30)
brain_seu_Pat90_all_sub <- RunUMAP(brain_seu_Pat90_all_sub, dims = 1:30)
DimPlot(brain_seu_Pat90_all_sub, group.by = "MajorCellType", label = TRUE, repel = TRUE)

# add new umap
vrEmbeddings(brain_vr_Pat90_all, type = "umap_filtered") <- Embeddings(brain_seu_Pat90_all_sub, reduction = "umap")
vrEmbeddingPlot(brain_vr_Pat90_all, embedding = "umap_filtered")

####
# annotate Hot tumors ####
####

####
## get immunescores ####
####

# get score per modules using Seurat
brain_vr_Pat90_all_assay1_seu <- VoltRon::as.Seurat(brain_vr_Pat90_all, cell.assay = "Xenium", type = "image")
brain_vr_Pat90_all_assay1_seu <- NormalizeData(brain_vr_Pat90_all_assay1_seu, scale.factor = 1000)
for(i in 1:length(score_markers)){
  cur_markers <- intersect(score_markers[[i]], SeuratObject::Features(brain_vr_Pat90_all_assay1_seu))
  brain_vr_Pat90_all_assay1_seu <- AddModuleScore(brain_vr_Pat90_all_assay1_seu, features = list(cur_markers), name = names(score_markers)[i], ctrl = 5)  
}

# transfer module scores to voltron objects
temp <- rep(NA, length(vrSpatialPoints(brain_vr_Pat90_all)))
names(temp) <- vrSpatialPoints(brain_vr_Pat90_all)
temp[names(brain_vr_Pat90_all_assay1_seu$ImmuneScores1)] <- brain_vr_Pat90_all_assay1_seu$ImmuneScores1
brain_vr_Pat90_all$ImmuneScores <- temp

# visualize score
vrSpatialFeaturePlot(brain_vr_Pat90_all, assay = "Xenium", features = c("ImmuneScores"),          
                     plot.segments = FALSE, pt.size = 1, n.tile = 300, alpha = 1)

####
## annotate hot tumor ####  
####

brain_vr_Pat90_all <- annotateSpatialData(brain_vr_Pat90_all, assay = "Assay1", use.image = TRUE, channel = "H&E")
brain_vr_Pat90_all <- annotateSpatialData(brain_vr_Pat90_all, assay = "Assay2", use.image = TRUE, channel = "H&E")
brain_vr_Pat90_all$annotation <- ifelse(grepl("hot", brain_vr_Pat90_all$annotation.1), brain_vr_Pat90_all$annotation.1, brain_vr_Pat90_all$annotation)

####
# save as AnnData ####
####

# interactive TissUUmaps tool
as.AnnData(brain_vr_Pat90_all, file = "../data/AnnDataData/Pat_90.h5ad", flip_coordinates = TRUE, 
           assay = "Xenium", method = "anndata", name = "main", channel = "H&E")
library(VoltRon)
library(Seurat)
library(harmony)
library(dplyr)
library(ComplexHeatmap)
library(xlsx)

####
# Import data ####
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

# get image data and register
brain_image_Pat24 <- importImageData(image = "../data/images/H&E/TMA.TIFF/Xenium-12.tif", tile.size = 10)
xen_reg <- registerSpatialData(object_list = list(brain_vr_Pat24, brain_image_Pat24))
vrImages(brain_vr_Pat24[["Assay3"]], name = "DAPI", channel = "H&E") <- vrImages(xen_reg$registered_spat[[2]], name = "image_1_reg")

####
# clustering ####
####

brain_vr_Pat24 <- subset(brain_vr, samples = "Pat24_19_1b")

# clustering with Seurat
brain_seu_Pat24 <- as.Seurat(brain_vr_Pat24, cell.assay = "Xenium", type = "image")
brain_seu_Pat24 <- NormalizeData(brain_seu_Pat24, scale.factor = 100)

# variable features
brain_seu_Pat24 <- FindVariableFeatures(brain_seu_Pat24, selection.method = "vst", nfeatures = 150)

# dimensionality reduction
brain_seu_Pat24 <- ScaleData(brain_seu_Pat24)

brain_seu_Pat24 <- RunPCA(brain_seu_Pat24, features = VariableFeatures(brain_seu_Pat24), npcs = 15)
brain_seu_Pat24 <- RunUMAP(brain_seu_Pat24, dims = 1:15)
DimPlot(brain_seu_Pat24, reduction = "umap", group.by = "Sample")

# clustering
brain_seu_Pat24 <- FindNeighbors(brain_seu_Pat24, dims = 1:15)
brain_seu_Pat24 <- FindClusters(brain_seu_Pat24, resolution = c(0.6,0.8,1.0,1.1,1.2))

# visualize
g_list <- list()
g_list[["Xenium_snn_res.0.6"]] <- DimPlot(brain_seu_Pat24, reduction = "umap", group.by = "Xenium_snn_res.0.6", label = T)
g_list[["Xenium_snn_res.0.8"]] <- DimPlot(brain_seu_Pat24, reduction = "umap", group.by = "Xenium_snn_res.0.8", label = T)
g_list[["Xenium_snn_res.1"]] <- DimPlot(brain_seu_Pat24, reduction = "umap", group.by = "Xenium_snn_res.1", label = T)
g_list[["Xenium_snn_res.1.1"]] <- DimPlot(brain_seu_Pat24, reduction = "umap", group.by = "Xenium_snn_res.1.1", label = T)
g_list[["Xenium_snn_res.1.2"]] <- DimPlot(brain_seu_Pat24, reduction = "umap", group.by = "Xenium_snn_res.1.2", label = T)
ggarrange(plotlist = g_list, ncol = 3, nrow = 2)

# res 0.8 looks good
brain_vr_Pat24$cluster <- as.character(brain_seu_Pat24$Xenium_snn_res.0.8)

####
# markers analysis ####
####

####
## check res 0.8 ####
####

# find top markers
Idents(brain_seu_Pat24) <- "Xenium_snn_res.0.8"
markers <- FindAllMarkers(brain_seu_Pat24)
top_markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.8, p_val_adj < 0.05, pct.1 > 0.5)

####
# annotation ####
####

# select a top_markers from marker analysis section
top_markers_table <- top_markers_table %>%
  group_by(cluster) %>%
  mutate(temp = gsub(" ", "", paste(Cell.type.label, collapse = ","))) %>%
  mutate(FinalAnnotation = names(table(strsplit(temp, split = ",")[[1]]))[which.max(table(strsplit(temp, split = ",")[[1]]))]) %>%
  select(-temp)  
top_markers_celltype <- top_markers_table %>%
  select(cluster, FinalAnnotation) %>%
  distinct() %>%
  group_by(FinalAnnotation) %>%
  mutate(FinalAnnotation_distinct = paste(FinalAnnotation, 1:n(), sep = "_"))
top_markers_table$FinalAnnotation <- top_markers_celltype$FinalAnnotation_distinct[match(top_markers_table$cluster, top_markers_celltype$cluster)]

# update umap
vrEmbeddings(brain_vr_Pat24, type = "umap", overwrite = TRUE) <- Embeddings(brain_seu_Pat24, reduction = "umap")

# mark cell types in VoltRon and Seurat object
brain_vr_Pat24$CellType <- top_markers_celltype$FinalAnnotation[match(brain_vr_Pat24$cluster, top_markers_celltype$cluster)]
brain_seu_Pat24$CellType <- top_markers_celltype$FinalAnnotation[match(brain_vr_Pat24$cluster, top_markers_celltype$cluster)]

####
## Update Annotation ####
####

# oldnames "Tumorcells_2" "Tumorcells_4" "Tumorcells_1" "Tumorcells_3"
oldnames <- unique(brain_vr_Pat24$CellType[grepl("Tumorcells", brain_vr_Pat24$CellType)])
newnames <- c("Tumorcells_2 (MKI67/TOP2A/CENPF/GAS2L3+)", "Tumorcells_4 (ABCB5/MET+)", "Tumorcells_1 (ABCB5/MET/CENPF/CDH1+)", "Tumorcells_3 (ABCB5/MET+)")

brain_vr_Pat24$CellType_TumorAnnotation <- brain_vr_Pat24$CellType
brain_seu_Pat24$CellType_TumorAnnotation <- brain_seu_Pat24$CellType
for(i in 1:length(oldnames)){
  brain_vr_Pat24$CellType_TumorAnnotation[brain_vr_Pat24$CellType_TumorAnnotation == oldnames[i]] <- newnames[i]
  brain_seu_Pat24$CellType_TumorAnnotation[brain_seu_Pat24$CellType_TumorAnnotation == oldnames[i]] <- newnames[i]
}

####
## Update annotation (Major Cell Type) ####
####

temp <- brain_vr_Pat24$CellType_TumorAnnotation
temp[grep("Endo_1", temp)] <- "Endothelial Cells"
temp[grepl("Immunecells", temp)] <- "Immune Cells"
temp[grepl("Oligo_1", temp)] <- "Oligodendrocytes"
temp[grep("VLMC", temp)] <- "VLMC/CAF"
temp[grep("TAMs", temp)] <- "TAMs"
temp[temp == "Tumorcells_2 (MKI67/TOP2A/CENPF/GAS2L3+)"] <- "Proliferating Tumor Cells 1"
temp[grepl("Tumorcells_1", temp)] <- "Tumor Cells 1"
temp[grepl("Tumorcells_3", temp)] <- "Tumor Cells 2"
temp[grepl("Tumorcells_4", temp)] <- "Tumor Cells 3"
temp[grepl("Neurons_1", temp)] <- "undefined"
temp[grepl("Neurons_2", temp)] <- "Neurons"
temp[grepl("Neurons_3", temp)] <- "Neurons"
temp[grepl("Astro", temp)] <- "Astrocytes"
brain_vr_Pat24$MajorCellType <- temp
# saveRDS(brain_vr_Pat24, file = "../data/VoltRonData/brain_vr_4600_annotated.rds")

####
# annotate Hot tumors ####
####

####
## get immunescores ####
####

# get score per modules using Seurat
brain_vr_Pat24_assay1_seu <- VoltRon::as.Seurat(brain_vr_Pat24, cell.assay = "Xenium", type = "image")
brain_vr_Pat24_assay1_seu <- NormalizeData(brain_vr_Pat24_assay1_seu, scale.factor = 1000)
for(i in 1:length(score_markers)){
  cur_markers <- intersect(score_markers[[i]], SeuratObject::Features(brain_vr_Pat24_assay1_seu))
  brain_vr_Pat24_assay1_seu <- AddModuleScore(brain_vr_Pat24_assay1_seu, features = list(cur_markers), name = names(score_markers)[i], ctrl = 5)  
}

# transfer module scores to voltron objects
temp <- rep(NA, length(vrSpatialPoints(brain_vr_Pat24)))
names(temp) <- vrSpatialPoints(brain_vr_Pat24)
temp[names(brain_vr_Pat24_assay1_seu$ImmuneScores1)] <- brain_vr_Pat24_assay1_seu$ImmuneScores1
brain_vr_Pat24$ImmuneScores <- temp

# visualize score
vrSpatialFeaturePlot(brain_vr_Pat24, assay = "Xenium", features = c("ImmuneScores"),          
                     plot.segments = FALSE, pt.size = 1, n.tile = 300, alpha = 1)

####
## annotate hot tumor ####  
####

brain_vr_Pat24 <- annotateSpatialData(brain_vr_Pat24, use.image = TRUE, channel = "H&E")
vrSpatialPlot(brain_vr_Pat24, group.by = "annotation", n.tile = 300)
# saveRDS(brain_vr_Pat24, file = "../data/VoltRonData/brain_vr_Pat24_annotated.rds")

####
# save as AnnData ####
####

# convert to anndata
as.AnnData(brain_vr_Pat24, file = "../data/AnnDataData/brain_anndata_Pat24.h5ad", flip_coordinates = TRUE, 
           name = "DAPI", channel = "H&E", assay = "Xenium")


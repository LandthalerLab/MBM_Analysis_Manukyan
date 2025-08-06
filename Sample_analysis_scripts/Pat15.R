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
brain_image_Pat15A_A <- importImageData(image.path = "../data/images/H&E/OME.TIFF/Xenium-3.tif", tile.size = 10)
brain_image_Pat15A_B <- importImageData(image.path = "../data/images/H&E/OME.TIFF/Xenium-16.tif", tile.size = 10)

# register and save images
xen_reg1 <- registerSpatialData(object_list = list(subset(brain_vr_Pat15, samples = "Pat15_14a_A"),
                                                   brain_image_Pat15A_A),
                                keypoints = readRDS(file = "../data/AuxiliaryData/Pat15_14a_A_keypoints.rds"))
# saveRDS(xen_reg1$keypoints, file = "../data/AuxiliaryData/Pat15_14a_A_keypoints.rds")
vrImages(brain_vr_Pat15[["Assay5"]], name = "DAPI", channel = "H&E") <- vrImages(xen_reg1$registered_spat[[2]], name = "image_1_reg")

xen_reg2 <- registerSpatialData(object_list = list(subset(brain_vr_Pat15, samples = "Pat15_14_a_B"),
                                                   brain_image_Pat15A_B),
                                keypoints = readRDS(file = "../data/AuxiliaryData/Pat15_14_a_B_keypoints.rds"))
# saveRDS(xen_reg2$keypoints, file = "../data/AuxiliaryData/Pat15_14_a_B_keypoints.rds")
vrImages(brain_vr_Pat15[["Assay4"]], name = "DAPI", channel = "H&E") <- vrImages(xen_reg2$registered_spat[[2]], name = "image_1_reg")

####
# clustering ####
####

# read registered VoltRon objects
brain_vr_Pat15 <- subset(brain_vr, assays = c("Assay4", "Assay5"))

# get clusters from original invidiual seurat objects
brain_vr_Pat15$OriginalClusters <- brain_vr_Pat15$Sample
brain_vr_Pat15$OriginalClusters[grepl("Assay4$",vrSpatialPoints(brain_vr_Pat15))] <- brain_seu_list$Pat15_14_a_B$Clusters
brain_vr_Pat15$OriginalClusters[grepl("Assay5$",vrSpatialPoints(brain_vr_Pat15))] <- brain_seu_list$Pat15_14a_A$Clusters

# clustering with Seurat
brain_seu_Pat15 <- as.Seurat(brain_vr_Pat15, cell.assay = "Xenium", type = "image")
brain_seu_Pat15 <- NormalizeData(brain_seu_Pat15, scale.factor = 100)

# dimensionality reduction
brain_seu_Pat15 <- ScaleData(brain_seu_Pat15)
brain_seu_Pat15 <- RunPCA(brain_seu_Pat15, features = Features(brain_seu_Pat15), npcs = 20)
brain_seu_Pat15 <- RunUMAP(brain_seu_Pat15, dims = 1:20)
DimPlot(brain_seu_Pat15, reduction = "umap", group.by = "Sample")

# clustering
brain_seu_Pat15 <- FindNeighbors(brain_seu_Pat15, dims = 1:20)
brain_seu_Pat15 <- FindClusters(brain_seu_Pat15, resolution = c(0.6,0.8,1.0,1.1,1.2))

# visualize
g_list <- list()
g_list[["Xenium_snn_res.0.6"]] <- DimPlot(brain_seu_Pat15, reduction = "umap", group.by = "Xenium_snn_res.0.6", label = T)
g_list[["Xenium_snn_res.0.8"]] <- DimPlot(brain_seu_Pat15, reduction = "umap", group.by = "Xenium_snn_res.0.8", label = T)
g_list[["Xenium_snn_res.1"]] <- DimPlot(brain_seu_Pat15, reduction = "umap", group.by = "Xenium_snn_res.1", label = T)
g_list[["Xenium_snn_res.1.1"]] <- DimPlot(brain_seu_Pat15, reduction = "umap", group.by = "Xenium_snn_res.1.1", label = T)
g_list[["Xenium_snn_res.1.2"]] <- DimPlot(brain_seu_Pat15, reduction = "umap", group.by = "Xenium_snn_res.1.2", label = T)
ggarrange(plotlist = g_list, ncol = 3, nrow = 2)

# res 0.8 looks good
brain_vr_Pat15$cluster <- as.character(brain_seu_Pat15$Xenium_snn_res.1)

####
# marker analysis ####
####

####
## check resolution 1.0 ####
####

# marker analysis
brain_seu_Pat15 <- as.Seurat(brain_vr_Pat15, type = "image", cell.assay = "Xenium")
Idents(brain_seu_Pat15) <- "Xenium_snn_res.1"
markers <- FindAllMarkers(brain_seu_Pat15)
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
top_markers_celltype <- top_markers %>%
  select(cluster, FinalAnnotation) %>%
  distinct() %>%
  group_by(FinalAnnotation) %>%
  mutate(FinalAnnotation = paste(FinalAnnotation, 1:n(), sep = "_"))

####
## Manual annotation, ignore Final ####
####

# cluster 0 half CD163 half APOE ?
# cluster 1 Tumorcells_1
# cluster 2 Melanocyic tumor (GPNMB+)
# cluster 3 Oligodendrocytes (EFHD1+)
# cluster 4 TAMs
# cluster 5 Tumorcells_2
# cluster 6 endothelial cells 3 (ACTA2++ and ?)
# cluster 7 Immunecells_1
# cluster 8 Immunecells_2
# cluster 9 TOP2A+ MKI67+ cells
# cluster 10 VLMC (DCN+)
# cluster 11 Tumorcells_3
# cluster 12 endothelialcells_1
# cluster 13 endothelialcells_2
# cluster 14 tumorcells_3 (NGFR/PTPRZ1+)
# cluster 15 Oligodendrocytes (CHI3L1+ HHATL+)
# cluster 16 Astrocytes (CHI3L1+ GFAP+)
# cluster 17 Tumorcells_4 (MET++)
# cluster 18 ISG15+ OAS3+ cells

# get final annotations
annotations <- c("half CD163 half APOE ?", "Tumorcells_1", "Melanocyic tumor (GPNMB+)", "Oligodendrocytes (EFHD1+)", "TAMs", "Tumorcells_2",
                 "endothelialcells_3 (ACTA2++ and ?)", "Immunecells_1", "Immunecells_2", "TOP2A+ MKI67+ cells", "VLMC (DCN+)", "Tumorcells_3",
                 "endothelialcells_1", "endothelialcells_2", "tumorcells_3 (NGFR/PTPRZ1+)", "Oligodendrocytes (CHI3L1+ HHATL+)", "Astrocytes (CHI3L1+ GFAP+)",
                 "Tumorcells_4 (MET++)", "ISG15+ OAS3+ cells")

# mark cell types in marker table
top_markers$FinalAnnotation <- annotations[match(as.numeric(as.character(top_markers$cluster))+1, 1:length(annotations))]

# mark cell types in VoltRon and Seurat object
brain_vr_Pat15$CellType <- annotations[match(as.numeric(as.character(brain_vr_Pat15$cluster))+1, 1:length(annotations))]
brain_seu_Pat15$CellType <- annotations[match(as.numeric(as.character(brain_vr_Pat15$cluster))+1, 1:length(annotations))]

####
## Update Annotation ####
####

oldnames <- annotations[grepl("Tumorcells|tumorcells", annotations)]
newnames <- c("Tumorcells_1 (BZW2/GAS2L3+)", "SPP1+ CADM1+", "Tumorcells_2 (LOX+) ", "Tumorcells_3 (GAS2L3+)" , "Tumorcells_4 (MET/BZW2/GAS2L3+)")

brain_vr_Pat15$CellType_TumorAnnotation <- brain_vr_Pat15$CellType
brain_seu_Pat15$CellType_TumorAnnotation <- brain_seu_Pat15$CellType
for(i in 1:length(oldnames)){
  brain_vr_Pat15$CellType_TumorAnnotation[brain_vr_Pat15$CellType_TumorAnnotation == oldnames[i]] <- newnames[i]
  brain_seu_Pat15$CellType_TumorAnnotation[brain_seu_Pat15$CellType_TumorAnnotation == oldnames[i]] <- newnames[i]
}

# update umap
vrEmbeddings(brain_vr_Pat15, type = "umap", overwrite = TRUE) <- Embeddings(brain_seu_Pat15, reduction = "umap")

####
## Update annotation (Major Cell Type) ####
####

temp <- brain_vr_Pat15$CellType_TumorAnnotation
temp[grep("endothelialcells", temp)] <- "Endothelial Cells"
temp[grepl("Immunecells", temp)] <- "Immune Cells"
temp[temp == "ISG15+ OAS3+ cells"] <- "IFN Cells"
temp[temp == "Astrocytes (CHI3L1+ GFAP+)"] <- "Astrocytes"
temp[temp == "TOP2A+ MKI67+ cells"] <- "Proliferating Tumor Cells 1"
temp[grep("VLMC", temp)] <- "VLMC/CAF"
temp[temp == "Melanocyic tumor (GPNMB+)"] <- "Melanocyic tumor (GPNMB+ NES+)"
temp[grepl("Tumorcells_1", temp)] <- "Tumor Cells 1"
temp[grepl("Tumorcells_2", temp)] <- "Tumor Cells 2"
temp[grepl("Tumorcells_3", temp)] <- "Tumor Cells 3"
temp[grepl("Tumorcells_4", temp)] <- "Tumor Cells 4"
temp[grepl("SPP1+", temp)] <- "Tumor Cells 5"
temp[temp == "half CD163 half APOE ?"] <- "Low Count Cells"
temp[temp == "Melanocyic tumor (GPNMB+ NES+)"] <- "Tumor Cells 6"
temp[temp == "Oligodendrocytes (CHI3L1+ HHATL+)"] <- "Tumor Cells 7"
brain_vr_Pat15$MajorCellType <- temp
# saveRDS(brain_vr_Pat15, file = "../data/VoltRonData/brain_vr_Pat15_annotated.rds")

####
# annotate Hot tumors ####
####

####
## get immunescores ####
####

# get score per modules using Seurat
brain_vr_Pat15_assay1_seu <- VoltRon::as.Seurat(brain_vr_Pat15, cell.assay = "Xenium", type = "image")
brain_vr_Pat15_assay1_seu <- NormalizeData(brain_vr_Pat15_assay1_seu, scale.factor = 1000)
for(i in 1:length(score_markers)){
  cur_markers <- intersect(score_markers[[i]], SeuratObject::Features(brain_vr_Pat15_assay1_seu))
  brain_vr_Pat15_assay1_seu <- AddModuleScore(brain_vr_Pat15_assay1_seu, features = list(cur_markers), name = names(score_markers)[i], ctrl = 5)  
}

# transfer module scores to voltron objects
temp <- rep(NA, length(vrSpatialPoints(brain_vr_Pat15)))
names(temp) <- vrSpatialPoints(brain_vr_Pat15)
temp[names(brain_vr_Pat15_assay1_seu$ImmuneScores1)] <- brain_vr_Pat15_assay1_seu$ImmuneScores1
brain_vr_Pat15$ImmuneScores <- temp

# visualize score
vrSpatialFeaturePlot(brain_vr_Pat15, assay = "Xenium", features = c("ImmuneScores"),          
                     plot.segments = FALSE, pt.size = 1, n.tile = 300, alpha = 1)

####
## annotate hot tumor ####  
####

brain_vr_Pat15 <- annotateSpatialData(brain_vr_Pat15, assay = "Assay4", use.image = TRUE, channel = "H&E")
brain_vr_Pat15 <- annotateSpatialData(brain_vr_Pat15, assay = "Assay5", use.image = TRUE, channel = "H&E")
brain_vr_Pat15$annotation <- apply(cbind(brain_vr_Pat15$annotation, brain_vr_Pat15$annotation.1), 1, function(x){
  ifelse("hot tumor" %in% x, "hot tumor", "undefined")
})
vrSpatialPlot(brain_vr_Pat15, group.by = "annotation", n.tile = 300)
# saveRDS(brain_vr_Pat15, file = "../data/VoltRonData/brain_vr_Pat15_annotated.rds")

####
# save as AnnData ####
####

# interactive TissUUmaps tool
as.AnnData(brain_vr_Pat15, file = "../data/AnnDataData/brain_anndata_Pat15.h5ad", flip_coordinates = TRUE, 
           assay = "Xenium", method = "anndata", name = "DAPI", channel = "H&E")
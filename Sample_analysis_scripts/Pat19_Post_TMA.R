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

# score markers (Interferon, Immune etc.)
score_markers <- read.xlsx("../data/Supplementary Information/score_markers.xlsx", sheetName = "Sheet1")
score_markers <- as.list(score_markers)
score_markers <- lapply(score_markers, function(x) x[!is.na(x)]) 

####
# overlay images ####
####

# get clustered data
brain_vr_Pat19_Post_TMA <- subset(brain_vr, assays = c("Assay10", "Assay11"))

# get image data and register
brain_image_Pat19_Post_TMAB  <- importImageData(image = "../data/images/H&E/TMA.TIFF/Xenium-6.tif", tile.size = 10)

# register and save images, use 2022 vs 1992
brain_vr_Pat19_Post_TMA_B <- subset(brain_vr_Pat19_Post_TMA, samples = "Pat19_Post_TMA_15_1_B")
xen_reg1 <- registerSpatialData(object_list = list(brain_vr_Pat19_Post_TMA_B, brain_image_Pat19_Post_TMAB))
vrImages(brain_vr_Pat19_Post_TMA[["Assay10"]], name = "DAPI", channel = "H&E") <- vrImages(xen_reg1$registered_spat[[2]], name = "main_reg")

# get image data and register
brain_image_Pat19_Post_TMAA  <- importImageData(image = "../data/images/H&E/TMA.TIFF/Xenium-10.tif", tile.size = 10)

# register and save images, use directly, dont change anything
brain_vr_Pat19_Post_TMA_A <- subset(brain_vr_Pat19_Post_TMA, samples = "Pat19_Post_TMA_15_1_A")
xen_reg1 <- registerSpatialData(object_list = list(brain_vr_Pat19_Post_TMA_A, brain_image_Pat19_Post_TMAA))
vrImages(brain_vr_Pat19_Post_TMA[["Assay11"]], name = "DAPI", channel = "H&E") <- vrImages(xen_reg1$registered_spat[[2]], name = "main_reg")

####
# clustering ####
####

# clustering with Seurat
brain_seu_Pat19_Post_TMA <- VoltRon::as.Seurat(brain_vr_Pat19_Post_TMA, cell.assay = "Xenium", type = "image")
brain_seu_Pat19_Post_TMA <- NormalizeData(brain_seu_Pat19_Post_TMA, scale.factor = 100)

# dimensionality reduction
brain_seu_Pat19_Post_TMA <- ScaleData(brain_seu_Pat19_Post_TMA)
brain_seu_Pat19_Post_TMA <- RunPCA(brain_seu_Pat19_Post_TMA, features = Features(brain_seu_Pat19_Post_TMA), npcs = 20)
brain_seu_Pat19_Post_TMA <- RunUMAP(brain_seu_Pat19_Post_TMA, dims = 1:20)
DimPlot(brain_seu_Pat19_Post_TMA, reduction = "umap", group.by = "Sample", split.by = "Sample")

# clustering
brain_seu_Pat19_Post_TMA <- FindNeighbors(brain_seu_Pat19_Post_TMA, dims = 1:20)
brain_seu_Pat19_Post_TMA <- FindClusters(brain_seu_Pat19_Post_TMA, resolution = c(0.6,0.8,1.0,1.1,1.2))

# visualize
g_list <- list()
g_list[["Xenium_snn_res.0.6"]] <- DimPlot(brain_seu_Pat19_Post_TMA, reduction = "umap", group.by = "Xenium_snn_res.0.6", label = T)
g_list[["Xenium_snn_res.0.8"]] <- DimPlot(brain_seu_Pat19_Post_TMA, reduction = "umap", group.by = "Xenium_snn_res.0.8", label = T)
g_list[["Xenium_snn_res.1"]] <- DimPlot(brain_seu_Pat19_Post_TMA, reduction = "umap", group.by = "Xenium_snn_res.1", label = T)
g_list[["Xenium_snn_res.1.1"]] <- DimPlot(brain_seu_Pat19_Post_TMA, reduction = "umap", group.by = "Xenium_snn_res.1.1", label = T)
g_list[["Xenium_snn_res.1.2"]] <- DimPlot(brain_seu_Pat19_Post_TMA, reduction = "umap", group.by = "Xenium_snn_res.1.2", label = T)
ggarrange(plotlist = g_list, ncol = 3, nrow = 2)

# res 1 looks good
brain_vr_Pat19_Post_TMA$cluster <- as.character(brain_seu_Pat19_Post_TMA$Xenium_snn_res.1)

####
# marker analysis ####
####

####
# subclustering 12 ####
####

Idents(brain_seu_Pat19_Post_TMA) <- "Xenium_snn_res.1"
brain_seu_Pat19_Post_TMA <- FindSubCluster(brain_seu_Pat19_Post_TMA, cluster = "12", graph.name = "Xenium_snn", 
                                      subcluster.name = "subcluster", resolution = 0.4)

####
## marker analysis ####
####

Idents(brain_seu_Pat19_Post_TMA_sub) <- "subcluster"
markers <- FindAllMarkers(brain_seu_Pat19_Post_TMA_sub)
markers <- markers %>% left_join(xenium_markers, by = c("gene" = "Genes"))
top_markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.8, p_val_adj < 0.05, pct.1 > 0.5)

####
## update clusters ####
####

# update cluster indices based on xenium.res.1
tmp <- brain_seu_Pat19_Post_TMA$subcluster
tmp[tmp=="12_1"] <- "17"
tmp[tmp=="12_3"] <- "17"
tmp[tmp=="12_0"] <- "12"
tmp[tmp=="12_2"] <- "12"
brain_vr_Pat19_Post_TMA$cluster <- as.character(tmp)

####
# annotation ####
####

####
## Manual annotation ####
# NOTE: This is based on some manual subclustering, see above
####

# cluster 0 Tumorcells_1 (PDGFD+)
# cluster 1 Tumorcells_2 (ABCB5/MET/POSTN+)
# cluster 2 IFI6+ SERPINA3+
# cluster 3 T cells 1
# cluster 4 TAMs_1
# cluster 5 Tumorcells_3 (IGFBP3/IGFBP5/LOX+)
# cluster 6 TAMs_2
# cluster 7 SOX4+ tumorcells?
# cluster 8 Tumorcells_4 (TOP2A/CENPF+)
# cluster 9 T cells 2
# cluster 10 VLMC_1 (DCN+)
# cluster 11 Tumorcells_5 (LOX+)
# cluster 12 VLMC_2 (DCN+)
# cluster 13 ISG15+ IFI6+ IFI44L+ cells
# cluster 14 endothelialcells_1
# cluster 15 undefined
# cluster 16 Tumorcells_6 (LIF+)
# cluster 17 undefined 

# get final annotations
annotations <- c("Tumorcells_1 (PDGFD+)", "Tumorcells_2 (ABCB5/MET/POSTN+)", "IFI6+ SERPINA3+", "T cells 1", "TAMs_1", "Tumorcells_3 (IGFBP3/IGFBP5/LOX+)", 
                 "TAMs_2", "SOX4+ tumorcells?", "Tumorcells_4 (TOP2A/CENPF+)", "T cells 2", "VLMC_1 (DCN+)", "Tumorcells_5 (LOX+)", "VLMC_2 (DCN+)", 
                 "ISG15+ IFI6+ IFI44L+ cells", "endothelialcells_1", "undefined", "Tumorcells_6 (LIF+)", "undefined")

# mark cell types in marker table
top_markers$FinalAnnotation <- annotations[match(as.numeric(as.character(top_markers$cluster))+1, 1:length(annotations))]

# mark cell types in VoltRon and Seurat object
brain_vr_Pat19_Post_TMA$CellType <- annotations[match(as.numeric(as.character(brain_vr_Pat19_Post_TMA$cluster))+1, 1:length(annotations))]
brain_seu_Pat19_Post_TMA$CellType <- annotations[match(as.numeric(as.character(brain_vr_Pat19_Post_TMA$cluster))+1, 1:length(annotations))]

# update umap
vrEmbeddings(brain_vr_Pat19_Post_TMA, type = "umap", overwrite = TRUE) <- Embeddings(brain_seu_Pat19_Post_TMA, reduction = "umap")

####
## Update annotation (Major Cell Type) ####
####

temp <- brain_vr_Pat19_Post_TMA$CellType
temp[grep("endothelialcells", temp)] <- "Endothelial Cells"
temp[grepl("T cell", temp)] <- "Immune Cells"
temp[temp == "ISG15+ IFI6+ IFI44L+ cells"] <- "IFN Cells"
temp[grepl("TOP2A", temp)] <- "Proliferating Tumor Cells 1"
temp[grep("VLMC", temp)] <- "VLMC/CAF"
temp[grep("TAMs", temp)] <- "TAMs"
temp[grepl("Tumorcells_1", temp)] <- "Tumor Cells 1"
temp[grepl("Tumorcells_2", temp)] <- "Tumor Cells 2"
temp[grepl("Tumorcells_3", temp)] <- "Tumor Cells 3"
temp[temp == "IFI6+ SERPINA3+"] <- "Tumor Cells 4"
temp[grepl("Tumorcells_5", temp)] <- "Tumor Cells 5"
temp[grepl("Tumorcells_6", temp)] <- "Tumor Cells 6"
temp[temp == "SOX4+ tumorcells?"] <- "Tumor Cells 7"
brain_vr_Pat19_Post_TMA$MajorCellType <- temp
# saveRDS(brain_vr_Pat19_Post_TMA, file = "../data/VoltRonData/brain_vr_Pat19_Post_TMA_annotated.rds")

####
# save as AnnData ####
####

# interactive TissUUmaps tool
as.AnnData(brain_vr_Pat19_Post_TMA, file = "../data/AnnDataData/brain_anndata_Pat19_Post_TMA.h5ad", flip_coordinates = TRUE, 
           assay = "Xenium", method = "anndata", name = "DAPI", channel = "H&E")

####
# annotate Hot tumors ####
####

####
## get immunescores ####
####

# get score per modules using Seurat
brain_vr_Pat19_Post_TMA_assay1_seu <- VoltRon::as.Seurat(brain_vr_Pat19_Post_TMA, cell.assay = "Xenium", type = "image")
brain_vr_Pat19_Post_TMA_assay1_seu <- NormalizeData(brain_vr_Pat19_Post_TMA_assay1_seu, scale.factor = 1000)
for(i in 1:length(score_markers)){
  cur_markers <- intersect(score_markers[[i]], SeuratObject::Features(brain_vr_Pat19_Post_TMA_assay1_seu))
  brain_vr_Pat19_Post_TMA_assay1_seu <- AddModuleScore(brain_vr_Pat19_Post_TMA_assay1_seu, features = list(cur_markers), name = names(score_markers)[i], ctrl = 5)  
}

# transfer module scores to voltron objects
temp <- rep(NA, length(vrSpatialPoints(brain_vr_Pat19_Post_TMA)))
names(temp) <- vrSpatialPoints(brain_vr_Pat19_Post_TMA)
temp[names(brain_vr_Pat19_Post_TMA_assay1_seu$ImmuneScores1)] <- brain_vr_Pat19_Post_TMA_assay1_seu$ImmuneScores1
brain_vr_Pat19_Post_TMA$ImmuneScores <- temp

# visualize score
vrSpatialFeaturePlot(brain_vr_Pat19_Post_TMA, assay = "Xenium", features = c("ImmuneScores"),          
                     plot.segments = FALSE, pt.size = 1, n.tile = 300, alpha = 1)

####
## annotate hot tumor ####  
####

brain_vr_Pat19_Post_TMA <- annotateSpatialData(brain_vr_Pat19_Post_TMA, assay = "Assay10", use.image = TRUE, channel = "H&E")
brain_vr_Pat19_Post_TMA <- annotateSpatialData(brain_vr_Pat19_Post_TMA, assay = "Assay11", use.image = TRUE, channel = "H&E")
brain_vr_Pat19_Post_TMA$annotation <- apply(cbind(brain_vr_Pat19_Post_TMA$annotation, brain_vr_Pat19_Post_TMA$annotation.1), 1, function(x){
  ifelse("hot tumor" %in% x, "hot tumor", "undefined")
})
vrSpatialPlot(brain_vr_Pat19_Post_TMA, group.by = "annotation", n.tile = 300)
# saveRDS(brain_vr_Pat19_Post_TMA, file = "../data/VoltRonData/brain_vr_Pat19_Post_TMA_annotated.rds")
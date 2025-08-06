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
Pat19_Post_section <- importXenium("../../../data/RadkeRedmerSpatial_Brain_05112023/20240124__160424__Landthaler_human_brain_24012024/output-XETG00046__0020704__Region_1__20240124__161116//",
                              resolution_level = 3, overwrite_resolution = TRUE, sample_name = "Pat19_Post")
Pat19_Pre_section <- importXenium("../../../data/RadkeRedmerSpatial_Brain_05112023/20240124__160424__Landthaler_human_brain_24012024/output-XETG00046__0020704__Region_2__20240124__161116///",
                              resolution_level = 3, overwrite_resolution = TRUE, sample_name = "Pat19_Pre")
brain_all <- merge(Pat19_Post_section, Pat19_Pre_section)

####
# overlay images ####
####

# width 7374 vs 7392
imgdata <- importImageData("../data/Images/H&E/sections/MBM_HE_1.2.24/Xenium ID0020704_sec2.tif")
xen_reg <- registerSpatialData(object_list = list(Pat19_Post_section, imgdata))
vrImages(Pat19_Post_section[["Assay1"]], name = "main", channel = "H&E") <- vrImages(xen_reg$registered_spat[[2]], name = "main_reg")

# width 4778 vs 4730
imgdata2 <- importImageData("../data/Images/H&E/sections/MBM_HE_1.2.24/Xenium ID0020704_sec1.tif")
xen_reg <- registerSpatialData(object_list = list(Pat19_Pre_section, imgdata2))
vrImages(Pat19_Pre_section[["Assay1"]], name = "main", channel = "H&E") <- vrImages(xen_reg$registered_spat[[2]], name = "main_reg")

####
# clustering ####
####

# clustering with Seurat
brain_seu_Pat19_Post_Pat19_Pre <- as.Seurat(brain_all, cell.assay = "Xenium", type = "image")
brain_seu_Pat19_Post_Pat19_Pre <- subset(brain_seu_Pat19_Post_Pat19_Pre, subset = nCount_Xenium > 5)
brain_seu_Pat19_Post_Pat19_Pre <- NormalizeData(brain_seu_Pat19_Post_Pat19_Pre, scale.factor = 100)

# dimensionality reduction
brain_seu_Pat19_Post_Pat19_Pre <- ScaleData(brain_seu_Pat19_Post_Pat19_Pre)
brain_seu_Pat19_Post_Pat19_Pre <- RunPCA(brain_seu_Pat19_Post_Pat19_Pre, features = Features(brain_seu_Pat19_Post_Pat19_Pre), npcs = 30)
brain_seu_Pat19_Post_Pat19_Pre <- RunUMAP(brain_seu_Pat19_Post_Pat19_Pre, dims = 1:20)

# clustering
brain_seu_Pat19_Post_Pat19_Pre <- FindNeighbors(brain_seu_Pat19_Post_Pat19_Pre, dims = 1:20)
brain_seu_Pat19_Post_Pat19_Pre <- FindClusters(brain_seu_Pat19_Post_Pat19_Pre, resolution = c(0.2, 0.3, 0.4))
brain_seu_Pat19_Post_Pat19_Pre <- FindClusters(brain_seu_Pat19_Post_Pat19_Pre, resolution = c(0.5, 0.7))
brain_seu_Pat19_Post_Pat19_Pre <- FindClusters(brain_seu_Pat19_Post_Pat19_Pre, resolution = c(0.9, 1.0))
brain_seu_Pat19_Post_Pat19_Pre <- FindClusters(brain_seu_Pat19_Post_Pat19_Pre, resolution = c(1.1, 1.2))

# visualize
g_list <- list()
g_list[["Xenium_snn_res.0.2"]] <- DimPlot(brain_seu_Pat19_Post_Pat19_Pre, reduction = "umap", group.by = "Xenium_snn_res.0.2", label = T) + NoLegend()
g_list[["Xenium_snn_res.0.3"]] <- DimPlot(brain_seu_Pat19_Post_Pat19_Pre, reduction = "umap", group.by = "Xenium_snn_res.0.3", label = T) + NoLegend()
g_list[["Xenium_snn_res.0.4"]] <- DimPlot(brain_seu_Pat19_Post_Pat19_Pre, reduction = "umap", group.by = "Xenium_snn_res.0.4", label = T) + NoLegend()
g_list[["Xenium_snn_res.0.5"]] <- DimPlot(brain_seu_Pat19_Post_Pat19_Pre, reduction = "umap", group.by = "Xenium_snn_res.0.5", label = T) + NoLegend()
g_list[["Xenium_snn_res.0.7"]] <- DimPlot(brain_seu_Pat19_Post_Pat19_Pre, reduction = "umap", group.by = "Xenium_snn_res.0.7", label = T) + NoLegend()
g_list[["Xenium_snn_res.0.9"]] <- DimPlot(brain_seu_Pat19_Post_Pat19_Pre, reduction = "umap", group.by = "Xenium_snn_res.0.9", label = T) + NoLegend()
g_list[["Xenium_snn_res.1"]] <- DimPlot(brain_seu_Pat19_Post_Pat19_Pre, reduction = "umap", group.by = "Xenium_snn_res.1", label = T) + NoLegend()
g_list[["Xenium_snn_res.1.1"]] <- DimPlot(brain_seu_Pat19_Post_Pat19_Pre, reduction = "umap", group.by = "Xenium_snn_res.1.1", label = T) + NoLegend()
g_list[["Xenium_snn_res.1.2"]] <- DimPlot(brain_seu_Pat19_Post_Pat19_Pre, reduction = "umap", group.by = "Xenium_snn_res.1.2", label = T) + NoLegend()
ggarrange(plotlist = g_list, ncol = 5, nrow = 2)

brain_seu_Pat19_Post_Pat19_Pre_filter <- brain_seu_Pat19_Post_Pat19_Pre[,!brain_seu_Pat19_Post_Pat19_Pre$Xenium_snn_res.0.9 %in% c("27","28","29","30")]

# get clusters to VoltRon
Clusters <- brain_seu_Pat19_Post_Pat19_Pre_filter$Xenium_snn_res.0.9[brain_seu_Pat19_Post_Pat19_Pre_filter$Sample == "Pat19_Post"]
Clusters <- as.character(Clusters[vrSpatialPoints(Pat19_Post_section)])
Clusters[is.na(Clusters)] <- "undefined"
Pat19_Post_section$Clusters <- Clusters
Clusters <- brain_seu_Pat19_Post_Pat19_Pre_filter$Xenium_snn_res.0.9[brain_seu_Pat19_Post_Pat19_Pre_filter$Sample == "Pat19_Pre"]
Clusters <- as.character(Clusters[vrSpatialPoints(brain_all, assay = "Assay2")])
Clusters[is.na(Clusters)] <- "undefined"
Pat19_Pre_section$Clusters <- Clusters

# get embeddings to VoltRon
embed <- Embeddings(brain_seu_Pat19_Post_Pat19_Pre_filter, reduction = "umap")
overlapping_cells <- rownames(embed)[rownames(embed) %in% vrSpatialPoints(Pat19_Post_section)]
overlapping_cells <- vrSpatialPoints(Pat19_Post_section)[vrSpatialPoints(Pat19_Post_section) %in% overlapping_cells]
vrEmbeddings(Pat19_Post_section, type = "umap") <- Embeddings(brain_seu_Pat19_Post_Pat19_Pre_filter, reduction = "umap")[overlapping_cells,]
embed <- Embeddings(brain_seu_Pat19_Post_Pat19_Pre_filter, reduction = "umap")
overlapping_cells <- rownames(embed)[rownames(embed) %in% gsub("Assay1", "Assay2", vrSpatialPoints(Pat19_Pre_section))]
overlapping_cells <- gsub("Assay1", "Assay2", vrSpatialPoints(Pat19_Pre_section))[gsub("Assay1", "Assay2", vrSpatialPoints(Pat19_Pre_section)) %in% overlapping_cells]
embed <- Embeddings(brain_seu_Pat19_Post_Pat19_Pre_filter, reduction = "umap")[overlapping_cells,]
rownames(embed) <- gsub("Assay2", "Assay1", rownames(embed))
vrEmbeddings(Pat19_Pre_section, type = "umap", overwrite = TRUE) <- embed

# visualize embeddings
brain_all <- merge(Pat19_Post_section, Pat19_Pre_section)
vrEmbeddingPlot(brain_all, embedding = "umap", group.by = "CellType_TumorAnnotation")

####
# Marker Analysis ####
####

# Xenium_snn_res.0.9 looks good
brain_seu_Pat19_Post_Pat19_Pre_filter <- brain_seu_Pat19_Post_Pat19_Pre[,!brain_seu_Pat19_Post_Pat19_Pre$Xenium_snn_res.0.9 %in% c("27","28","29","30")]
Idents(brain_seu_Pat19_Post_Pat19_Pre_filter) <- "Xenium_snn_res.0.9"
markers <- FindAllMarkers(brain_seu_Pat19_Post_Pat19_Pre_filter)
top_markers <- markers %>% filter(avg_log2FC > 0.5, p_val_adj < 0.05) %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 20)
  unique_markers <- top_markers %>%
  dplyr::filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
  group_by(gene) %>% mutate(n_gene = n()) %>%
  filter(n_gene < 2) %>% select(-n_gene)

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

# cluster 0 Tumorcells_1 (MET+/SOX6+)
# cluster 1 PDGFD+ NES+
# cluster 2 PTPRZ1+ (partly OLIG1+)
# cluster 3 Tumorcells_2 (ABCB5+)
# cluster 4 Tumorcells_3 (SOX4+)
# cluster 5 TAMs_2 (C1QC+)
# cluster 6 Immunecells_1
# cluster 7 7_?
# cluster 8 VLMC (DCN+??)
# cluster 9 Tumorcells_4 (TOP2A+)
# cluster 10 Tumorcells_5 (TOP2A+)
# cluster 11 Tumorcells_6 (IGFBP3+/IGFBP5+)
# cluster 12 Immunecells_2
# cluster 13 mix_13 (cluster again, partly endo)
# cluster 14 TAMs_3 (RNASET2+)
# cluster 15 endothelialcells
# cluster 16 TAMs_1 (C1QC+)
# cluster 17 mix_17 (EFHD1+/TRPC6+)
# cluster 18 18_?
# cluster 19 ISG15+ OAS3+ MX1+ MX2+ cells (subcluster?)
# cluster 20 Tumorcells_7 (IDH1+)
# cluster 21 Astrocyte (CADM1+)
# cluster 22 22_?
# cluster 23 23_mix (HILDPA+/CCND3+)
# cluster 24 Tumorcells_8 (TOP2A+)
# cluster 25 25_? (NCAM1+)
# cluster 26 VLMC (DCN+)

annotations <- c("Tumorcells_1 (MET+/SOX6+)", 
                 "PDGFD+ NES+", 
                 "PTPRZ1+ (partly OLIG1+)", 
                 "Tumorcells_2 (ABCB5+)", 
                 "Tumorcells_3 (SOX4+)", 
                 "TAMs_2 (C1QC+)",
                 "Immunecells_1", 
                 "7_?", 
                 "VLMC (DCN+??)", 
                 "Tumorcells_4 (TOP2A+)", 
                 "Tumorcells_5 (TOP2A+)",
                 "Tumorcells_6 (IGFBP3+/IGFBP5+)",
                 "Immunecells_2", 
                 "mix_13 (cluster again, partly endo)", 
                 "TAMs_3 (RNASET2+)",
                 "endothelialcells", 
                 "TAMs_1 (C1QC+)", 
                 "mix_17 (EFHD1+/TRPC6+)", 
                 "18_?", 
                 "ISG15+ OAS3+ MX1+ MX2+ cells (subcluster?)", 
                 "Tumorcells_7 (IDH1+)", 
                 "Astrocyte (CADM1+)", 
                 "22_?", 
                 "23_mix (HILDPA+/CCND3+)", 
                 "Tumorcells_8 (TOP2A+)", 
                 "25_? (NCAM1+)", 
                 "VLMC (DCN+)")

# mark cell types in marker table
top_markers$FinalAnnotation <- annotations[match(as.numeric(top_markers$cluster)+1, 1:length(annotations))]

# mark cell types in VoltRon object
Pat19_Post_section$CellType <- annotations[match(as.numeric(Pat19_Post_section$Clusters)+1, 1:length(annotations))]
Pat19_Pre_section$CellType <- annotations[match(as.numeric(Pat19_Pre_section$Clusters)+1, 1:length(annotations))]

# mark cell types in Seurat object
brain_seu_Pat19_Post_Pat19_Pre$CellType <- annotations[match(as.numeric(as.character(brain_seu_Pat19_Post_Pat19_Pre$Xenium_snn_res.0.9))+1, 1:length(annotations))]

####
## Update Annotation ####
####

# update tumor annotation
top_markers <- read.xlsx("../data/MarkerData/brain_seu_Pat19_Post_Pat19_Pre_marker_analysis.xlsx", sheetName = "results")
Pat19_Post_section <- readRDS(file = "../data/VoltRonData/Pat19_Post_section_annotated.rds")
Pat19_Pre_section <- readRDS(file = "../data/VoltRonData/Pat19_Pre_section_annotated.rds")
brain_seu_Pat19_Post_Pat19_Pre <- readRDS(file = "../data/SeuratData/brain_seu_Pat19_Post_Pat19_Pre_annotated.rds")
  
oldnames <- annotations[grepl("Tumorcells", annotations)]
newnames <- c("Tumorcells_1 (MET+)", "Tumorcells_2 (ABCB5/POSTN+)", "SOX4+ cells (tumor?)", 
              "Tumorcells_4 (MKI67/TOP2A/CENPF+)", "Tumorcells_5 (TOP2A/CENPF+)", "Tumorcells_6 (LOX+)", 
              "Tumorcells_7 (DDR2+)", "Tumorcells_8 (TOP2A/CENPF+)")
Pat19_Post_section$CellType_TumorAnnotation <- Pat19_Post_section$CellType_withimmuncells
Pat19_Pre_section$CellType_TumorAnnotation <- Pat19_Pre_section$CellType_withimmuncells
brain_seu_Pat19_Post_Pat19_Pre$CellType_TumorAnnotation <- brain_seu_Pat19_Post_Pat19_Pre$CellType_withimmuncells
top_markers$UpdatedAnnotation <- top_markers$FinalAnnotation
for(i in 1:length(oldnames)){
  Pat19_Post_section$CellType_TumorAnnotation[Pat19_Post_section$CellType_TumorAnnotation == oldnames[i]] <- newnames[i]
  Pat19_Pre_section$CellType_TumorAnnotation[Pat19_Pre_section$CellType_TumorAnnotation == oldnames[i]] <- newnames[i]
  brain_seu_Pat19_Post_Pat19_Pre$CellType_TumorAnnotation[brain_seu_Pat19_Post_Pat19_Pre$CellType_TumorAnnotation == oldnames[i]] <- newnames[i]
  top_markers$UpdatedAnnotation[top_markers$UpdatedAnnotation == oldnames[i]] <- newnames[i]
}

####
## Update annotation (Major Cell Type) ####
####

Pat19_Post_section <- readRDS(file = "../data/VoltRonData/Pat19_Post_section_annotated.rds")
Pat19_Pre_section <- readRDS(file = "../data/VoltRonData/Pat19_Pre_section_annotated.rds")
brain_seu_Pat19_Post_Pat19_Pre <- readRDS(file = "../data/SeuratData/brain_seu_Pat19_Post_Pat19_Pre_annotated.rds")

temp <- brain_seu_Pat19_Post_Pat19_Pre$CellType_TumorAnnotation
temp_subcluster <- brain_seu_Pat19_Post_Pat19_Pre$`subcluster_7_?`
temp[grepl("VLMC", temp)] <- "VLMC/CAF"
temp[grepl("endothelial", temp)] <- "Endothelial Cells"
temp[grepl("TAM", temp)] <- "TAMs"
# tumor cells
temp[temp== "SOX4+ cells (tumor?)"] <- "Tumor Cells 1"
temp[temp== "PTPRZ1+ (partly OLIG1+)"] <- "Tumor Cells 2"
temp[temp== "Tumorcells_7 (DDR2+)"] <- "Tumor Cells 3"
temp[temp== "Tumorcells_1 (MET+)"] <- "Tumor Cells 4"
temp[temp== "Tumorcells_6 (LOX+)"] <- "Tumor Cells 5"
temp[temp == "22_?"] <- "Tumor Cells 6"
temp[temp == "Tumorcells_5 (TOP2A/CENPF+)"] <- "Proliferating Tumor Cells 1"
temp[temp == "Tumorcells_4 (MKI67/TOP2A/CENPF+)"] <- "Proliferating Tumor Cells 2"
temp[temp == "Tumorcells_8 (TOP2A/CENPF+)"] <- "Proliferating Tumor Cells 3"
temp[temp == "PDGFD+ NES+"] <- "NES+ Cells"
temp[temp == "Tumorcells_2 (ABCB5/POSTN+)"] <- "NES+ Cells"
# immune cells
temp[grepl("T cells", temp)] <- "Immune Cells"
temp[temp == "Immunecells_1_11"] <- "Immune Cells"
temp[temp == "Immunecells_2_13?"] <- "Immune Cells"
temp[grep("IL7R", temp)] <- "Immune Cells"
temp[grepl("T cells", temp)] <- "Immune Cells"
temp[grepl("CTLA4\\+ CCR8\\+ cells", temp)] <- "Immune Cells"
temp[grepl("CYTIP\\+ cells", temp)] <- "Immune Cells"
temp[grepl("CD3G\\+ TRAC\\+ cells", temp)] <- "Immune Cells"
temp[grepl("CD8A\\+ cells", temp)] <- "Immune Cells"
temp[grepl("CCL4\\+ cells", temp)] <- "Immune Cells"
temp[grepl("CYTIP\\+ cells", temp)] <- "Immune Cells"
temp[grepl("CCL5\\+ cells", temp)] <- "Immune Cells"
temp[grepl("CD48\\+ cells", temp)] <- "Immune Cells"
temp[grepl("CD48\\+ cells", temp)] <- "Immune Cells"
temp[grepl("CXCR4\\+ cells", temp)] <- "Immune Cells"
temp[temp == "CD8A cytotoxic (CCL5+ NKG7+ CCL4+)"] <- "Immune Cells"
# ifn cells
temp[temp == "ISG15+ OAS3+ MX1+ MX2+ cells (subcluster?)"] <- "IFN Cells"
# low count
temp[temp == "PTPRC++ cells"] <- "Low Count Cells"
temp[temp == "CD2+ cells"] <- "Low Count Cells"
temp[temp == "mix_17 (EFHD1+/TRPC6+)"] <- "Low Count Cells"
temp[temp == "23_mix (HILDPA+/CCND3+)"] <- "Low Count Cells"
temp[temp == "18_?"] <- "Low Count Cells"
# undefined
temp[temp == "mix_13 (cluster again, partly endo)"] <- "undefined"
temp[temp == "7_?"] <- "undefined"
temp[temp_subcluster == "7_?_2"] <- "Astrocytes"
temp[temp == "CD69+ TCF7+ cells (mix?)"] <- "undefined"
temp[temp == "25_? (NCAM1+)"] <- "undefined"
brain_seu_Pat19_Post_Pat19_Pre$MajorCellType <- temp
saveRDS(brain_seu_Pat19_Post_Pat19_Pre, file = "../data/SeuratData/brain_seu_Pat19_Post_Pat19_Pre_annotated.rds")

Pat19_Post_section$MajorCellType <- rep(NA, length(vrSpatialPoints(Pat19_Post_section)))
Pat19_Pre_section$MajorCellType <- rep(NA, length(vrSpatialPoints(Pat19_Pre_section)))
temp <- setNames(Pat19_Post_section$MajorCellType, vrSpatialPoints(Pat19_Post_section))
Pat19_Post_section$MajorCellType  <- brain_seu_Pat19_Post_Pat19_Pre$MajorCellType[vrSpatialPoints(Pat19_Post_section)]
Pat19_Pre_section$MajorCellType  <- brain_seu_Pat19_Post_Pat19_Pre$MajorCellType[gsub("Assay1", "Assay2", vrSpatialPoints(Pat19_Pre_section))]
# saveRDS(Pat19_Post_section, file = "../data/VoltRonData/Pat19_Post_section_annotated.rds")
# saveRDS(Pat19_Pre_section, file = "../data/VoltRonData/Pat19_Pre_section_annotated.rds")

####
# Make new filtered UMAP ####
####

# subset for umap
brain_seu_Pat19_Post_Pat19_Pre_sub <- brain_seu_Pat19_Post_Pat19_Pre[,!brain_seu_Pat19_Post_Pat19_Pre$MajorCellType %in% c("undefined", "Low Count Cells", "NES+ Cells")]
brain_seu_Pat19_Post_Pat19_Pre_sub <- NormalizeData(brain_seu_Pat19_Post_Pat19_Pre_sub, scale.factor = 1000)
brain_seu_Pat19_Post_Pat19_Pre_sub <- ScaleData(brain_seu_Pat19_Post_Pat19_Pre_sub)
brain_seu_Pat19_Post_Pat19_Pre_sub <- RunPCA(brain_seu_Pat19_Post_Pat19_Pre_sub, features = Features(brain_seu_Pat19_Post_Pat19_Pre_sub), npcs = 30)
brain_seu_Pat19_Post_Pat19_Pre_sub <- RunUMAP(brain_seu_Pat19_Post_Pat19_Pre_sub, dims = 1:20)
DimPlot(brain_seu_Pat19_Post_Pat19_Pre_sub, group.by = "MajorCellType", label = TRUE, repel = TRUE)

# insert embedding
embed <- Embeddings(brain_seu_Pat19_Post_Pat19_Pre_sub, reduction = "umap")
vrEmbeddings(Pat19_Post_section, type = "umap_filtered") <- embed[intersect(vrSpatialPoints(Pat19_Post_section),rownames(embed)),]
ind <- intersect(rownames(embed), gsub("Assay1", "Assay2", vrSpatialPoints(Pat19_Pre_section)))
embed2 <- embed[ind,]
rownames(embed2) <- gsub("Assay2", "Assay1", rownames(embed2))
vrEmbeddings(Pat19_Pre_section, type = "umap_filtered") <- embed2

####
# Subclustering ####
####

####
## Subclustering 7_? ####
####

Idents(brain_seu_Pat19_Post_Pat19_Pre) <- "CellType_TumorAnnotation"
brain_seu_Pat19_Post_Pat19_Pre <- FindSubCluster(brain_seu_Pat19_Post_Pat19_Pre, cluster = "7_?", subcluster.name = "subcluster_7_?", 
                                        resolution = 0.2, graph.name = "Xenium_snn")
Idents(brain_seu_Pat19_Post_Pat19_Pre) <- "subcluster_7_?"
m1 <- FindMarkers(object = brain_seu_Pat19_Post_Pat19_Pre, ident.1 = "7_?_2")
# 7_?_2 and 7_?_1 are GFAP+
# saveRDS(brain_seu_Pat19_Post_Pat19_Pre, file = "../data/SeuratData/brain_seu_Pat19_Post_Pat19_Pre_annotated.rds")

Pat19_Post_section$`subcluster_7_?` <- rep(NA, length(vrSpatialPoints(Pat19_Post_section)))
Pat19_Pre_section$`subcluster_7_?` <- rep(NA, length(vrSpatialPoints(Pat19_Pre_section)))
temp <- setNames(Pat19_Post_section$`subcluster_7_?`, vrSpatialPoints(Pat19_Post_section))
Pat19_Post_section$`subcluster_7_?`  <- brain_seu_Pat19_Post_Pat19_Pre$`subcluster_7_?`[vrSpatialPoints(Pat19_Post_section)]
Pat19_Pre_section$`subcluster_7_?`  <- brain_seu_Pat19_Post_Pat19_Pre$`subcluster_7_?`[gsub("Assay1", "Assay2", vrSpatialPoints(Pat19_Pre_section))]

####
# Individual Marker Analysis ####
####

# check markers
Pat19_Post_section <- readRDS(file = "../data/VoltRonData/Pat19_Post_section_annotated.rds")
Pat19_Pre_section <- readRDS(file = "../data/VoltRonData/Pat19_Pre_section_annotated.rds")

# check TAM3 cluster
Pat19_Post_section_subset <- subset(Pat19_Post_section, interactive = TRUE)
Pat19_Post_section_subset <- Pat19_Post_section_subset$subsets[[1]]
g1 <- vrSpatialPlot(Pat19_Post_section_subset, group.by = "CellType_TumorAnnotation", group.ids = "TAMs_3 (RNASET2+)", 
                    plot.segments = TRUE, background = c("black")) + NoLegend()
g2 <- vrSpatialFeaturePlot(Pat19_Post_section_subset, features = c("CD68", "IDH1"), 
                           plot.segments = TRUE, background = c("black"), collapse = FALSE)
g1 | (g2[[1]]/g2[[2]])

####
# CellType Heatmap ####
####

# get data
Pat19_Post_section <- readRDS(file = "../data/VoltRonData/Pat19_Post_section_annotated.rds")
Pat19_Pre_section <- readRDS(file = "../data/VoltRonData/Pat19_Pre_section_annotated.rds")
brain_vr <- merge(Pat19_Post_section, Pat19_Pre_section)
  
# heatmap markers
defining_markers <- read.xlsx("../data/MarkerData/brain_seu_Pat19_Post_Pat19_Pre_supercelltype_topmarkers.xlsx", sheetIndex = 1)

# cell type
celltypes <- unique(brain_vr$MajorCellType)
celltypes <- celltypes[!(celltypes %in% c("Low Count Cells", "undefined", "Astrocyte (CADM1+)", "NES+ Cells"))]
brain_vr <- subset(brain_vr, subset = MajorCellType %in% celltypes)

# change celltype names
temp <- brain_vr$MajorCellType
temp[grepl("Tumor", brain_vr$MajorCellType)] <- "Tumor Cells"
temp[grepl("Proliferating", brain_vr$MajorCellType)] <- "Proliferating Tumor Cells"
brain_vr$SuperCellType <- temp

# normalize and aggregate
brain_vr <- normalizeData(brain_vr, sizefactor = 1000)
brain_vr_data <- vrData(brain_vr, norm = TRUE)
brain_vr_data_agg <- aggregate(t(brain_vr_data), list(brain_vr$Sample, brain_vr$SuperCellType), mean)
sample_heatmap <- brain_vr_data_agg$Group.1
celltype_heatmap <- brain_vr_data_agg$Group.2
brain_vr_data_agg <- brain_vr_data_agg[,-1*c(1,2)]
# brain_vr_data_agg <- t(brain_vr_data_agg)

# scale data
rownames_brain_vr_data_agg <- rownames(brain_vr_data_agg)
brain_vr_data_agg <- apply(brain_vr_data_agg, 2, scale)
rownames(brain_vr_data_agg) <- rownames_brain_vr_data_agg

# get markers
brain_vr_data_agg <- brain_vr_data_agg[,defining_markers$gene]

# visualize
rownames(brain_vr_data_agg) <- paste0(celltype_heatmap, "_", sample_heatmap)
g1 <- Heatmap(brain_vr_data_agg, row_names_rot = 0, column_title_rot = 45, column_names_rot = 45,
              cluster_row_slices = FALSE, cluster_column_slices = FALSE, cluster_rows = FALSE,
              column_split = defining_markers$cluster, 
              show_row_names = TRUE, column_names_side = "bottom")
draw(g1, padding = unit(c(20, 35, 2, 2), "mm"))


####
# annotate Hot tumors ####
####

Pat19_Post_section <- readRDS(file = "../data/VoltRonData/Pat19_Post_section_annotated.rds")

####
## get immunescores ####
####

# get score per modules using Seurat
Pat19_Post_section_assay1_seu <- VoltRon::as.Seurat(Pat19_Post_section, cell.assay = "Xenium", type = "image")
Pat19_Post_section_assay1_seu <- NormalizeData(Pat19_Post_section_assay1_seu, scale.factor = 1000)
for(i in 1:length(score_markers)){
  cur_markers <- intersect(score_markers[[i]], SeuratObject::Features(Pat19_Post_section_assay1_seu))
  Pat19_Post_section_assay1_seu <- AddModuleScore(Pat19_Post_section_assay1_seu, features = list(cur_markers), name = names(score_markers)[i], ctrl = 5)  
}

# transfer module scores to voltron objects
temp <- rep(NA, length(vrSpatialPoints(Pat19_Post_section)))
names(temp) <- vrSpatialPoints(Pat19_Post_section)
temp[names(Pat19_Post_section_assay1_seu$ImmuneScores1)] <- Pat19_Post_section_assay1_seu$ImmuneScores1
Pat19_Post_section$ImmuneScores <- temp

# visualize score
vrSpatialFeaturePlot(Pat19_Post_section, assay = "Xenium", features = c("ImmuneScores"),          
                     plot.segments = FALSE, pt.size = 1, n.tile = 300, alpha = 1)

####
## annotate hot tumor ####  
####

Pat19_Post_section <- annotateSpatialData(Pat19_Post_section, use.image = TRUE, channel = "H&E")
vrSpatialPlot(Pat19_Post_section, group.by = "annotation", n.tile = 300)
# saveRDS(Pat19_Post_section, file = "../data/VoltRonData/Pat19_Post_section_annotated.rds")

####
# save as AnnData ####
####

# interactive TissUUmaps tool
Pat19_Post_Pat19_Pre_section <- merge(Pat19_Post_section, Pat19_Pre_section)
as.AnnData(Pat19_Post_Pat19_Pre_section, file = "../data/AnnDataData/Pat_19.h5ad", flip_coordinates = TRUE, 
           assay = "Xenium", method = "anndata", name = "main", channel = "H&E")
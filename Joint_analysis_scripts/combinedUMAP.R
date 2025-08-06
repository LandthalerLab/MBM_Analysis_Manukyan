# library(VoltRon)
library(Seurat)
library(harmony)
library(dplyr)
library(patchwork)
library(xlsx)
library(ComplexHeatmap)

####
# Import Spatial Data ####
####

# BZW2 expressing TMAs
list_of_samples <- list()
list_of_samples[[1]] <- readRDS("../../data/VoltRonData/brain_vr_Pat23_annotated.rds")
list_of_samples[[8]] <- readRDS("../../data/VoltRonData/brain_vr_Pat90_TMA_annotated.rds")
list_of_samples[[9]] <- readRDS("../../data/VoltRonData/brain_vr_Pat90_MBM_TMA_annotated.rds")
list_of_samples[[3]] <- readRDS("../../data/VoltRonData/brain_vr_Pat15_annotated.rds")
list_of_samples[[6]] <- readRDS("../../data/VoltRonData/brain_vr_Pat19_Post_TMA_annotated.rds")
list_of_samples[[2]] <- readRDS("../../data/VoltRonData/brain_vr_Pat24_annotated.rds")
list_of_samples[[5]] <- readRDS("../../data/VoltRonData/brain_vr_Pat19_Pre_TMA_annotated.rds")
list_of_samples[[10]] <- readRDS("../../data/VoltRonData/brain_vr_Pat3_annotated.rds")
list_of_samples[[4]] <- readRDS("../../data/VoltRonData/brain_vr_CNS1_annotated.rds")
list_of_samples[[11]] <- readRDS("../../data/VoltRonData/brain_vr_CNS2_annotated.rds")
list_of_samples[[7]] <- readRDS("../../data/VoltRonData/brain_vr_Pat6_annotated.rds")

# Essential tumor markers
defining_markers <- read.xlsx("../data/Supplementary Information/Final_defining markers.xlsx", sheetIndex = 1)

# xenium annotations 
xenium_markers <- read.xlsx("../data/Supplementary Information/xenium_brain_annotation.xlsx", sheetIndex = 1)

####
# Process Samples ####
####

# merge samples
brain_vr <- merge(list_of_samples[[1]], list_of_samples[-1])

# filter out some cell types
# Low count, undefined and nes cells are noisy
# IFN cells are noisy
celltypes <- unique(brain_vr$MajorCellType)
celltypes <- celltypes[!(celltypes %in% c("Low Count Cells", "undefined", "NES+ cells", "IFN Cells", "8?", "Astrocytes_1", "Oligodendrocytes (EFHD1+)", "IGFBP3+ cells"))]
brain_vr <- subset(brain_vr, subset = MajorCellType %in% celltypes)

# take out glis3+ cells
`%notin%` <- Negate("%in%")
brain_vr <- subset(brain_vr, subset = CellType_TumorAnnotation %notin% "Astrocytes (GLIS3+?)")

# filter out low counts
brain_vr <- subset(brain_vr, subset = Count > 40)

# change celltype names
temp <- brain_vr$MajorCellType
samples <- brain_vr$Sample
temp[grepl("Tumor", brain_vr$MajorCellType)] <- "Tumor Cells"
temp[grepl("Proliferating", brain_vr$MajorCellType)] <- "Proliferating Tumor Cells"
temp[grepl("Tumor", brain_vr$MajorCellType) & grepl("N3575", samples)] <- paste0(temp[grepl("Tumor", brain_vr$MajorCellType) & grepl("N3575", samples)], " (NGFR+)")
temp[grepl("Astrocytes \\(Tumor\\?\\)", brain_vr$MajorCellType)] <- "Astrocytes (tumor?)"
brain_vr$SuperCellType <- temp

####
# Analyze all samples together ####
####

# visualize unintegrated umap
brain_vr <- normalizeData(brain_vr, method = "LogNorm", sizefactor = 1000)
brain_vr <- getPCA(brain_vr, features = vrFeatures(brain_vr), dims = 20, type = "pca_combined")
brain_vr <- getUMAP(brain_vr, dims = 1:30, data.type = "pca_combined", umap.key = "umap_combined")
vrEmbeddingPlot(brain_vr, embedding = "umap_combined", group.by = "Sample")

####
# Harmony integration ####
####

# harmony based data integration
my_pca <- vrEmbeddings(brain_vr, type = "pca_combined")
my_metadata <- Metadata(brain_vr)
my_harmony_embeddings <- RunHarmony(my_pca, my_metadata, "Sample")
vrEmbeddings(brain_vr, type = "pca_harmony") <- my_harmony_embeddings
brain_vr <- getUMAP(brain_vr, dims = 1:20, data.type = "pca_harmony", umap.key = "umap_harmony")
vrEmbeddingPlot(brain_vr, embedding = "umap_harmony", group.by = "Sample")
vrEmbeddingPlot(brain_vr, embedding = "umap_harmony", group.by = "SuperCellType")

####
# save as AnnData ####
####

# merge cells 
brain_vr <- merge(list_of_samples[[1]], list_of_samples[-1])

# change celltype names
temp <- brain_vr$MajorCellType
samples <- brain_vr$Sample
temp[grepl("Tumor", brain_vr$MajorCellType)] <- "Tumor Cells"
temp[grepl("Tumor", brain_vr$MajorCellType) & grepl("N3575", samples)] <- "Tumor Cells (NGFR+)"
temp[grepl("Proliferating", brain_vr$MajorCellType)] <- "Proliferating Tumor Cells"
temp[grepl("Astrocytes \\(Tumor\\?\\)", brain_vr$MajorCellType)] <- "Astrocytes (tumor?)"
brain_vr$SuperCellType <- temp

# clean metadata
metadata <- brain_vr@metadata@cell
metadata <- metadata[,c("Count", "Sample", "FeatureCount", "MajorCellType", "annotation", "ImmuneScores", "SuperCellType")]
brain_vr@metadata@cell <- metadata

# clean embeddings
embed_names <- vrEmbeddingNames(brain_vr)
embed_names <- embed_names[!embed_names %in% c("pca", "umap")]
for(emb in embed_names)
  vrEmbeddings(brain_vr, type = emb, overwrite = TRUE) <- NULL

# update sample metadata from excel
sample_orig_names <- read.xlsx("../data/Metadata/Patients_Nomenclature.xlsx", sheetIndex = 1)
sample_orig_names <- sample_orig_names[-14,]
rownames(sample_orig_names) <- sample_orig_names$ID
sample_orig_names$Pat.ID <- gsub(" ", "_", sample_orig_names$Pat.ID)
sample_orig_names$Pat.ID <- gsub("\\/", "_", sample_orig_names$Pat.ID)
sample_orig_names$Type <- gsub(" ", "", sample_orig_names$Type)

# interactive TissUUmaps tool
as.AnnData(brain_vr, file = "../data/AnnDataData/TMA.h5ad", flip_coordinates = TRUE, 
           assay = "Xenium", method = "anndata", name = "DAPI", channel = "H&E")
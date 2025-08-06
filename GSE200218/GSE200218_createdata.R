# library(VoltRon)
library(Seurat)
library(Giotto)
library(harmony)
library(dplyr)
library(Matrix)

####
# Import Data ####
####

# genes
gene_names <- read.csv("../../../data/RadkeRedmerSpatial_Brain_05112023/GSE200218/GSE200218_sc_sn_gene_names.csv.gz")
gene_names <- gene_names[,1]

# metadata
metadata <- read.csv("../../../data/RadkeRedmerSpatial_Brain_05112023/GSE200218/GSE200218_sc_sn_metadata.csv", row.names = 1)

# data matrix
datax <- readMM("../../../data/RadkeRedmerSpatial_Brain_05112023/GSE200218/GSE200218_sc_sn_counts.mtx")
rownames(datax) <- make.unique(gene_names)
colnames(datax) <- rownames(metadata)

# integrated data
integ_data <- read.csv("../../../data/RadkeRedmerSpatial_Brain_05112023/GSE200218/GSE200218_sc_sn_integrated_data.csv", row.names = 1)
integ_data <- methods::as(as.matrix(integ_data), 'CsparseMatrix')
colnames(integ_data) <- rownames(metadata)

####
# Create Seurat Object ####
####

seu <- CreateSeuratObject(datax, meta.data = metadata)
seu[["integrated"]] <- CreateAssayObject(integ_data)
# saveRDS(seu, file = "../../../data/RadkeRedmerSpatial_Brain_05112023/GSE200218/seu_integrated.rds")

####
# Run analysis ####
####

DefaultAssay(seu) <- "integrated"
seu <- ScaleData(seu)
seu <- RunPCA(seu, features = Features(seu), npcs = 60)
seu <- RunUMAP(seu, dims = 1:50)
saveRDS(seu, file = "../../../data/RadkeRedmerSpatial_Brain_05112023/GSE200218/seu_integrated.rds")

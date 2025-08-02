library(VoltRon)
library(Seurat)
library(harmony)
library(dplyr)
library(Matrix)
library(xlsx)
library(Seurat)
library(Giotto)
library(harmony)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(tidyverse)
library(magrittr)
library(scCustomize)
library(here)
library(glue)
library(officer)
library(rvg)
library(viridis)
library(dplyr)
library(patchwork)

####
# Import Data ####
####

scMBM <- readRDS ("...directory/....rds")

# score markers (Interferon, Immune etc.)
score_markers <- read.table("...directory/.....txt",sep="\t",h=T,stringsAsFactors = T)
score_markers <- as.list(score_markers)
score_markers <- lapply(score_markers, function(x) x[!is.na(x)]) 

###
# Visualization ####
####

# visualize 
DimPlot_scCustom(scMBM, reduction = "umap", group.by = "malignant", figure_plot = TRUE, pt.size=0.01, raster = F)
FeaturePlot(scMBM, reduction = "umap", features = "gene")

# normalize 
seu <- NormalizeData(seu, assay = "RNA")

# get score per modules using Seurat
for(i in 1:length(score_markers)){
  cur_markers <- intersect(score_markers[[i]], SeuratObject::Features(scMBM))
  scMBM <- AddModuleScore(scMBM, features = list(cur_markers), name = names(score_markers)[i], ctrl = 5, slot = "data", assay = "RNA")  
}

# write 
# write.table(scMBM@meta.data[["name"]],file='scMBM_.....txt', sep=",")

# visualize 
FeaturePlot_scCustom(seurat_object = scMBM, features = c("feature"),figure_plot = TRUE,                        ,
                     colors_use = viridis_plasma_dark_high, na_color = "lightgray")+NoLegend()
DimPlot_scCustom(seurat_object = scMBM,group.by = "cell_type_fine", raster=F)

# reduces resolution to match subsequent analysis 
scMBM <- FindClusters(scMBM, resolution = 0.1) 

# visualize 
DimPlot(, reduction = "umap",label = T)

# feature scatter plot
FeatureScatter(
  scMBM,
  feature1="#1",
  feature2="#2",
  cells = NULL,
  shuffle = FALSE,
  seed = 1,
  group.by = NULL,
  split.by = "malignant",
  cols = c("blue","red"),
  pt.size = 1,
  shape.by = NULL,
  span = NULL,
  smooth = FALSE,
  combine = TRUE,
  slot = "data",
  plot.cor = TRUE,
  ncol = NULL,
  raster = FALSE,
  raster.dpi = c(512, 512),
  jitter = TRUE,
  log = FALSE
)

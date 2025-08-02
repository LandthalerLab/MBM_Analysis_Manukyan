library(VoltRon)
library(Seurat)
library(harmony)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(RColorBrewer)
library(scCustomize)
library(ggplot2)
library(scales)
library(Seurat)
library(dplyr)
library(cowplot)
library(devtools)
library(magrittr)
library(viridis)

####
## Import Spatial Data ####
####

####TMA samples
brain_vr <- readRDS("C:/Users/tredm/Desktop/Spatial manuscript_Project 1_final/Figure 4/Final/Corrected UMAP/brain_vr_combined_umap.rds")

####
## Visualize ####
####

# cell types
celltypes <- unique(brain_vr$SuperCellType)
set.seed(1)
colors <- sample(hue_pal()(length(celltypes)))
names(colors) <- celltypes
g_main <- vrEmbeddingPlot(brain_vr, embedding = "umap_harmony", group.by = "SuperCellType",n.tile = 500, colors = colors)+NoLegend()
# ggsave("Combined_TMA_UMAP_Donor_Annotations_final.pdf", device = "pdf", width = 20, height = 15, dpi = 700)

###Color Code
colors <- DiscretePalette_scCustomize(num_colors = 11, palette = "varibow", shuffle_pal = TRUE)
names(colors) <- unique(brain_vr$SuperCellType)

###Super cell types
vrEmbeddingPlot(brain_vr,assay = "Xenium", embedding = "umap_harmony", group.by = "SuperCellType",n.tile = 600, colors = colors)+NoLegend()

###ImmuneScore
vrEmbeddingFeaturePlot(brain_vr, assay = "Xenium", features = c("ImmuneScores"), combine.features = F, n.tile = 500, embedding = "umap_harmony",pt.size = 2)+ 
  scale_fill_viridis(option="plasma")

###Markers
vrEmbeddingFeaturePlot(brain_vr, assay = "Xenium", features = c("MET"), combine.features = F, n.tile = 500, embedding = "umap_harmony",pt.size = 2)+ 
  scale_fill_viridis(option="plasma")

###Visuaization of Donors/Tumor samples
temp <- brain_vr$Sample
temp <- sapply(temp, function(x) strsplit(x, split = "_")[[1]][1])
brain_vr$Donor <- temp
donors <- unique(brain_vr$Donor)
set.seed(1)
colors <- DiscretePalette_scCustomize(num_colors = 11, palette = "varibow", shuffle_pal = TRUE, length(celltypes))
names(colors) <- donors
vrEmbeddingPlot(brain_vr, embedding = "umap_harmony", group.by = "Donor",pt.size=2, label=F,font.size=3, colors = colors)+NoLegend()

##### SCORE CALCULATION/VISUALIZATION OF TMA SAMAPLES

##ImmuneScore (ESTIMATE) --> used for Figure 4
score_markers <-c("HLA-B","CD74","HLA-DRA","TAP1","HLA-DPB1","CTSS","PTPRC","CD48","LGALS9","FCER1G","IL7R","ITGB2","CD69","GZMB",
                  "CD2","TGFB1","CORO1A","CCL5","CCR7","CD52","KLRB1","NKG7","GNLY")

# get score per modules using Seurat
brain_vr_assay1_seu <- VoltRon::as.Seurat(brain_vr, cell.assay = "Xenium", type = "image")
brain_vr_assay1_seu <- NormalizeData(brain_vr_assay1_seu, scale.factor = 1000)
cur_markers <- intersect(score_markers, SeuratObject::Features(brain_vr_assay1_seu))
brain_vr_assay1_seu <- AddModuleScore(brain_vr_assay1_seu, features = list(cur_markers), name = "ImmuneScore", ctrl = 5)

# transfer module scores to voltron objects
temp <- rep(NA, length(vrSpatialPoints(brain_vr)))
names(temp) <- vrSpatialPoints(brain_vr)
temp[names(brain_vr_assay1_seu$ImmuneScore1)] <- brain_vr_assay1_seu$ImmuneScore1
brain_vr$ImmuneScore <- temp

# visualize score

vrSpatialFeaturePlot(brain_vr, assay = "Assay12", features = "ImmuneScores", background.color = "white",plot.segments = T,pt.size=1)+
  scale_fill_viridis(option="plasma")

vrEmbeddingFeaturePlot(brain_vr, assay = "Xenium", features = c("ImmuneScore"), combine.features = F, n.tile = 500, embedding = "umap_harmony",pt.size = 2)+ 
  scale_fill_viridis(option="plasma")
# ggsave("MajorUMAP_Cellular_ImmuneScore_ntile500.tiff", units="in", width=4, height=3, dpi=3000,background="transparent") 




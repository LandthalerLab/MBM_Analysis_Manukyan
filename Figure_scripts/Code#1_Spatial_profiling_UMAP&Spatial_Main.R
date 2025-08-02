library(VoltRon)
library(Seurat)
library(Giotto)
library(harmony)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(here)
library(glue)
library(officer)
library(rvg)
library(viridis) 
library(devtools)
library(ggrepel)
library(RColorBrewer)
library(tidyverse)
library(magrittr)
library(scCustomize)
library(Polychrome)

####

# get TMA samples 
Xenium_file <- readRDS("brain_vr_T170_annotated.rds")
Xenium_file <- readRDS("brain_vr_H27660_annotated.rds")
Xenium_file <- readRDS("brain_vr_H15340_annotated.rds")

# data
brain_nonTMA_all <- readRDS("brain_nonTMA_all.rds")

## Color palette
colors <- DiscretePalette_scCustomize(num_colors = 20, palette = "varibow", shuffle_pal = T)
names(colors) <- unique(Xenium_file$MajorCellType)

#Visualize
vrEmbeddingPlot(Xenium_file, assay = c("Xenium"), colors=colors, group.by = "MajorCellType", embedding = "umap" ,pt.size = 1.5, label = F, font.size = 1.5)

vrEmbeddingFeaturePlot(Xenium_file, assay = c("Xenium"),n.tile = 600, embedding = "umap",features = c("NGFR"),norm = F,pt.size = 1,combine.features = F)+ 
    scale_fill_viridis(option = "plasma",limits = c(0,5))
  
vrSpatialPlot(Xenium_file, assay = "Xenium", group.by = "IGFBP3",scale.image = FALSE,colors=colors, background.color = "white",plot.segments = T,pt.size=10)

vrSpatialFeaturePlot(Xenium_file, assay = "Assay6",features = c("LOX"),combine.features = F, background.color = "white",plot.segments = TRUE, alpha = 0.6)+
  scale_fill_viridis(option="plasma",limits = c(0,2.5))


######################Single cell presentation of selected areas!!!

# here enter the associated assay number found in SampleMetadata(object)
object_subset <- subset(Xenium_file, assay = "Xenium") 
object_subset2 <- subset(object_subset, interactive = TRUE)
object_subset <- object_subset2$subsets[[1]]

colors <- DiscretePalette_scCustomize(num_colors = (20), palette = "varibow", shuffle_pal = TRUE)#Du musst hier exact die Anzahl Farben angeben, die den CellTypes oder Clusetern entspricht!! Also vorher feststellen
names(colors) <- unique(object_subset$MajorCellType)
vrSpatialPlot(object_subset, assay = "Xenium", group.by = "MajorCellType",pt.size = 0.1, background.color = "white",colors=colors,
                 plot.segments = T,n.tile=600, scale.image = FALSE)+NoLegend()

vrSpatialFeaturePlot(object_subset, assay = "Xenium", features = c("SOX4"),combine.features = FALSE, background.color = "white",plot.segments = F,n.tile=600, scale.image = FALSE)+
scale_fill_viridis(option="plasma")

vrSpatialFeaturePlot(object_subset, features = c("TAP1","SOX4"),combine.features = TRUE, background.color = "white",plot.segments = TRUE,scale.image = FALSE)
                   
#####Marker Analysis of tumor specimen or zoomed areas

brain_seu <- VoltRon::as.Seurat(object_subset,cell.assay = "Xenium", type = "image")
brain_seu <- subset(brain_seu, subset = Count > 10)
brain_seu <- NormalizeData(brain_seu, scale.factor = 1000)

Idents(brain_seu) <- "MajorCellType"
markers <- FindAllMarkers(brain_seu)

topmarkers <- markers %>%
  filter(pct.1 > 0.5, avg_log2FC > 0.5, p_val_adj < 0.01) %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 20)

#group_by("MajorCellType") %>%
#top_n(n = 10, wt = avg_log2FC)
# write.table(as.data.frame(topmarkers),file='File Name_MajorCellTypes_Top20Markers.txt', sep=",")

library(ComplexHeatmap)
marker_features <- unique(topmarkers$gene)
vrHeatmapPlot(object_subset,assay = "Xenium", features = marker_features, group.by = "MajorCellType", 
              show_row_names = TRUE, font.size = 10)

##ImmuneScore (ESTIMATE) --> used for Figure 4
score_markers <-c("HLA-B","CD74","HLA-DRA","TAP1","HLA-DPB1","CTSS","PTPRC","CD48","LGALS9","FCER1G","IL7R","ITGB2","CD69","GZMB",
                  "CD2","TGFB1","CORO1A","CCL5","CCR7","CD52","KLRB1","NKG7","GNLY")

# get score per modules using Seurat
object_subset_assay1_seu <- VoltRon::as.Seurat(object_subset, cell.assay = "Xenium", type = "image")
object_subset_assay1_seu <- NormalizeData(object_subset_assay1_seu, scale.factor = 1000)
cur_markers <- intersect(score_markers, SeuratObject::Features(object_subset_assay1_seu))
object_subset_assay1_seu <- AddModuleScore(object_subset_assay1_seu, features = list(cur_markers), name = "ImmuneScore", ctrl = 5)

# transfer module scores to voltron objects
temp <- rep(NA, length(vrSpatialPoints(object_subset)))
names(temp) <- vrSpatialPoints(object_subset)
temp[names(object_subset_assay1_seu$ImmuneScore1)] <- object_subset_assay1_seu$ImmuneScore1
object_subset$ImmuneScore <- temp

# visualize score
n.tile = 500
vrEmbeddingFeaturePlot(brain_vr, assay = "Xenium", features = c("ImmuneScore"), combine.features = F, n.tile = 500, embedding = "umap_harmony",pt.size = 2)+ 
  scale_fill_viridis(option="plasma")
























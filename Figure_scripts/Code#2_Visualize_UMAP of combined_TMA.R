####
####
# Description: This code produces a UMAP combining of all TMA samples (Supplementary figure 5)
####
####

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
# UMAP of all TMAs ####
####

colors <- DiscretePalette_scCustomize(num_colors = 4, palette = "varibow", shuffle_pal = TRUE)
names(colors) <- unique(brain_vr$Sample)

# data
brain_vr <- readRDS("F:/MBM_TMA profiling/Xenium analyses NEU/VoltRonData/Voltron files/brain_nonTMA_all.rds")
brain_vr <- readRDS("F:/MBM_TMA profiling/Xenium analyses NEU/brain_all_merged.rds")

## visualize samples separately on UMAP
g_main <- vrEmbeddingPlot(brain_vr, embedding = "umap",group.by = "Sample", colors=colors) + geom_point(size = 20)
g_list <- list()
for(i in 1:length(brain_vr@sample.metadata$Sample)){
  g_list[[i]] <- vrEmbeddingPlot(brain_vr, assay = rownames(brain_vr@sample.metadata)[i],
                                 embedding = "umap",pt.size=0.5, group.by = "Sample") +
    scale_color_manual(values = colors[i]) + NoLegend() +
    labs(title = brain_vr@sample.metadata$Sample[i])
}
g_all <- ggarrange(plotlist = g_list, ncol = 4, nrow = 4)
g_main | g_all
# ggsave("Combined_non-TMA_UMAP_non-annoated_final.pdf", device = "pdf", width = 25, height = 10, dpi = 700)







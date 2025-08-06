library(VoltRon)
library(Seurat)
library(harmony)
library(dplyr)
library(ComplexHeatmap)

####
# Get VoltRon objects ####
####

# import (change this path to where the xenium output is)
brain <- importXenium("../../data/RadkeRedmerSpatial_Brain_05112023/output-XETG00046__0011815__Region_1__20231020__132520/", resolution_level = 3, overwrite_resolution = TRUE)
brain <- modulateImage(brain, brightness = 200)

# label samples
samples <- c("H15340-23_6", "H27660_23_3", "N4600_19_1b", "N5609_14_a_B",
             "N5609_14a_A", "N3575_19_1A_C", "N3575_19_1A_B", "N3575_19_1A_A",
             "T170_21_13", "N3864_15_1_B", "N3864_15_1_A", "N5627_15",
             "N3401_15", "H25296_21_1C", "H26352_23_D_B", "H26352_23_D_A")

# get image subsets
image_subsets <- readRDS("data/AuxiliaryData/image_subset_list.rds")[[1]]

# separate samples
brain_list <- list()
for(i in 1:length(image_subsets)){
  brain_list[[i]] <- subset(brain, image = image_subsets[[i]])
  brain_list[[i]]$Sample <- samples[i]
}

# merge samples
brain_vr <- merge(brain_list[[1]], brain_list[-1])

# Visualize a section and its raw counts
vrSpatialFeaturePlot(brain_vr, assay = "Assay9", features = "Count")

####
# Analyze all samples together ####
####

## filter ####
brain_vr <- subset(brain_vr, subset = Count > 5)
brain_vr$FeatureCount <- apply(vrData(brain_vr, norm = FALSE),2,function(x) sum(x>0))
vrScatterPlot(brain_vr, feature.1 = "Count", feature.2 = "FeatureCount")
brain_vr <- normalizeData(brain_vr, method = "LogNorm", sizefactor = 100)

## get variable features ####
brain_vr <- getPCA(brain_vr, features = vrFeatures(brain_vr), dims = 30)
brain_vr <- getUMAP(brain_vr, dims = 1:30)

# merged
brain_vr <- readRDS("data/VoltRonData/brain_all_merged.rds")
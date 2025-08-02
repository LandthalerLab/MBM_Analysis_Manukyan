library(VoltRon)
library(Seurat)
library(Giotto)
library(harmony)
library(dplyr)
library(ComplexHeatmap)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(viridis)
#library(xlsx)

####
# TAP1 vs BZW2 in H25296 and H26352 nonTMA ####
####

# get data
Xenium_file <- readRDS(file = "...directory/.....rds")

####
## visualize ####
####

# visualize TAP1 and BZW2 areas in Sample Pat 90, primary mucosal melanoma and MBM (Figure 2 and Figure S2)
temp <- Xenium_file$MajorCellType
temp[temp %in% c("Tumor Cells 5", "Tumor Cells 6", "Tumor Cells 7")] <- "Tumor Cells (BZW2+)"
temp[temp %in% c("Proliferating Tumor Cells 2")] <- "Tumor Cells (BZW2/TOP2A+)"
temp[temp %in% c("Tumor Cells 4")] <- "Tumor Cells (TAP1+)"
Xenium_file$MajorCellType <- temp
vrSpatialPlot(Xenium_file, assay = "Assay1", group.by = "MajorCellType", pt.size=0.01, 
              group.ids = c("Tumor Cells (BZW2+)", "Tumor Cells (BZW2/TOP2A+)", "Tumor Cells (TAP1+)", "Immune Cells"), 
              background = c("main", "H&E"), alpha = 1, 
              colors = list(`Immune Cells` = "yellow", `Tumor Cells (BZW2/TOP2A+)` = "blue", `Tumor Cells (BZW2+)` = "lightblue", `Tumor Cells (TAP1+)` = "red"))
# ggsave("Spatialplot_BZW2_TAP1_areas.pdf", device = "pdf", width = 20, height = 16, dpi = 700)

####
# Compare BZW2 and other genes across hot and cold ####
####

####
## Import Spatial Data ####
####

# BZW2 expressing TMAs
list_of_samples <- list()
list_of_samples[[1]] <- readRDS("...directory/.....rds")
list_of_samples[[2]] <- readRDS("...directory/.....rds")

####
## get module genes ####
####

score_markers <- readRDS("...directory/.....txt",sep="\t",h=T,stringsAsFactors = T)
score_markers <- as.list(score_markers)
score_markers <- lapply(score_markers, function(x) x[!is.na(x)])

####
## all cells ####
####

####
### aggregate samples ####
####

# aggregate hot and cold tumor
all_agg_data <- NULL
for(samp in list_of_samples){
  
  # sample
  print(unique(samp$Sample))  
  
  # get data
  brain_vr_data <- normalizeData(samp, sizefactor = 1000)
  datax <- vrData(brain_vr_data, norm = FALSE)
  metadata <- Metadata(brain_vr_data)
  
  # hot and cold tumor
  if(!"annotation" %in% colnames(metadata))
    metadata$annotation <- "undefined"
    metadata$annotation[grepl("hot tumor",  metadata$annotation)] <- c("hottumor")
    metadata$annotation[!grepl("hottumor",  metadata$annotation)] <- c("coldtumor")
  
  # aggregate data
  agg_data <- aggregate(t(datax), list(metadata$annotation, metadata$Sample), sum)
  rownames(agg_data) <- paste(agg_data[,1],agg_data[,2], sep = "_")
  agg_data <- t(agg_data[,-1*c(1,2)])
  all_agg_data <- cbind(all_agg_data, agg_data)
}

# normalize for count/size depth, counts per million (CPM)
agg_data_norm <- sweep(all_agg_data, 2, colSums(all_agg_data), FUN = "/")
agg_data_norm <- agg_data_norm*1000000

# create condition variable and sample
condition <- sapply(colnames(agg_data_norm), function(x) strsplit(x, split = "_")[[1]][1])
sample  <- sapply(colnames(agg_data_norm), function(x) {
  temp <- strsplit(x, split = "_")[[1]]
  paste(temp[2:length(temp)], collapse = "_")
})

####
### aggregate scores ####
####

score_agg_data_norm <- NULL
for(i in 1:length(score_markers)){
  cur_markers <- intersect(score_markers[[i]], rownames(agg_data_norm))
  print(cur_markers)
  cur_data <- agg_data_norm[cur_markers,]
  cur_data <- colSums(cur_data)
  score_agg_data_norm <- rbind(score_agg_data_norm, cur_data)
}
rownames(score_agg_data_norm) <- names(score_markers)
# write.table(score_agg_data_norm,file='Spatial_molecular_program_scores_selected_final.txt', sep=",")

####
### compare markers across hot and cold ####
####

markers <- c("TAP1")
agg_data_norm_vis <- t(agg_data_norm[markers,,drop = FALSE])
sample_markers <- "Scores"
score_agg_data_norm_vis <- t(score_agg_data_norm[sample_markers,,drop = FALSE])
agg_data_norm_vis <- data.frame(id = rownames(agg_data_norm_vis), agg_data_norm_vis, score_agg_data_norm_vis, condition = condition, sample = sample)
agg_data_norm_vis <- melt(agg_data_norm_vis)
ggplot(agg_data_norm_vis, aes(x = condition, y = value, color = condition)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point() + 
  stat_compare_means(method = "t.test", label.y = 200,size=3)+
  geom_line(aes(group = sample),
            color = "black", linetype="dashed") + 
  scale_color_manual(values = c("blue", "red")) +
  facet_wrap(~variable, scales = "free_y")

####
## tumor cells only ####
####

####
### aggregate samples ####
####

# aggregate hot and cold tumor
all_agg_data <- NULL
for(samp in list_of_samples){
  
  # sample
  print(unique(samp$Sample))  
  
  # get data
  brain_vr_data <- normalizeData(samp, sizefactor = 1000)
  datax <- vrData(brain_vr_data, norm = FALSE)
  metadata <- Metadata(brain_vr_data)
  
  # hot and cold tumor
  if(!"annotation" %in% colnames(metadata))
    metadata$annotation <- "undefined"
  metadata$annotation[grepl("hot tumor",  metadata$annotation)] <- c("hottumor")
  metadata$annotation[!grepl("hottumor",  metadata$annotation)] <- c("coldtumor")
  
  # separate Tumors
  datax <- datax[,grepl("^Tumor", metadata$MajorCellType) | grepl("IFN Cells", metadata$MajorCellType)]
  metadata <- metadata[grepl("^Tumor", metadata$MajorCellType) | grepl("IFN Cells", metadata$MajorCellType),]
  
  # aggregate data
  agg_data <- aggregate(t(datax), list(metadata$annotation, metadata$Sample), sum)
  rownames(agg_data) <- paste(agg_data[,1],agg_data[,2], sep = "_")
  agg_data <- t(agg_data[,-1*c(1,2)])
  all_agg_data <- cbind(all_agg_data, agg_data)
}


# normalize for count/size depth, counts per million (CPM)
agg_data_norm <- sweep(all_agg_data, 2, colSums(all_agg_data), FUN = "/")
agg_data_norm <- agg_data_norm*1000000

# create condition variable and sample
condition <- sapply(colnames(agg_data_norm), function(x) strsplit(x, split = "_")[[1]][1])
sample  <- sapply(colnames(agg_data_norm), function(x) {
  temp <- strsplit(x, split = "_")[[1]]
  paste(temp[2:length(temp)], collapse = "_")
})

####
### compare markers across hot and cold ####
####

agg_data_norm_vis <- agg_data_norm[c("MX1", "IGFBP5", "IFNB1", "STAT1", "BZW2", "GJA1", "HILPDA", "CTLA4"),]
agg_data_norm_vis <- data.frame(id = colnames(agg_data_norm_vis), melt(t(agg_data_norm_vis)), condition = condition, sample = sample)
colnames(agg_data_norm_vis) <- c("id", "id2", "gene", "value", "condition", "sample")
agg_data_norm_vis <- na.omit(agg_data_norm_vis)
ggplot(agg_data_norm_vis, aes(x = condition, y = value, color = condition)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point() +
  stat_compare_means(method = "t.test", label.y = 1000,size=3)+
  geom_line(aes(group = sample),
            color = "black", linetype="dashed") + 
  scale_color_manual(values = c("blue", "red")) + 
  facet_wrap(~gene, ncol = 4, scales = "free_y") + 
  ylab("Normalized Expression") + 
  xlab("Region")
# ggsave(filename = "figure/hot_vs_cold_expression_barplot.pdf", plot = last_plot(), device = "pdf", width = 11, height = 8)




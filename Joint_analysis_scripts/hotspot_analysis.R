# library(VoltRon)
library(Seurat)
library(harmony)
library(dplyr)
library(ComplexHeatmap)
library(reshape2)
library(xlsx)
library(DESeq2)
library(patchwork)
library(Rfast)

####
# Compare BZW2 and other genes across hot and cold ####
####

####
## Import Spatial Data ####
####

# BZW2 expressing TMAs
list_of_samples <- list()
list_of_samples[[1]] <- readRDS("../../data/VoltRonData/brain_vr_Pat23_annotated.rds")
list_of_samples[[2]] <- readRDS("../../data/VoltRonData/brain_vr_Pat90_TMA_annotated.rds")
list_of_samples[[3]] <- readRDS("../../data/VoltRonData/brain_vr_Pat90_MBM_TMA_annotated.rds")
list_of_samples[[4]] <- readRDS("../../data/VoltRonData/Pat90_all_section_annotated.rds")

####
## get module genes ####
####

score_markers <- read.xlsx("../../../data/Supplementary Information/score_markers.xlsx", sheetName = "Sheet1")
score_markers <- as.list(score_markers)
score_markers <- lapply(score_markers, function(x) x[!is.na(x)])

####
## compute additional scores ####
####

# compute all other than ImmuneScores
score_markers <- score_markers[names(score_markers) != "ImmuneScores"]
for(i in 1:length(list_of_samples)){
  samp <- list_of_samples[[i]]
  print(unique(samp$Sample))
  samp_seu <- VoltRon::as.Seurat(samp, cell.assay = "Xenium", type = "image")
  samp_seu <- NormalizeData(samp_seu, scale.factor = 1000)
  for(j in 1:length(score_markers)){
    cur_markers <- intersect(score_markers[[i]], SeuratObject::Features(samp_seu))
    samp_seu <- AddModuleScore(samp_seu, features = list(cur_markers), name = names(score_markers)[j], ctrl = 5)  
  } 
  
  # transfer module scores to voltron objects
  meta_seu <- samp_seu@meta.data
  colnames(meta_seu) <- gsub("1$", "", colnames(meta_seu))
  meta_seu <- meta_seu[,names(score_markers), drop = FALSE]
  for(col in colnames(meta_seu)){
    samp <- VoltRon::addMetadata(samp, type = "cell", assay = "Xenium", value = meta_seu[,col], label = col)
  }
  
  list_of_samples[[i]] <- samp
}

####
## find automated hot spots ####
####

####
### N3575 ####
####

# calculate hot spots for immune score
list_of_samples[[1]] <- getSpatialNeighbors(list_of_samples[[1]], radius = 40, method = "radius")
list_of_samples[[1]] <- getHotSpotAnalysis(list_of_samples[[1]], features = "ImmuneScores", graph.type = "radius", alpha.value = 0.05)

# visualize
g1 <- vrSpatialFeaturePlot(list_of_samples[[1]], features = c("BZW2"), plot.segments = TRUE, alpha = 1, ncol = 1)
g2 <- vrSpatialFeaturePlot(list_of_samples[[1]], features = c("ImmuneScores"), plot.segments = TRUE, alpha = 1, ncol = 1)
g3 <- vrSpatialPlot(list_of_samples[[1]], group.by = c("ImmuneScores_hotspot_flag"), plot.segments = TRUE, alpha = 1, ncol = 1)
g4 <- vrSpatialPlot(list_of_samples[[1]], group.by = c("annotation"), plot.segments = TRUE, alpha = 1, ncol = 1)
g1 | g2 | g3 | g4

####
### H25296 ####
####

# calculate hot spots for immune score
list_of_samples[[2]] <- getSpatialNeighbors(list_of_samples[[2]], radius = 40, method = "radius")
list_of_samples[[2]] <- getHotSpotAnalysis(list_of_samples[[2]], features = "ImmuneScores", graph.type = "radius", alpha.value = 0.05)

# visualize
g1 <- vrSpatialFeaturePlot(list_of_samples[[2]], features = c("BZW2"), plot.segments = TRUE, alpha = 1, ncol = 1)
g2 <- vrSpatialFeaturePlot(list_of_samples[[2]], features = c("ImmuneScores"), plot.segments = TRUE, alpha = 1, ncol = 1)
g3 <- vrSpatialPlot(list_of_samples[[2]], group.by = c("ImmuneScores_hotspot_flag"), plot.segments = TRUE, alpha = 1, ncol = 1)
g4 <- vrSpatialPlot(list_of_samples[[2]], group.by = c("annotation"), plot.segments = TRUE, alpha = 1, ncol = 1)
(g1 | g2) / (g3 | g4)

####
### H26352 ####
####

# calculate hot spots for immune score
list_of_samples[[3]] <- getSpatialNeighbors(list_of_samples[[3]], assay = "Assay15", radius = 50, method = "radius")
list_of_samples[[3]] <- getHotSpotAnalysis(list_of_samples[[3]], assay = "Assay15", features = "ImmuneScores", graph.type = "radius", alpha.value = 0.05)

# visualize
g1 <- vrSpatialFeaturePlot(list_of_samples[[3]], assay = "Assay15", features = c("BZW2"), plot.segments = TRUE, alpha = 1, ncol = 1)
g2 <- vrSpatialFeaturePlot(list_of_samples[[3]], assay = "Assay15", features = c("ImmuneScores"), plot.segments = TRUE, alpha = 1, ncol = 1)
g3 <- vrSpatialPlot(list_of_samples[[3]], assay = "Assay15", group.by = c("ImmuneScores_hotspot_flag"), plot.segments = TRUE, alpha = 1, ncol = 1)
g4 <- vrSpatialPlot(list_of_samples[[3]], assay = "Assay15", group.by = c("annotation"), plot.segments = TRUE, alpha = 1, ncol = 1)
(g1 | g2) / (g3 | g4)


####
### H2526 ####
####

# calculate hot spots for immune score
list_of_samples[[4]] <- getSpatialNeighbors(list_of_samples[[4]], radius = 50, method = "radius")
list_of_samples[[4]] <- getHotSpotAnalysis(list_of_samples[[4]], features = "ImmuneScores", graph.type = "radius")

# visualize
g1 <- vrSpatialFeaturePlot(list_of_samples[[4]],  features = c("BZW2"), plot.segments = FALSE, alpha = 1, ncol = 1, n.tile = 300)
g2 <- vrSpatialFeaturePlot(list_of_samples[[4]],  features = c("ImmuneScores"), plot.segments = FALSE, alpha = 1, ncol = 1, n.tile = 300)
g3 <- vrSpatialPlot(list_of_samples[[4]], group.by = c("ImmuneScores_hotspot_flag"), plot.segments = FALSE, alpha = 1, ncol = 1, n.tile = 300)
g4 <- vrSpatialPlot(list_of_samples[[4]], group.by = c("annotation"), plot.segments = FALSE, alpha = 1, ncol = 1, n.tile = 300)
(g1 | g2 | g3 | g4)

####
### merge samples ####
####

# merge data
list_of_samples[[5]] <- list_of_samples[[4]]
list_of_samples[[5]] <- subset(list_of_samples[[5]], assays = "Xenium")
brain_vr <- merge(list_of_samples[[1]], list_of_samples[c(2,3,5)])

# change celltype names
temp <- brain_vr$MajorCellType
samples <- brain_vr$Sample
temp[grepl("Tumor", brain_vr$MajorCellType)] <- "Tumor Cells"
temp[grepl("Proliferating", brain_vr$MajorCellType)] <- "Proliferating Tumor Cells"
temp[grepl("Tumor", brain_vr$MajorCellType) & grepl("N3575", samples)] <- paste0(temp[grepl("Tumor", brain_vr$MajorCellType) & grepl("N3575", samples)], " (NGFR+)")
# temp[grepl("Tumor", brain_vr$MajorCellType) & grepl("N3575", samples)] <- "Tumor Cells (NGFR+)"
temp[grepl("Astrocytes \\(Tumor\\?\\)", brain_vr$MajorCellType)] <- "Astrocytes (tumor?)"
brain_vr$SuperCellType <- temp

# tag hot tumor
tmp <- brain_vr$ImmuneScores_hotspot_flag
tmp[is.na(tmp)] <- "cold"
tmp <- ifelse(tmp == "hot","hot","cold")
brain_vr$Tag_HotTumorCell <- ifelse(tmp == "hot" & (grepl("^Tumor", brain_vr$MajorCellType) | grepl("IFN Cells", brain_vr$MajorCellType)),
                                "Hot_TumorCell", 
                                ifelse((grepl("^Tumor", brain_vr$MajorCellType) | grepl("IFN Cells", brain_vr$MajorCellType)),
                                       "Tumor", "others"))
table(brain_vr$Sample,brain_vr$Tag_HotTumorCell)

# subset tumor
spatialpoints <- vrSpatialPoints(brain_vr)
brain_vr_subset <- subset(brain_vr, spatialpoints = spatialpoints[grepl("Tumor", brain_vr$SuperCellType)])

# table
tmp <- brain_vr_subset$ImmuneScores_hotspot_flag
tmp[is.na(tmp)] <- "cold"
brain_vr_subset$Sample_hot <- paste0(brain_vr_subset$Sample, "_", ifelse(tmp == "hot","hot","cold"))
table(brain_vr_subset$SuperCellType, brain_vr_subset$Sample_hot)

####
## visualize tumor, hot tumors and other cells ####
####

# TMA
assays <- c("Assay1", "Assay2", "Assay3", "Assay7", "Assay9")
g1 <- vrSpatialFeaturePlot(brain_vr, assay = assays, features = c("BZW2"), plot.segments = TRUE, alpha = 1, ncol = 1)
g2 <- vrSpatialFeaturePlot(brain_vr, assay = assays, features = c("ImmuneScores"), plot.segments = TRUE, alpha = 1, ncol = 1)
g3 <- vrSpatialPlot(brain_vr, assay = assays, group.by = c("ImmuneScores_hotspot_flag"), plot.segments = TRUE, alpha = 1, ncol = 1)
g5 <- vrSpatialPlot(brain_vr, assay = assays, group.by = c("Tag_HotTumorCell"), plot.segments = TRUE, alpha = 1, ncol = 1,
                    background.color = "black", colors = list(Hot_TumorCell = "blue", others = "grey", Tumor = "yellow"))
ggpubr::ggarrange(plotlist = list(g1, g2, g3, g5), nrow = 1, ncol = 4, widths = c(1,1,1.5,1))

# Non-TMA
assays <- c("Assay12", "Assay13")
g1 <- vrSpatialFeaturePlot(brain_vr, assay = assays, features = c("BZW2"), alpha = 1, ncol = 1, n.tile = 300)
g2 <- vrSpatialFeaturePlot(brain_vr, assay = assays, features = c("ImmuneScores"), alpha = 1, ncol = 1, n.tile = 300)
g3 <- vrSpatialPlot(brain_vr, assay = assays, group.by = c("ImmuneScores_hotspot_flag"), alpha = 1, ncol = 1, n.tile = 300)
g5 <- vrSpatialPlot(brain_vr, assay = assays, group.by = c("Tag_HotTumorCell"), alpha = 1, ncol = 1, , n.tile = 300,
                    background.color = "black", colors = list(Hot_TumorCell = "blue", others = "grey", Tumor = "yellow"))
ggpubr::ggarrange(plotlist = list(g1, g2, g3, g5), nrow = 1, ncol = 4, widths = c(1,1,1.5,1))

# TMA
vrSpatialFeaturePlot(brain_vr, assay = "Assay3", features = c("ImmuneScores"), plot.segments = TRUE, alpha = 1, ncol = 1) + 
  scale_fill_viridis_c()

####
## aggregate samples for testing ####
####

####
### all cells ####
####

# aggregate hot and cold tumor
all_agg_data <- NULL
for(samp in list_of_samples[c(1,2,3,5)]){
  
  # sample
  print(unique(samp$Sample))  
  
  # get data
  datax <- vrData(brain_vr, norm = FALSE)
  metadata <- Metadata(brain_vr)
  
  # hot and cold tumor
  tmp <- metadata$ImmuneScores_hotspot_flag
  tmp[is.na(tmp)] <- "cold"
  metadata$ImmuneScores_hotspot_flag <- ifelse(tmp == "hot","hot","cold")
  
  # aggregate data
  agg_data <- aggregate(t(datax), list(metadata$ImmuneScores_hotspot_flag, metadata$Sample), sum)
  rownames(agg_data) <- paste(agg_data[,1],agg_data[,2], sep = "_")
  agg_data <- t(agg_data[,-1*c(1,2)])
  all_agg_data <- cbind(all_agg_data, agg_data)
}
write.table(all_agg_data, file = "data/aggregated_rawcounts_allcells.txt", quote = FALSE, row.names = TRUE)
all_agg_data <- read.table("data/aggregated_rawcounts_allcells.txt")

# normalize for count/size depth, counts per million (CPM)
agg_data_norm <- sweep(all_agg_data, 2, colSums(all_agg_data), FUN = "/")
agg_data_norm <- log2(agg_data_norm*10000)

# create condition variable and sample
condition <- sapply(colnames(agg_data_norm), function(x) strsplit(x, split = "_")[[1]][1])
sample  <- sapply(colnames(agg_data_norm), function(x) {
  temp <- strsplit(x, split = "_")[[1]]
  paste(temp[2:length(temp)], collapse = "_")
})

####
### tumor cells ####
####

# aggregate hot and cold tumor
all_agg_data_tumor <- NULL
for(samp in list_of_samples[c(1,2,3,5)]){
  
  # sample
  print(unique(samp$Sample))  
  
  # get data
  datax <- vrData(brain_vr, norm = FALSE)
  metadata <- Metadata(brain_vr)
  
  # hot and cold tumor
  tmp <- metadata$ImmuneScores_hotspot_flag
  tmp[is.na(tmp)] <- "cold"
  metadata$ImmuneScores_hotspot_flag <- ifelse(tmp == "hot","hot","cold")
  
  # separate Tumors
  datax <- datax[,grepl("^Tumor", metadata$MajorCellType) | grepl("IFN Cells", metadata$MajorCellType)]
  metadata <- metadata[grepl("^Tumor", metadata$MajorCellType) | grepl("IFN Cells", metadata$MajorCellType),]
  
  # aggregate data
  agg_data <- aggregate(t(datax), list(metadata$ImmuneScores_hotspot_flag, metadata$Sample), sum)
  rownames(agg_data) <- paste(agg_data[,1],agg_data[,2], sep = "_")
  agg_data <- t(agg_data[,-1*c(1,2)])
  all_agg_data_tumor <- cbind(all_agg_data_tumor, agg_data)
}

# save
write.table(all_agg_data_tumor, file = "data/aggregated_rawcounts_tumorcells.txt", quote = FALSE, row.names = TRUE)
all_agg_data_tumor <- read.table("data/aggregated_rawcounts_tumorcells.txt")

# normalize for count/size depth, counts per million (CPM)
agg_data_tumor_norm <- sweep(all_agg_data_tumor, 2, colSums(all_agg_data_tumor), FUN = "/")
agg_data_tumor_norm <- log2(agg_data_tumor_norm*10000)

####
## compare markers across hot and cold ####
####

####
## set markers for tumor or all cells ####
####

# comparison list
agg_data_tumor_norm_vis <- agg_data_tumor_norm[c("TAP1", "BZW2", "SOX4", "PDCD1LG2","BCAN","IGFBP3","SOX6"),]    
agg_data_norm_vis <- agg_data_norm[c("HLA-DRA", "PSENEN", "MX1", "LGALS9","NKG7","CCL5","CD8A","ITGB2","ITGB7","C1QC","MET"),]

####
## visualize hot vs cold tumor ####
####

agg_data_norm_vis <- rbind(agg_data_tumor_norm_vis,agg_data_norm_vis)
agg_data_norm_vis <- data.frame(id = colnames(agg_data_norm_vis), reshape2::melt(t(agg_data_norm_vis)), condition = condition, sample = sample)
colnames(agg_data_norm_vis) <- c("id", "id2", "gene", "value", "condition", "sample")
agg_data_norm_vis <- na.omit(agg_data_norm_vis)
ggplot(agg_data_norm_vis, aes(x = condition, y = value, color = condition)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point() +
  geom_line(aes(group = sample),
            color = "black", linetype="dashed") + 
  scale_color_manual(values = c("blue", "red", "grey")) + 
  facet_wrap(~gene, ncol = 8, scales = "free_y") + 
  ylab("Normalized Expression") + 
  xlab("Region")
# ggsave(filename = "figure/hot_vs_cold_expression_barplot.pdf", plot = last_plot(), device = "pdf", width = 11, height = 8)

# violin plot
brain_vr_subset <- subset(brain_vr, Tag_HotTumorCell %in% c("Hot_TumorCell", "Tumor"))
samples <- unique(brain_vr_subset$Sample)
g_list <- list()
for(i in 1:length(samples)){
  print(samples[i])
  brain_vr_temp <- subset(brain_vr_subset, samples = samples[i])
  brain_vr_temp <- subset(brain_vr_temp, Tag_HotTumorCell %in% c("Hot_TumorCell", "Tumor"))
  g_list[[samples[i]]] <- vrViolinPlot(brain_vr_temp, group.by = "Tag_HotTumorCell", features = c("BZW2", "TAP1"), plot.points = TRUE, pt.size = 0.4)
}
ggpubr::ggarrange(plotlist = g_list, nrow = 2, ncol = 4, labels = samples)

####
## DE analysis across hot and cold ####
####

####
### pseudo-bulk level ####
####

# get regular deseq2 results
# all cells
dds <- DESeqDataSetFromMatrix(countData = all_agg_data,
                              DataFrame(condition, sample),
                              design= ~ condition + sample)
dds <- DESeq(dds)
results_data <- as.data.frame(results(dds, contrast = c("condition", "hot", "cold")))

# tumor
dds <- DESeqDataSetFromMatrix(countData = all_agg_data_tumor,
                              DataFrame(condition, sample),
                              design= ~ condition + sample)
dds <- DESeq(dds)
results_data_tumor <- as.data.frame(results(dds, contrast = c("condition", "hot", "cold")))

####
### single cell level ####
####

# in single cell resolution
brain_vr_subset <- subset(brain_vr, Tag_HotTumorCell %in% c("Hot_TumorCell", "Tumor"))
samples <- unique(brain_vr_subset$Sample)
allmarkers <- NULL
for(i in 1:length(samples)){
  print(samples[i])
  brain_vr_temp <- subset(brain_vr_subset, samples = samples[i])
  brain_vr_temp_seu <- CreateSeuratObject(vrData(brain_vr_temp, norm = FALSE), meta.data = Metadata(brain_vr_temp))
  brain_vr_temp_seu <- Seurat::NormalizeData(brain_vr_temp_seu, scale.factor = 1000)
  Idents(brain_vr_temp_seu) <- "Tag_HotTumorCell"
  if(length(unique(Idents(brain_vr_temp_seu))) > 1){
    markers <- FindAllMarkers(brain_vr_temp_seu, logfc.threshold = 0, min.pct = 0)
    allmarkers <- rbind(allmarkers,
                        data.frame(markers, sample = samples[i])) 
  }
}
allmarkers <- allmarkers[allmarkers$cluster == "Hot_TumorCell",]

####
## Proximity Analysis ####
####

# marker analysis
brain_vr_temp <- subset(brain_vr, assays = "Assay12")
brain_vr_temp_seu <- CreateSeuratObject(vrData(brain_vr_temp, norm = FALSE), meta.data = Metadata(brain_vr_temp))
brain_vr_temp_seu <- Seurat::NormalizeData(brain_vr_temp_seu, scale.factor = 1000)
Idents(brain_vr_temp_seu) <- "MajorCellType"
markers <- FindAllMarkers(brain_vr_temp_seu, logfc.threshold = 0, min.pct = 0)

# cells
brain_vr_subset <- subset(brain_vr, assays = "Assay12")
brain_vr_subset$hotcoldmarker <- ifelse(brain_vr_subset$MajorCellType %in% c("Proliferating Tumor Cells 2", "Tumor Cells 5", "Tumor Cells 6", "Tumor Cells 7"),
  "BZW2+ Tumor Cells",
  ifelse(brain_vr_subset$MajorCellType %in% c("Tumor Cells 2", "Tumor Cells 4"),
         "TAP1+ Tumor Cells",
         ifelse(brain_vr_subset$MajorCellType == "Immune Cells", "Immune Cells", "Other Cells")))

# get coords
coords <- vrCoordinates(brain_vr, assay = "Assay12")
coords_TAP1 <- coords[brain_vr_subset$hotcoldmarker == "TAP1+ Tumor Cells",]
coords_BZW2 <- coords[brain_vr_subset$hotcoldmarker == "BZW2+ Tumor Cells",]
coords_Immune <- coords[brain_vr_subset$hotcoldmarker == "Immune Cells",]

# compute distance
coords_TAP1 <- na.omit(coords_TAP1)
coords_BZW2 <- na.omit(coords_BZW2)
coords_Immune <- na.omit(coords_Immune)
dist1_TAP1 <- Rfast::dista(coords_TAP1, coords_Immune)
dist1_TAP1_min <- apply(dist1_TAP1, 1, min)
dist1_BZW2_1 <- Rfast::dista(coords_BZW2[1:80000,], coords_Immune)
dist1_BZW2_1_min <- apply(dist1_BZW2_1, 1, min)
dist1_BZW2_2 <- Rfast::dista(coords_BZW2[80001:nrow(coords_BZW2),], coords_Immune)
dist1_BZW2_2_min <- apply(dist1_BZW2_2, 1, min)
dist1_BZW2_min <- c(dist1_BZW2_1_min, dist1_BZW2_2_min)

# plot distances
plotdata <- data.frame(dist = c(dist1_TAP1_min, dist1_BZW2_min),
                       group = c(rep("TAP1+ Tumor Cells", length(dist1_TAP1_min)),
                                 rep("BZW2+ Tumor Cells", length(dist1_BZW2_min))))
write.table(plotdata, file = "data/proximity_analysis_BZW2_TAP1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
plotdata <- read.table("data/proximity_analysis_BZW2_TAP1.txt", header = TRUE, sep = "\t")
plotdata$dist <- plotdata$dist * 0.4250 # resolution=2 pixel to micron ratio
plotdata$dist <- log2(plotdata$dist)
g1 <- ggplot(data = plotdata, aes(x = group, y = dist, fill = group)) + 
  geom_violin() + 
  labs(title = "Distance to Nearest Immune Cell")

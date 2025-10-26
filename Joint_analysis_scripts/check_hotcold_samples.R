# library(VoltRon)
library(Seurat)
library(harmony)
library(dplyr)
library(ComplexHeatmap)
library(xlsx)

####
# Import Spatial Data ####
####

list_of_samples <- list()
list_of_samples[[1]] <- readRDS("../../data/VoltRonData/brain_vr_Pat23_annotated.rds")
list_of_samples[[2]] <- readRDS("../../data/VoltRonData/brain_vr_Pat24_annotated.rds")
list_of_samples[[3]] <- readRDS("../../data/VoltRonData/brain_vr_Pat15_annotated.rds")
list_of_samples[[4]] <- readRDS("../../data/VoltRonData/brain_vr_CNS1_annotated.rds")
list_of_samples[[5]] <- readRDS("../../data/VoltRonData/brain_vr_Pat19_Pre_TMA_annotated.rds")
list_of_samples[[6]] <- readRDS("../../data/VoltRonData/brain_vr_Pat19_Post_TMA_annotated.rds")
list_of_samples[[7]] <- readRDS("../../data/VoltRonData/brain_vr_Pat6_annotated.rds")
list_of_samples[[8]] <- readRDS("../../data/VoltRonData/brain_vr_Pat90_TMA_annotated.rds")
list_of_samples[[9]] <- readRDS("../../data/VoltRonData/brain_vr_Pat90_MBM_TMA_annotated.rds")
list_of_samples[[10]] <- readRDS("../../data/VoltRonData/brain_vr_Pat3_annotated.rds")
list_of_samples[[11]] <- readRDS("../../data/VoltRonData/brain_vr_CNS2_annotated.rds")
list_of_samples[[12]] <- readRDS("../../data/VoltRonData/Pat19_Pre_section_annotated.rds")
list_of_samples[[13]] <- readRDS("../../data/VoltRonData/Pat19_Post_section_annotated.rds")
list_of_samples[[14]] <- readRDS("../../data/VoltRonData/Pat90_all_section_annotated.rds")

####
# get module genes ####
####

score_markers <- read.xlsx("../../data/Supplementary Information/score_markers.xlsx", sheetName = "Sheet1")
score_markers <- as.list(score_markers)
score_markers <- lapply(score_markers, function(x) x[!is.na(x)])

####
# aggregate scores ####
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
write.csv(agg_data_norm, file = "data/hotcoldtumor_aggregate_norm.csv", quote = FALSE, row.names = TRUE)
agg_data_norm <- read.csv("data/hotcoldtumor_aggregate_norm.csv", row.names = 1)

# create condition variable and sample
condition <- sapply(colnames(agg_data_norm), function(x) strsplit(x, split = "_")[[1]][1])
sample  <- sapply(colnames(agg_data_norm), function(x) {
  temp <- strsplit(x, split = "_")[[1]]
  paste(temp[2:length(temp)], collapse = "_")
})

# score data
score_agg_data_norm <- NULL
for(i in 1:length(score_markers)){
  cur_markers <- intersect(score_markers[[i]], rownames(agg_data_norm))
  print(cur_markers)
  cur_data <- agg_data_norm[cur_markers,]
  cur_data <- colSums(cur_data)
  score_agg_data_norm <- rbind(score_agg_data_norm, cur_data)
}

####
## compare BZW2 across samples ####
####

# aggregate samples
agg_data_sample <- aggregate(t(all_agg_data), list(sample), sum)
rownames(agg_data_sample) <- agg_data_sample[,1]
agg_data_sample <- t(agg_data_sample[,-1])

# normalize for count/size depth, counts per million (CPM)
agg_data_sample_norm <- sweep(agg_data_sample, 2, colSums(agg_data_sample), FUN = "/")
agg_data_sample_norm <- agg_data_sample_norm*1000000

# BZW2 expression across samples
agg_data_sample_norm_vis <- agg_data_sample_norm["BZW2",]
agg_data_sample_norm_vis <- data.frame(id = names(agg_data_sample_norm_vis), value = melt(agg_data_sample_norm_vis))
agg_data_sample_norm_vis$level <- ifelse(agg_data_sample_norm_vis$value < 2500, "low (BZW2)", 
                                  ifelse(agg_data_sample_norm_vis$value < 6000, "medium (BZW2)", "high (BZW2)"))
ggplot(agg_data_sample_norm_vis, aes(x = id, y = value, fill = level)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust = 1))+
  geom_hline(yintercept = 2500) +
  geom_hline(yintercept = 6000) + 
  ylab("Normalized BZW2 Expression") + 
  xlab("Sample")
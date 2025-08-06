# library(VoltRon)
library(Seurat)
library(harmony)
library(dplyr)
library(Matrix)
library(xlsx)
library(ggplot2)
library(ComplexHeatmap)

####
# Import Data ####
####

seu <- readRDS("../../../../data/RadkeRedmerSpatial_Brain_05112023/GSE200218/seu_integrated.rds")

# score markers (Interferon, Immune etc.)
score_markers <- read.xlsx("../../data/Supplementary Information/score_markers.xlsx", sheetName = "Sheet1")
score_markers <- as.list(score_markers)
score_markers <- lapply(score_markers, function(x) x[!is.na(x)]) 

####
# Visualization ####
####

DimPlot(seu, reduction = "umap", group.by = "cell_type_fine")
FeaturePlot(seu, reduction = "umap", features = "ITGB7")

####
# Interferon scoring ####
####

# normalize 
seu <- NormalizeData(seu, assay = "RNA")

# get score per modules using Seurat
for(i in 1:length(score_markers)){
  cur_markers <- intersect(score_markers[[i]], SeuratObject::Features(seu))
  seu <- AddModuleScore(seu, features = list(cur_markers), name = names(score_markers)[i], ctrl = 5, slot = "data", assay = "RNA")  
}

g1 <- DimPlot(seu, reduction = "umap", group.by = "cell_type_fine", label = TRUE) + NoLegend()
g2 <- FeaturePlot(seu, reduction = "umap", features = "ImmuneScores1")

####
# Tumor cell subclustering ####
####

# compare features
DefaultAssay(seu) <- "RNA"
Seurat::FeatureScatter(seu, feature1 = "ImmuneScores1", feature2 = "BZW2", slot = "data")
g1 <- FeaturePlot(seu, features = c("ImmuneScores1"))
g2 <- FeaturePlot(seu, features = c("BZW2"))
VlnPlot(seu, features = c("TAP1", "BZW2"), group.by = "cell_type_fine")

# subset tumor cells 
seu_tumor_integrated <- subset(seu, subset = cell_type_fine == "Tumor cells")
DefaultAssay(seu_tumor_integrated) <- "integrated"
seu_tumor_integrated <- ScaleData(seu_tumor_integrated)
seu_tumor_integrated <- RunPCA(seu_tumor_integrated, npcs = 30, features = Features(seu_tumor_integrated))
ElbowPlot(seu_tumor_integrated)
seu_tumor_integrated <- RunUMAP(seu_tumor_integrated, dims = 1:20)
seu_tumor_integrated <- NormalizeData(seu_tumor_integrated, assay = "RNA")
g1 <- FeaturePlot(seu_tumor_integrated, features = c("TAP1"), pt.size = 1.5, slot = "data", order = TRUE)
g2 <- FeaturePlot(seu_tumor_integrated, features = c("BZW2"), pt.size = 1.5, slot = "data", order = TRUE)
g1 | g2

# clustering
seu_tumor_integrated <- FindNeighbors(seu_tumor_integrated, dims = 1:10)
seu_tumor_integrated <- FindClusters(seu_tumor_integrated, resolution = c(0.5, 0.6, 0.7, 0.8, 0.9))
DimPlot(seu_tumor_integrated, reduction = "umap", group.by = "integrated_snn_res.0.5", label = TRUE) + NoLegend()
saveRDS(seu_tumor_integrated, file = "seu_tumor_integrated.rds")

# visualize
g1 <- DimPlot(seu_tumor_integrated, reduction = "umap", group.by = "integrated_snn_res.0.5", label = TRUE) + NoLegend()
g2 <- VlnPlot(seu_tumor_integrated, features = c("TAP1"), slot = "data", group.by = "integrated_snn_res.0.5")
g3 <- VlnPlot(seu_tumor_integrated, features = c("BZW2"), slot = "data", group.by = "integrated_snn_res.0.5")
# g1 | g2 | g3
g4 <- FeaturePlot(seu_tumor_integrated, features = c("TAP1"), pt.size = 1.5, slot = "data", order = TRUE)
g5 <- FeaturePlot(seu_tumor_integrated, features = c("BZW2"), pt.size = 1.5, slot = "data", order = TRUE)
(g1 | g2 | g3)/(g4 | g5)

g4 <- FeaturePlot(seu_tumor_integrated, features = c("TAP1"), pt.size = 1.5, slot = "counts", order = TRUE)
g5 <- FeaturePlot(seu_tumor_integrated, features = c("BZW2"), pt.size = 1.5, slot = "counts", order = TRUE)
(g4 | g5)

# marker analysis
Idents(seu_tumor_integrated) <- "integrated_snn_res.0.5"
DefaultAssay(seu_tumor_integrated) <- "RNA"
markers <- FindAllMarkers(seu_tumor_integrated, min.pct = 0, logfc.threshold = 0)
write.table(markers, file = "tumorcells_integrated_snn_res.0.5_markers.tsv", sep = "\t", quote = FALSE)

####
# Tumor cell subclustering (no integration) ####
####

seu_tumor <- subset(seu, subset = cell_type_fine == "Tumor cells")
DefaultAssay(seu_tumor) <- "RNA"
seu_tumor <- NormalizeData(seu_tumor)
seu_tumor <- FindVariableFeatures(seu_tumor)
seu_tumor <- ScaleData(seu_tumor)
seu_tumor <- RunPCA(seu_tumor, npcs = 30)
ElbowPlot(seu_tumor)
seu_tumor <- RunUMAP(seu_tumor, dims = 1:30)
DimPlot(seu_tumor, reduction = "umap", group.by = "ID", label = FALSE)

# visualize
g1 <- FeaturePlot(seu_tumor, features = c("TAP1"), pt.size = 1.5, slot = "data", order = TRUE)
g2 <- FeaturePlot(seu_tumor, features = c("BZW2"), pt.size = 1.5, slot = "data", order = TRUE)
g1 | g2


####
## thresholding ####
####

# data
datax <- GetAssayData(seu_tumor, assay = "RNA", layer = "data")

# visualize markers
plot_data <- as.data.frame(t(datax[c("BZW2", "TAP1"),]))
plot_data <- reshape2::melt(plot_data)
ggplot() + 
  geom_point(mapping = aes(x = variable, y = value), data = plot_data, position = position_jitter())

# thresholds
thresholds <- apply(plot_data, 2, function(x) {
  temp <- x[x!=0]
  quantile(temp, probs = 0.75)
})

####
### investigate thresholding ####
####

plot_data2 <- plot_data
plot_data2[,1] <- plot_data2[,1] > thresholds[1]
plot_data2[,2] <- plot_data2[,2] > thresholds[2]
table(plot_data2$BZW2, plot_data2$TAP1)

plot_data2 <- plot_data
plot_data2[,1] <- cut(plot_data2[,1], breaks = 10)
plot_data2[,2] <- cut(plot_data2[,2], breaks = 10)
tablex <- as.matrix(rbind(table(plot_data2$BZW2, plot_data2$TAP1)))
ComplexHeatmap::Heatmap(log1p(tablex), cluster_columns = FALSE, cluster_rows = FALSE)

####
### flag TAP1+ and BZW2+ regions ####
####

seu_tumor$hotcoldregions <- "undefined"
seu_tumor$hotcoldregions[plot_data2[,1] & plot_data2[,2]] <- "double"
seu_tumor$hotcoldregions[!plot_data2[,1] & plot_data2[,2]] <- "TAP1"
seu_tumor$hotcoldregions[plot_data2[,1] & !plot_data2[,2]] <- "BZW2"

# marker analysis
seu_tumor_test <- subset(seu_tumor, subset = hotcoldregions!="double")
Idents(seu_tumor_test) <- "hotcoldregions"
DefaultAssay(seu_tumor_test) <- "RNA"
markers <- FindAllMarkers(seu_tumor_test, min.pct = 0.3, logfc.threshold = 0)
w
### Based on https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf ###
library( "EnhancedVolcano" )
library( "GenomicFeatures" )
library( "Rsamtools" )
library( "GenomicAlignments" )
library( "DESeq2")
library( "BiocParallel" )
library( "AnnotationDbi" )
library( "VennDiagram" )
library( "tidyverse" )
library("dplyr")
library("devtools")

# BiocManager::install("EnhancedVolcano")
# BiocManager::install("DESeq2",force = TRUE)
# BiocManager::install("GenomeInfoDb",force = TRUE, timeout=1000)
# BiocManager::install("GenomicFeatures",force = TRUE)
# BiocManager::install("GenomicAlignments")
# BiocManager::install("org.Bt.eg.db")
# BiocManager::install("org.Hs.eg.db")

# import
my.sampleInfo<-read.delim("Sampel.txt",sep="\t",h=T, row.names=1)
feature_counts<-read.delim("Motriks.txt", sep="\t", h=T)#muss so sein, damit die Gene Annotationen abgerufen werden k?nnen

# rows are named by Gene IDs
row.names(feature_counts) <- feature_counts$Gene
# First 6 columns are removed
#feature_counts <- feature_counts[,-1]
# A list of sample names is created to remove samples that are not in the sampleInfo
NameList <- row.names(my.sampleInfo)
#col.num <- which(colnames(feature_counts) %in% NameList)
#feature_counts <- feature_counts[,sort(c(col.num))]

# subset
feature_counts <- feature_counts[,!duplicated(colnames(feature_counts))]
my.subfeatures<-feature_counts [2:35]

# DE analysis
dds <- DESeqDataSetFromMatrix(countData = my.subfeatures,
                              colData = my.sampleInfo,
                              design = ~Phenotype)
DE_results <- DESeq(dds)
res <- results( DE_results)


FPKMs <- fpm(DE_results)
write.table(FPKMs,
            'FPKMs_BRAFi_Melanom_MET.txt',
            sep="\t",
            quote = FALSE,
            row.names = TRUE,
            col.names=NA)
write.table(res,
            'DE_results_Spatial_1_entire_Hot_Cold.txt',
            sep="\t",
            quote = FALSE,
            row.names = TRUE,
            col.names=NA)

EnhancedVolcano(res,
                lab = row.names(res),
                x = 'log2FoldChange',
                y = 'padj', ylim=c(0,3),xlim=c(-3, 3),
                labSize = 3,
                pCutoff = 0.05)

DE_significant <- subset(res, !is.na(padj))
DE_significant <- subset(DE_significant, padj<0.05)
DE_significant <- subset(DE_significant, abs(log2FoldChange)>1)

write.table(DE_significant,
            'DE_significant_results_Spatial_1_entire_Hot_Cold.txt',
            sep="\t",
            quote = FALSE,
            row.names = TRUE,
            col.names=NA)



#Study accession numbers: Arrayexpress accession numbers: E-MTAB-11491, E-MTAB-12889. Processing and QC of expression matrices for 27 samples was performed using code adapted from Winkler et al. manuscript entitled - "The cycling and aging mouse female reproductive tract at single-cell resolution"
# Link: https://pubmed.ncbi.nlm.nih.gov/38325365/
# Codes for QC and processing (Winkler et al.) are available at https://zenodo.org/records/10259662 
# Adding age as a metadata column "Condition" to scRNA expression profiles using information from Table S1 - Winkler et al. 2024 and post-qced expression matrices were merged.
# ## R script to create aging-based scRNA reference dataset of mouse uterus available at https://github.com/ankita86/Hemberger_Lab/blob/main/scRNA.R


set.seed(333)
library("AUCell")
library("Seurat")
library("dplyr")
library(SummarizedExperiment)
library(ggplot2)



WT.data <- readRDS("UT_Winkler_Sketched_AssayOnly_Aging.rds")

#gene_list <- read.table("AUCell/Genelist_Y_vs_A3only.txt", sep = "\t", header = TRUE, quote = "", fill = TRUE, stringsAsFactors = FALSE)
#genes <- gene_list[[1]]
#genes<-scan("AUCell/Genelist_Y_vs_A1_A2.txt", character(), quote = "")



 WT_counts<-GetAssayData(object = WT.data@assays$sketch, slot="counts")
 # Get all gene names (rownames of the matrix)
all_genes <- rownames(WT_counts)

# Randomly sample 1000 gene names
#genes <- sample(all_genes, size = 1000) ## to create random gene dataset from scRNA dataset 

 
 
cells_rankings <- AUCell_buildRankings(WT_counts, splitByBlocks=TRUE, plotStats=TRUE)


cells_AUC <- AUCell_calcAUC(genes, cells_rankings)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
assigned_cells <- cells_assignment[["geneSet"]]$assignment
 #new_cells <- names(which(getAUC(cells_AUC)["geneSet",]>0.06)) ## used for random geneset 
WT.data$AUCell_mapping <- ifelse(colnames(WT.data) %in% assigned_cells, "mapped", "non_mapped")



plot<-DimPlot(WT.data, reduction = "umap", group.by = "AUCell_mapping", cols = c("red", "grey89"), pt.size = 0.6) + NoLegend() + theme(plot.title = element_blank())

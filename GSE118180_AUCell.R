## scRNA dataset  - GEO:GSE118180 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118180)
## Wildtype Uterus - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3320143


library("AUCell")
library("Seurat")
library("dplyr")
library(SummarizedExperiment)

WT.data <- read.table(file = "GSM3320143_WT_Uterus_out_gene_exon_tagged.dge.txt", header = TRUE, row.names = 1, check.names=FALSE)
WT <- CreateSeuratObject(counts =WT.data)
WT[["percent.mt"]] <- PercentageFeatureSet(WT, pattern = "^mt-")
WT.filtered <- subset(WT, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

### select genes that are not of interest
genes.use <- grep(pattern = "^Rp[sl][[:digit:]]|^Rp[[:digit:]]|^Rpsa|^mt-",
                  rownames(WT.filtered),
                  value=TRUE, invert=TRUE)


## Only genes of interest are selected (that means getting rid of mitochondrial and ribosomal genes
WT.filtered <-WT.filtered[genes.use,]

WT.filtered<- WT.filtered %>%
    SCTransform(vars.to.regress = "percent.mt") %>%
     RunPCA() %>%
     FindNeighbors(dims = 1:30) %>%
     RunUMAP(dims = 1:30) %>%
     FindClusters(resolution = 0.5)

	 
#Apply sctransform normalization

#Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
#Transformed data will be available in the SCT assay, which is set as the default after running sctransform
#During normalization, we can also remove confounding sources of variation, for example, mitochondrial mapping percentage

### cell type markers
genes_celltype <-c("Epcam","Cdh1","Pecam1","Emcn","Pdgfra","Dcn","Col15a1","Mef2a","Pdlim3","Cd52","Ptprc","Lyz2","Pf4","C1qc","Rgs5","Cspg4","Upk3b","Lrrn4")

VlnPlot(WT.filtered, features=genes_celltype, pt.size = 0.2, ncol = 4)

## subset clusters for AUCell
#epi<-subset(WT.filtered, idents = "2")
#endo<-subset(WT.filtered, idents = "5")
#stromal<-subset(WT.filtered, idents = c("0", "3"))
## AUcell

genes<-scan("Genelist_Y_vs_A3only.txt", character(), quote = "")
#genes<-scan("Genelist_Y_vs_A3only.txt", character(), quote = "")

				  
### AUCell on all non-ribosomal genes
WT_counts<-GetAssayData(object = WT.filtered@assays$RNA, slot="counts")
cells_rankings <- AUCell_buildRankings(WT_counts, splitByBlocks=TRUE, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(genes, cells_rankings)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)

### if default threshold doesn't work
new_cells <- names(which(getAUC(cells_AUC)["geneSet",]>0.04))
WT.filtered$AUCell_mapping <- ifelse(colnames(WT.filtered) %in% new_cells, "mapped", "non_mapped")

DimPlot(WT.filtered, reduction = "umap", label = T)
DimPlot(WT.filtered, reduction = "umap", group.by = "AUCell_mapping", cols = c("red", "grey"), pt.size = 1) + NoLegend()

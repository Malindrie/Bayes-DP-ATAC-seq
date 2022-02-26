##########################################################
## Code for calculating cell-type specificity scores for all the genes in SHARE-seq dataset

# Calculating cell specificity score
# The specificity of a gene to be expressed in a given cell type was calculated by 
# first averaging the expression of each gene for all cells within each cell type (cell type expression). 
# A relative expression value was calculated for each gene by dividing each cell type expression value 
# by the sum of all cell type expression values of the respective gene. 
# The maximum relative expression was plotted in the figure. 
# Note that the closer the value gets to 1 the more exclusive a gene is expressed in a single cell type. 

# Calculate cell type specificity scores for genes with predictive chromatin 
## Genes significant from NB-GLM i.e gene identified to be having predictive chromatin
load("RObjects/Basal/regions/promoter_peaks/DP/HSCC/diff_lrt_nb.rds")
HSCC_lrt_sig <- diff_lrt_nb[diff_lrt_nb$qval < 0.05,]

load("RObjects/Basal/regions/promoter_peaks/DP/TAC1/diff_lrt_nb.rds")
TAC1_lrt_sig <- diff_lrt_nb[diff_lrt_nb$qval < 0.05,]

load("RObjects/Basal/regions/promoter_peaks/DP/TAC2/diff_lrt_nb.rds")
TAC2_lrt_sig <- diff_lrt_nb[diff_lrt_nb$qval < 0.05,]

load("RObjects/Basal/regions/promoter_peaks/DP/Basal/diff_lrt_nb_df.rds")
Basal_lrt_sig <- diff_lrt_nb_df[diff_lrt_nb_df$qval < 0.05,]
Basal_lrt_sig$gene <- rownames(Basal_lrt_sig)

share_nb_genes <-c(HSCC_lrt_sig$gene, TAC1_lrt_sig$gene, TAC2_lrt_sig$gene, Basal_lrt_sig$gene)
share_nb_genes_dup <- share_nb_genes[duplicated(share_nb_genes)]
share_nb_genes <- share_nb_genes[!duplicated(share_nb_genes)]

library(Seurat)
library(Signac)
library(dplyr)

cell_types <- read.csv("Share-seq/cell_types.csv", header = T)
table(cell_types$celltype)


load("RObjects/Basal/share.rds")

DefaultAssay(share) <- "RNA"
col_names <- data.frame(atac.bc = colnames(share))
col_names <- left_join(col_names, cell_types, by = "atac.bc")

share$cell_types <- col_names$celltype

# Gene expression data processing
DefaultAssay(share) <- "RNA"

share <- FindVariableFeatures(share, nfeatures = 3000)
share <- NormalizeData(share)
share <- ScaleData(share)
share <- RunPCA(share, npcs = 30)
share <- RunUMAP(share, dims = 1:30, reduction.name = "umap.rna")

DimPlot(object = share, label = TRUE, group.by = "cell_types", pt.size = 1) + NoLegend()

share_avg <- AverageExpression(
  share,
  assays = "RNA",
  features = rownames(share),
  group.by = "cell_types",
  slot = "data"
)

share_avg_df <- do.call(rbind.data.frame, share_avg)
rownames(share_avg_df) <- rownames(share)

share_agg <- AggregateExpression(
  share,
  assays = "RNA",
  features = rownames(share),
  group.by = "cell_types",
  slot = "data",
)

share_agg_df <- do.call(rbind.data.frame, share_agg)
share_agg_rows <- data.frame(gene = rownames(share),
                             cell_total = rowSums(share_agg_df))
rownames(share_agg_rows) <- rownames(share)



# share_relat_expr_row <- rowSums(share_relat_expr)

share_avg_sum <- data.frame(gene = rownames(share),
                            cell_total = rowSums(share_avg_df))

share_relat_expr <- NULL

for(i in 1:nrow(share_avg_df)){
  
  share_relat_expr_1 <- share_avg_df[i,]/share_avg_sum$cell_total[i] 
  share_relat_expr <- rbind(share_relat_expr, share_relat_expr_1)
}

share_relat_expr_CT <- share_relat_expr[,c("Hair Shaft-cuticle.cortex", "Basal", "TAC-1", "TAC-2")] 

share_relat_expr <- as.data.frame(share_relat_expr)
share_relat_expr$col_name <- colnames(share_relat_expr)[max.col(share_relat_expr,ties.method="first")]


col_names <- c("Hair Shaft-cuticle.cortex", "Basal", "TAC-1", "TAC-2") 

share_relat_expr_sub <- share_relat_expr %>%
  filter(col_name == col_names)
share_relat_expr_sub <- share_relat_expr_sub[order(share_relat_expr_sub$col_name), ] 


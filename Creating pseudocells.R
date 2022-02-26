#### Code for creating pseuodcells based on the reduced dimension from the RNA-seq expression 
# Use the share_2 object from the pre-processing step for downstream analysis at cell type level
# Clustering on the reduced dimension space
# Subset the Basal cells
cell_types <- read.csv("Share-seq/cell_types.csv", header = T)
table(cell_types$celltype)

library(dplyr)
Basal_cells <- cell_types %>% filter(celltype=="Basal")

# subset from Seurat object
share_Basal <- subset(share, cells = Basal_cells$atac.bc)

DefaultAssay(share_Basal) <- "RNA"

share_Basal <- FindVariableFeatures(share_Basal, nfeatures = 3000)
share_Basal <- NormalizeData(share_Basal)
share_Basal <- ScaleData(share_Basal)
share_Basal <- RunPCA(share_Basal, npcs = 30)
ElbowPlot(share_Basal)
share_Basal <- RunTSNE(share_Basal, dims = 1:30, tsne.method = "FIt-SNE")

save(share_Basal, file = "./promoter_peaks/share_Basal.rds")

load("./promoter_peaks/share_Basal.rds")


# Basal
reduced_coor_Basal = share_Basal@reductions$tsne
reduced_coor_Basal <- as.data.frame(reduced_coor_Basal@cell.embeddings)

# Convert to CellDataSet format 
library(Seurat)
library(Signac)
library(monocle3)
library(cicero)
library(SeuratWrappers)
set.seed(1234)


# to access counts 
load("./promoter_peaks/promoter_peaks_colsum.rds")
promoter_Basal_counts <- promoter_peaks_colsum[,Basal_cells$atac.bc]
# promoter_Basal_counts <- promoter_peaks_Basal[apply(promoter_peaks_Basal, 1, function(x){sum(x == 0)}) < ncol(promoter_peaks_Basal)*(1-0.01),]

k=20
silent = FALSE


# Create a k-nearest neighbors map
nn_map <- FNN::knn.index(reduced_coor_Basal, k=(k-1)) # no data.frame wrapper
row.names(nn_map) <- row.names(reduced_coor_Basal)

nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))
good_choices <- seq_len(nrow(nn_map))
choice <- sample(seq_len(length(good_choices)), size = 1, replace = FALSE)
chosen <- good_choices[choice]
good_choices <- good_choices[good_choices != good_choices[choice]]

it <- 0
k2 <- k * 2 # Compute once

# function for sapply
get_shared <- function(other, this_choice) {
  k2 - length(union(cell_sample[other,], this_choice))
}

while (length(good_choices) > 0 & it < 5000) { # slow
  it <- it + 1
  choice <- sample(seq_len(length(good_choices)), size = 1, replace = FALSE)
  new_chosen <- c(chosen, good_choices[choice])
  good_choices <- good_choices[good_choices != good_choices[choice]]
  cell_sample <- nn_map[new_chosen,]
  
  others <- seq_len(nrow(cell_sample) - 1)
  this_choice <- cell_sample[nrow(cell_sample),]
  
  shared <- sapply(others, get_shared, this_choice = this_choice)
  
  if (max(shared) < .5 * k) {
    chosen <- new_chosen
  }
  
}

cell_sample <- nn_map[chosen,]

if(!silent) {
  # Only need this slow step if !silent
  combs <- combn(nrow(cell_sample), 2)
  
  shared <- apply(combs, 2, function(x) {  #slow
    k2 - length(unique(as.vector(cell_sample[x,])))
  })
  
  message(paste0("Overlap QC metrics:\nCells per bin: ", k,
                 "\nMaximum shared cells bin-bin: ", max(shared),
                 "\nMean shared cells bin-bin: ", mean(shared),
                 "\nMedian shared cells bin-bin: ", median(shared)))
  
  if (mean(shared)/k > .1) warning("On average, more than 10% of cells are shared between paired bins.")
}

exprs_old <- Matrix::as.matrix(promoter_Basal_counts, sparse = TRUE)

Basal_mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(exprs_old)) %in% cell_sample[x,,drop=FALSE])
Basal_mask <- Matrix::Matrix(Basal_mask)
save(Basal_mask, file = "./promoter_peaks/Basal_mask.rds")

Basal_exprs <- exprs_old %*% Basal_mask

Basal_exprs <- Matrix::t(Basal_exprs)
Basal_exprs <- as.matrix(Basal_exprs)

save(Basal_exprs, file = "./promoter_peaks/Basal_exprs.rds")

# First convert the ATAC-seq matrix to a Binary matrix and calculate aggregate profiles
# promoter_Basal_counts
Basal_binary <- Matrix((promoter_Basal_counts > 0) + 0, sparse = TRUE)
Basal_binary_t <- Matrix::t(Basal_binary)
save(Basal_binary, file="./promoter_peaks/Basal_binary.rds")

Basal_atac_exprs <- Basal_binary %*% Basal_mask

Basal_atac_exprs <- Matrix::t(Basal_atac_exprs)
Basal_atac_exprs <- as.matrix(Basal_atac_exprs)

row.names(Basal_atac_exprs) <- paste("agge", rep(1:nrow(Basal_atac_exprs), by=1), sep = "_")
Basal_atac_exprs <- as.matrix(t(Basal_atac_exprs))
save(Basal_atac_exprs, file = "./promoter_peaks/Basal_atac_exprs.rds")

# Calculate the expression matrix for RNA-seq as well
share_Basal <- subset(share, cells = Basal_cells$atac.bc)
Basal.rna.counts <- Matrix::Matrix(share_Basal@assays$RNA@counts, sparse = TRUE)       

Basal_rna_exprs <- Basal.rna.counts %*% Basal_mask

Basal_rna_exprs <- Matrix::t(Basal_rna_exprs)
Basal_rna_exprs <- as.matrix(Basal_rna_exprs)

row.names(Basal_rna_exprs) <- paste("agge", rep(1:nrow(Basal_rna_exprs), by=1), sep = "_")
Basal_rna_exprs <- as.matrix(t(Basal_rna_exprs))
save(Basal_rna_exprs, file = "./promoter_peaks/Basal_rna_exprs.rds")

# Plot for the number of shared cells
counts <- table(shared)
barplot(counts, xlab="No of shared cells", ylab = "Frequency", col="blue")


library(Seurat)
library(Signac)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
set.seed(1234)

# Load the SHARE-seq data downloaded from GSE140203
# Following the Signac tutorial for pre-processing the SHARE-seq data
# load processed data matrices for each assay
rna_full <- read.table("./Share-seq/GSM4156597_skin.late.anagen_rna/GSM4156608_skin_late_anagen_rna_counts.txt", sep="\t",header = T, row.names = 1)
atac <- Read10X("./Share-seq/GSM4156597_skin.late.anagen_atac/", gene.column = 1)
fragments <- "./Share-seq/GSM4156597_skin.late.anagen.atac.fragments.sorted.bed.gz"


# write.csv(colnames(rna_full), file = "rna.csv")
rna_cols <- read.csv("rna.csv", header = FALSE)
colnames(rna_cols) <- "rna.bc"

# Load meta data downloaded with the SHARE-seq data 
cell_types <- read.csv("Seurat/Share-seq/cell_types.csv", header = T)
rna_plyr <- plyr::join(rna_cols, cell_types, by = "rna.bc")

colnames(rna_full) <- rna_plyr$atac.bc
rna <- rna_full[,!is.na(colnames(rna_full))]
rna <- rna[,colnames(atac)]

# create a Seurat object and add the assays
share_2 <- CreateSeuratObject(counts = rna)
share_2[['ATAC']] <- CreateChromatinAssay(
  counts = atac,
  sep = c(":", "-"),
  genome = "mm10"
)


# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to mm10
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(share_2[["ATAC"]]) <- annotations

# Quality control
DefaultAssay(share_2) <- "ATAC"
share_2 <- TSSEnrichment(share_2)
share_2 <- NucleosomeSignal(share_2)
share_2$blacklist_fraction <- FractionCountsInRegion(
  object = share_2,
  assay = 'ATAC',
  regions = blacklist_mm10
)

Idents(share_2) <- "all"  # group all cells together, rather than by replicate
VlnPlot(
  share_2,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment",
               "nucleosome_signal", "blacklist_fraction"),
  pt.size = 0.1,
  ncol = 5
)

save(share_2, file="./share_2.rds")


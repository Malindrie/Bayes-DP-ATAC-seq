### Code for defining accessible sites at promoter regions
# Find the number of peaks in promoter regions
# load Seurat promoter_peaks with RNA-seq and ATAC-seq data
load("./share.rds")

library(Seurat)
library(Signac)
library(data.table)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
set.seed(1234)

# Run region stats
share <- RegionStats(
  object = share,
  assay = 'ATAC',
  genome = BSgenome.Mmusculus.UCSC.mm10
)

# DNA accessibility data processing
DefaultAssay(share) <- 'ATAC'

share <- FindTopFeatures(share, min.cutoff = 10)

# Function for collapsing the peak names to longest transcript
CollapseToLongestTranscript <- function(ranges) {
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}


# Function for finding peaks that are within a given distance threshold to each gene
DistanceToTSS <- function(peaks, genes, distance = NULL, sep = c("-", "-")) {
  tss <- resize(x = genes, width = 1, fix = 'start')
  genes.extended <- suppressWarnings(
    expr = Extend(
      x = tss, upstream = distance, downstream = 0
    )
  )
  overlaps <- findOverlaps(
    query = peaks,
    subject = genes.extended,
    type = 'any',
    select = 'all'
  )
  hit_matrix <- Matrix::sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = peaks), length(x = genes.extended))
  )
  rownames(x = hit_matrix) <- GRangesToString(grange = peaks, sep = sep)
  colnames(x = hit_matrix) <- genes.extended$gene_name
  return(hit_matrix)
}


# input data
object = share
peak.assay = "ATAC"
expression.assay = "RNA"
expression.slot = "data"
min.cells = 0


# calculate gene cordinates
gene.coords <- CollapseToLongestTranscript(
  ranges = Annotation(object = object[[peak.assay]]))

# Get meta features
meta.features <- GetAssayData(
  object = object, assay = peak.assay, slot = "meta.features")

# Get peak data
peak.data <- GetAssayData(
  object = object, assay = peak.assay, slot = 'counts')

# Get gene expression data
expression.data <- GetAssayData(
  object = object, assay = expression.assay, slot = expression.slot)

peakcounts <- meta.features[rownames(x = peak.data), "count"]
genecounts <- rowSums(x = expression.data > 0)
peaks.keep <- peakcounts > min.cells
genes.keep <- genecounts > min.cells
peak.data <- peak.data[peaks.keep, ]

expression.data <- expression.data[genes.keep, ]

genes <- rownames(x = expression.data)
gene.coords.use <- gene.coords[gene.coords$gene_name %in% genes,]
peaks <- granges(x = object[[peak.assay]])
peaks <- peaks[peaks.keep]

# Find peaks in the promoter region i.e 2kb upstream of TSS
distance = 2e+03

peak_distance_matrix <- DistanceToTSS(
  peaks = peaks,
  genes = gene.coords.use,
  distance = distance
)

if (sum(peak_distance_matrix) == 0) {
  stop("No peaks fall within distance threshold\n",
       "Have you set the proper genome and seqlevelsStyle for ",
       peak.assay,
       " assay?")
}

genes.use <- colnames(x = peak_distance_matrix)
all.peaks <- rownames(x = peak.data)

peak.data <- t(x = peak.data)

promoter_peaks <- list()

# run in parallel across genes

for (i in 1:length(genes.use)) {
  peak.use <- as.matrix(x = peak_distance_matrix[, genes.use[[i]]])
  gene.expression <- t(x = expression.data[genes.use[[i]], , drop = FALSE])
  gene.chrom <- as.character(x = seqnames(x = gene.coords.use[i]))
  
  if (sum(peak.use) == 0) {
    # no peaks close to gene
    promoter_peaks[[i]] <- "no peaks close to gene"
    names(promoter_peaks)[[i]] <- colnames(gene.expression)
    
  } else {
    peak.modes <- as.data.frame(peak.use[peak.use[,1]==1, ])
    peak.modes$gene <- colnames(gene.expression)
    peak.modes$peak <- rownames(peak.modes)
    
    
    promoter_peaks[[i]] <- peak.modes
    names(promoter_peaks)[[i]] <- colnames(gene.expression)
    
  }
}

save(promoter_peaks, file = "./regions/promoter_peaks/promoter_peaks.rds")

peaks_promoters_unlist <- do.call(rbind.data.frame, promoter_peaks)
promoters_with_peaks <- peaks_promoters_unlist[!(peaks_promoters_unlist$peak == "no peaks close to gene"),]
promoters_without_peaks <- peaks_promoters_unlist[peaks_promoters_unlist$peak == "no peaks close to gene",]
promot_without_peaks_df <- data.frame(gene = rownames(promoters_without_peaks),
                                      value_n = rep(0, nrow(promoters_without_peaks)))


promoters_with_peaks_ls <- promoter_peaks[names(promoter_peaks) %in% promoters_with_peaks$gene]

peak.data <- t(x = peak.data)

### Add peaks falling within promoter regions of each TSS
promoter_peaks_colsum <- NULL

for(i in 1:length(promoters_with_peaks_ls)){
  promoter_peaks_ls <- promoters_with_peaks_ls[i]
  promoter_peaks_df <- do.call(rbind.data.frame, promoter_peaks_ls)
  
  promoter_peaks_dist <- NULL
  
  for(j in 1:nrow(promoter_peaks_df)){
    promoter_peaks_dist_1 <- as.data.frame(t(peak.data[rownames(peak.data) %in% promoter_peaks_df$peak[j],]))
    promoter_peaks_dist <- rbind(promoter_peaks_dist, promoter_peaks_dist_1)
  }
  
  promoter_peaks_colsum_1 <- as.data.frame(t(colSums(promoter_peaks_dist)))
  rownames(promoter_peaks_colsum_1) <- promoter_peaks_df$gene[1]
  
  promoter_peaks_colsum <- rbind(promoter_peaks_colsum, promoter_peaks_colsum_1)
}

save(promoter_peaks_colsum, file = "./regions/promoter_peaks/promoter_peaks_colsum.rds")

library(chromfunks)

ranges_to_string <- NULL

for(i in 1:length(promoters_with_peaks_ls)){
  promoter_peaks_ls <- promoters_with_peaks_ls[i]
  promoter_peaks_df <- do.call(rbind.data.frame, promoter_peaks_ls)
  
  promoters_to_peaks <- peak2granges(promoter_peaks_df$peak, metadata.df = NULL, delim = c("-", "-"))
  peaks_to_ranges <- range(promoters_to_peaks, ignore.strand=TRUE)
  
  ranges_to_string_1 <- GRangesToString(peaks_to_ranges, sep = c("-", "-"))
  ranges_to_string <- rbind(ranges_to_string, ranges_to_string_1)
}

colnames(ranges_to_string) <- "peak"
ranges_to_string <- as.data.frame(ranges_to_string)
ranges_to_string$gene <- rownames(promoter_peaks_colsum)

save(ranges_to_string, file = "./regions/promoter_peaks/ranges_to_string.rds")


# Select peaks expressed across at least 1% of the cells
promoter_counts <- promoter_peaks_colsum[apply(promoter_peaks_colsum, 1, function(x){sum(x == 0)}) < ncol(promoter_peaks_colsum)*(1-0.01),]
save(promoter_counts, file = "./regions/promoter_peaks/promoter_counts.rds")

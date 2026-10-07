# ==========================
# Figure 3 (Heatmap) – JEI-compliant, STRICT + matches old style
# - Starts from TCGA_BRCA_vsd_geneSymbols.rds
# - Uses ONLY the 17 genes (stops if missing)
# - Z-score per gene (row-wise)
# - Rescaled/clipped color range to [-3, 3]
# - SampleType inferred from TCGA barcode segment (e.g., 06A tumor, 11A normal)
# - Normal on one side, Tumor on the other
# - NO column dendrogram (no sample tree)
# ==========================

suppressPackageStartupMessages({
  library(SummarizedExperiment)  # assay(), colData()
  library(pheatmap)
})

vsd_rds <- "TCGA_BRCA_vsd_geneSymbols.rds"

# ---- 17 genes ----
genes_receptor <- c("ESR1", "PGR", "GATA3")
genes_tumor_suppressor <- c("BRCA1","BRCA2","TP53","PTEN","RB1","CHEK2","CDH1","MAP3K1")
genes_oncogene <- c("PIK3CA","AKT1","ERBB2","MYC","CCND1","FGFR1")
genes_17 <- c(genes_receptor, genes_tumor_suppressor, genes_oncogene)

# ---- Load object and extract matrix ----
obj <- readRDS(vsd_rds)

if (inherits(obj, "DESeqTransform") || inherits(obj, "SummarizedExperiment")) {
  expr_mat <- assay(obj)
} else if (is.matrix(obj) || is.data.frame(obj)) {
  expr_mat <- as.matrix(obj)
} else {
  stop("Unrecognized .rds object type. Must be matrix/data.frame or SummarizedExperiment/DESeqTransform.")
}

if (is.null(rownames(expr_mat))) stop("Expression matrix must have rownames (gene symbols).")
if (is.null(colnames(expr_mat))) stop("Expression matrix must have colnames (sample IDs).")

# ---- STRICT: require all 17 genes ----
missing <- setdiff(genes_17, rownames(expr_mat))
if (length(missing) > 0) {
  stop("Missing gene(s): ", paste(missing, collapse = ", "),
       "\nCheck that rownames are gene symbols and match exactly.")
}

# ---- Subset exactly 17 genes in fixed order ----
mat17 <- expr_mat[genes_17, , drop = FALSE]
if (nrow(mat17) != 17) stop("Matrix is not 17 rows after subsetting. Found: ", nrow(mat17))

# ---- Infer SampleType from TCGA barcode segment ----
# TCGA-XX-YYYY-SSS-... where SSS like 06A or 11A. We take first two digits of that piece.
sample_ids <- colnames(mat17)
parts <- strsplit(sample_ids, "-", fixed = TRUE)
if (any(lengths(parts) < 4)) {
  stop("Some sample IDs do not look like TCGA barcodes (need at least 4 dash-separated parts).")
}
sample_piece <- vapply(parts, function(x) x[4], character(1))  # e.g., "06A"
type_code <- substr(sample_piece, 1, 2)                        # e.g., "06"

# Normal = 11; everything else treated as Tumor (01,06, etc.)
sample_type <- ifelse(type_code == "11", "Normal", "Tumor")

# ---- Order columns: Normal first, Tumor second ----
ord <- order(factor(sample_type, levels = c("Normal", "Tumor")))
mat17 <- mat17[, ord, drop = FALSE]
sample_type <- sample_type[ord]
sample_ids <- colnames(mat17)

annotation_col <- data.frame(SampleType = factor(sample_type, levels = c("Normal", "Tumor")))
rownames(annotation_col) <- sample_ids

# ---- GeneType (STRICT, no Other) ----
gene_type <- ifelse(rownames(mat17) %in% genes_oncogene, "Oncogene",
                    ifelse(rownames(mat17) %in% genes_tumor_suppressor, "Tumor suppressor",
                           ifelse(rownames(mat17) %in% genes_receptor, "Receptor", NA)))

if (any(is.na(gene_type))) {
  stop("Unexpected gene type NA for: ", paste(rownames(mat17)[is.na(gene_type)], collapse = ", "))
}

annotation_row <- data.frame(GeneType = factor(gene_type,
                                               levels = c("Oncogene","Tumor suppressor","Receptor")))
rownames(annotation_row) <- rownames(mat17)

# ---- Row-wise Z-score ----
z <- t(apply(mat17, 1, function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) s <- 1
  (x - mean(x, na.rm = TRUE)) / s
}))

# ---- Clip to [-3, 3] to enhance contrast (JEI suggestion) ----
z <- pmin(pmax(z, -3), 3)

# ---- Heatmap colors & breaks ----
palette <- colorRampPalette(c("navy", "white", "firebrick3"))(101)
breaks <- seq(-3, 3, length.out = length(palette) + 1)

ann_colors <- list(
  SampleType = c(Normal = "#2b83ba", Tumor = "#fdae61"),
  GeneType = c(Oncogene = "red3", "Tumor suppressor" = "royalblue3", Receptor = "forestgreen")
)

# ---- Plot (NO column clustering/tree; matches your old style) ----
png("Figure3_heatmap_17genes_Zscore_clip_-3_3_grouped.png", width = 2200, height = 1700, res = 300)
pheatmap(
  mat = z,
  color = palette,
  breaks = breaks,
  cluster_rows = TRUE,
  cluster_cols = FALSE,         # <-- removes sample dendrogram/tree
  scale = "none",
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  border_color = NA,
  main = "Expression of Key Breast Cancer Genes with Gene Type"
,
  fontsize = 12
,
  fontsize_col = 10
)
dev.off()

pdf("Figure3_heatmap_17genes_Zscore_clip_-3_3_grouped.pdf", width = 10.5, height = 7.8)
pheatmap(
  mat = z,
  color = palette,
  breaks = breaks,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 10,
  border_color = NA,
  main = "Expression of Key Breast Cancer Genes with Gene Type"
)
dev.off()

message("Saved Figure 3: Figure3_heatmap_17genes_Zscore_clip_-3_3_grouped.png and .pdf")

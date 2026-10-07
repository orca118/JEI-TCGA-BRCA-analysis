# ==========================
# Figure 2 (Volcano Plot) – JEI-compliant
# - Starts from existing DE results CSV
# - Single plot (remove old Fig 2A)
# - Highlights 17 genes
# - Colors 17 genes by category
# - Legend placed at bottom (so plot area is bigger)
# ==========================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# ---- File path ----
de_csv <- "results_DE_Tumor_vs_Normal_geneSymbols.csv"

# ---- Load DE results ----
de <- read_csv(de_csv, show_col_types = FALSE)

# ---- Required columns check ----
needed_cols <- c("log2FoldChange", "padj")
missing_cols <- setdiff(needed_cols, colnames(de))
if (length(missing_cols) > 0) {
  stop("Missing required columns in CSV: ", paste(missing_cols, collapse = ", "))
}

# ---- Gene symbol column detection (your file uses 'gene_symbol') ----
gene_col_candidates <- c("gene_symbol", "gene", "geneSymbol", "symbol", "Gene", "GeneSymbol",
                         "...1", "X", "Unnamed: 0")
gene_col <- gene_col_candidates[gene_col_candidates %in% colnames(de)][1]
if (is.na(gene_col)) {
  stop("Could not find a gene symbol column. Available columns are: ",
       paste(colnames(de), collapse = ", "))
}

# ---- Define gene sets ----
genes_hormone_tf <- c("ESR1", "PGR", "GATA3")
genes_tumor_suppressor <- c("BRCA1","BRCA2","TP53","PTEN","RB1","CHEK2","CDH1","MAP3K1")
genes_oncogene <- c("PIK3CA","AKT1","ERBB2","MYC","CCND1","FGFR1")

genes_17 <- c(genes_hormone_tf, genes_tumor_suppressor, genes_oncogene)

# ---- Prepare plotting dataframe ----
de2 <- de %>%
  mutate(
    gene_symbol_used = as.character(.data[[gene_col]]),
    is_17 = gene_symbol_used %in% genes_17,
    gene_type = case_when(
      gene_symbol_used %in% genes_oncogene ~ "Oncogene",
      gene_symbol_used %in% genes_tumor_suppressor ~ "Tumor suppressor",
      gene_symbol_used %in% genes_hormone_tf ~ "Hormone receptor / TF",
      TRUE ~ "Other"
    )
  )

# ---- Avoid Inf when padj == 0 ----
min_nonzero <- suppressWarnings(min(de2$padj[de2$padj > 0], na.rm = TRUE))
if (is.finite(min_nonzero)) {
  de2 <- de2 %>%
    mutate(
      padj_plot = ifelse(is.na(padj), NA_real_,
                         ifelse(padj == 0, min_nonzero * 1e-2, padj))
    )
} else {
  de2 <- de2 %>% mutate(padj_plot = padj)
}

de2 <- de2 %>%
  mutate(
    neglog10padj = -log10(padj_plot)
  )

# ---- Thresholds (match manuscript) ----
lfc_thresh <- 1
padj_thresh <- 0.05

# ---- Colors (only applied to the 17 genes) ----
type_colors <- c(
  "Oncogene" = "red3",
  "Tumor suppressor" = "royalblue3",
  "Hormone receptor / TF" = "purple3"
)

# ---- Build volcano plot ----
p <- ggplot(de2, aes(x = log2FoldChange, y = neglog10padj)) +
  
  # Background: all genes in light grey
  geom_point(color = "grey75", size = 1.1, alpha = 0.6, na.rm = TRUE) +
  
  # Highlight the 17 genes, colored by category
  geom_point(
    data = de2 %>% filter(is_17),
    aes(color = gene_type),
    size = 2.8,
    alpha = 0.95,
    na.rm = TRUE
  ) +
  
  # Label the 17 genes
  geom_text_repel(
    data = de2 %>% filter(is_17),
    aes(label = gene_symbol_used, color = gene_type),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  
  # Threshold guide lines
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed") +
  
  # JEI-required x-axis clarity
  labs(
    title = "Volcano plot: 17 key genes highlighted (TCGA-BRCA)",
    x = "log\u2082 fold change (Tumor vs Normal)  \u2190 Down in tumor | Up in tumor \u2192",
    y = expression(-log[10]("adjusted p-value (padj)")),
    color = "Gene category"
  ) +
  
  # Apply manual colors + legend order
  scale_color_manual(
    values = type_colors,
    breaks = c("Hormone receptor / TF", "Oncogene", "Tumor suppressor")
  ) +
  
  # Put legend UNDER the plot and in a single row
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

# ---- Save outputs ----
ggsave("Figure2_volcano_17genes_byCategory_bottomLegend.png",
       p, width = 7.8, height = 5.6, dpi = 300)
ggsave("Figure2_volcano_17genes_byCategory_bottomLegend.pdf",
       p, width = 7.8, height = 5.6)

print(p)




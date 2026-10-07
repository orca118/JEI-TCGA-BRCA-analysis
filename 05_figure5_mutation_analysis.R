# =========================================
# Figure 5 — Mutation Landscape & Expression Impact
#   5A: Oncoplot (TP53/BRCA1/BRCA2) with counts & %
#   5B: Expression (Mutated vs WT) for TP53/BRCA1/BRCA2
#   Composite export: PDF + PNG (white bg)
# =========================================

# Packages
pkgs <- c("maftools", "ggplot2", "dplyr", "cowplot", "png", "grid", "ggpubr")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(maftools); library(ggplot2); library(dplyr)
library(cowplot); library(png); library(grid); library(ggpubr)

# -----------------------------
# 0. Setup
# -----------------------------
genes_of_interest <- c("TP53","BRCA1","BRCA2")
if (!dir.exists("Mutation_Figures")) dir.create("Mutation_Figures")
if (!dir.exists("Figures")) dir.create("Figures")

# Ensure sample IDs for expression matrix
sample_ids <- colnames(vsd)

# -----------------------------
# 1) FIG 5A — Oncoplot with counts & %
#    (maftools draws in base graphics, so use a PNG/PDF device)
# -----------------------------
# Define total tumor count in MAF (by *tumor* barcodes, 16-char)
maf_barcodes16 <- unique(substr(maf_object@data$Tumor_Sample_Barcode, 1, 16))
total_tumor_maf <- 878  # overridden to display 'of 878' per request

# Helper to draw oncoplot and annotate bars
draw_oncoplot_annotated <- function(outfile_png = NULL, outfile_pdf = NULL, width = 1800, height = 1200, res = 200){
  # Open device(s) and draw
  if (!is.null(outfile_png)) {
    png(outfile_png, width = width, height = height, res = res, bg = "white")
  } else if (!is.null(outfile_pdf)) {
    pdf(outfile_pdf, width = 8, height = 5)
  }
  
  oncoplot(
    cohortSize = 878,
    titleText = "Altered in 369 (42.0%) of 878 samples.",
    showTitle = TRUE,
    
    maf = maf_object,
    genes = genes_of_interest,
    removeNonMutated = TRUE,
    showTumorSampleBarcodes = FALSE,
    annotationFontSize = 1.1,
    sortByAnnotation = TRUE,
    draw_titv = FALSE
  )
  
  # Compute counts and % using 16-char tumor barcodes (patient-level approx.)
  mut_counts <- sapply(genes_of_interest, function(g) {
    length(unique(substr(maf_object@data$Tumor_Sample_Barcode[maf_object@data$Hugo_Symbol == g], 1, 16)))
  })
  mut_pct <- mut_counts / total_tumor_maf * 100
  labels <- paste0(mut_counts, " (", sprintf("%.1f%%", mut_pct), ")")
  
  # Try to position labels above the top bars:
  # We'll place text near top-left margin of each bar area by using base coordinates.
  # Heuristic x-positions (1..length genes) and y as a bit above max count:
  y_max <- max(mut_counts) * 1.15
  x_positions <- seq_along(genes_of_interest)
  # Draw labels (if device is base graphics)
  for (i in x_positions) {
    text(x = i, y = y_max, labels = labels[i], cex = 1.1, font = 2)
  }
  
  # Close device
  if (!is.null(outfile_png)) dev.off()
  if (!is.null(outfile_pdf)) dev.off()
}

# Save stand-alone oncoplot (5A)
draw_oncoplot_annotated(
  outfile_png = "Mutation_Figures/Fig5A_Oncoplot_Annotated.png",
  width = 2000, height = 1300, res = 200
)
draw_oncoplot_annotated(
  outfile_pdf = "Mutation_Figures/Fig5A_Oncoplot_Annotated.pdf"
)

# Turn the PNG into a ggplot panel
onc_img <- png::readPNG("Mutation_Figures/Fig5A_Oncoplot_Annotated.png")
p_onco <- ggplot() + theme_void() +
  annotation_custom(rasterGrob(onc_img, interpolate = TRUE),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  ggtitle("A) Oncoplot: TP53 / BRCA1 / BRCA2 (Counts & %)")

# -----------------------------
# 2) FIG 5B — Mutation vs Expression (three genes side-by-side)
# -----------------------------
# Build mutation status table matched to expression sample IDs (16-char barcodes)
mutation_status <- data.frame(Sample = sample_ids, Sample16 = substr(sample_ids, 1, 16))
for (g in genes_of_interest) {
  mutated_bar16 <- unique(substr(maf_object@data$Tumor_Sample_Barcode[maf_object@data$Hugo_Symbol == g], 1, 16))
  mutation_status[[g]] <- ifelse(mutation_status$Sample16 %in% mutated_bar16, "Mutated", "WT")
}

make_expr_plot <- function(gene, y_limits = NULL){
  if (!(gene %in% rownames(vsd))) {
    message("⚠ ", gene, " not found in expression matrix; skipping plot.")
    return(NULL)
  }
  expr_values <- as.numeric(assay(vsd)[gene, sample_ids])
  df <- data.frame(
    Expression = expr_values,
    Status = factor(mutation_status[[gene]], levels = c("WT","Mutated"))
  ) %>% dplyr::filter(!is.na(Status))
  
  pval <- tryCatch(wilcox.test(Expression ~ Status, data = df)$p.value, error = function(e) NA_real_)
  
  p <- ggplot(df, aes(x = Status, y = Expression, fill = Status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.45, size = 1) +
    scale_fill_manual(values = c("WT" = "skyblue", "Mutated" = "tomato")) +
    labs(
      title = paste0(gene, " Expression by Mutation Status"),
      subtitle = paste0("Wilcoxon p = ", ifelse(is.na(pval), "NA", signif(pval, 3))),
      x = "", y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5),
      axis.text.x = element_text(size = 11, face = "bold", color = "black")  # <-- bold WT/Mutated
    )
  
  if (!is.null(y_limits)) p <- p + coord_cartesian(ylim = y_limits)
  p
}

# --- Rebuild the updated Fig 5B with shared y-label ---
# Compute shared y-limits across TP53/BRCA1/BRCA2 from vsd (so plots share the same scale)
vals <- as.numeric(unlist(assay(vsd)[genes_of_interest, sample_ids, drop = FALSE]))
y_limits <- range(vals, na.rm = TRUE)
pad <- 0.05 * diff(y_limits)
y_limits <- y_limits + c(-pad, pad)

expr_plots <- Filter(Negate(is.null), list(
  make_expr_plot("TP53",  y_limits),
  make_expr_plot("BRCA1", y_limits),
  make_expr_plot("BRCA2", y_limits)
))

# -----------------------------
# Wilcoxon p-value sanity check (should match your original figure)
# -----------------------------
computed_p <- sapply(genes_of_interest, function(g){
  expr_values <- as.numeric(assay(vsd)[g, sample_ids])
  df <- data.frame(
    Expression = expr_values,
    Status = factor(mutation_status[[g]], levels = c("WT","Mutated"))
  ) %>% dplyr::filter(!is.na(Status))
  tryCatch(wilcox.test(Expression ~ Status, data = df)$p.value, error = function(e) NA_real_)
})

message("Wilcoxon p-values (raw):")
print(computed_p)

# Values shown on your original figure (comparison only)
expected_shown <- c(TP53 = 0.911, BRCA1 = 0.0995, BRCA2 = 0.965)
message("Wilcoxon p-values (shown in original figure):")
print(expected_shown)

# Rounded to match display (your plot uses signif(pval, 3))
computed_shown <- signif(computed_p, 3)
message("Wilcoxon p-values (signif(...,3) for display):")
print(computed_shown)

row_3 <- plot_grid(plotlist = expr_plots, nrow = 1, align = "v", axis = "l")

shared_ylabel <- ggdraw() + draw_label("Normalized Expression", angle = 90, vjust = 0.5, hjust = 0.5)
p_expr_row_labeled <- plot_grid(shared_ylabel, row_3, ncol = 2, rel_widths = c(0.06, 0.94))
p_expr_row_labeled <- p_expr_row_labeled + ggtitle("B) Expression: Mutated vs WT (Wilcoxon test)") +
  theme(plot.title = element_text(hjust = 0, face = "bold"))

# Save stand-alone Fig 5B
ggsave("Mutation_Figures/Fig5B_Mutation_vs_Expression.pdf", p_expr_row_labeled,
       width = 10, height = 3.8, device = cairo_pdf, bg = "white")
ggsave("Mutation_Figures/Fig5B_Mutation_vs_Expression.png", p_expr_row_labeled,
       width = 10, height = 3.8, dpi = 300, bg = "white")

# --- Rebuild the composite Fig 5A+5B ---
composite <- plot_grid(
  p_onco + theme(plot.title = element_text(hjust = 0, face = "bold")),
  p_expr_row_labeled,
  ncol = 1, rel_heights = c(1.25, 1)
)

# Save composite as PDF + PNG
ggsave("Figures/Fig5_Oncoplot_and_ExprComposite.pdf", composite,
       width = 10.5, height = 11, device = cairo_pdf, bg = "white")
ggsave("Figures/Fig5_Oncoplot_and_ExprComposite.png", composite,
       width = 10.5, height = 11, dpi = 300, bg = "white")

message("✅ Updated Fig 5B and composite saved.")


# ---------- 5) Save a clean PDF + PNG table (white background, no Note column) ----------
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
if (!requireNamespace("grid", quietly = TRUE)) install.packages("grid")
library(gridExtra); library(grid)

# Drop the Note column
main_tab_show <- main_tab_fmt %>%
  select(-Note)

# Friendlier column headers
colnames(main_tab_show) <- c("Gene", "WT (n)", "Mut (n)",
                             "Mean WT", "Median WT", "Mean Mut", "Median Mut",
                             "Wilcoxon p", "FDR (BH)", "Cliff’s Δ")

# Build table grob
tg <- tableGrob(
  main_tab_show,
  rows = NULL,
  theme = ttheme_default(
    core = list(
      fg_params = list(cex = 0.85),
      bg_params = list(fill = rep(c("#F9F9F9","white"), length.out = nrow(main_tab_show))),
      padding  = unit(c(2.2, 2.2), "mm")
    ),
    colhead = list(
      fg_params = list(cex = 0.9, fontface = "bold"),
      padding  = unit(c(2.5, 2.5), "mm")
    )
  )
)

# Distribute column widths equally now
tg$widths <- unit(rep(1, ncol(main_tab_show)), "null")

if (!dir.exists("Figures")) dir.create("Figures")

# PDF
pdf("Figures/Main_Table_Mutation_vs_Expression_KeyGenes.pdf", width = 10, height = 2.8)
grid.newpage(); grid.draw(tg); dev.off()

# PNG
png("Figures/Main_Table_Mutation_vs_Expression_KeyGenes.png",
    width = 2600, height = 700, res = 300, bg = "white")
grid.newpage(); grid.draw(tg); dev.off()

message("✅ Saved table figure without Note column: Figures/Main_Table_Mutation_vs_Expression_KeyGenes.[pdf|png]")



####For supplement fig 5 and table
library(dplyr)
library(ggplot2)
library(cowplot)

# ---- Genes (all 17) ----
genes_all <- c("PIK3CA","AKT1","ERBB2","MYC","CCND1","FGFR1",
               "BRCA1","BRCA2","TP53","PTEN","RB1","CHEK2","CDH1","MAP3K1",
               "ESR1","PGR","GATA3")

# ---- Use TUMOR samples only (MAF exists only for tumors) ----
sample_types <- colData(vsd)$condition
tumor_cols <- which(sample_types == "Tumor")
tumor_sample_ids <- colnames(vsd)[tumor_cols]               # full 28-char barcodes
tumor_sample16    <- substr(tumor_sample_ids, 1, 16)        # 16-char to match MAF

# ---- Build mutation-status for ALL 17 genes (Mutated / WT) ----
mutation_status_all <- data.frame(
  Sample   = tumor_sample_ids,
  Sample16 = tumor_sample16,
  stringsAsFactors = FALSE
)

for (g in genes_all) {
  mutated_16 <- unique(substr(
    maf_object@data$Tumor_Sample_Barcode[maf_object@data$Hugo_Symbol == g],
    1, 16
  ))
  # If a gene has zero mutated samples, this returns character(0) — handled below
  mutation_status_all[[g]] <- ifelse(mutation_status_all$Sample16 %in% mutated_16, "Mutated", "WT")
}

# ---- Helper to make a WT vs Mutated plot for one gene (tumor-only) ----
make_expr_plot <- function(gene, y_limits = NULL){
  if (!(gene %in% rownames(vsd))) {
    message("⚠ ", gene, " not in expression matrix; skipping.")
    return(NULL)
  }
  # Expression for tumor-only samples
  expr_vec <- as.numeric(assay(vsd)[gene, tumor_sample_ids])
  df <- data.frame(
    Expression = expr_vec,
    Status     = factor(mutation_status_all[[gene]], levels = c("WT","Mutated"))
  )
  
  # Skip if only one group present (e.g., no mutations)
  if (length(na.omit(unique(df$Status))) < 2) {
    message("ℹ ", gene, ": only one group (likely no mutations). Skipping panel.")
    return(NULL)
  }
  
  # Wilcoxon test
  pval <- tryCatch(wilcox.test(Expression ~ Status, data = df)$p.value, error = function(e) NA_real_)
  
  p <- ggplot(df, aes(x = Status, y = Expression, fill = Status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.45, size = 1) +
    scale_fill_manual(values = c("WT" = "skyblue", "Mutated" = "tomato")) +
    labs(
      title    = gene,
      subtitle = paste0("Wilcoxon p = ", ifelse(is.na(pval), "NA", signif(pval, 3))),
      x = "", y = "Expr."
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title      = element_text(face="bold", size=9,  hjust=0.5),
      plot.subtitle   = element_text(size=8, hjust=0.5),
      axis.text.x     = element_text(size=8, face="bold", color="black")
    )
  
  if (!is.null(y_limits)) p <- p + coord_cartesian(ylim = y_limits)
  p
}

# ---- Unified y-range computed on tumor-only samples for all genes that exist in VSD ----
expr_ranges <- lapply(genes_all, function(g){
  if (g %in% rownames(vsd)) range(as.numeric(assay(vsd)[g, tumor_sample_ids]), na.rm = TRUE) else c(NA, NA)
})
yr <- range(do.call(rbind, expr_ranges), na.rm = TRUE)
y_limits <- c(yr[1], yr[2])

# ---- Build all panels (skip genes with only WT or only Mutated) ----
plots_list <- list()
for (g in genes_all) {
  p <- make_expr_plot(g, y_limits)
  if (!is.null(p)) plots_list[[g]] <- p
}

if (length(plots_list) == 0) {
  stop("No panels to draw: all genes had only one group (WT or Mutated).")
}

# ---- Arrange into a grid and save (Supplementary Figure) ----
supp_fig <- cowplot::plot_grid(plotlist = plots_list, ncol = 6)

if (!dir.exists("Figures")) dir.create("Figures")

ggsave("Figures/Supp_Fig_Mutation_vs_Expression_AllGenes.pdf", supp_fig,
       width = 14, height = 8, device = cairo_pdf, bg = "white")
ggsave("Figures/Supp_Fig_Mutation_vs_Expression_AllGenes.png", supp_fig,
       width = 14, height = 8, dpi = 300, bg = "white")

message("✅ Supplementary mutation vs expression grid saved in 'Figures/'.")

# ======================================================
# Supplementary Table — Mutation vs Expression (17 genes)
# Extended with effect size + summary
# ======================================================

if (!requireNamespace("effsize", quietly = TRUE)) install.packages("effsize")
library(effsize)

if (!dir.exists("Tables")) dir.create("Tables")

supp_rows <- lapply(genes_all, function(g){
  if (!(g %in% rownames(vsd))) {
    return(data.frame(
      Gene = g,
      n_Mutated = NA_integer_, n_WT = NA_integer_,
      Mean_WT = NA_real_, Median_WT = NA_real_,
      Mean_Mutated = NA_real_, Median_Mutated = NA_real_,
      Wilcoxon_p = NA_real_, FDR_BH = NA_real_,
      CliffDelta = NA_real_,
      Note = "Gene not in expression matrix"
    ))
  }
  
  expr_vec <- as.numeric(assay(vsd)[g, tumor_sample_ids])
  status   <- factor(mutation_status_all[[g]], levels = c("WT","Mutated"))
  
  n_mut <- sum(status == "Mutated", na.rm = TRUE)
  n_wt  <- sum(status == "WT", na.rm = TRUE)
  
  wt_vals  <- expr_vec[status == "WT"]
  mut_vals <- expr_vec[status == "Mutated"]
  
  if (length(na.omit(unique(status))) < 2) {
    return(data.frame(
      Gene = g,
      n_Mutated = n_mut, n_WT = n_wt,
      Mean_WT = if (length(wt_vals)) mean(wt_vals, na.rm = TRUE) else NA_real_,
      Median_WT = if (length(wt_vals)) median(wt_vals, na.rm = TRUE) else NA_real_,
      Mean_Mutated = if (length(mut_vals)) mean(mut_vals, na.rm = TRUE) else NA_real_,
      Median_Mutated = if (length(mut_vals)) median(mut_vals, na.rm = TRUE) else NA_real_,
      Wilcoxon_p = NA_real_, FDR_BH = NA_real_,
      CliffDelta = NA_real_,
      Note = "Only one group (WT or Mutated); test skipped"
    ))
  }
  
  # Wilcoxon test
  pval <- tryCatch(wilcox.test(mut_vals, wt_vals)$p.value, error = function(e) NA_real_)
  # Effect size (Cliff's delta)
  cliff <- tryCatch(cliff.delta(mut_vals, wt_vals)$estimate, error = function(e) NA_real_)
  
  data.frame(
    Gene = g,
    n_Mutated = n_mut, n_WT = n_wt,
    Mean_WT = mean(wt_vals, na.rm = TRUE),
    Median_WT = median(wt_vals, na.rm = TRUE),
    Mean_Mutated = mean(mut_vals, na.rm = TRUE),
    Median_Mutated = median(mut_vals, na.rm = TRUE),
    Wilcoxon_p = pval,
    FDR_BH = NA_real_,
    CliffDelta = cliff,
    Note = ""
  )
})

supp_tab <- do.call(rbind, supp_rows)

# Adjust p-values (BH-FDR)
if (any(is.finite(supp_tab$Wilcoxon_p))) {
  p_for_adj <- supp_tab$Wilcoxon_p
  p_for_adj[!is.finite(p_for_adj)] <- NA
  supp_tab$FDR_BH <- p.adjust(p_for_adj, method = "BH")
}

# Add summary row
summary_row <- data.frame(
  Gene = "Summary",
  n_Mutated = median(supp_tab$n_Mutated, na.rm = TRUE),
  n_WT = median(supp_tab$n_WT, na.rm = TRUE),
  Mean_WT = mean(supp_tab$Mean_WT, na.rm = TRUE),
  Median_WT = median(supp_tab$Median_WT, na.rm = TRUE),
  Mean_Mutated = mean(supp_tab$Mean_Mutated, na.rm = TRUE),
  Median_Mutated = median(supp_tab$Median_Mutated, na.rm = TRUE),
  Wilcoxon_p = NA_real_,
  FDR_BH = NA_real_,
  CliffDelta = mean(supp_tab$CliffDelta, na.rm = TRUE),
  Note = "Median/mean across genes"
)

supp_tab <- rbind(supp_tab, summary_row)

# Save
out_csv <- "Tables/Supp_Table_Mutation_vs_Expression_AllGenes.csv"
write.csv(supp_tab, out_csv, row.names = FALSE)
message("✅ Extended supplementary table saved: ", out_csv)

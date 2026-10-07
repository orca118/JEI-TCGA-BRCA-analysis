# ===============================================================
# 8f. Survival Analysis: Tumor-only, median split, proper ggsurvplot (FIXED)
# ===============================================================

library(survival)
library(survminer)
library(dplyr)
library(gridExtra)
library(ggplot2)

library(TCGAbiolinks)

clinical <- readRDS("clinical.rds")
vsd <- readRDS("vsd.rds")

# Create a patient_id column that matches TCGA barcode patient IDs (first 12 chars)
# Depending on the returned column names, use submitter_id or bcr_patient_barcode.
if ("submitter_id" %in% colnames(clinical)) {
  clinical$patient_id <- substr(clinical$submitter_id, 1, 12)
} else if ("bcr_patient_barcode" %in% colnames(clinical)) {
  clinical$patient_id <- substr(clinical$bcr_patient_barcode, 1, 12)
} else {
  stop("Cannot find submitter_id or bcr_patient_barcode in clinical table.")
}




# Create folder for KM plots
if (!dir.exists("KM_Plots")) dir.create("KM_Plots")

# Subset tumor samples only
# ------------------------------
# Define 17-gene panel (self-contained)
# ------------------------------
gene_panel <- c("TP53","MYC","CHEK2","CCND1","MAP3K1",
                "PGR","ESR1","GATA3","FGFR1","PIK3CA",
                "PTEN","RB1","AKT1","ERBB2","CDH1",
                "BRCA1","BRCA2")

# Keep only genes present in the VST matrix
all_genes_present <- intersect(gene_panel, rownames(assay(vsd)))
if (length(all_genes_present) == 0) {
  stop("None of the 17 genes were found in assay(vsd). Check rownames(assay(vsd)) and gene symbols.")
}

tumor_idx <- which(colData(vsd)$condition == "Tumor")
vsd_sub <- assay(vsd)[all_genes_present, tumor_idx, drop = FALSE]

# Match clinical data
tumor_barcodes <- colnames(vsd_sub)
tumor_patients <- substr(tumor_barcodes, 1, 12)

clinical_sub <- clinical[clinical$patient_id %in% tumor_patients, ]

# Build OS variables (prefer death, else last follow-up)
clinical_sub$OS_time <- as.numeric(ifelse(
  is.na(clinical_sub$days_to_death),
  clinical_sub$days_to_last_follow_up,
  clinical_sub$days_to_death
))
clinical_sub$OS_event <- ifelse(clinical_sub$vital_status == "Dead", 1,
                          ifelse(clinical_sub$vital_status == "Alive", 0, NA))

# Filter out missing OS data
valid_idx <- which(!is.na(clinical_sub$OS_time) & !is.na(clinical_sub$OS_event))
clinical_sub <- clinical_sub[valid_idx, ]
# Reorder expression columns to align with clinical patient order
match_idx <- match(clinical_sub$patient_id, substr(colnames(vsd_sub), 1, 12))
vsd_sub <- vsd_sub[, match_idx, drop = FALSE]

# Optionally filter very low-expression genes (robust)
keep_genes <- rowSums(is.finite(vsd_sub)) >= (0.8 * ncol(vsd_sub))
vsd_sub <- vsd_sub[keep_genes, , drop = FALSE]

# Initialize results table
surv_results <- data.frame(
  Gene     = rownames(vsd_sub),
  HR       = NA_real_,
  lower_CI = NA_real_,
  upper_CI = NA_real_,
  pvalue   = NA_real_
)

# For combined PDF (main plots only)
km_plots <- list()

# ------- Loop over genes with robust printing into devices -------
for (gene in rownames(vsd_sub)) {
  tryCatch({
    gene_expr <- as.numeric(vsd_sub[gene, ])
    med <- median(gene_expr, na.rm = TRUE)
    surv_group <- ifelse(gene_expr > med, "High", "Low")

    keep <- is.finite(gene_expr) & !is.na(surv_group) &
            !is.na(clinical_sub$OS_time) & !is.na(clinical_sub$OS_event)

    gene_df <- data.frame(
      OS_time  = clinical_sub$OS_time[keep],
      OS_event = clinical_sub$OS_event[keep],
      Group    = factor(surv_group[keep], levels = c("Low","High"))
    )

    # Require at least 10 per group
    if (length(unique(gene_df$Group)) < 2 || any(table(gene_df$Group) < 10)) {
      message("⚠ Skipping ", gene, ": not enough samples per group")
      next
    }

    # Cox model
    fit_cox <- coxph(Surv(OS_time, OS_event) ~ Group, data = gene_df)
    s <- summary(fit_cox)

    surv_results[surv_results$Gene == gene, c("HR","lower_CI","upper_CI","pvalue")] <-
      c(s$coefficients[1, "exp(coef)"],
        s$conf.int[1, "lower .95"],
        s$conf.int[1, "upper .95"],
        s$coefficients[1, "Pr(>|z|)"])

    # KM fit
    fit_km <- survfit(Surv(OS_time, OS_event) ~ Group, data = gene_df)

    # Build ggsurvplot object (includes risk table)
    ggs <- ggsurvplot(
      fit_km,
      data       = gene_df,
      pval       = TRUE,
      risk.table = TRUE,
      title      = paste("Survival by", gene, "expression"),
      legend.title = "Expression",
      legend.labs  = c("Low", "High")
    )

    # --- Save INDIVIDUAL PDF: must **print(ggs)** to draw both plot & table ---
    pdf(file.path("KM_Plots", paste0("KM_", gene, ".pdf")), width = 6, height = 5)
    print(ggs)                # <- critical: prints the arranged ggsurvplot object
    dev.off()

    # Also save PNG (hi-res)
    png(file.path("KM_Plots", paste0("KM_", gene, ".png")), width = 2000, height = 1600, res = 300)
    print(ggs)
    dev.off()

    # For combined PDF (main plot only)
    km_plots[[gene]] <- ggs$plot

  }, error = function(e) {
    message("❌ Error for ", gene, ": ", conditionMessage(e))
  })
}

# ------- Save combined KM plots PDF if any collected -------
if (length(km_plots) > 0) {
  pdf(file.path("KM_Plots", "KM_AllGenes.pdf"), width = 10, height = 12)
  do.call(grid.arrange, c(km_plots, ncol = 2))
  dev.off()

  # Optional: combined PNG
  png(file.path("KM_Plots", "KM_AllGenes.png"), width = 2400, height = 2800, res = 300)
  do.call(grid.arrange, c(km_plots, ncol = 2))
  dev.off()
} else {
  message("ℹ No KM plots were generated for combination (insufficient groups for all genes).")
}

# ------- Save survival summary table -------
write.csv(surv_results, file.path("KM_Plots", "Survival_Analysis_SelectedGenes.csv"), row.names = FALSE)
message("✅ Survival analysis completed. PDFs/PNGs and CSV saved in KM_Plots/")


# ===============================================================

library(ggplot2)
library(dplyr)

# Load survival results
surv_results <- read.csv("KM_Plots/Survival_Analysis_SelectedGenes.csv")

# Filter out genes with NA HR
surv_results <- surv_results %>% filter(!is.na(HR))

# Create significance category: Significant, Borderline, Not Significant
surv_results <- surv_results %>%
  mutate(Significance = case_when(
    pvalue < 0.05 ~ "Significant",
    pvalue >= 0.05 & pvalue < 0.1 ~ "Borderline",
    TRUE ~ "Not Significant"
  ))


# Order genes for plotting and add gene-type grouping
# ------------------------------
# Add gene functional categories (for Figure 4D readability)
# ------------------------------
gene_type_map <- data.frame(
  Gene = c(
    # Oncogenes
    "PIK3CA","AKT1","ERBB2","MYC","CCND1","FGFR1",
    # Tumor suppressors
    "BRCA1","BRCA2","TP53","PTEN","RB1","CHEK2","CDH1","MAP3K1",
    # Hormone receptor–related
    "ESR1","PGR","GATA3"
  ),
  GeneType = c(
    rep("Oncogene", 6),
    rep("Tumor suppressor", 8),
    rep("Hormone receptor–related", 3)
  ),
  stringsAsFactors = FALSE
)

surv_results <- dplyr::left_join(surv_results, gene_type_map, by = "Gene")
surv_results$GeneType <- factor(
  surv_results$GeneType,
  levels = c("Oncogene", "Tumor suppressor", "Hormone receptor–related")
)

# Order genes by functional category, then by HR (optional: keeps plot tidy)
surv_results <- surv_results %>%
  arrange(GeneType, desc(HR))

# Set factor levels so the first rows appear at the top within facets
surv_results$Gene <- factor(surv_results$Gene, levels = rev(surv_results$Gene))

# Create forest plot (genes grouped by functional category)
p <- ggplot(
  surv_results,
  aes(x = HR, y = Gene)
) +
  geom_pointrange(aes(xmin = lower_CI, xmax = upper_CI, color = Significance), size = 0.9) +
  geom_vline(xintercept = 1, linetype = 2, color = "red") +
  facet_grid(GeneType ~ ., scales = "free_y", space = "free_y", switch = "y") +
  xlab("Hazard Ratio (HR)") + ylab("") +
  ggtitle("Overall survival associations for selected genes (TCGA-BRCA)") +
  scale_color_manual(values = c(
    "Significant" = "steelblue",
    "Borderline" = "skyblue",
    "Not Significant" = "gray"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    panel.spacing.y = unit(0.6, "lines")
  )

# Save the plot as PNG (journal format)
ggsave("KM_Plots/Forest_HR_AllGenes_Grouped.png", plot = p, width = 7.5, height = 8.5, dpi = 300, bg = "white")
message("✅ Forest plot saved as KM_Plots/Forest_HR_AllGenes_Grouped.png")
message("✅ Forest plot saved as KM_Plots/Forest_HR_AllGenes_Borderline.pdf")

# ================================
# Build Figure 4 composite (KM x3 + Forest)
# ================================
packages <- c("ggplot2", "cowplot", "png", "grid")
for (p in packages) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(ggplot2); library(cowplot); library(png); library(grid)

# Paths to KM and forest images (edit names if yours differ)
km_dir <- "KM_Plots"
fig_dir <- "Figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir)

km_files <- c(
  file.path(km_dir, "KM_BRCA1.png"),
  file.path(km_dir, "KM_BRCA2.png"),
  file.path(km_dir, "KM_TP53.png")
)
forest_file <- file.path(km_dir, "Forest_HR_AllGenes_Grouped.png")

# Check existence
missing <- c(km_files[!file.exists(km_files)], forest_file[!file.exists(forest_file)])
if (length(missing) > 0) {
  stop("Missing required image(s):\n", paste(missing, collapse = "\n"))
}

# Helper: read PNG into ggplot-friendly object
as_gg <- function(path) {
  img <- png::readPNG(path)
  ggplot() + theme_void() +
    annotation_custom(rasterGrob(img, interpolate = TRUE), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
}

p_brca1 <- as_gg(km_files[1]) + ggtitle("A) BRCA1: Overall Survival")
p_brca2 <- as_gg(km_files[2]) + ggtitle("B) BRCA2: Overall Survival")
p_tp53  <- as_gg(km_files[3]) + ggtitle("C) TP53: Overall Survival")
p_forest <- as_gg(forest_file) + ggtitle("D) Cox Forest Plot: 17-Gene Panel")

# Titles style
for (p in list(p_brca1, p_brca2, p_tp53, p_forest)) {
  p$labels$title <- p$labels$title # no-op; just clarifies intent
}
title_theme <- theme(plot.title = element_text(hjust = 0, face = "bold", size = 12))
p_brca1 <- p_brca1 + title_theme
p_brca2 <- p_brca2 + title_theme
p_tp53  <- p_tp53  + title_theme
p_forest <- p_forest + title_theme

# Assemble: top row (3 KMs), bottom row (forest full width)
top_row <- plot_grid(p_brca1, p_brca2, p_tp53, ncol = 3, rel_widths = c(1,1,1))
layout <- plot_grid(top_row, p_forest, ncol = 1, rel_heights = c(1, 1.2))

# Save PDF + PNG (white background)
ggsave(file.path(fig_dir, "Fig4_Composite.pdf"), plot = layout, width = 12, height = 10, device = cairo_pdf, bg = "white")
ggsave(file.path(fig_dir, "Fig4_Composite.png"), plot = layout, width = 12, height = 10, dpi = 300, bg = "white")

message("✅ Saved composite to Figures/Fig4_Composite.pdf and .png")











#####check the ci, hr ratio for these 17 genes, 

# ===============================================================
# Summarize survival results for all 17 genes
# ===============================================================

# Load survival results CSV
surv_results <- read.csv("KM_Plots/Survival_Analysis_SelectedGenes.csv")

# Ensure expected columns exist
if (!all(c("Gene", "HR", "lower_CI", "upper_CI", "pvalue") %in% colnames(surv_results))) {
  stop("❌ CSV missing required columns: Gene, HR, lower_CI, upper_CI, pvalue")
}

# Loop through all rows (genes)
for (i in seq_len(nrow(surv_results))) {
  gene <- surv_results$Gene[i]
  HR <- surv_results$HR[i]
  lower_CI <- surv_results$lower_CI[i]
  upper_CI <- surv_results$upper_CI[i]
  pval <- surv_results$pvalue[i]
  
  cat("\n", gene, " Survival Analysis Results\n", sep = "")
  cat(strrep("-", nchar(gene) + 28), "\n", sep = "")
  
  if (is.na(HR)) {
    cat("No valid survival analysis results (insufficient data).\n")
    next
  }
  
  cat("Hazard Ratio (HR): ", round(HR, 3), "\n")
  cat("95% CI: ", round(lower_CI, 3), " - ", round(upper_CI, 3), "\n")
  cat("p-value: ", signif(pval, 3), "\n\n")
  
  # Interpretation
  if (HR > 1) {
    cat("Interpretation: Higher ", gene, " expression is associated with a HIGHER risk of death", sep = "")
  } else {
    cat("Interpretation: Higher ", gene, " expression is associated with a LOWER risk of death", sep = "")
  }
  
  if (pval < 0.05) {
    cat(" (statistically significant).\n")
  } else if (pval < 0.1) {
    cat(" (borderline, suggestive but not significant).\n")
  } else {
    cat(" (not statistically significant).\n")
  }
}

cat("\n✅ Completed survival interpretation for all genes.\n")

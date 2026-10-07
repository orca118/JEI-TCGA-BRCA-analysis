# ===============================================================
# CLEAN SELF-CONTAINED SURVIVAL ANALYSIS SCRIPT (FIXED)
# ===============================================================

library(survival)
library(survminer)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(TCGAbiolinks)
library(DESeq2)

# ------------------------------
# Load saved expression object
# ------------------------------
vsd <- readRDS("TCGA_BRCA_vsd_geneSymbols.rds")

# ------------------------------
# Pull clinical data
# ------------------------------
clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

if ("submitter_id" %in% colnames(clinical)) {
  clinical$patient_id <- substr(clinical$submitter_id, 1, 12)
} else if ("bcr_patient_barcode" %in% colnames(clinical)) {
  clinical$patient_id <- substr(clinical$bcr_patient_barcode, 1, 12)
} else {
  stop("Cannot find submitter_id or bcr_patient_barcode in clinical table.")
}

# ------------------------------
# Gene panel
# ------------------------------
gene_panel <- c("TP53","MYC","CHEK2","CCND1","MAP3K1",
                "PGR","ESR1","GATA3","FGFR1","PIK3CA",
                "PTEN","RB1","AKT1","ERBB2","CDH1",
                "BRCA1","BRCA2")

# ------------------------------
# Subset tumor samples
# ------------------------------
all_genes_present <- intersect(gene_panel, rownames(assay(vsd)))

tumor_idx <- which(colData(vsd)$condition == "Tumor")
vsd_sub <- assay(vsd)[all_genes_present, tumor_idx, drop = FALSE]

# ------------------------------
# Match clinical data (FIXED ORDER)
# ------------------------------
tumor_barcodes <- colnames(vsd_sub)
tumor_patients <- substr(tumor_barcodes, 1, 12)

clinical_sub <- clinical[clinical$patient_id %in% tumor_patients, ]

# ------------------------------
# Build survival variables
# ------------------------------
clinical_sub$OS_time <- as.numeric(ifelse(
  is.na(clinical_sub$days_to_death),
  clinical_sub$days_to_last_follow_up,
  clinical_sub$days_to_death
))

clinical_sub$OS_event <- ifelse(clinical_sub$vital_status == "Dead", 1,
                          ifelse(clinical_sub$vital_status == "Alive", 0, NA))

# ------------------------------
# Save confirmation
# ------------------------------
print("Script loaded successfully and ready for analysis.")

#!/usr/bin/env Rscript
# =============================================================================
# Cross-Species Medulloblastoma Classification
# Random Forest in PCA space - unified cross-species classification
# =============================================================================

library(Biobase)
library(GEOquery)
library(dplyr)
library(tidyr)
library(ggplot2)
library(randomForest)

cat("=== CROSS-SPECIES MB CLASSIFICATION (Random Forest) ===\n\n")

# ============================================================================
# 1. LOAD DATA
# ============================================================================
cat("1. Loading data...\n")

# --- Human (Cavalli et al., Cancer Cell 2017) ---
human_list <- readRDS("data/human/GSE85217_raw.rds")
expr_human <- exprs(human_list[[1]])
pdata_human <- pData(human_list[[1]])
subgroups <- pdata_human[["subgroup:ch1"]]

# Convert Ensembl IDs to gene symbols
ortho_human <- readRDS("data/orthology/human_ensembl_to_symbol.rds")
symbol_map <- setNames(ortho_human$hgnc_symbol, ortho_human$ensembl_gene_id)
rownames(expr_human) <- symbol_map[gsub("_at$", "", rownames(expr_human))]
expr_human <- expr_human[!is.na(rownames(expr_human)) & rownames(expr_human) != "", ]
expr_agg <- aggregate(expr_human, by = list(Gene = rownames(expr_human)), FUN = mean)
rownames(expr_agg) <- expr_agg$Gene
expr_human <- as.matrix(expr_agg[, -1])

# --- Zebrafish (Casey et al., Cell Reports 2024) ---
zf_counts <- read.table("data/zebrafish/GSE242897/GSE242897_19521raw_counts.txt.gz",
                        header = TRUE, row.names = 1, sep = "\t")
zmb_cols <- grep("PTCH1_P53MUT", colnames(zf_counts), value = TRUE)
zf_expr <- log2(as.matrix(zf_counts[, zmb_cols]) + 1)

ortho_zf <- readRDS("data/orthology/zebrafish_to_human.rds")
zf_to_human <- setNames(ortho_zf$hsapiens_homolog_associated_gene_name, ortho_zf$ensembl_gene_id)
rownames(zf_expr) <- zf_to_human[rownames(zf_expr)]
zf_expr <- zf_expr[!is.na(rownames(zf_expr)) & rownames(zf_expr) != "", ]
zf_agg <- aggregate(zf_expr, by = list(Gene = rownames(zf_expr)), FUN = mean)
rownames(zf_agg) <- zf_agg$Gene
zf_expr <- as.matrix(zf_agg[, -1])

# --- Chicken/HDMB03 xenografts (This study) - Xenografts only ---
chicken <- read.table("data/chicken/BULK_Chicken_expression.tsv",
                      header = TRUE, row.names = 1, sep = "\t")
# Xenograft samples (parenchyma) AND pre-graft HDMB03 culture samples (
# compares the cell line before and after grafting)
xeno_cols <- grep("parenchyma", colnames(chicken), value = TRUE)
chicken_expr <- log2(as.matrix(chicken[, xeno_cols]) + 1)
rownames(chicken_expr) <- toupper(rownames(chicken_expr))
cult_cols <- grep("HDMB03", colnames(chicken), value = TRUE)
culture_expr <- log2(as.matrix(chicken[, cult_cols]) + 1)
rownames(culture_expr) <- toupper(rownames(culture_expr))

# --- Mouse (GSE65888: Kawauchi, Gibson) ---
gse65888 <- getGEO(filename = "data/mouse/GSE65888-GPL1261_series_matrix.txt.gz", getGPL = FALSE)
expr_mouse <- exprs(gse65888)
pdata_mouse <- pData(gse65888)

gpl <- getGEO(filename = "data/mouse/GPL1261.soft.gz")
gpl_table <- Table(gpl)
probe_to_symbol <- setNames(gpl_table[["Gene Symbol"]], gpl_table$ID)
rownames(expr_mouse) <- probe_to_symbol[rownames(expr_mouse)]
expr_mouse <- expr_mouse[!is.na(rownames(expr_mouse)) & rownames(expr_mouse) != "", ]

# Convert mouse to human symbols
ortho_mouse <- readRDS("data/orthology/mouse_to_human.rds")
mouse_to_human <- setNames(ortho_mouse$human_symbol, ortho_mouse$mouse_symbol)
rownames(expr_mouse) <- mouse_to_human[rownames(expr_mouse)]
expr_mouse <- expr_mouse[!is.na(rownames(expr_mouse)) & rownames(expr_mouse) != "", ]
mouse_agg <- aggregate(expr_mouse, by = list(Gene = rownames(expr_mouse)), FUN = mean)
rownames(mouse_agg) <- mouse_agg$Gene
expr_mouse <- as.matrix(mouse_agg[, -1])

# Annotate mouse models
mouse_info <- data.frame(
  sample_id = colnames(expr_mouse),
  title = pdata_mouse$title,
  stringsAsFactors = FALSE
)
mouse_info$model <- case_when(
  grepl("Ptch", mouse_info$title) & grepl("Cdkn2c", mouse_info$title) ~ "SHH_Ptch",
  grepl("Nestin-Cre", mouse_info$title) & grepl("p53", mouse_info$title) ~ "SHH_Nestin",
  grepl("Blbp-MYC|Gad2-MYC|Ptf1a-MYC", mouse_info$title, ignore.case = TRUE) ~ "G3_IUE",
  grepl("Atoh1-CreER|Prom1-CreER|MycRFP|Myc-RFP", mouse_info$title, ignore.case = TRUE) ~ "G3_Ortho",
  TRUE ~ NA_character_
)
mouse_info <- mouse_info[!is.na(mouse_info$model), ]
expr_mouse <- expr_mouse[, mouse_info$sample_id]

# --- Mouse GTML (GSE36594: Swartling - MYCN-driven MB) ---
gse36594 <- getGEO(filename = "data/mouse_gtml/GSE36594_series_matrix.txt.gz", getGPL = FALSE)
expr_gtml <- exprs(gse36594)
pdata_gtml <- pData(gse36594)

# Get GTML tumor samples (32 tumors, exclude cell line)
gtml_idx <- grep("GTML Medulloblastoma Tumor", pdata_gtml$title)
expr_gtml <- expr_gtml[, gtml_idx]
pdata_gtml <- pdata_gtml[gtml_idx, ]

# For Mouse Exon array, use gene symbol from characteristics
# The rownames are probeset IDs - we need to convert to gene symbols
# Load GPL6096 annotation
gpl6096 <- getGEO(filename = "data/mouse_gtml/GPL6096.soft.gz")
gpl6096_table <- Table(gpl6096)
# This array has "gene_assignment" column with gene symbols
gene_assign <- gpl6096_table[["gene_assignment"]]
names(gene_assign) <- gpl6096_table$ID
# Extract first gene symbol from assignment string
extract_symbol <- function(x) {
  if (is.na(x) || x == "---") return(NA)
  parts <- strsplit(as.character(x), " // ")[[1]]
  if (length(parts) >= 2) return(parts[2])
  return(NA)
}
probe_to_symbol_gtml <- sapply(gene_assign, extract_symbol)
rownames(expr_gtml) <- probe_to_symbol_gtml[rownames(expr_gtml)]
expr_gtml <- expr_gtml[!is.na(rownames(expr_gtml)) & rownames(expr_gtml) != "", ]

# Aggregate by gene
gtml_agg <- aggregate(expr_gtml, by = list(Gene = rownames(expr_gtml)), FUN = mean)
rownames(gtml_agg) <- gtml_agg$Gene
expr_gtml <- as.matrix(gtml_agg[, -1])

# Convert mouse to human symbols (reuse orthology)
rownames(expr_gtml) <- mouse_to_human[rownames(expr_gtml)]
expr_gtml <- expr_gtml[!is.na(rownames(expr_gtml)) & rownames(expr_gtml) != "", ]
gtml_agg2 <- aggregate(expr_gtml, by = list(Gene = rownames(expr_gtml)), FUN = mean)
rownames(gtml_agg2) <- gtml_agg2$Gene
expr_gtml <- as.matrix(gtml_agg2[, -1])

# --- Mouse GTML_2023 and GMYC (GSE139240: Swartling/Mainwaring - Nature Comm 2023) ---
# RNA-seq counts with mouse Ensembl IDs
gse139240 <- read.table("data/mouse_gmyc/GSE139240_Counts_matrix.txt.gz",
                        header = TRUE, row.names = 1, sep = "\t")

# Split into GTML_2023 and GMYC conditions
gtml2023_cols <- c("GTML_1", "GTML_2", "GTML_3", "GTML_4", "GTML_5", "GTML_6", "GTML_7", "GTML_8")
gmyc_cols <- c("GMYC_1", "GMYC_2", "GMYC_3", "GMYC_4", "GMYC_5", "GMYC_6")

expr_gtml2023 <- as.matrix(gse139240[, gtml2023_cols])
expr_gmyc <- as.matrix(gse139240[, gmyc_cols])

# Log2 transform (RNA-seq counts)
expr_gtml2023 <- log2(expr_gtml2023 + 1)
expr_gmyc <- log2(expr_gmyc + 1)

# Convert mouse Ensembl IDs to symbols using biomaRt (if available) or simple mapping
# First, clean Ensembl IDs (remove version number)
rownames(expr_gtml2023) <- gsub("\\..*", "", rownames(expr_gtml2023))
rownames(expr_gmyc) <- gsub("\\..*", "", rownames(expr_gmyc))

# Load or create mouse Ensembl to symbol mapping
mouse_ens_file <- "data/orthology/mouse_ensembl_to_symbol.rds"
if (file.exists(mouse_ens_file)) {
  mouse_ens <- readRDS(mouse_ens_file)
} else {
  # Create mapping using biomaRt
  cat("   Creating mouse Ensembl to symbol mapping...\n")
  library(biomaRt)
  mart <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror = "useast")
  mouse_ens <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
                     filters = "ensembl_gene_id",
                     values = unique(c(rownames(expr_gtml2023), rownames(expr_gmyc))),
                     mart = mart)
  saveRDS(mouse_ens, mouse_ens_file)
}
mouse_ens_to_symbol <- setNames(mouse_ens$mgi_symbol, mouse_ens$ensembl_gene_id)

# Convert to mouse symbols
rownames(expr_gtml2023) <- mouse_ens_to_symbol[rownames(expr_gtml2023)]
rownames(expr_gmyc) <- mouse_ens_to_symbol[rownames(expr_gmyc)]

# Remove NA and empty rownames
expr_gtml2023 <- expr_gtml2023[!is.na(rownames(expr_gtml2023)) & rownames(expr_gtml2023) != "", ]
expr_gmyc <- expr_gmyc[!is.na(rownames(expr_gmyc)) & rownames(expr_gmyc) != "", ]

# Aggregate by gene symbol
gtml2023_agg <- aggregate(expr_gtml2023, by = list(Gene = rownames(expr_gtml2023)), FUN = mean)
rownames(gtml2023_agg) <- gtml2023_agg$Gene
expr_gtml2023 <- as.matrix(gtml2023_agg[, -1])

gmyc_agg <- aggregate(expr_gmyc, by = list(Gene = rownames(expr_gmyc)), FUN = mean)
rownames(gmyc_agg) <- gmyc_agg$Gene
expr_gmyc <- as.matrix(gmyc_agg[, -1])

# Convert mouse to human symbols (using same orthology as above)
rownames(expr_gtml2023) <- mouse_to_human[rownames(expr_gtml2023)]
rownames(expr_gmyc) <- mouse_to_human[rownames(expr_gmyc)]

expr_gtml2023 <- expr_gtml2023[!is.na(rownames(expr_gtml2023)) & rownames(expr_gtml2023) != "", ]
expr_gmyc <- expr_gmyc[!is.na(rownames(expr_gmyc)) & rownames(expr_gmyc) != "", ]

# Final aggregation
gtml2023_agg2 <- aggregate(expr_gtml2023, by = list(Gene = rownames(expr_gtml2023)), FUN = mean)
rownames(gtml2023_agg2) <- gtml2023_agg2$Gene
expr_gtml2023 <- as.matrix(gtml2023_agg2[, -1])

gmyc_agg2 <- aggregate(expr_gmyc, by = list(Gene = rownames(expr_gmyc)), FUN = mean)
rownames(gmyc_agg2) <- gmyc_agg2$Gene
expr_gmyc <- as.matrix(gmyc_agg2[, -1])

cat("   Human:", ncol(expr_human), "| Zebrafish:", ncol(zf_expr),
    "| Chicken:", ncol(chicken_expr), "| HDMB03 culture:", ncol(culture_expr),
    "| Mouse GSE65888:", ncol(expr_mouse),
    "| GTML:", ncol(expr_gtml), "| GTML_2023:", ncol(expr_gtml2023),
    "| GMYC:", ncol(expr_gmyc), "\n")

# ============================================================================
# 2. GENE SELECTION & NORMALIZATION
# ============================================================================
cat("2. Selecting genes and normalizing...\n")

# Top 2500 variable genes in human
N_GENES <- 2500
gene_vars <- apply(expr_human, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:N_GENES])

# Find common genes across all species (including all mouse Group 3 models)
common_genes <- Reduce(intersect, list(top_genes, rownames(zf_expr),
                                        rownames(chicken_expr), rownames(expr_mouse),
                                        rownames(expr_gtml), rownames(expr_gtml2023),
                                        rownames(expr_gmyc)))
cat("   Common genes:", length(common_genes), "\n")

# Z-score using human reference parameters
human_means <- rowMeans(expr_human[common_genes, ])
human_sds <- apply(expr_human[common_genes, ], 1, sd)
human_sds[human_sds == 0] <- 1

normalize_to_human <- function(expr, genes, means, sds) {
  avail <- intersect(genes, rownames(expr))
  expr_z <- (expr[avail, ] - means[avail]) / sds[avail]
  expr_z[!is.finite(expr_z)] <- 0
  return(t(expr_z))
}

X_human <- normalize_to_human(expr_human, common_genes, human_means, human_sds)
X_zf <- normalize_to_human(zf_expr, common_genes, human_means, human_sds)
X_ch <- normalize_to_human(chicken_expr, common_genes, human_means, human_sds)
X_cult <- normalize_to_human(culture_expr, common_genes, human_means, human_sds)
X_mouse <- normalize_to_human(expr_mouse, common_genes, human_means, human_sds)
X_gtml <- normalize_to_human(expr_gtml, common_genes, human_means, human_sds)
X_gtml2023 <- normalize_to_human(expr_gtml2023, common_genes, human_means, human_sds)
X_gmyc <- normalize_to_human(expr_gmyc, common_genes, human_means, human_sds)

# ============================================================================
# 3. PCA & PROJECTION
# ============================================================================
cat("3. Running PCA on human data...\n")

pca <- prcomp(X_human, center = FALSE, scale. = FALSE)
var_exp <- summary(pca)$importance[2, 1:5] * 100
cat("   Variance explained (PC1-5):", paste(round(var_exp, 1), collapse = ", "), "%\n")

N_PCS <- 20
human_pcs <- pca$x[, 1:N_PCS]
zf_pcs <- X_zf %*% pca$rotation[, 1:N_PCS]
ch_pcs <- X_ch %*% pca$rotation[, 1:N_PCS]
cult_pcs <- X_cult %*% pca$rotation[, 1:N_PCS]
mouse_pcs <- X_mouse %*% pca$rotation[, 1:N_PCS]
gtml_pcs <- X_gtml %*% pca$rotation[, 1:N_PCS]
gtml2023_pcs <- X_gtml2023 %*% pca$rotation[, 1:N_PCS]
gmyc_pcs <- X_gmyc %*% pca$rotation[, 1:N_PCS]

# ============================================================================
# 4. TRAIN RANDOM FOREST
# ============================================================================
cat("4. Training Random Forest on PCA coordinates...\n")

set.seed(42)
train_idx <- sample(seq_len(nrow(human_pcs)), size = round(0.8 * nrow(human_pcs)))
test_idx <- setdiff(seq_len(nrow(human_pcs)), train_idx)

X_train <- human_pcs[train_idx, ]
y_train <- factor(subgroups[train_idx])
X_test <- human_pcs[test_idx, ]
y_test <- factor(subgroups[test_idx])

cat("   Train:", nrow(X_train), "| Test:", nrow(X_test), "\n")

set.seed(42)
rf_model <- randomForest(x = X_train, y = y_train, ntree = 1000, importance = TRUE)
cat("   OOB error rate:", round(rf_model$err.rate[1000, "OOB"] * 100, 1), "%\n")

# ============================================================================
# 5. PREDICT ALL SAMPLES
# ============================================================================
cat("5. Classifying all samples...\n")

# Human test set
pred_test <- predict(rf_model, X_test)
probs_test <- predict(rf_model, X_test, type = "prob")
human_acc <- mean(pred_test == y_test) * 100
cat("   Human test accuracy:", round(human_acc, 1), "%\n")

# Zebrafish
pred_zf <- predict(rf_model, zf_pcs)
probs_zf <- predict(rf_model, zf_pcs, type = "prob")
cat("   Zebrafish predictions:", as.character(pred_zf), "\n")

# Chicken
pred_ch <- predict(rf_model, ch_pcs)
probs_ch <- predict(rf_model, ch_pcs, type = "prob")
cat("   Chicken predictions:", as.character(pred_ch), "\n")

# HDMB03 culture (pre-graft)
pred_cult <- predict(rf_model, cult_pcs)
probs_cult <- predict(rf_model, cult_pcs, type = "prob")
cat("   HDMB03 culture predictions:", as.character(pred_cult), "\n")

# Mouse GSE65888
pred_mouse <- predict(rf_model, mouse_pcs)
probs_mouse <- predict(rf_model, mouse_pcs, type = "prob")

# Mouse GTML (Swartling 2012)
pred_gtml <- predict(rf_model, gtml_pcs)
probs_gtml <- predict(rf_model, gtml_pcs, type = "prob")
cat("   GTML predictions:", table(pred_gtml), "\n")

# Mouse GTML_2023 (Swartling/Mainwaring 2023)
pred_gtml2023 <- predict(rf_model, gtml2023_pcs)
probs_gtml2023 <- predict(rf_model, gtml2023_pcs, type = "prob")
cat("   GTML_2023 predictions:", table(pred_gtml2023), "\n")

# Mouse GMYC (Swartling/Mainwaring 2023)
pred_gmyc <- predict(rf_model, gmyc_pcs)
probs_gmyc <- predict(rf_model, gmyc_pcs, type = "prob")
cat("   GMYC predictions:", table(pred_gmyc), "\n")

# ============================================================================
# 6. COMPILE RESULTS
# ============================================================================
# IMPORTANT: GTML and GMYC are ALL Group 3 models (MYCN/MYC-driven, SHH-independent)
# Force actual = "Group3" for all these samples regardless of classifier prediction
all_results <- rbind(
  data.frame(sample_id = rownames(X_test), species = "Human",
             model = as.character(y_test), actual = as.character(y_test),
             predicted = as.character(pred_test), probs_test, stringsAsFactors = FALSE),
  data.frame(sample_id = rownames(zf_pcs), species = "Zebrafish", model = "ZMB",
             actual = "SHH", predicted = as.character(pred_zf),
             probs_zf, stringsAsFactors = FALSE),
  data.frame(sample_id = rownames(ch_pcs), species = "Chicken",
             model = "Xenograft",
             actual = "Group3", predicted = as.character(pred_ch),
             probs_ch, stringsAsFactors = FALSE),
  data.frame(sample_id = rownames(cult_pcs), species = "Chicken_culture",
             model = "Culture",
             actual = "Group3", predicted = as.character(pred_cult),
             probs_cult, stringsAsFactors = FALSE),
  data.frame(sample_id = rownames(mouse_pcs), species = "Mouse",
             model = mouse_info$model,
             actual = case_when(grepl("SHH", mouse_info$model) ~ "SHH",
                                grepl("WNT", mouse_info$model) ~ "WNT",
                                grepl("G3", mouse_info$model) ~ "Group3"),
             predicted = as.character(pred_mouse), probs_mouse, stringsAsFactors = FALSE),
  # GTML (Swartling 2012) - FORCED Group3
  data.frame(sample_id = rownames(gtml_pcs), species = "Mouse_GTML",
             model = "GTML",
             actual = "Group3", predicted = as.character(pred_gtml),
             probs_gtml, stringsAsFactors = FALSE),
  # GTML_2023 (Swartling/Mainwaring 2023) - FORCED Group3
  data.frame(sample_id = rownames(gtml2023_pcs), species = "Mouse_GTML2023",
             model = "GTML_2023",
             actual = "Group3", predicted = as.character(pred_gtml2023),
             probs_gtml2023, stringsAsFactors = FALSE),
  # GMYC (Swartling/Mainwaring 2023) - FORCED Group3
  data.frame(sample_id = rownames(gmyc_pcs), species = "Mouse_GMYC",
             model = "GMYC",
             actual = "Group3", predicted = as.character(pred_gmyc),
             probs_gmyc, stringsAsFactors = FALSE)
)

# ============================================================================
# 7. CREATE FIGURE
# ============================================================================
cat("6. Creating figure...\n")

plot_data <- all_results %>%
  dplyr::select(sample_id, species, model, actual, predicted, WNT, SHH, Group3, Group4) %>%
  pivot_longer(cols = c(WNT, SHH, Group3, Group4),
               names_to = "Subgroup", values_to = "Probability")

plot_data$Subgroup <- factor(plot_data$Subgroup, levels = c("WNT", "SHH", "Group3", "Group4"))

# Add study reference labels
# Note: Mouse models from GSE65888 (Kawauchi), GSE36594 (Swartling 2012), GSE139240 (Swartling 2023)
plot_data$label <- case_when(
  plot_data$species == "Human" & plot_data$model == "WNT" ~ "Cavalli et al.\nCancer Cell 2017\nWNT",
  plot_data$species == "Human" & plot_data$model == "SHH" ~ "Cavalli et al.\nCancer Cell 2017\nSHH",
  plot_data$species == "Human" & plot_data$model == "Group3" ~ "Cavalli et al.\nCancer Cell 2017\nGroup3",
  plot_data$species == "Human" & plot_data$model == "Group4" ~ "Cavalli et al.\nCancer Cell 2017\nGroup4",
  plot_data$species == "Zebrafish" ~ "Casey et al.\nCell Rep 2024\nZMB",
  plot_data$species == "Chicken_culture" ~ "This study\nHDMB03 culture",
  plot_data$species == "Chicken" ~ "This study\nXenograft",
  plot_data$species == "Mouse" & plot_data$model == "SHH_Ptch" ~ "Kawauchi et al.\nCancer Cell 2012\nSHH_Ptch",
  plot_data$species == "Mouse" & plot_data$model == "SHH_Nestin" ~ "Kawauchi et al.\nOncogene 2017\nSHH_Nestin",
  plot_data$species == "Mouse" & plot_data$model == "G3_IUE" ~ "Kawauchi et al.\nOncogene 2017\nG3_IUE",
  plot_data$species == "Mouse" & plot_data$model == "G3_Ortho" ~ "Kawauchi et al.\nCancer Cell 2012\nG3_Ortho",
  plot_data$species == "Mouse_GTML" ~ "Swartling et al.\nCancer Cell 2012\nGTML (n=32)",
  plot_data$species == "Mouse_GTML2023" ~ "Mainwaring et al.\nNat Comm 2023\nGTML (n=8)",
  plot_data$species == "Mouse_GMYC" ~ "Mainwaring et al.\nNat Comm 2023\nGMYC (n=6)",
  TRUE ~ NA_character_
)

order_vec <- c(
 "Cavalli et al.\nCancer Cell 2017\nWNT", "Cavalli et al.\nCancer Cell 2017\nSHH",
 "Cavalli et al.\nCancer Cell 2017\nGroup3", "Cavalli et al.\nCancer Cell 2017\nGroup4",
 "Casey et al.\nCell Rep 2024\nZMB", "This study\nHDMB03 culture", "This study\nXenograft",
 "Kawauchi et al.\nCancer Cell 2012\nSHH_Ptch",
 "Kawauchi et al.\nOncogene 2017\nSHH_Nestin", "Kawauchi et al.\nOncogene 2017\nG3_IUE",
 "Kawauchi et al.\nCancer Cell 2012\nG3_Ortho",
 "Swartling et al.\nCancer Cell 2012\nGTML (n=32)",
 "Mainwaring et al.\nNat Comm 2023\nGTML (n=8)",
 "Mainwaring et al.\nNat Comm 2023\nGMYC (n=6)"
)
plot_data$label <- factor(plot_data$label, levels = order_vec)
plot_data <- plot_data[!is.na(plot_data$label), ]

colors <- c("WNT" = "#3498db", "SHH" = "#e74c3c", "Group3" = "#2ecc71", "Group4" = "#9b59b6")

p <- ggplot(plot_data, aes(x = label, y = Probability, fill = Subgroup)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7, alpha = 0.7, outlier.shape = NA) +
  geom_point(aes(color = Subgroup),
             position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.02, dodge.width = 0.8),
             size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors, guide = "none") +
  geom_vline(xintercept = c(4.5, 7.5), linetype = "dashed", alpha = 0.3) +
  labs(title = "Cross-Species MB Classification (Random Forest in PCA Space)",
       subtitle = sprintf("RF (ntree=1000) on %d PCs from %d genes | Human test accuracy: %.1f%%",
                          N_PCS, length(common_genes), human_acc),
       x = "", y = "RF Probability") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        plot.title = element_text(face = "bold"),
        legend.position = "right")

# Create output directories if needed
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

ggsave("figures/cross_species_RF_PCA.png", p, width = 18, height = 7, dpi = 150)
cat("   Saved: figures/cross_species_RF_PCA.png\n")

# Save results
write.csv(all_results, "results/rf_pca_classification_results.csv", row.names = FALSE)
cat("   Saved: results/rf_pca_classification_results.csv\n")

# ============================================================================
# 8. TERNARY PLOT (SHH / G3 / G4)
# ============================================================================
cat("7. Creating ternary plot (SHH/G3/G4)...\n")

# Load ggtern
if (!require(ggtern, quietly = TRUE)) {
  install.packages("ggtern")
  library(ggtern)
}

# Prepare data for ternary plot
ternary_data <- all_results %>%
  dplyr::select(sample_id, species, model, SHH, Group3, Group4) %>%
  mutate(
    # Normalize to sum to 1 (exclude WNT)
    total = SHH + Group3 + Group4,
    SHH_norm = SHH / total,
    G3_norm = Group3 / total,
    G4_norm = Group4 / total
  )

# Create condition labels (one per condition)
ternary_data$condition <- case_when(
  ternary_data$species == "Human" & ternary_data$model == "SHH" ~ "Human SHH",
  ternary_data$species == "Human" & ternary_data$model == "Group3" ~ "Human G3",
  ternary_data$species == "Human" & ternary_data$model == "Group4" ~ "Human G4",
  ternary_data$species == "Human" & ternary_data$model == "WNT" ~ "Human WNT",
  ternary_data$species == "Zebrafish" ~ "Zebrafish (Casey)",
  ternary_data$species == "Chicken" ~ "HDMB03 Xenograft",
  ternary_data$species == "Chicken_culture" ~ "HDMB03 culture",
  ternary_data$species == "Mouse" & ternary_data$model == "SHH_Ptch" ~ "Mouse SHH_Ptch",
  ternary_data$species == "Mouse" & ternary_data$model == "SHH_Nestin" ~ "Mouse SHH_Nestin",
  ternary_data$species == "Mouse" & ternary_data$model == "G3_IUE" ~ "Mouse G3_IUE",
  ternary_data$species == "Mouse" & ternary_data$model == "G3_Ortho" ~ "Mouse G3_Ortho",
  ternary_data$species == "Mouse_GTML" ~ "Mouse GTML",
  ternary_data$species == "Mouse_GTML2023" ~ "Mouse GTML_2023",
  ternary_data$species == "Mouse_GMYC" ~ "Mouse GMYC",
  TRUE ~ "Other"
)

# Filter out Human WNT for this plot (focus on SHH/G3/G4)
ternary_data <- ternary_data %>% filter(condition != "Human WNT")

# Define colors for conditions
condition_colors <- c(
  "Human SHH" = "#e74c3c",
  "Human G3" = "#2ecc71",
  "Human G4" = "#9b59b6",
  "Zebrafish (Casey)" = "#e67e22",
  "HDMB03 Xenograft" = "#d35400",
  "HDMB03 culture" = "#f39c12",
  "Mouse SHH_Ptch" = "#c0392b",
  "Mouse SHH_Nestin" = "#922b21",
  "Mouse G3_IUE" = "#27ae60",
  "Mouse G3_Ortho" = "#1e8449",
  "Mouse GTML" = "#3498db",
  "Mouse GTML_2023" = "#2980b9",
  "Mouse GMYC" = "#1abc9c"
)

# Calculate centroids for each condition
centroids <- ternary_data %>%
  group_by(condition) %>%
  summarise(
    SHH_norm = mean(SHH_norm),
    G3_norm = mean(G3_norm),
    G4_norm = mean(G4_norm),
    .groups = "drop"
  )

# Define label positions manually for each condition (hand-tuned)
# Format: condition, label_x (SHH), label_y (G3), label_z (G4)
label_positions <- tribble(
  ~condition, ~label_x, ~label_y, ~label_z, ~short_label,
  "Human SHH",        0.85, 0.05, 0.10, "SHH",
  "Human G3",         0.10, 0.85, 0.05, "G3",
  "Human G4",         0.05, 0.10, 0.85, "G4",
  "Zebrafish (Casey)", 0.70, 0.05, 0.25, "ZF",
  "HDMB03 Xenograft", 0.20, 0.65, 0.15, "Xenograft",
  "HDMB03 culture",   0.15, 0.30, 0.55, "Culture",
  "Mouse SHH_Ptch",   0.55, 0.35, 0.10, "SHH_Ptch",
  "Mouse SHH_Nestin", 0.70, 0.15, 0.15, "SHH_Nestin",
  "Mouse G3_IUE",     0.25, 0.60, 0.15, "G3_IUE",
  "Mouse G3_Ortho",   0.35, 0.45, 0.20, "G3_Ortho",
  "Mouse GTML",       0.40, 0.50, 0.10, "GTML",
  "Mouse GTML_2023",  0.50, 0.40, 0.10, "GTML_2023",
  "Mouse GMYC",       0.30, 0.55, 0.15, "GMYC"
)

# Merge with centroids
centroids_labels <- centroids %>%
  left_join(label_positions, by = "condition")

# Create ternary plot with centroids and labels
p_tern <- ggtern(ternary_data, aes(x = SHH_norm, y = G3_norm, z = G4_norm)) +
  # Individual points (alpha 50%)
  geom_point(aes(color = condition), size = 2, alpha = 0.5) +
  # Centroids as larger squares with black border
  geom_point(data = centroids, aes(fill = condition),
             shape = 22, size = 6, alpha = 1, color = "black", stroke = 1) +
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values = condition_colors) +
  labs(
    title = "Cross-Species MB Classification - Ternary Plot",
    subtitle = "Normalized RF probabilities (WNT excluded) | Squares = centroids",
    x = "P(SHH)",
    y = "P(Group3)",
    z = "P(Group4)",
    color = "Condition",
    fill = "Condition"
  ) +
  theme_bw() +
  theme_showarrows() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)),
         fill = guide_legend(override.aes = list(size = 4)))

# Add annotations for each centroid using ggtern's annotate
for (i in 1:nrow(centroids_labels)) {
  row <- centroids_labels[i, ]
  # Arrow from label to centroid
  p_tern <- p_tern +
    annotate("segment",
             x = row$label_x, y = row$label_y, z = row$label_z,
             xend = row$SHH_norm, yend = row$G3_norm, zend = row$G4_norm,
             arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
             color = "gray30", linewidth = 0.3) +
    annotate("text",
             x = row$label_x, y = row$label_y, z = row$label_z,
             label = row$short_label,
             size = 2.5, fontface = "bold")
}

ggsave("figures/cross_species_ternary.png", p_tern, width = 10, height = 8, dpi = 150)
cat("   Saved: figures/cross_species_ternary.png\n")

cat("\n=== DONE ===\n")

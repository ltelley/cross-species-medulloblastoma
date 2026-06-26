#!/usr/bin/env Rscript
# =============================================================================
# Cross-Species MB PCA - TRUE PROJECTION onto Human PCA space
#
# 1. Train PCA on HUMAN data only (763 samples)
# 2. Project animal data onto that PCA space using the same rotation
# This ensures animal samples are positioned relative to human subgroups
# =============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(biomaRt)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

message("\n========================================")
message(">>> Cross-Species MB PCA (Projection onto Human)")
message("========================================")

# =============================================================================
# Load data (using cached files from previous runs)
# =============================================================================

# --- Human ---
message("\n>>> Loading Human data...")
gse_human <- readRDS("data/human/GSE85217_raw.rds")[[1]]
expr_human_raw <- exprs(gse_human)
pheno_human <- pData(gse_human)
pheno_human$subgroup <- gsub(".*: ", "", pheno_human$`subgroup:ch1`)

ens2sym <- readRDS("data/orthology/human_ensembl_to_symbol.rds")
ens2sym_clean <- ens2sym %>% filter(hgnc_symbol != "") %>% distinct(ensembl_gene_id, .keep_all = TRUE)
gene_map <- setNames(ens2sym_clean$hgnc_symbol, ens2sym_clean$ensembl_gene_id)

ensembl_ids <- gsub("_at$", "", rownames(expr_human_raw))
mapped_genes <- gene_map[ensembl_ids]
valid_idx <- !is.na(mapped_genes) & mapped_genes != ""
expr_human <- expr_human_raw[valid_idx, ]
rownames(expr_human) <- toupper(mapped_genes[valid_idx])

if (any(duplicated(rownames(expr_human)))) {
  iqrs <- apply(expr_human, 1, IQR, na.rm = TRUE)
  df <- data.frame(gene = rownames(expr_human), iqr = iqrs, idx = 1:nrow(expr_human))
  df <- df %>% group_by(gene) %>% slice_max(iqr, n = 1, with_ties = FALSE) %>% ungroup()
  expr_human <- expr_human[df$idx, ]
  rownames(expr_human) <- df$gene
}

# --- Mouse ---
message(">>> Loading Mouse data...")
gse_mouse <- getGEO(filename = "data/mouse/GSE33199_series_matrix.txt.gz")
expr_mouse_raw <- exprs(gse_mouse)
pheno_mouse <- pData(gse_mouse)

gpl_mouse <- getGEO("GPL1261", destdir = "data/mouse/")
gpl_table <- Table(gpl_mouse)
probe2gene_mouse <- setNames(as.character(gpl_table$`Gene Symbol`), gpl_table$ID)
probe2gene_mouse <- probe2gene_mouse[!is.na(probe2gene_mouse) & probe2gene_mouse != ""]
common_probes <- intersect(rownames(expr_mouse_raw), names(probe2gene_mouse))
expr_mouse_tmp <- expr_mouse_raw[common_probes, ]
genes_mouse <- probe2gene_mouse[common_probes]
iqrs <- apply(expr_mouse_tmp, 1, IQR, na.rm = TRUE)
df <- data.frame(probe = common_probes, gene = genes_mouse, iqr = iqrs)
best <- df %>% group_by(gene) %>% slice_max(iqr, n = 1, with_ties = FALSE) %>% ungroup()
expr_mouse <- expr_mouse_tmp[best$probe, ]
rownames(expr_mouse) <- toupper(best$gene)

sample_source <- as.character(pheno_mouse$source_name_ch1)
mouse_subgroup <- case_when(
  grepl("Myc tumor|cmyc", sample_source, ignore.case = TRUE) ~ "G3",
  grepl("Mycn tumor", sample_source, ignore.case = TRUE) ~ "G3",
  grepl("Shh-type|Ptc", sample_source, ignore.case = TRUE) ~ "SHH",
  grepl("Wnt-type", sample_source, ignore.case = TRUE) ~ "WNT",
  grepl("control", sample_source, ignore.case = TRUE) ~ "Normal",
  TRUE ~ "Unknown"
)

# --- Zebrafish ---
message(">>> Loading Zebrafish data...")
zf_counts <- read.delim(gzfile("data/zebrafish/GSE242897/GSE242897_19521raw_counts.txt.gz"), row.names = 1)
zf_ortho <- readRDS("data/orthology/zebrafish_to_human.rds")
zf_casey <- zf_ortho %>%
  filter(hsapiens_homolog_perc_id > 30, hsapiens_homolog_associated_gene_name != "") %>%
  group_by(ensembl_gene_id) %>% slice_max(hsapiens_homolog_perc_id, n = 1, with_ties = FALSE) %>% ungroup()
ortho_map_zf <- setNames(toupper(zf_casey$hsapiens_homolog_associated_gene_name), zf_casey$ensembl_gene_id)

common_genes_zf <- intersect(rownames(zf_counts), names(ortho_map_zf))
zf_log2 <- log2(as.matrix(zf_counts[common_genes_zf, ]) + 1)
rownames(zf_log2) <- ortho_map_zf[common_genes_zf]
if (any(duplicated(rownames(zf_log2)))) {
  rmeans <- rowMeans(zf_log2, na.rm = TRUE)
  df <- data.frame(gene = rownames(zf_log2), rmean = rmeans, idx = 1:nrow(zf_log2))
  df <- df %>% arrange(desc(rmean)) %>% distinct(gene, .keep_all = TRUE)
  zf_log2 <- zf_log2[df$idx, ]
  rownames(zf_log2) <- df$gene
}
tumor_cols_zf <- grep("PTCH1", colnames(zf_log2), value = TRUE)
expr_zebrafish <- zf_log2[, tumor_cols_zf]

# --- Chicken ---
message(">>> Loading Chicken data...")
expr_chicken <- as.matrix(read.delim("data/chicken/BULK_Chicken_expression.tsv", row.names = 1, check.names = FALSE))
meta_chicken <- read.delim("data/chicken/BULK_Chicken_metadata.tsv")
tumor_cols_ck <- meta_chicken$sample_id[meta_chicken$condition == "Tumor"]
expr_chicken <- expr_chicken[, tumor_cols_ck, drop = FALSE]
rownames(expr_chicken) <- toupper(rownames(expr_chicken))


# =============================================================================
# COMMON GENES
# =============================================================================
message("\n>>> Finding common genes...")

gene_lists <- list(
  Human = rownames(expr_human),
  Mouse = rownames(expr_mouse),
  Zebrafish = rownames(expr_zebrafish),
  Chicken = rownames(expr_chicken)
)
common_genes <- Reduce(intersect, gene_lists)
message("  Common genes: ", length(common_genes))

expr_human <- expr_human[common_genes, ]
expr_mouse <- expr_mouse[common_genes, ]
expr_zebrafish <- expr_zebrafish[common_genes, ]
expr_chicken <- expr_chicken[common_genes, ]


# =============================================================================
# TRAIN PCA ON HUMAN DATA
# =============================================================================
message("\n>>> Training PCA on Human data...")

# Z-score human genes
human_means <- rowMeans(expr_human)
human_sds <- apply(expr_human, 1, sd)
expr_human_z <- (expr_human - human_means) / human_sds
expr_human_z[is.na(expr_human_z)] <- 0

# Run PCA on human
pca_human <- prcomp(t(expr_human_z), center = TRUE, scale. = TRUE)
var_pct <- summary(pca_human)$importance[2, ] * 100

message("  Human PCA trained")
message("  PC1: ", round(var_pct[1], 1), "% variance")
message("  PC2: ", round(var_pct[2], 1), "% variance")


# =============================================================================
# PROJECT ANIMAL DATA ONTO HUMAN PCA SPACE
# =============================================================================
message("\n>>> Projecting animal data...")

# Function to project data onto PCA space
project_onto_pca <- function(new_data, pca_object, train_means, train_sds) {
  # Z-score using HUMAN gene means/sds
  new_z <- (new_data - train_means) / train_sds
  new_z[is.na(new_z)] <- 0

  # Center using PCA center and scale
  new_centered <- t(new_z) - pca_object$center
  new_scaled <- sweep(new_centered, 2, pca_object$scale, "/")

  # Project onto PC space
  projected <- new_scaled %*% pca_object$rotation

  return(projected)
}

# Project each animal dataset
mouse_proj <- project_onto_pca(expr_mouse, pca_human, human_means, human_sds)
zebrafish_proj <- project_onto_pca(expr_zebrafish, pca_human, human_means, human_sds)
chicken_proj <- project_onto_pca(expr_chicken, pca_human, human_means, human_sds)


# =============================================================================
# BUILD COMBINED RESULTS
# =============================================================================
message("\n>>> Building results...")

# Human PCA coordinates
pca_human_df <- data.frame(
  sample_id = rownames(pca_human$x),
  PC1 = pca_human$x[, 1],
  PC2 = pca_human$x[, 2],
  PC3 = pca_human$x[, 3],
  species = "Human",
  subgroup = pheno_human$subgroup,
  stringsAsFactors = FALSE
)

# Mouse
pca_mouse_df <- data.frame(
  sample_id = rownames(mouse_proj),
  PC1 = mouse_proj[, 1],
  PC2 = mouse_proj[, 2],
  PC3 = mouse_proj[, 3],
  species = "Mouse",
  subgroup = mouse_subgroup,
  stringsAsFactors = FALSE
)

# Zebrafish
pca_zebrafish_df <- data.frame(
  sample_id = rownames(zebrafish_proj),
  PC1 = zebrafish_proj[, 1],
  PC2 = zebrafish_proj[, 2],
  PC3 = zebrafish_proj[, 3],
  species = "Zebrafish",
  subgroup = "SHH",
  stringsAsFactors = FALSE
)

# Chicken
pca_chicken_df <- data.frame(
  sample_id = rownames(chicken_proj),
  PC1 = chicken_proj[, 1],
  PC2 = chicken_proj[, 2],
  PC3 = chicken_proj[, 3],
  species = "Chicken",
  subgroup = "G3",
  stringsAsFactors = FALSE
)

# Combine
pca_df <- rbind(pca_human_df, pca_mouse_df, pca_zebrafish_df, pca_chicken_df)


# =============================================================================
# VISUALIZATION
# =============================================================================
message("\n>>> Creating figures...")

sg_colors <- c(WNT = "#FFD700", SHH = "#E31A1C",
               Group3 = "#33A02C", G3 = "#33A02C",
               Group4 = "#1F78B4", G4 = "#1F78B4",
               Normal = "grey60", Unknown = "grey30")

sp_shapes <- c(Human = 16, Mouse = 15, Zebrafish = 17, Chicken = 8)

# Main PCA plot
p1 <- ggplot() +
  geom_point(data = pca_df %>% filter(species == "Human"),
             aes(x = PC1, y = PC2, color = subgroup),
             shape = 16, size = 1.2, alpha = 0.35) +
  geom_point(data = pca_df %>% filter(species != "Human"),
             aes(x = PC1, y = PC2, color = subgroup, shape = species),
             size = 3.5, alpha = 0.9, stroke = 0.8) +
  scale_color_manual(values = sg_colors, name = "Subgroup") +
  scale_shape_manual(values = sp_shapes, name = "Species") +
  labs(
    title = "Cross-Species Medulloblastoma PCA",
    subtitle = paste0("Animals projected onto human PCA | ", length(common_genes), " genes"),
    x = paste0("PC1 (", round(var_pct[1], 1), "%)"),
    y = paste0("PC2 (", round(var_pct[2], 1), "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave("figures/PCA_cross_species_MB.pdf", p1, width = 12, height = 9)
ggsave("figures/PCA_cross_species_MB.png", p1, width = 12, height = 9, dpi = 300)

# PC1 vs PC3
p2 <- ggplot() +
  geom_point(data = pca_df %>% filter(species == "Human"),
             aes(x = PC1, y = PC3, color = subgroup),
             shape = 16, size = 1.2, alpha = 0.35) +
  geom_point(data = pca_df %>% filter(species != "Human"),
             aes(x = PC1, y = PC3, color = subgroup, shape = species),
             size = 3.5, alpha = 0.9, stroke = 0.8) +
  scale_color_manual(values = sg_colors) +
  scale_shape_manual(values = sp_shapes) +
  labs(
    title = "Cross-Species Medulloblastoma PCA",
    x = paste0("PC1 (", round(var_pct[1], 1), "%)"),
    y = paste0("PC3 (", round(var_pct[3], 1), "%)")
  ) +
  theme_minimal(base_size = 14)

ggsave("figures/PCA_cross_species_PC1vPC3.pdf", p2, width = 12, height = 9)

# Save results
write.csv(pca_df, "results/pca_coordinates.csv", row.names = FALSE)


# =============================================================================
# CORRELATION WITH HUMAN CENTROIDS
# =============================================================================
message(">>> Computing correlations...")

# Combine all z-scored expression
expr_all_z <- cbind(
  (expr_human - human_means) / human_sds,
  (expr_mouse - human_means) / human_sds,
  (expr_zebrafish - human_means) / human_sds,
  (expr_chicken - human_means) / human_sds
)
expr_all_z[is.na(expr_all_z)] <- 0

human_sgs <- c("WNT", "SHH", "Group3", "Group4")
centroids <- sapply(human_sgs, function(sg) {
  ids <- pca_df$sample_id[pca_df$species == "Human" & pca_df$subgroup == sg]
  rowMeans(expr_all_z[, ids, drop = FALSE], na.rm = TRUE)
})

animal_ids <- pca_df$sample_id[pca_df$species != "Human"]
animal_expr <- expr_all_z[, animal_ids, drop = FALSE]
cor_mat <- cor(animal_expr, centroids, use = "pairwise.complete.obs")

anno_row <- pca_df %>%
  filter(species != "Human") %>%
  select(sample_id, species, subgroup) %>%
  as.data.frame()
rownames(anno_row) <- anno_row$sample_id
anno_row$sample_id <- NULL

pdf("figures/correlation_heatmap.pdf", width = 10, height = max(8, nrow(cor_mat) * 0.3))
pheatmap(cor_mat, cluster_cols = FALSE,
         annotation_row = anno_row,
         main = "Correlation with Human MB Subgroup Centroids",
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         breaks = seq(-1, 1, length.out = 101),
         fontsize_row = 7)
dev.off()

write.csv(cor_mat, "results/correlation_with_human_centroids.csv")


# =============================================================================
# SUMMARY
# =============================================================================
message("\n========================================")
message(">>> RESULTS SUMMARY")
message("========================================")

cat("\n=== PCA CENTROIDS ===\n")
centroids_pca <- pca_df %>%
  group_by(species, subgroup) %>%
  summarise(n = n(), PC1 = round(mean(PC1), 1), PC2 = round(mean(PC2), 1), .groups = 'drop')
print(centroids_pca)

cat("\n=== CORRELATION WITH HUMAN SUBGROUPS (mean by model) ===\n")
cor_summary <- pca_df %>%
  filter(species != "Human") %>%
  select(sample_id, species, subgroup) %>%
  left_join(rownames_to_column(as.data.frame(cor_mat), "sample_id")) %>%
  group_by(species, subgroup) %>%
  summarise(across(all_of(human_sgs), ~round(mean(.), 3)), .groups = 'drop')
print(cor_summary)

cat("\n=== VALIDATION ===\n")
# Check if expected clustering is observed
mouse_g3_cor <- mean(cor_mat[rownames(cor_mat) %in% pca_df$sample_id[pca_df$species == "Mouse" & pca_df$subgroup == "G3"], "Group3"])
mouse_shh_cor <- mean(cor_mat[rownames(cor_mat) %in% pca_df$sample_id[pca_df$species == "Mouse" & pca_df$subgroup == "SHH"], "SHH"])
zf_shh_cor <- mean(cor_mat[rownames(cor_mat) %in% pca_df$sample_id[pca_df$species == "Zebrafish"], "SHH"])
chicken_g3_cor <- mean(cor_mat[rownames(cor_mat) %in% pca_df$sample_id[pca_df$species == "Chicken"], "Group3"])

cat("Mouse G3 vs Human Group3:", round(mouse_g3_cor, 3), ifelse(mouse_g3_cor > 0.1, "OK", "LOW"), "\n")
cat("Mouse SHH vs Human SHH:", round(mouse_shh_cor, 3), ifelse(mouse_shh_cor > 0.1, "OK", "LOW"), "\n")
cat("Zebrafish SHH vs Human SHH:", round(zf_shh_cor, 3), ifelse(zf_shh_cor > 0.1, "OK", "LOW"), "\n")
cat("Chicken G3 vs Human Group3:", round(chicken_g3_cor, 3), ifelse(chicken_g3_cor > 0.1, "OK", "LOW"), "\n")

message("\n========================================")
message(">>> FILES SAVED")
message("========================================")
message("  figures/PCA_cross_species_MB.pdf/png")
message("  figures/PCA_cross_species_PC1vPC3.pdf")
message("  figures/correlation_heatmap.pdf")
message("  results/pca_coordinates.csv")
message("  results/correlation_with_human_centroids.csv")
message("========================================\n")

#!/usr/bin/env Rscript
# Random Forest in PCA space - unified cross-species classification
# Train RF on human PCA coordinates, apply to projected animal data

library(Biobase)
library(GEOquery)
library(dplyr)
library(tidyr)
library(ggplot2)
library(randomForest)

cat("=== RF IN PCA SPACE - CROSS-SPECIES CLASSIFICATION ===\n\n")

# ========== LOAD ALL DATA ==========
cat("Loading data...\n")

# Human
human_list <- readRDS("data/human/GSE85217_raw.rds")
expr_human <- exprs(human_list[[1]])
pdata_human <- pData(human_list[[1]])
subgroups <- pdata_human[["subgroup:ch1"]]

ortho_human <- readRDS("data/orthology/human_ensembl_to_symbol.rds")
symbol_map <- setNames(ortho_human$hgnc_symbol, ortho_human$ensembl_gene_id)
rownames(expr_human) <- symbol_map[gsub("_at$", "", rownames(expr_human))]
expr_human <- expr_human[!is.na(rownames(expr_human)) &
                           rownames(expr_human) != "", ]
expr_agg <- aggregate(expr_human,
                      by = list(Gene = rownames(expr_human)), FUN = mean)
rownames(expr_agg) <- expr_agg$Gene
expr_human <- as.matrix(expr_agg[, -1])

# Zebrafish
zf_counts <- read.table(
  "data/zebrafish/GSE242897/GSE242897_19521raw_counts.txt.gz",
  header = TRUE, row.names = 1, sep = "\t")
zmb_cols <- grep("PTCH1_P53MUT", colnames(zf_counts), value = TRUE)
zf_expr <- log2(as.matrix(zf_counts[, zmb_cols]) + 1)
ortho_zf <- readRDS("data/orthology/zebrafish_to_human.rds")
zf_to_human <- setNames(ortho_zf$hsapiens_homolog_associated_gene_name,
                        ortho_zf$ensembl_gene_id)
rownames(zf_expr) <- zf_to_human[rownames(zf_expr)]
zf_expr <- zf_expr[!is.na(rownames(zf_expr)) & rownames(zf_expr) != "", ]
zf_agg <- aggregate(zf_expr, by = list(Gene = rownames(zf_expr)), FUN = mean)
rownames(zf_agg) <- zf_agg$Gene
zf_expr <- as.matrix(zf_agg[, -1])

# Chicken
chicken <- read.table("data/chicken/BULK_Chicken_expression.tsv",
                      header = TRUE, row.names = 1, sep = "\t")
chicken_expr <- log2(as.matrix(chicken) + 1)
rownames(chicken_expr) <- toupper(rownames(chicken_expr))

# Mouse
gse65888 <- getGEO(filename = "data/mouse/GSE65888-GPL1261_series_matrix.txt.gz",
                   getGPL = FALSE)
expr_mouse <- exprs(gse65888)
pdata_mouse <- pData(gse65888)

gpl <- getGEO(filename = "data/mouse/GPL1261.soft.gz")
gpl_table <- Table(gpl)
probe_to_symbol <- setNames(gpl_table[["Gene Symbol"]], gpl_table$ID)
rownames(expr_mouse) <- probe_to_symbol[rownames(expr_mouse)]
expr_mouse <- expr_mouse[!is.na(rownames(expr_mouse)) &
                           rownames(expr_mouse) != "", ]

ortho_mouse <- readRDS("data/orthology/mouse_to_human.rds")
mouse_to_human <- setNames(ortho_mouse$human_symbol, ortho_mouse$mouse_symbol)
rownames(expr_mouse) <- mouse_to_human[rownames(expr_mouse)]
expr_mouse <- expr_mouse[!is.na(rownames(expr_mouse)) &
                           rownames(expr_mouse) != "", ]
mouse_agg <- aggregate(expr_mouse,
                       by = list(Gene = rownames(expr_mouse)), FUN = mean)
rownames(mouse_agg) <- mouse_agg$Gene
expr_mouse <- as.matrix(mouse_agg[, -1])

mouse_info <- data.frame(
  sample_id = colnames(expr_mouse),
  title = pdata_mouse$title,
  stringsAsFactors = FALSE
)
mouse_info$model <- case_when(
  grepl("Ptch", mouse_info$title) &
    grepl("Cdkn2c", mouse_info$title) ~ "SHH_Ptch",
  grepl("Nestin-Cre", mouse_info$title) &
    grepl("p53", mouse_info$title) ~ "SHH_Nestin",
  grepl("Ctnnb1", mouse_info$title, ignore.case = TRUE) ~ "WNT",
  grepl("Blbp-MYC|Gad2-MYC|Ptf1a-MYC", mouse_info$title,
        ignore.case = TRUE) ~ "G3_IUE",
  grepl("Atoh1-CreER|Prom1-CreER|MycRFP|Myc-RFP", mouse_info$title,
        ignore.case = TRUE) ~ "G3_Ortho",
  TRUE ~ NA_character_
)
mouse_info <- mouse_info[!is.na(mouse_info$model), ]
expr_mouse <- expr_mouse[, mouse_info$sample_id]

cat("  Human:", ncol(expr_human), "| Zebrafish:", ncol(zf_expr),
    "| Chicken:", ncol(chicken_expr), "| Mouse:", ncol(expr_mouse), "\n")

# ========== SELECT GENES AND NORMALIZE ==========
cat("\nSelecting genes and normalizing...\n")

# Top variable genes
N_GENES <- 2500
gene_vars <- apply(expr_human, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:N_GENES])

# Find common genes across all species
common_genes <- Reduce(intersect, list(
  top_genes,
  rownames(zf_expr),
  rownames(chicken_expr),
  rownames(expr_mouse)
))
cat("  Common genes:", length(common_genes), "\n")

# Z-score normalize using human parameters (per gene)
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
X_mouse <- normalize_to_human(expr_mouse, common_genes, human_means, human_sds)

# ========== PCA ON HUMAN ==========
cat("\nRunning PCA on human data...\n")
pca <- prcomp(X_human, center = FALSE, scale. = FALSE)
var_exp <- summary(pca)$importance[2, 1:10] * 100
cat("  Variance explained (PC1-5):",
    paste(round(var_exp[1:5], 1), collapse = ", "), "%\n")

# Use top 20 PCs for classification
N_PCS <- 20
human_pcs <- pca$x[, 1:N_PCS]

# Project animal data
zf_pcs <- X_zf %*% pca$rotation[, 1:N_PCS]
ch_pcs <- X_ch %*% pca$rotation[, 1:N_PCS]
mouse_pcs <- X_mouse %*% pca$rotation[, 1:N_PCS]

# ========== TRAIN/TEST SPLIT ==========
cat("\nSplitting human data...\n")
set.seed(42)
train_idx <- sample(seq_len(nrow(human_pcs)), size = round(0.8 * nrow(human_pcs)))
test_idx <- setdiff(seq_len(nrow(human_pcs)), train_idx)

X_train <- human_pcs[train_idx, ]
y_train <- factor(subgroups[train_idx])
X_test <- human_pcs[test_idx, ]
y_test <- factor(subgroups[test_idx])

cat("  Train:", nrow(X_train), "| Test:", nrow(X_test), "\n")

# ========== TRAIN RANDOM FOREST ==========
cat("\nTraining Random Forest on PCA coordinates...\n")
set.seed(42)
rf_model <- randomForest(
  x = X_train,
  y = y_train,
  ntree = 1000,
  importance = TRUE
)

cat("  OOB error rate:", round(rf_model$err.rate[1000, "OOB"] * 100, 1), "%\n")

# ========== EVALUATE ==========
cat("\nEvaluating...\n")

# Human test set
pred_test <- predict(rf_model, X_test)
probs_test <- predict(rf_model, X_test, type = "prob")
human_acc <- mean(pred_test == y_test) * 100
cat("  Human test accuracy:", round(human_acc, 1), "%\n")

# Zebrafish
pred_zf <- predict(rf_model, zf_pcs)
probs_zf <- predict(rf_model, zf_pcs, type = "prob")
cat("  Zebrafish predictions:", as.character(pred_zf), "\n")

# Chicken
pred_ch <- predict(rf_model, ch_pcs)
probs_ch <- predict(rf_model, ch_pcs, type = "prob")
cat("  Chicken predictions:", as.character(pred_ch), "\n")

# Mouse
pred_mouse <- predict(rf_model, mouse_pcs)
probs_mouse <- predict(rf_model, mouse_pcs, type = "prob")

# ========== COMPILE RESULTS ==========
all_results <- rbind(
  data.frame(
    sample_id = rownames(X_test),
    species = "Human", model = as.character(y_test),
    actual = as.character(y_test), predicted = as.character(pred_test),
    probs_test, stringsAsFactors = FALSE
  ),
  data.frame(
    sample_id = rownames(zf_pcs),
    species = "Zebrafish", model = "ZMB",
    actual = "SHH", predicted = as.character(pred_zf),
    probs_zf, stringsAsFactors = FALSE
  ),
  data.frame(
    sample_id = rownames(ch_pcs),
    species = "Chicken",
    model = ifelse(grepl("parenchyma", rownames(ch_pcs)), "Xenograft", "Culture"),
    actual = "Group3", predicted = as.character(pred_ch),
    probs_ch, stringsAsFactors = FALSE
  ),
  data.frame(
    sample_id = rownames(mouse_pcs),
    species = "Mouse", model = mouse_info$model,
    actual = case_when(
      grepl("SHH", mouse_info$model) ~ "SHH",
      grepl("WNT", mouse_info$model) ~ "WNT",
      grepl("G3", mouse_info$model) ~ "Group3"
    ),
    predicted = as.character(pred_mouse),
    probs_mouse, stringsAsFactors = FALSE
  )
)

# ========== ACCURACY ==========
cat("\n=== RF-PCA CLASSIFICATION ACCURACY ===\n")
accuracy <- all_results %>%
  group_by(species, model) %>%
  summarise(
    n = n(),
    correct = sum(predicted == actual),
    accuracy = round(mean(predicted == actual) * 100, 1),
    .groups = "drop"
  )
print(as.data.frame(accuracy), row.names = FALSE)

cat("\n=== BY SPECIES ===\n")
species_acc <- all_results %>%
  group_by(species) %>%
  summarise(
    n = n(), correct = sum(predicted == actual),
    accuracy = round(mean(predicted == actual) * 100, 1),
    .groups = "drop"
  )
print(as.data.frame(species_acc), row.names = FALSE)

# ========== CREATE FIGURE (JITTER) ==========
cat("\nCreating jitter figure...\n")

plot_data <- all_results %>%
  select(sample_id, species, model, actual, predicted, WNT, SHH, Group3, Group4) %>%
  pivot_longer(cols = c(WNT, SHH, Group3, Group4),
               names_to = "Subgroup", values_to = "Probability")

plot_data$Subgroup <- factor(plot_data$Subgroup,
                              levels = c("WNT", "SHH", "Group3", "Group4"))

# Add study reference labels
plot_data$label <- case_when(
  plot_data$species == "Human" & plot_data$model == "WNT" ~
    "Cavalli et al.\nCancer Cell 2017\nWNT",
  plot_data$species == "Human" & plot_data$model == "SHH" ~
    "Cavalli et al.\nCancer Cell 2017\nSHH",
  plot_data$species == "Human" & plot_data$model == "Group3" ~
    "Cavalli et al.\nCancer Cell 2017\nGroup3",
  plot_data$species == "Human" & plot_data$model == "Group4" ~
    "Cavalli et al.\nCancer Cell 2017\nGroup4",
  plot_data$species == "Zebrafish" ~
    "Casey et al.\nCell Rep 2024\nZMB",
  plot_data$species == "Chicken" & plot_data$model == "Culture" ~
    "This study\nCulture",
  plot_data$species == "Chicken" & plot_data$model == "Xenograft" ~
    "This study\nXenograft",
  plot_data$species == "Mouse" & plot_data$model == "WNT" ~
    "Gibson et al.\nNature 2010\nWNT",
  plot_data$species == "Mouse" & plot_data$model == "SHH_Ptch" ~
    "Uziel et al.\nGenes Dev 2005\nSHH_Ptch",
  plot_data$species == "Mouse" & plot_data$model == "SHH_Nestin" ~
    "Kawauchi et al.\nOncogene 2017\nSHH_Nestin",
  plot_data$species == "Mouse" & plot_data$model == "G3_IUE" ~
    "Kawauchi et al.\nOncogene 2017\nG3_IUE",
  plot_data$species == "Mouse" & plot_data$model == "G3_Ortho" ~
    "Kawauchi et al.\nCancer Cell 2012\nG3_Ortho",
  TRUE ~ NA_character_
)

order_vec <- c(
  "Cavalli et al.\nCancer Cell 2017\nWNT",
  "Cavalli et al.\nCancer Cell 2017\nSHH",
  "Cavalli et al.\nCancer Cell 2017\nGroup3",
  "Cavalli et al.\nCancer Cell 2017\nGroup4",
  "Casey et al.\nCell Rep 2024\nZMB",
  "This study\nCulture",
  "This study\nXenograft",
  "Gibson et al.\nNature 2010\nWNT",
  "Uziel et al.\nGenes Dev 2005\nSHH_Ptch",
  "Kawauchi et al.\nOncogene 2017\nSHH_Nestin",
  "Kawauchi et al.\nOncogene 2017\nG3_IUE",
  "Kawauchi et al.\nCancer Cell 2012\nG3_Ortho"
)
plot_data$label <- factor(plot_data$label, levels = order_vec)
plot_data <- plot_data[!is.na(plot_data$label), ]

colors <- c("WNT" = "#3498db", "SHH" = "#e74c3c",
            "Group3" = "#2ecc71", "Group4" = "#9b59b6")

# Boxplot with transparent points overlay
p <- ggplot(plot_data, aes(x = label, y = Probability, fill = Subgroup)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7,
               alpha = 0.7, outlier.shape = NA) +
  geom_point(aes(color = Subgroup),
             position = position_jitterdodge(jitter.width = 0.1,
                                              jitter.height = 0.02,
                                              dodge.width = 0.8),
             size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors, guide = "none") +
  geom_vline(xintercept = c(4.5, 5.5, 7.5), linetype = "dashed", alpha = 0.3) +
  labs(title = "Cross-Species MB Classification (RF in PCA Space)",
       subtitle = sprintf("RF trained on %d PCs from %d genes",
                          N_PCS, length(common_genes)),
       x = "", y = "RF Probability") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        plot.title = element_text(face = "bold"),
        legend.position = "right")

ggsave("figures/cross_species_RF_PCA.png", p,
       width = 16, height = 7, dpi = 150)
cat("Saved: figures/cross_species_RF_PCA.png\n")

# Save
write.csv(all_results, "results/rf_pca_classification_results.csv",
          row.names = FALSE)
cat("Saved: results/rf_pca_classification_results.csv\n")

cat("\n=== DONE ===\n")

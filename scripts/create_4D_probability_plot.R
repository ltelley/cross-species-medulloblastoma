#!/usr/bin/env Rscript
# 4D probability plot with confidence ellipses
# Each dimension = probability of being WNT, SHH, Group3, Group4

library(dplyr)
library(ggplot2)
library(patchwork)

cat("=== 4D PROBABILITY PLOT ===\n\n")

# Load results
results <- read.csv("results/rf_pca_classification_results.csv")

# Create condition labels
results$condition <- paste0(results$species, "_", results$model)

# Simplify condition names for legend
results$condition <- case_when(
  results$species == "Human" ~ paste0("Human_", results$model),
  results$species == "Zebrafish" ~ "Zebrafish_ZMB",
  results$species == "Chicken" ~ paste0("Chicken_", results$model),
  results$species == "Mouse" ~ paste0("Mouse_", results$model),
  TRUE ~ results$condition
)

cat("Conditions:\n")
print(table(results$condition))

# Define colors per condition
condition_colors <- c(
  # Human subgroups
  "Human_WNT" = "#3498db",
  "Human_SHH" = "#e74c3c",
  "Human_Group3" = "#2ecc71",
  "Human_Group4" = "#9b59b6",
  # Zebrafish
  "Zebrafish_ZMB" = "#e67e22",
  # Chicken
  "Chicken_Culture" = "#1abc9c",
  "Chicken_Xenograft" = "#16a085",
  # Mouse
  "Mouse_WNT" = "#5dade2",
  "Mouse_SHH_Ptch" = "#ec7063",
  "Mouse_SHH_Nestin" = "#f1948a",
  "Mouse_G3_IUE" = "#58d68d",
  "Mouse_G3_Ortho" = "#82e0aa"
)

# Create a 2x2 grid with specific pairs and ellipses
cat("\nCreating 2x2 grid with ellipses...\n")

# Function to create scatter with ellipse
make_scatter_ellipse <- function(data, x_var, y_var, title = "", show_legend = FALSE) {
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]],
                        color = condition, fill = condition)) +
    geom_point(alpha = 0.6, size = 2) +
    stat_ellipse(level = 0.8, geom = "polygon", alpha = 0.1) +
    scale_color_manual(values = condition_colors) +
    scale_fill_manual(values = condition_colors) +
    labs(title = title, x = paste0("P(", x_var, ")"), y = paste0("P(", y_var, ")")) +
    theme_bw()

  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(legend.position = "right",
                   legend.text = element_text(size = 7),
                   legend.title = element_text(size = 9))
  }
  return(p)
}

# Create 4 key plots (last one with legend)
p1 <- make_scatter_ellipse(results, "SHH", "Group3", "SHH vs Group3")
p2 <- make_scatter_ellipse(results, "SHH", "Group4", "SHH vs Group4")
p3 <- make_scatter_ellipse(results, "Group3", "Group4", "Group3 vs Group4")
p4 <- make_scatter_ellipse(results, "WNT", "SHH", "WNT vs SHH", show_legend = TRUE)

combined <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Cross-Species MB Classification in 4D Probability Space",
    subtitle = "RF probabilities with 80% confidence ellipses"
  )

ggsave("figures/cross_species_4D_grid.png", combined, width = 14, height = 10, dpi = 150)
cat("Saved: figures/cross_species_4D_grid.png\n")

# Create a single comprehensive plot: SHH vs Group3 (the most informative)
cat("\nCreating main SHH vs Group3 plot...\n")

p_main <- ggplot(results, aes(x = SHH, y = Group3, color = condition, fill = condition)) +
  geom_point(alpha = 0.7, size = 2.5) +
  stat_ellipse(level = 0.8, geom = "polygon", alpha = 0.15, linewidth = 0.8) +
  scale_color_manual(values = condition_colors, name = "Condition") +
  scale_fill_manual(values = condition_colors, name = "Condition") +
  labs(
    title = "Cross-Species MB Classification",
    subtitle = "RF probability space with 80% confidence ellipses",
    x = "P(SHH)",
    y = "P(Group3)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  ) +
  coord_fixed()

ggsave("figures/cross_species_SHH_vs_G3.png", p_main, width = 10, height = 8, dpi = 150)
cat("Saved: figures/cross_species_SHH_vs_G3.png\n")

cat("\n=== DONE ===\n")

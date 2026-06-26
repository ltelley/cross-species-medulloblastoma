#!/usr/bin/env Rscript
# Radar/Spider plot - 4D visualization of RF probabilities
# Each axis = probability of WNT, SHH, Group3, Group4

library(dplyr)
library(tidyr)
library(ggplot2)

cat("=== RADAR PLOT - 4D PROBABILITY VISUALIZATION ===\n\n")

# Load results
results <- read.csv("results/rf_pca_classification_results.csv")

# Create condition labels
results$condition <- case_when(
  results$species == "Human" ~ paste0("Human_", results$model),
  results$species == "Zebrafish" ~ "Zebrafish_ZMB",
  results$species == "Chicken" ~ paste0("Chicken_", results$model),
  results$species == "Mouse" ~ paste0("Mouse_", results$model),
  TRUE ~ "Other"
)

# Use individual samples (not mean per condition)
# Add sample_id to make each sample unique
results$sample_label <- paste0(results$condition, "_", seq_len(nrow(results)))

cat("Samples per condition:\n")
print(table(results$condition))

# Reshape for radar plot - each sample is a separate polygon
radar_data <- results %>%
  select(sample_id, condition, WNT, SHH, Group3, Group4) %>%
  pivot_longer(cols = c(WNT, SHH, Group3, Group4),
               names_to = "Subgroup",
               values_to = "Probability")

# Set factor order for the axes (clockwise: WNT, SHH, G3, G4)
radar_data$Subgroup <- factor(radar_data$Subgroup,
                               levels = c("WNT", "SHH", "Group3", "Group4"))

# Define colors
condition_colors <- c(
  "Human_WNT" = "#3498db",
  "Human_SHH" = "#e74c3c",
  "Human_Group3" = "#2ecc71",
  "Human_Group4" = "#9b59b6",
  "Zebrafish_ZMB" = "#e67e22",
  "Chicken_Culture" = "#1abc9c",
  "Chicken_Xenograft" = "#16a085",
  "Mouse_WNT" = "#5dade2",
  "Mouse_SHH_Ptch" = "#ec7063",
  "Mouse_SHH_Nestin" = "#f1948a",
  "Mouse_G3_IUE" = "#58d68d",
  "Mouse_G3_Ortho" = "#82e0aa"
)

# Order conditions logically
condition_order <- c(
  "Human_WNT", "Human_SHH", "Human_Group3", "Human_Group4",
  "Zebrafish_ZMB",
  "Chicken_Culture", "Chicken_Xenograft",
  "Mouse_WNT", "Mouse_SHH_Ptch", "Mouse_SHH_Nestin",
  "Mouse_G3_IUE", "Mouse_G3_Ortho"
)
radar_data$condition <- factor(radar_data$condition, levels = condition_order)

# Create faceted radar plot - one per condition
cat("\nCreating faceted radar plot...\n")

p_radar <- ggplot(radar_data, aes(x = Subgroup, y = Probability,
                                   color = condition)) +
  geom_point(size = 3, alpha = 0.7) +
  coord_polar() +
  scale_color_manual(values = condition_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
  facet_wrap(~condition, ncol = 4) +
  labs(title = "Cross-Species MB Classification - 4D Probability Space",
       subtitle = "RF probability per condition") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey85")
  )

ggsave("figures/cross_species_radar.png", p_radar, width = 14, height = 10, dpi = 150)
cat("Saved: figures/cross_species_radar.png\n")

# Alternative: Faceted radar by species
cat("\nCreating faceted radar plot...\n")

# Add species info
radar_data$species <- case_when(
  grepl("Human", radar_data$condition) ~ "Human",
  grepl("Zebrafish", radar_data$condition) ~ "Zebrafish",
  grepl("Chicken", radar_data$condition) ~ "Chicken",
  grepl("Mouse", radar_data$condition) ~ "Mouse"
)

# Filter to animal models only for cleaner view
animal_data <- radar_data %>% filter(species != "Human")

p_animal <- ggplot(animal_data, aes(x = Subgroup, y = Probability,
                                     group = condition, color = condition)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  coord_polar() +
  scale_color_manual(values = condition_colors, name = "Model") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  facet_wrap(~species, ncol = 3) +
  labs(title = "Animal Model MB Classification - 4D Probability Space",
       subtitle = "Radar plot: mean RF probability per model") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 10, face = "bold"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey80")
  )

ggsave("figures/cross_species_radar_animals.png", p_animal, width = 14, height = 6, dpi = 150)
cat("Saved: figures/cross_species_radar_animals.png\n")

cat("\n=== DONE ===\n")

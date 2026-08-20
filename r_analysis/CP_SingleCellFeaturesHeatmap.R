# Heatmap of all morphology features in single_cell_features_and_clusters.csv
# Image-level 10–90% trimmed means; same colors / Time annotation as the
# main all-times analysis. Does not modify existing notebooks.

library(tidyverse)
library(pheatmap)
source("astro_analysis_helpers.R")

DATA_PATH <- "C:/Users/talku/OneDrive/Habib Lab/CellProfiler_Results/single_cell_features_and_clusters.csv"
OUTPUT_DIR <- "C:/Users/talku/OneDrive/Documents/HabibLab/output/single_cell_features/heatmaps"
TIME_ORDER <- c("ctrl", "4h", "24h", "72h", "7d")
TRIM_LOWER <- 0.10
TRIM_UPPER <- 0.90

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cells <- read_csv(DATA_PATH, show_col_types = FALSE)

meta_cols <- c(
  "sample", "original_cell_id", "compact_cell_id",
  "measurement_backend", "cluster", "crop_size_pixels"
)
feature_cols <- setdiff(names(cells), meta_cols)
feature_cols <- feature_cols[vapply(cells[feature_cols], is.numeric, logical(1))]

cells <- cells %>%
  mutate(
    Time = factor(
      str_extract(sample, "^(ctrl|4h|24h|72h|7d)"),
      levels = TIME_ORDER
    ),
    PhotoID = str_extract(sample, "[0-9]+$"),
    ImageNumber = as.integer(factor(sample))
  )

if (anyNA(cells$Time)) {
  stop(
    "Could not parse Time from sample names: ",
    paste(unique(cells$sample[is.na(cells$Time)]), collapse = ", ")
  )
}

image_means <- summarise_trimmed_means(
  data = cells,
  features = feature_cols,
  lower = TRIM_LOWER,
  upper = TRIM_UPPER
) %>%
  left_join(
    cells %>%
      group_by(ImageNumber) %>%
      summarise(
        PhotoID = first(PhotoID),
        Time = first(Time),
        .groups = "drop"
      ),
    by = "ImageNumber"
  ) %>%
  arrange(Time, PhotoID) %>%
  mutate(ImageLabel = paste0("Photo", PhotoID, "_", Time))

mean_matrix <- image_means %>%
  select(all_of(feature_cols)) %>%
  as.matrix()
rownames(mean_matrix) <- image_means$ImageLabel

TIME_COLORS <- get_time_colors(levels(image_means$Time))
row_annotation <- data.frame(
  Time = image_means$Time,
  row.names = rownames(mean_matrix)
)
ann_colors <- list(Time = TIME_COLORS[levels(row_annotation$Time)])

scaled <- scale_features(mean_matrix)

plot_image_heatmap(
  scaled_mat = scaled,
  row_ann = row_annotation,
  ann_colors = ann_colors,
  title = paste0(
    "All features — single_cell_features_and_clusters\n",
    "image-level ", TRIM_LOWER * 100, "–", TRIM_UPPER * 100,
    "% trimmed means, z-scored per feature"
  ),
  file = file.path(OUTPUT_DIR, "all_features.png"),
  fontsize_row = 10,
  fontsize_col = 10
)

cat("Images:", nrow(image_means), "\n")
cat("Times:\n")
print(table(image_means$Time))
cat("Features:", paste(feature_cols, collapse = ", "), "\n")
cat("Saved:", file.path(OUTPUT_DIR, "all_features.png"), "\n")

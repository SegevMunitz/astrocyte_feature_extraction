# Intensity-only boxplots from AstroResultsAstrocytes_with_nuclei.csv
# Image-level 10–90% trimmed means per Time (same convention as morphology notebooks).

library(tidyverse)
source("astro_analysis_helpers.R")

DATA_PATH <- "C:/Users/talku/OneDrive/Habib Lab/CellProfiler_Results/AstroResultsAstrocytes_with_nuclei.csv"
OUTPUT_DIR <- "C:/Users/talku/OneDrive/Documents/HabibLab/output/intensity_boxplots"

TIME_ORDER <- c("ctrl", "4h", "24h", "72h", "7d")
TRIM_LOWER <- 0.10
TRIM_UPPER <- 0.90
FEATURES_PER_PAGE <- 8

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cells_data <- load_astro_cells(DATA_PATH)
intensity_cols <- names(cells_data)[
  str_detect(names(cells_data), "^Intensity_") &
    !str_detect(names(cells_data), "Center_[XY]$")
]
intensity_cols <- intensity_cols[
  vapply(cells_data[intensity_cols], is.numeric, logical(1))
]
if (length(intensity_cols) == 0) {
  stop("No Intensity_* columns found in: ", DATA_PATH)
}

image_metadata <- build_image_metadata(cells_data, TIME_ORDER)
TIME_COLORS <- get_time_colors(levels(image_metadata$Time))

image_means <- build_image_means(
  cells_data = cells_data,
  cell_features = intensity_cols,
  skeleton_features = character(),
  image_metadata = image_metadata,
  use_trimmed = TRUE,
  trim_lower = TRIM_LOWER,
  trim_upper = TRIM_UPPER
)

# Rank by Kruskal–Wallis effect size so the most time-varying intensity features come first.
intensity_stats <- compute_time_feature_stats(
  image_means = image_means,
  feature_cols = intensity_cols,
  use_trimmed = TRUE,
  trim_lower = TRIM_LOWER,
  trim_upper = TRIM_UPPER
) %>%
  arrange(desc(epsilon_squared))

write_csv(intensity_stats, file.path(OUTPUT_DIR, "intensity_time_stats.csv"))

ranked_features <- intensity_stats$Feature
subtitle <- paste0(
  "Each dot = one image; 10–90% trimmed mean per image; ordered by Kruskal–Wallis ε²"
)

pages <- split(
  ranked_features,
  ceiling(seq_along(ranked_features) / FEATURES_PER_PAGE)
)

for (page_index in seq_along(pages)) {
  page_features <- pages[[page_index]]
  start_i <- (page_index - 1) * FEATURES_PER_PAGE + 1
  end_i <- start_i + length(page_features) - 1
  outfile <- file.path(
    OUTPUT_DIR,
    sprintf("intensity_boxplots_%02d_%02d.png", start_i, end_i)
  )
  plot_time_boxplots(
    image_means = image_means,
    top_features = page_features,
    time_colors = TIME_COLORS,
    title = sprintf("Intensity features %d–%d", start_i, end_i),
    subtitle = subtitle,
    file = outfile,
    width = 12,
    height = max(6, 2.2 * ceiling(length(page_features) / 2))
  )
  message("Wrote ", outfile)
}

# Dedicated panel for IntegratedIntensity features (log10), usually the most interpretable.
integrated <- ranked_features[str_detect(ranked_features, "IntegratedIntensity")]
if (length(integrated) > 0) {
  log_means <- image_means
  for (feature in integrated) {
    x <- log_means[[feature]]
    if (any(x <= 0, na.rm = TRUE)) {
      offset <- abs(min(x, na.rm = TRUE)) + 1
      log_means[[feature]] <- log10(x + offset)
    } else {
      log_means[[feature]] <- log10(x)
    }
  }
  outfile <- file.path(OUTPUT_DIR, "intensity_boxplots_integrated_log10.png")
  plot_time_boxplots(
    image_means = log_means,
    top_features = integrated,
    time_colors = TIME_COLORS,
    title = "Integrated intensity (log10)",
    subtitle = "Image-level 10–90% trimmed means",
    file = outfile,
    width = 12,
    height = max(5, 2.4 * ceiling(length(integrated) / 2))
  )
  message("Wrote ", outfile)
}

message("Intensity columns: ", length(intensity_cols))
message("Images: ", nrow(image_means))
message("Outputs: ", OUTPUT_DIR)

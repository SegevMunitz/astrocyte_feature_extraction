# Shared helpers for astrocyte morphology analysis notebooks.

coalesce_cols <- function(df, pattern) {
  cols <- names(df)[stringr::str_detect(names(df), pattern)]
  if (length(cols) == 0) {
    return(rep(NA_character_, nrow(df)))
  }
  purrr::reduce(df[cols], ~ dplyr::coalesce(as.character(.x), as.character(.y)))
}

short_feature_names <- function(cols) {
  cols %>%
    stringr::str_replace("^AreaShape_", "") %>%
    stringr::str_replace("^Intensity_", "Int_") %>%
    stringr::str_replace("^ObjectSkeleton_", "Skel_")
}

make_output_folders <- function(root, subfolders) {
  paths <- stats::setNames(file.path(root, subfolders), subfolders)
  invisible(lapply(paths, dir.create, showWarnings = FALSE, recursive = TRUE))
  as.list(paths)
}

get_time_colors <- function(times) {
  palette <- c(
    "ctrl" = "#999999",
    "4h" = "#E41A1C",
    "24h" = "#377EB8",
    "72h" = "#4DAF4A",
    "7d" = "#984EA3"
  )
  missing <- setdiff(times, names(palette))
  if (length(missing) > 0) {
    extra <- stats::setNames(
      grDevices::colorRampPalette(c("#AAAAAA", "#555555"))(length(missing)),
      missing
    )
    palette <- c(palette, extra)
  }
  palette[times]
}

HEATMAP_COLORS <- grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
COSINE_COLORS <- HEATMAP_COLORS

scale_features <- function(mat) {
  scaled <- scale(mat)
  scaled[is.nan(scaled)] <- 0
  colnames(scaled) <- short_feature_names(colnames(scaled))
  scaled
}

cosine_similarity_matrix <- function(mat) {
  X <- scale(mat)
  X[is.nan(X)] <- 0
  norms <- sqrt(colSums(X^2))
  norms[norms == 0] <- 1
  X_norm <- sweep(X, 2, norms, "/")
  sim <- t(X_norm) %*% X_norm
  rownames(sim) <- colnames(sim) <- short_feature_names(colnames(mat))
  sim
}

summarise_trimmed_means <- function(data, features, lower = 0.10, upper = 0.90) {
  data %>%
    dplyr::group_by(ImageNumber) %>%
    dplyr::group_modify(function(grp, ...) {
      feat_means <- purrr::map_dbl(features, function(feat) {
        x <- grp[[feat]]
        qs <- stats::quantile(x, c(lower, upper), na.rm = TRUE)
        x_trim <- x[x >= qs[1] & x <= qs[2]]
        if (length(x_trim) == 0) NA_real_ else mean(x_trim, na.rm = TRUE)
      })
      names(feat_means) <- features

      cell_counts <- purrr::map_int(features, function(feat) {
        x <- grp[[feat]]
        qs <- stats::quantile(x, c(lower, upper), na.rm = TRUE)
        sum(x >= qs[1] & x <= qs[2], na.rm = TRUE)
      })

      tibble::tibble(
        !!!feat_means,
        CellCount_trimmed = min(cell_counts)
      )
    }) %>%
    dplyr::ungroup()
}

METADATA_COLS <- c(
  "ImageNumber", "ObjectNumber", "PhotoID", "Time",
  "CellCount", "CellCount_trimmed", "ImageLabel"
)

feature_matrix <- function(image_means, feature_cols) {
  keep <- setdiff(feature_cols, METADATA_COLS)
  keep <- intersect(keep, names(image_means))
  mat <- image_means %>%
    dplyr::ungroup() %>%
    dplyr::select(dplyr::all_of(keep)) %>%
    as.matrix()
  rownames(mat) <- image_means$ImageLabel
  if (any(colnames(mat) %in% METADATA_COLS)) {
    stop("Metadata columns leaked into the feature matrix: ",
         paste(intersect(colnames(mat), METADATA_COLS), collapse = ", "))
  }
  mat
}

log_transform_intensity_columns <- function(mat, pattern = "IntegratedIntensity") {
  cols <- colnames(mat)[stringr::str_detect(colnames(mat), pattern)]
  for (col in cols) {
    x <- mat[, col]
    if (any(x <= 0, na.rm = TRUE)) {
      offset <- abs(min(x, na.rm = TRUE)) + 1
      mat[, col] <- log10(x + offset)
    } else {
      mat[, col] <- log10(x)
    }
    colnames(mat)[colnames(mat) == col] <- paste0(col, "_log10")
  }
  mat
}

summarise_image_skeleton <- function(data, skeleton_features) {
  if (length(skeleton_features) == 0) {
    return(tibble::tibble(ImageNumber = integer()))
  }

  data %>%
    dplyr::group_by(ImageNumber) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(skeleton_features),
        ~ mean(.x, na.rm = TRUE),
        .names = "{.col}"
      ),
      .groups = "drop"
    )
}

load_astro_cells <- function(path) {
  cells_data <- readr::read_csv(path, skip = 1, show_col_types = FALSE)
  cells_data$PhotoID <- coalesce_cols(cells_data, "^Metadata_PhotoID")
  cells_data$Time <- coalesce_cols(cells_data, "^Metadata_Time")
  cells_data
}

identify_feature_columns <- function(cells_data) {
  cell_feature_cols <- names(cells_data)[
    stringr::str_detect(names(cells_data), "^(AreaShape_|Intensity_)") &
      !stringr::str_detect(names(cells_data), "Center_[XY]$")
  ]
  cell_feature_cols <- cell_feature_cols[
    vapply(cells_data[cell_feature_cols], is.numeric, logical(1))
  ]

  skeleton_feature_cols <- names(cells_data)[
    stringr::str_detect(names(cells_data), "^ObjectSkeleton_")
  ]
  skeleton_feature_cols <- skeleton_feature_cols[
    vapply(cells_data[skeleton_feature_cols], is.numeric, logical(1))
  ]

  feature_cols <- setdiff(
    c(cell_feature_cols, skeleton_feature_cols),
    METADATA_COLS
  )
  feature_cols_core <- cell_feature_cols[
    !stringr::str_detect(cell_feature_cols, "Zernike") &
      !stringr::str_detect(cell_feature_cols, "HuMoment") &
      !stringr::str_detect(cell_feature_cols, "InertiaTensor")
  ]
  feature_cols_core <- setdiff(
    c(feature_cols_core, skeleton_feature_cols),
    METADATA_COLS
  )

  list(
    cell = cell_feature_cols,
    skeleton = skeleton_feature_cols,
    all = feature_cols,
    core = feature_cols_core
  )
}

build_image_metadata <- function(cells_data, time_order) {
  extra_times <- setdiff(unique(cells_data$Time), time_order)
  full_order <- c(time_order, extra_times)

  cells_data %>%
    dplyr::group_by(ImageNumber) %>%
    dplyr::summarise(
      PhotoID = dplyr::first(PhotoID),
      Time = dplyr::first(Time),
      CellCount = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Time = factor(Time, levels = intersect(full_order, unique(Time))),
      ImageLabel = paste0("Photo", PhotoID, "_", Time)
    )
}

build_image_means <- function(
  cells_data,
  cell_features,
  skeleton_features,
  image_metadata,
  use_trimmed = TRUE,
  trim_lower = 0.10,
  trim_upper = 0.90
) {
  if (use_trimmed) {
    cell_means <- summarise_trimmed_means(
      data = cells_data,
      features = cell_features,
      lower = trim_lower,
      upper = trim_upper
    )
  } else {
    cell_means <- cells_data %>%
      dplyr::group_by(ImageNumber) %>%
      dplyr::summarise(
        dplyr::across(dplyr::all_of(cell_features), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      )
  }

  skeleton_means <- summarise_image_skeleton(cells_data, skeleton_features)

  image_means <- cell_means %>%
    dplyr::ungroup() %>%
    dplyr::left_join(skeleton_means, by = "ImageNumber") %>%
    dplyr::left_join(image_metadata, by = "ImageNumber") %>%
    dplyr::arrange(Time, PhotoID) %>%
    dplyr::ungroup()

  image_means
}

plot_image_heatmap <- function(
  scaled_mat,
  row_ann,
  ann_colors,
  title,
  file,
  fontsize_row = 10,
  fontsize_col = 9
) {
  grDevices::png(file, width = 3200, height = 1600, res = 150)
  pheatmap::pheatmap(
    scaled_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_method = "ward.D2",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    annotation_row = row_ann,
    annotation_colors = ann_colors,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    fontsize = 12,
    border_color = NA,
    main = title,
    color = HEATMAP_COLORS,
    breaks = seq(-2.5, 2.5, length.out = 101)
  )
  grDevices::dev.off()
}

plot_cosine_heatmap <- function(
  cosine_mat,
  title,
  file,
  fontsize = 9
) {
  cosine_dist <- stats::as.dist(1 - cosine_mat)
  row_clust <- stats::hclust(cosine_dist, method = "ward.D2")

  grDevices::png(file, width = 2600, height = 2600, res = 150)
  pheatmap::pheatmap(
    cosine_mat,
    cluster_rows = row_clust,
    cluster_cols = row_clust,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = fontsize,
    fontsize_row = fontsize,
    fontsize_col = fontsize,
    border_color = NA,
    main = title,
    color = COSINE_COLORS,
    breaks = seq(-1, 1, length.out = 101)
  )
  grDevices::dev.off()
}

compute_time_feature_stats <- function(image_means, feature_cols, use_trimmed, trim_lower, trim_upper) {
  n_groups <- nlevels(image_means$Time)
  n_images <- nrow(image_means)

  purrr::map_dfr(feature_cols, function(feat) {
    values <- image_means[[feat]]
    test <- stats::kruskal.test(values ~ image_means$Time)
    H <- unname(test$statistic)
    epsilon_sq <- (H - n_groups + 1) / (n_images - n_groups)

    tibble::tibble(
      Feature = feat,
      FeatureLabel = short_feature_names(feat),
      H_statistic = H,
      df = n_groups - 1,
      p_value = test$p.value,
      epsilon_squared = max(epsilon_sq, 0),
      mean_method = if (use_trimmed) {
        paste0(trim_lower * 100, "-", trim_upper * 100, "% trim + skeleton per cell")
      } else {
        "all cells + skeleton per cell"
      }
    )
  }) %>%
    dplyr::mutate(p_adj = stats::p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(dplyr::desc(epsilon_squared))
}

plot_time_feature_ranking <- function(
  time_feature_stats,
  time_colors,
  file,
  title,
  subtitle,
  n_features = 15
) {
  plot_data <- time_feature_stats %>%
    dplyr::slice_head(n = n_features) %>%
    dplyr::mutate(FeatureLabel = forcats::fct_reorder(FeatureLabel, epsilon_squared))

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = epsilon_squared, y = FeatureLabel)) +
    ggplot2::geom_col(fill = "#377EB8") +
    ggplot2::geom_vline(xintercept = 0.14, linetype = "dashed", color = "gray40") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Effect size (epsilon-squared, ε²)",
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 17),
      plot.subtitle = ggplot2::element_text(size = 12, color = "#444444"),
      axis.text = ggplot2::element_text(size = 13),
      axis.title = ggplot2::element_text(size = 14)
    )

  ggplot2::ggsave(
    file,
    plot = p,
    width = 10,
    height = max(7, 0.32 * n_features + 2),
    dpi = 150
  )
  p
}

plot_top_features_h_pvalue <- function(
  time_feature_stats,
  n_features = 25,
  file,
  title,
  subtitle
) {
  plot_data <- time_feature_stats %>%
    dplyr::slice_head(n = n_features) %>%
    dplyr::mutate(
      FeatureLabel = forcats::fct_reorder(.data$FeatureLabel, .data$epsilon_squared),
      neg_log10_p = -log10(.data$p_value)
    ) %>%
    tidyr::pivot_longer(
      c("H_statistic", "neg_log10_p"),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = factor(
        .data$metric,
        levels = c("H_statistic", "neg_log10_p"),
        labels = c("H statistic", "-log10(p-value)")
      )
    )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$value, y = .data$FeatureLabel)) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$metric),
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(~ metric, scales = "free_x", ncol = 2) +
    ggplot2::scale_fill_manual(values = c("#377EB8", "#E41A1C")) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 17),
      plot.subtitle = ggplot2::element_text(size = 12, color = "#444444"),
      axis.text = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 13),
      strip.text = ggplot2::element_text(face = "bold", size = 12)
    )

  ggplot2::ggsave(file, plot = p, width = 12, height = 10, dpi = 150)
  p
}

plot_time_boxplots <- function(
  image_means,
  top_features,
  time_colors,
  title,
  subtitle,
  file = NULL,
  width = 11,
  height = 8
) {
  label_levels <- short_feature_names(top_features)
  boxplot_data <- image_means %>%
    dplyr::select(Time, dplyr::all_of(top_features)) %>%
    tidyr::pivot_longer(-Time, names_to = "Feature", values_to = "Value") %>%
    dplyr::mutate(Feature = factor(short_feature_names(Feature), levels = label_levels))

  p <- ggplot2::ggplot(boxplot_data, ggplot2::aes(x = Time, y = Value, fill = Time)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
    ggplot2::facet_wrap(~ Feature, scales = "free_y", ncol = 2) +
    ggplot2::scale_fill_manual(values = time_colors) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Time",
      y = "Image-level value"
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", size = 17),
      plot.subtitle = ggplot2::element_text(size = 12),
      strip.text = ggplot2::element_text(face = "bold", size = 13),
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14)
    )

  if (!is.null(file)) {
    ggplot2::ggsave(file, plot = p, width = width, height = height, dpi = 150)
  }
  p
}

SKELETON_INTENSITY_FEATURES <- c(
  "Intensity_IntegratedIntensity_GFAP_Images",
  "ObjectSkeleton_NumberNonTrunkBranches_MorphologicalSkeleton",
  "ObjectSkeleton_TotalObjectSkeletonLength_MorphologicalSkeleton"
)

select_top_plus_extras <- function(ranked_stats, n, extras = SKELETON_INTENSITY_FEATURES) {
  top <- ranked_stats %>% dplyr::slice_head(n = n)
  extra_rows <- ranked_stats %>%
    dplyr::filter(.data$Feature %in% extras, !.data$Feature %in% top$Feature)
  c(top$Feature, extra_rows$Feature)
}

dedup_features_by_cosine <- function(
  ranked_stats,
  feature_mat,
  candidate_features,
  cosine_threshold = 0.85,
  force_include = character(0)
) {
  candidates <- ranked_stats %>%
    dplyr::filter(.data$Feature %in% candidate_features) %>%
    dplyr::arrange(dplyr::desc(.data$epsilon_squared)) %>%
    dplyr::pull(.data$Feature)

  candidates <- intersect(candidates, colnames(feature_mat))
  if (length(candidates) == 0) {
    return(intersect(force_include, colnames(feature_mat)))
  }

  cos_mat <- cosine_similarity_matrix(feature_mat[, candidates, drop = FALSE])
  selected <- character(0)

  for (feat in candidates) {
    if (length(selected) == 0) {
      selected <- feat
      next
    }
    feat_lab <- short_feature_names(feat)
    sel_labs <- short_feature_names(selected)
    max_sim <- max(abs(cos_mat[feat_lab, sel_labs, drop = FALSE]))
    if (max_sim < cosine_threshold) {
      selected <- c(selected, feat)
    }
  }

  extras <- setdiff(intersect(force_include, colnames(feature_mat)), selected)
  c(selected, extras)
}


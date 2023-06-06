# Enriched heatmap script, adapted from original script received from Sahu:
# Also parts adapted from: https://github.com/wmegchel/starrseq2020/blob/master/Figure1/Fig1D_STARRseq_EnrichedHeatmap.R

# Author: Konsta Karttunen
# 10.2.2022
# Updated 30.12.2022

# PARAMETERS

# Input:
# 1. List of bed files in a named GRangesList
# 2. List of .bw files for which the signal is plotted for
# 3. List of names for the signal files
# 4. A dataframe of parameters for the matrix normalization
# 5. Output path and file name (default output is .pdf)
# 6. Integer value of how many basepairs the input peaks are extended.

suppressPackageStartupMessages({
  library(BiocParallel)
  library(parallel)
  library(EnrichedHeatmap)
  library(rtracklayer)
  library(circlize)
  library(data.table)
  library(tidyverse)
  library(matrixStats)
})

# https://www.bioconductor.org/packages/release/bioc/vignettes/EnrichedHeatmap/inst/doc/EnrichedHeatmap.html
makeEnrichedHeatmap <- function(beds_list, bws_list, names_list, norm_params, output, extension, p) {
  # Splitting the enriched heatmap plot by the bed files:
  part <- unlist(lapply(seq_along(beds_list), function(i) {
    name <- names(beds_list)[i]
    len <- length(beds_list[[i]])
    c(rep(name, len))
  }))

  # Unlist does not work for some reason thus using lapply:
  tmp_targets <- GRanges(bind_rows(lapply(beds_list, as_tibble)))
  tmp_targets_extended <- resize(tmp_targets, fix = "center", width = extension * 2)

  bws <- mclapply(seq_along(bws_list), mc.cores = p, function(i) {
    import.bw(bws_list[i], format = "bigwig", selection = BigWigSelection(tmp_targets_extended))
  })

  normMatrices <- mclapply(seq_along(bws), mc.cores = p, function(i) {
    normalizeToMatrix(
      signal = bws[[i]],
      target = resize(tmp_targets, fix = "center", width = 1),
      value_column = norm_params[i, ]$value_column,
      w = norm_params[i, ]$window,
      keep = c(norm_params[i, ]$keep_low, norm_params[i, ]$keep_high),
      background = norm_params[i, ]$background,
      target_ratio = norm_params[i, ]$target_ratio,
      mean_mode = norm_params[i, ]$mean_mode,
      extend = extension,
      smooth = norm_params[i, ]$smooth
    )
  })

  # Heatmap script
  ht_list <- NULL
  for (i in seq_along(normMatrices)) {
    col_max <- max(round(quantile(matrixStats::rowMaxs(normMatrices[[i]]), probs = 0.75), 1))

    ht_list <- ht_list + EnrichedHeatmap(
      column_title_gp = gpar(fontsize = 10),
      mat = normMatrices[[i]],
      row_split = part,
      cluster_rows = T,
      col = colorRamp2(c(0, .5 * col_max, col_max), c(heat_cols_red[1], heat_cols_red[2], heat_cols_red[3])),
      width = unit(20, "mm"),
      use_raster = T,
      raster_quality = 10,
      raster_device = "CairoPNG",
      raster_by_magick = F,
      name = names_list[i],
      column_title = names_list[i],
      axis_name = c("-3kb", "0", "+3kb"),
      axis_name_gp = gpar(fontsize = 10),

      # Annotations for the top of the heatmap
      top_annotation = HeatmapAnnotation(
        enriched = anno_enriched(
          gp = gpar(col = cluster_pal),
          height = unit(1.25, "cm"),
          axis_param = list(
            side = "right",
            facing = "inside"
          ),
        )
      ),
      gap = unit(1, "mm"),

      # Parameters for the legend
      heatmap_legend_param = list(
        title = NULL,
        color_bar = "continuous",
        legend_direction = "horizontal",
        legend_width = unit(15, "mm"),
        grid_border = "#555555", at = c(0, col_max),
        labels = c("0", col_max)
      )
    )
  }

  ht_list <- Heatmap(part,
    cluster_rows = F,
    cluster_row_slices = F,
    col = structure(cluster_pal, names = names(beds_list)),
    name = "Cluster",
    show_row_names = F,
    show_heatmap_legend = F,
    width = unit(5, "mm"),
    use_raster = T,
    raster_quality = 10,
    raster_by_magick = F,
    raster_device = "CairoPNG"
  ) + ht_list

  pdf(output, width = 12, height = 8)
  draw(
    ht_list,
    split = part,
    heatmap_legend_side = "bottom",
    ht_gap = unit(1, "mm"),
    padding = unit(c(3, 3, 3, 3), "mm")
  )
  dev.off()
}

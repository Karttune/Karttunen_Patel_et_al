library(tidyverse)
library(parallel)
library(profileplyr)

deeptools_mat <- read.table("/data/kzkarttu/gradu/raw_data/GP5D/220914_deeptools_GP5D_WT_STARR_peaks_histone_modifications_input_normalized_matrix.out", skip = 2, header = T)
deeptools_mat_scaled <- scale(deeptools_mat)

clusters <- 2:9
cores <- 40

read.table("/data/kzkarttu/gradu/scripts/R_scripts/gradu/data/deeptools/220914_deeptools_GP5D_WT_STARR_peaks_histone_modifications_normalized_matrix.out", skip = 2, header = T)

# Using mclapply to parallerize the kmeans cluster calculations to n cores:  
kmeans_clusters <- mclapply(seq_along(clusters), mc.cores = cores, function(i) {
  kmeans(deeptools_mat_scaled, clusters[i], nstart = 50, iter.max = 15)
})

path <- paste0("data/kmeans/", today, "_GP5D_STARR_kmeans_input_normalized_scaled.rds")
saveRDS(kmeans_clusters, path)
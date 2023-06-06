library(bedtoolsr)
library(rstatix)
library(parallel)

chrom_sizes <- "annotations/hg38.canonical.chrom.sizes"
mask <- "annotations/BSgenome_hg38_mask.bed"

# Setting the seed to the same value for each run for reproducibility
seeds <- c(1000:1999)

shuffle_random_regions <- function(peak_path, n, rmsk, level = "subfamily", cores_n = 60) {
  cores <- cores_n

  if (level == "subfamily") {
    shuf_counts <- tibble(subf = sort(unique(rmsk$subf)), freq = 0)
  } else if (level == "class") {
    shuf_counts <- tibble(subf = sort(unique(rmsk$class)), freq = 0)
  } else if (level == "lineage") {
    shuf_counts <- tibble(subf = sort(unique(rmsk$lineage)), freq = 0)
  } else {
    print("Level has to be 'subfamily', 'class' or 'lineage'! (And there must be a corresponding column in the repeatmasker annotation)")
    break
  }

  out <- bind_rows(shuf_counts, mclapply(1:n, mc.cores = cores, function(i) {
    shuffled <- bt.shuffle(i = peak_path, g = chrom_sizes, excl = mask, chrom = T, seed = seeds[i])
    colnames(shuffled) <- narrowpeak_colnames[1:length(shuffled)]

    shuffled <- GRanges(shuffled)

    if (i %% 10 == 0) {
      print(paste0(i, " shuffles"))
    }

    ovl <- findOverlaps(shuffled, rmsk, ignore.strand = T)
    nonte_ovl <- data.frame(Var1 = "Non_TE", Freq = queryLength(ovl) - length(unique(queryHits(ovl))))

    # Making a table of how many TEs were overlapped and turning into a data frame:
    # If level is set as "subfamily", the overlaps are counted by subfamily, also accepts "family" or "class"

    if (level == "subfamily") {
      res <- as.data.frame(table(rmsk[subjectHits(ovl), ]$subf)) %>%
        bind_rows(., nonte_ovl)
    } else if (level == "class") {
      res <- as.data.frame(table(rmsk[subjectHits(ovl), ]$class)) %>%
        bind_rows(., nonte_ovl)
    } else if (level == "lineage") {
      res <- as.data.frame(table(rmsk[subjectHits(ovl), ]$lineage)) %>%
        bind_rows(., nonte_ovl)
    }

    if (nrow(res) != 0) {
      if (level == "subfamily") {
        res %>%
          dplyr::rename(subf = Var1) %>%
          group_by(subf) %>% # Evaluate following calls for each value in the rowname column
          dplyr::summarise(freq = sum(Freq)) # sum all non-grouping variables
      } else if (level == "class") {
        res %>%
          dplyr::rename(class = Var1) %>%
          group_by(class) %>%
          dplyr::summarise(freq = sum(Freq))
      } else if (level == "lineage") {
        res %>%
          dplyr::rename(lineage = Var1) %>%
          group_by(lineage) %>%
          dplyr::summarise(freq = sum(Freq))
      }
    }
  }))

  if (level == "subfamily") {
    final_out <- out %>%
      group_by(subf) %>%
      dplyr::summarise(freq = sum(freq)) %>%
      mutate(mean_freq = freq / n)
  } else if (level == "class") {
    final_out <- out %>%
      group_by(class) %>%
      dplyr::summarise(freq = sum(freq)) %>%
      mutate(mean_freq = freq / n)
  } else if (level == "lineage") {
    final_out <- out %>%
      group_by(lineage) %>%
      dplyr::summarise(freq = sum(freq)) %>%
      mutate(mean_freq = freq / n)
  }
  return(final_out)
}

process_binomial_testing <- function(shuf, ovl_counts, grp_factor, threshold = 0.01) {
  
  if(!grp_factor %in% c("subf", "class", "lineage")) {
    print("Grouping factor must to be 'subf', 'class' or 'lineage'! (And there must be a corresponding column in the repeatmasker annotation)")
    break
  }
  
  out <- shuf %>%
    mutate(exp = mean_freq / sum(mean_freq)) %>%
    left_join(., ovl_counts, by = grp_factor) %>%
    mutate(
      n = replace(n, is.na(n), 0),
      n_all = sum(n)
    ) %>%
    rowwise() %>%
    mutate(p = binom_test(n, n_all, exp, alternative = "greater")$p) %>%
    ungroup() %>%
    mutate(p_adj = p.adjust(p, method = "BH")) %>%
    filter(mean_freq != 0) %>% # To avoid division by zero, filtering out subfamilies that had no overlaps
    mutate(
      obs_exp = n / mean_freq,
      sym = case_when(
        p_adj <= 0.0001 ~ "****",
        p_adj > 0.0001 & p_adj <= 0.001 ~ "***",
        p_adj > 0.001 & p_adj <= 0.01 ~ "**",
        p_adj > 0.01 & p_adj <= 0.05 ~ "*",
        p_adj > 0.05 ~ "ns"
      )
    )

  if (grp_factor == "subf") {
    out %>%
      mutate(
        class = str_extract(subf, "LINE$|SINE$|LTR$|DNA$"),
        sig = factor(ifelse(p_adj < threshold, class, "insig"))
      )
  } else if (grp_factor == "lineage") {
    out %>%
      mutate(
        sig = factor(ifelse(p_adj < threshold, lineage, "insig"))
      )
  } else {
    out %>%
      mutate(
        sig = factor(ifelse(p_adj < threshold, class, "insig"))
      )
  }
}

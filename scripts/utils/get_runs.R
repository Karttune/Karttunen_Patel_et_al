library(tidyverse)
library(parallel)
library(GenomicRanges)

get_cgcalls <- function(regs.sub, reads.cpg) {
  cgcalls <- bind_rows(mclapply(seq(nrow(regs.sub)), mc.cores = 52, function(i) {
    print(i)
    reg <- regs.sub[i, ]
    callreg <- regs.sub[i, ]
    reginfo <- reg %>%
      dplyr::rename(regstart = start, regend = end) %>%
      dplyr::select(-chrom)
    # subset reads overlapping this region
    cg.reg <- reads.cpg %>%
      filter(chrom == reg$chrom, start <= reg$start, end >= reg$end)
    calls <- tibble()
    if (nrow(cg.reg) != 0) {
      # fix calls
      calls <- mbedByCall(cg.reg, region = callreg, verbose = F) %>%
        redo_mcall(1.5)
      # add info and label, get distance
      calls <- calls %>%
        bind_cols(reginfo[rep(1, nrow(calls)), ])
    }
    calls
  }))
  cgcalls <- cgcalls %>%
    mutate(
      center = (regstart + regend) / 2,
      distance = ifelse(strand == "-", center - start, start - center)
    )
  cgcalls
}

get_gcruns <- function(regs.sub, reads.gpc) {
  gcruns <- bind_rows(mclapply(mc.cores = 52, seq(nrow(regs.sub)), function(i) {
    reg <- regs.sub[i, ]
    # subset reads overlapping this region
    gc.reg <- reads.gpc %>%
      filter(chrom == reg$chrom, start <= reg$start, end >= reg$end)
    runs <- tibble()

    if (nrow(gc.reg) != 0) {
      callreg <- regs.sub[i, ]
      reginfo <- reg %>%
        dplyr::rename(regstart = start, regend = end) %>%
        dplyr::select(-chrom)
      # fix calls
      gccalls <- mbedByCall(gc.reg, verbose = F) %>%
        redo_mcall(1)

      # smooth
      calls.reg <- gccalls %>%
        group_by(qname)
      calls.list <- calls.reg %>%
        group_split(.keep = T)
      smooth.list <- lapply(calls.list, smoothCalls, reg = callreg)

      calls.smooth <- bind_rows(smooth.list)
      # add info and label, get distance
      runs <- getRuns_fast(calls.smooth)
      runs <- runs %>%
        bind_cols(reginfo[rep(1, nrow(runs)), ])
    }
    runs
  }))

  gcruns <- gcruns %>%
    mutate(
      center = (regstart + regend) / 2,
      # get distance from center
      start_dist = ifelse(strand == "-", center - start, start - center),
      end_dist = ifelse(strand == "-", center - end, end - center),
      # get abs values of distance from the center
      distance = case_when(
        start_dist <= 0 & end_dist >= 0 ~ 0,
        start_dist > 0 ~ start_dist,
        end_dist < 0 ~ abs(end_dist)
      ),
      acc = ifelse(values == 0, "Closed", "Open")
    )
  gcruns
}

get_heatmap_runs <- function(dists, gcruns) {
  
  runs.heat <- bind_rows(mclapply(mc.cores = 52, dists, function(d) {
    x <- gcruns[gcruns$start_dist <= d & gcruns$end_dist >= d & gcruns$width <= 510, ]
    y <- x %>%
      group_by(tfap2a, te, width, values) %>%
      dplyr::summarize(n = dplyr::n())
    # roll mean along width for this distance
    bind_rows(lapply(seq(0, 500), function(s) {
      y[which(y$width >= s & y$width <= s + 10), ] %>%
        group_by(tfap2a, te, values) %>%
        dplyr::summarize(n = sum(n)) %>%
        mutate(width = s)
    })) %>%
      mutate(d = d)
  }))
}


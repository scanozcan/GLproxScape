# =============================================================================
# caspex_extras.R
# Supplementary result plots & diagnostics for the CasPEX pipeline.
# Source AFTER caspex_analysis.R:
#
#   source("caspex_analysis.R")
#   source("caspex_extras.R")
#   result <- run_caspex(gene = "ATP7B", grnas = ..., data_files = ...)
#   extras <- run_caspex_extras(result, out_dir = "caspex_output/extras")
#
# Every top-level function takes the invisible list returned by run_caspex()
# and returns either a ggplot object, a named list of ggplots, or (for
# expensive computations) a results object that its matching plot_* function
# can render.
#
# Contents
#   A. Biological interpretation
#       plot_tf_one_pager, plot_tf_family_enrichment,
#       plot_event_density, plot_composite_vs_specificity,
#       plot_tf_cooccurrence, rank_binding_events
#   B. Statistical robustness
#       run_permutation_null, plot_permutation_null,
#       run_sigma_sensitivity, plot_sigma_sensitivity,
#       run_event_jackknife, plot_event_jackknife,
#       plot_nnls_residual
#   C. QC / data sanity
#       plot_volcano_per_region, plot_region_correlation,
#       plot_pval_histograms, plot_motif_vs_nnls
#   D. Coverage-aware diagnostics
#       plot_coverage_rescue_scatter,
#       run_covfloor_sensitivity, plot_covfloor_sensitivity,
#       plot_coverage_stack
#   Wrapper
#       run_caspex_extras
# =============================================================================

# Package-internal dependencies (`%||%`, `COLS`, `.safe_pdf`,
# `theme_caspex`, every plot helper from caspex_analysis.R, every
# ChIP-Atlas helper from caspex_chipatlas.R) are visible at parse time
# because all R/ files share the package namespace. The script-style
# library() preamble and the "verify caspex_analysis.R was sourced"
# guard that used to live here are no longer needed.

#' Write a multi-page sensitivity-sweep PDF.
#'
#' Internal helper used by B2 (sigma) and D2 (cov_floor) so neither
#' panel grid gets so dense that it's unreadable. Splits the long-form
#' sweep data.frame `df` (must contain a `tf` column) into chunks of
#' `tfs_per_page` TFs and prints one page per chunk by calling
#' `plot_fn(df_chunk)`.
#'
#' @param df Long-form sweep data.frame with a `tf` column.
#' @param plot_fn Plot function that takes a subset of `df` and returns a ggplot.
#' @param out_path Output PDF path.
#' @param tfs_per_page Number of TFs to render per page (default 12).
#' @param width PDF width in inches (default 12).
#' @param height PDF height in inches (default 8).
#' @return Invisibly, NULL.
#' @noRd
.save_sensitivity_paginated <- function(df, plot_fn, out_path,
                                         tfs_per_page = 12,
                                         width = 12, height = 8) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  tfs <- sort(unique(df$tf))
  if (length(tfs) == 0) return(invisible(NULL))
  pdf(out_path, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  for (i in seq(1L, length(tfs), by = tfs_per_page)) {
    chunk <- tfs[i:min(i + tfs_per_page - 1L, length(tfs))]
    print(plot_fn(df[df$tf %in% chunk, , drop = FALSE]))
  }
}

# =============================================================================
# A.1  Per-TF one-pager dashboard
# =============================================================================

#' One page summarising a single TF: per-region bars, spatial footprint,
#' binding-event track, and a compact stats block.
#'
#' @param result   The invisible list returned by run_caspex()
#' @param tf_name  TF symbol to summarise
#' @return A patchwork ggplot assembly
#' @export
plot_tf_one_pager <- function(result, tf_name) {
  if (!tf_name %in% result$long_data$protein)
    stop(tf_name, " not found in result$long_data")

  ld       <- result$long_data
  pos_map  <- result$pos_map
  motif_hits <- if (!is.null(result$motif_results[[tf_name]]))
    result$motif_results[[tf_name]]$hits else integer(0)

  # Top panel: deconvolution view (full signal + events + motif ticks).
  # Honour whichever binding-path mode produced `result` \u2014 otherwise the
  # one-pager would silently draw smoothed-NNLS events even when the main
  # run used coverage-aware scoring, and would disagree with result$binding_events.
  # Pull upstream/downstream from the result so the A1 x-axis matches
  # the actual window the run used (e.g. -750/+500 for Pizzolato/Gao),
  # not plot_binding_deconvolution's legacy default (-2500/+500).
  p_top <- plot_binding_deconvolution(
    tf_name, ld, pos_map, motif_hits,
    weight_mode      = result$weight_mode      %||% "mod_t",
    cov_floor        = result$cov_floor        %||% 0.05,
    kernel_sigma     = result$kernel_sigma     %||% 250,
    upstream         = result$upstream         %||% 2500,
    downstream       = result$downstream       %||% 500
  ) + labs(title = NULL, subtitle = NULL)

  # Middle panel: per-region bar of logFC with p-value stars
  df <- ld[ld$protein == tf_name & ld$region %in% names(pos_map), ]
  df$pos    <- as.numeric(pos_map[df$region])
  df <- df[order(df$pos), ]
  df$region <- factor(df$region, levels = df$region)
  df$stars  <- cut(df$pval,
                   breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf),
                   labels = c("****", "***", "**", "*", "ns"))
  lfc_range <- range(c(0, df$lfc), na.rm = TRUE)
  p_bars <- ggplot(df, aes(x = region, y = lfc)) +
    geom_col(aes(fill = pval <= 0.05), alpha = 0.8, width = 0.7) +
    geom_text(aes(y = ifelse(lfc >= 0, lfc + diff(lfc_range) * 0.02,
                                       lfc - diff(lfc_range) * 0.02),
                  label = stars,
                  vjust = ifelse(lfc >= 0, 0, 1)),
              size = 3.2, fontface = "bold") +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    scale_fill_manual(values = c(`TRUE` = COLS$high, `FALSE` = COLS$guide),
                      labels = c(`TRUE` = "p \u2264 0.05",
                                 `FALSE` = "p > 0.05"),
                      name = NULL) +
    labs(x = NULL, y = "logFC") +
    theme_caspex() +
    theme(legend.position = "bottom")

  # Bottom panel: stats block as text
  sp_row <- result$spatial_df[result$spatial_df$protein == tf_name, , drop = FALSE]
  if (nrow(sp_row) == 0) sp_row <- data.frame(centroid = NA, spread = NA,
                                              composite = NA, n_regions = NA)
  ev <- if (!is.null(result$binding_events))
    result$binding_events[result$binding_events$tf == tf_name, , drop = FALSE] else
      data.frame()
  pwm_id <- if (!is.null(result$motif_results[[tf_name]]))
    result$motif_results[[tf_name]]$pwm$id else "n/a"
  mode_lbl <- paste0("coverage-aware (cov_floor=",
                     result$cov_floor %||% 0.05, ")")
  stats_txt <- sprintf(
    "%s | centroid %s bp | spread %s bp | composite %s | n_regions %s\nmatrix %s | PWM hits %d | events %d | mode: %s",
    tf_name,
    format(sp_row$centroid[1], nsmall = 1),
    format(sp_row$spread[1], nsmall = 1),
    format(sp_row$composite[1], nsmall = 3),
    sp_row$n_regions[1],
    pwm_id, length(motif_hits), nrow(ev), mode_lbl
  )
  p_stats <- ggplot() +
    annotate("text", x = 0, y = 0.5, label = stats_txt,
             hjust = 0, vjust = 0.5, size = 3.4,
             family = "mono", color = COLS$neutral) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_void()

  (p_top / p_bars / p_stats) +
    plot_layout(heights = c(6, 3, 0.8)) +
    plot_annotation(
      title    = paste0(tf_name, " \u2014 CasPEX summary"),
      subtitle = paste0("Region enrichment \u00b7 spatial footprint \u00b7 binding-event calls"),
      theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                       plot.subtitle = element_text(color = "grey50"))
    )
}

# =============================================================================
# A.2  TF-family enrichment among the motif-scanned TFs
# =============================================================================

#' Barplot of TF-family representation among the top-ranked TFs.
#'
#' Queries JASPAR for each TF's "family" annotation (cached across calls
#' via a list in the returned plot's environment). If JASPAR family info
#' is unavailable, the bar is tagged "unknown".
#'
#' @param result The result object from run_caspex()
#' @param top_n  Number of top-composite TFs to include (default 30)
#' @return ggplot
#' @export
plot_tf_family_enrichment <- function(result, top_n = 30) {
  spatial <- head(result$spatial_df, top_n)
  tfs <- as.character(spatial$protein)

  fam <- vapply(tfs, function(tf) {
    tryCatch({
      u <- paste0("https://jaspar.elixir.no/api/v1/matrix/?tf_name=", tf,
                  "&collection=CORE&format=json&page_size=1")
      r <- httr::GET(u, httr::timeout(10))
      if (httr::status_code(r) != 200) return("unknown")
      js <- httr::content(r, "parsed", simplifyVector = FALSE)
      if (length(js$results) == 0) return("unknown")
      f <- js$results[[1]]$family
      if (is.null(f) || length(f) == 0) return("unknown")
      paste(unlist(f), collapse = "/")
    }, error = function(e) "unknown")
  }, character(1))

  df <- data.frame(tf = tfs, family = fam, stringsAsFactors = FALSE)
  fam_tab <- as.data.frame(table(df$family))
  names(fam_tab) <- c("family", "n")
  fam_tab <- fam_tab[order(fam_tab$n, decreasing = TRUE), ]
  fam_tab$family <- factor(fam_tab$family, levels = fam_tab$family)

  ggplot(fam_tab, aes(x = family, y = n)) +
    geom_col(fill = COLS$guide, alpha = 0.85) +
    geom_text(aes(label = n), vjust = -0.3, size = 3, fontface = "bold",
              color = COLS$neutral) +
    labs(x = "JASPAR TF family",
         y = paste0("# TFs in top-", top_n),
         title = paste0("TF family representation among top-", top_n,
                        " composite TFs"),
         subtitle = paste0(nrow(df), " TFs queried | ",
                           sum(df$family == "unknown"),
                           " with no family annotation")) +
    theme_caspex() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# =============================================================================
# A.3  TSS-relative event density
# =============================================================================

#' Histogram of predicted binding-event positions across all TFs.
#'
#' @param result   Result object
#' @param binwidth Histogram binwidth in bp (default 50)
#' @param weighted Weight each event by its NNLS weight (default TRUE)
#' @export
plot_event_density <- function(result, binwidth = 50, weighted = TRUE) {
  ev <- result$binding_events
  if (is.null(ev) || nrow(ev) == 0)
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No binding events called"))
  ev$w <- if (weighted) ev$weight else 1

  ggplot(ev, aes(x = position, weight = w)) +
    geom_histogram(binwidth = binwidth, fill = COLS$guide,
                   color = "white", alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = COLS$tss, linewidth = 0.7) +
    annotate("text", x = 0, y = 0, label = " TSS",
             color = COLS$tss, hjust = 0, vjust = -0.5, size = 3.2,
             fontface = "bold") +
    scale_x_continuous(labels = scales::comma) +
    labs(x = "Position (bp, TSS-relative)",
         y = if (weighted) "Summed event weight" else "# events",
         title = "TSS-relative distribution of predicted binding events",
         subtitle = sprintf("%d events across %d TFs | binwidth %d bp | %sweighted",
                            nrow(ev), length(unique(ev$tf)), binwidth,
                            if (weighted) "" else "un")) +
    theme_caspex()
}

# =============================================================================
# A.4  Composite vs specificity scatter
# =============================================================================

#' 2-D scatter of composite enrichment vs per-TF specificity score.
#'
#' Specificity is defined per TF as the maximum over regions of
#' \eqn{lfc_R - mean_{R' \neq R}(lfc_{R'})}.
#'
#' @param result Output of \code{\link{run_caspex}}.
#' @param top_label Number of TFs to label on the plot (top by composite
#'   score). Default 15.
#' @return A ggplot.
#' @export
plot_composite_vs_specificity <- function(result, top_label = 15) {
  ld <- result$long_data[result$long_data$protein %in% result$spatial_df$protein &
                            result$long_data$region %in% names(result$pos_map), ]

  spec <- vapply(unique(ld$protein), function(p) {
    sub <- ld[ld$protein == p, ]
    if (nrow(sub) < 2) return(NA_real_)
    mean_by_r <- tapply(sub$lfc, sub$region, mean, na.rm = TRUE)
    rs <- names(mean_by_r)
    vals <- vapply(rs, function(r) mean_by_r[r] - mean(mean_by_r[rs != r]),
                   numeric(1))
    max(vals, na.rm = TRUE)
  }, numeric(1))

  df <- data.frame(
    protein   = names(spec),
    specificity = as.numeric(spec),
    composite = result$spatial_df$composite[
      match(names(spec), result$spatial_df$protein)],
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$composite) & !is.na(df$specificity), ]
  df$score <- df$composite * pmax(df$specificity, 0)   # for labelling rank
  df <- df[order(df$score, decreasing = TRUE), ]
  df$label <- ifelse(seq_len(nrow(df)) <= top_label, df$protein, "")

  ggplot(df, aes(x = composite, y = specificity)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey60", linewidth = 0.3) +
    geom_point(aes(size = composite), color = COLS$guide,
               alpha = 0.55, shape = 16) +
    geom_text(aes(label = label), vjust = -0.9,
              size = 3, fontface = "bold", color = COLS$neutral) +
    scale_size_continuous(range = c(1.2, 7), guide = "none") +
    labs(x = "Composite enrichment (global strength)",
         y = "Max per-region specificity (lfc_R \u2212 mean lfc_other)",
         title = "Broadly strong vs region-specific TFs",
         subtitle = paste0("Top-right: strong + region-focal | ",
                           "Top-left: focal but weaker | ",
                           "Bottom-right: diffuse")) +
    theme_caspex()
}

# =============================================================================
# B.1  Permutation null for the composite score
# =============================================================================

#' Permutation null: shuffle region labels B times, recompute the spatial
#' model, and build an empirical null distribution of composite scores.
#'
#' Returns a list with the observed composite per TF and the per-TF
#' permutation p-value (= fraction of null scores >= observed).
#' @param result (see function body).
#' @param n_perm (see function body).
#' @param seed (see function body).
#' @param weight_mode (see function body).
#' @export
run_permutation_null <- function(result, n_perm = 500, seed = 42,
                                  weight_mode = NULL) {
  set.seed(seed)
  weight_mode <- weight_mode %||% result$weight_mode %||% "mod_t"
  ld <- result$long_data[result$long_data$isTF, ]
  ld <- ld[ld$region %in% names(result$pos_map), ]
  pos <- result$pos_map
  proteins <- unique(ld$protein)

  observed <- setNames(result$spatial_df$composite,
                        as.character(result$spatial_df$protein))

  message("Running ", n_perm, " permutations across ", length(proteins), " TFs...")
  nulls <- matrix(NA_real_, nrow = length(proteins), ncol = n_perm,
                  dimnames = list(proteins, NULL))
  regs <- names(pos)

  for (b in seq_len(n_perm)) {
    if (b %% max(1, floor(n_perm / 10)) == 0)
      message("  permutation ", b, " / ", n_perm)
    perm_pos <- setNames(sample(pos), regs)
    per_tf <- vapply(proteins, function(p) {
      sub <- ld[ld$protein == p, ]
      sub$pos <- perm_pos[as.character(sub$region)]
      sub <- sub[!is.na(sub$pos) & !is.na(sub$lfc), ]
      if (nrow(sub) < 2) return(NA_real_)
      # Mirror compute_spatial(): floor lfc BEFORE weighting, then floor w
      sub$lfc <- pmax(sub$lfc, 0)
      w <- pmax(compute_region_weight(sub, mode = weight_mode), 0)
      if (sum(w) < 1e-8) return(0)
      mean(w) * log1p(nrow(sub))
    }, numeric(1))
    nulls[, b] <- per_tf
  }

  pvals <- vapply(proteins, function(p) {
    obs <- observed[[p]]
    if (is.null(obs) || is.na(obs)) return(NA_real_)
    v <- nulls[p, ]
    v <- v[!is.na(v)]
    if (length(v) == 0) return(NA_real_)
    (sum(v >= obs) + 1) / (length(v) + 1)
  }, numeric(1))

  out <- data.frame(
    protein   = names(pvals),
    composite = observed[names(pvals)],
    perm_p    = pvals,
    null_mean = rowMeans(nulls, na.rm = TRUE),
    null_sd   = apply(nulls, 1, sd, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  out$perm_fdr <- p.adjust(out$perm_p, method = "BH")
  out <- out[order(out$perm_p), ]
  list(summary = out, nulls = nulls, n_perm = n_perm, weight_mode = weight_mode)
}

#' Plot the permutation null: per-TF observed composite overlaid on its
#' null distribution. Shows top_n most-significant TFs.
#' @param perm_result (see function body).
#' @param top_n (see function body).
#' @export
plot_permutation_null <- function(perm_result, top_n = 20) {
  s <- head(perm_result$summary, top_n)
  s$protein <- factor(s$protein, levels = rev(s$protein))
  s$sig <- s$perm_fdr <= 0.05

  ggplot(s, aes(y = protein)) +
    geom_errorbarh(aes(xmin = null_mean - null_sd,
                       xmax = null_mean + null_sd),
                   color = "grey70", height = 0) +
    geom_point(aes(x = null_mean), color = "grey50", size = 2, shape = 4) +
    geom_point(aes(x = composite, fill = sig),
               shape = 21, size = 3.5, color = "black", stroke = 0.4) +
    scale_fill_manual(values = c(`TRUE` = COLS$high, `FALSE` = COLS$guide),
                      labels = c(`TRUE` = "FDR \u2264 0.05",
                                 `FALSE` = "FDR > 0.05"),
                      name = NULL) +
    geom_text(aes(x = composite,
                  label = sprintf("p=%.3g, q=%.3g", perm_p, perm_fdr)),
              hjust = -0.1, size = 2.8, color = COLS$neutral) +
    labs(x = "Composite score (observed vs null mean \u00b1 1 sd)",
         y = NULL,
         title = paste0("Permutation null (B=", perm_result$n_perm,
                        ") for composite score"),
         subtitle = paste0("Top ", top_n, " TFs by permutation p-value | ",
                           "weight mode: ", perm_result$weight_mode)) +
    theme_caspex() +
    theme(legend.position = "bottom",
          panel.grid.major.y = element_blank())
}

# =============================================================================
# B.2  Kernel sigma sensitivity
# =============================================================================

#' For each TF, rerun the deconvolution at multiple kernel widths and
#' capture the top event position per sigma. Stable centroids = robust.
#'
#' Uses the coverage-aware deconvolution at every σ, matching the
#' decoder that produced `result$binding_events`.
#' @param result (see function body).
#' @param sigmas (see function body).
#' @param tfs (see function body).
#' @param weight_mode (see function body).
#' @param cov_floor (see function body).
#' @param edge_guard_frac (see function body).
#' @export
run_sigma_sensitivity <- function(result,
                                   sigmas = c(100, 200, 250, 300, 500),
                                   tfs    = NULL,
                                   weight_mode      = NULL,
                                   cov_floor        = NULL,
                                   edge_guard_frac  = NULL) {
  weight_mode      <- weight_mode      %||% result$weight_mode      %||% "mod_t"
  cov_floor        <- cov_floor        %||% result$cov_floor        %||% 0.05
  edge_guard_frac  <- edge_guard_frac  %||% result$edge_guard_frac  %||% 0.15
  if (is.null(tfs)) {
    # Every TF that produced JASPAR motif-anchored binding events.
    # Intersect with motif_results PLUS motif_results_extra so augment-
    # pool TFs (motif_scan_pool='spatial_all') are visible -- without
    # the extra set, FOXP4-type augment-pool TFs are silently filtered
    # out even though they have motif-anchored events.
    be <- result$binding_events
    motif_tfs_with_events <- if (!is.null(be$motif_based))
      unique(be$tf[be$motif_based %in% TRUE]) else unique(be$tf)
    all_motif_names <- union(names(result$motif_results),
                              names(result$motif_results_extra %||% list()))
    tfs <- intersect(motif_tfs_with_events, all_motif_names)
  }
  if (length(tfs) == 0) return(invisible(NULL))

  # Combined motif lookup: try motif_results first, fall back to
  # motif_results_extra for augment-pool TFs.
  .get_hits <- function(tf) {
    h <- result$motif_results[[tf]]$hits
    if (is.null(h)) h <- result$motif_results_extra[[tf]]$hits
    if (is.null(h)) integer(0) else h
  }

  message("Sigma sensitivity across ", length(sigmas),
          " kernels \u00d7 ", length(tfs), " TFs...")
  rows <- list()
  for (tf in tfs) {
    hits <- .get_hits(tf)
    for (sg in sigmas) {
      ev <- predict_binding_events_coverage_aware(
              tf, result$long_data, result$pos_map, hits,
              kernel_sigma = sg, weight_mode = weight_mode,
              cov_floor = cov_floor,
              edge_guard_frac = edge_guard_frac)
      if (nrow(ev) == 0) next
      ev <- ev[order(ev$weight, decreasing = TRUE), ]
      top_n <- min(3, nrow(ev))
      rows[[length(rows) + 1]] <- data.frame(
        tf = tf, sigma = sg,
        rank = seq_len(top_n),
        position = ev$position[seq_len(top_n)],
        weight = ev$weight[seq_len(top_n)],
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) return(invisible(NULL))
  do.call(rbind, rows)
}

#' Plot results of run_sigma_sensitivity()
#'
#' Two-panel summary: per-sigma event count, and per-event spread vs.
#' kernel sigma.
#'
#' @return A patchwork ggplot assembly.
#' @param sigma_result (see function body).
#' @export
plot_sigma_sensitivity <- function(sigma_result) {
  if (is.null(sigma_result) || nrow(sigma_result) == 0)
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No sigma-sensitivity data"))

  ggplot(sigma_result, aes(x = sigma, y = position,
                            group = interaction(tf, rank),
                            color = factor(rank))) +
    geom_line(alpha = 0.5) +
    geom_point(aes(size = weight), alpha = 0.85) +
    facet_wrap(~ tf, ncol = 3, scales = "free_y") +
    scale_color_manual(values = c(`1` = COLS$high, `2` = COLS$mid,
                                   `3` = COLS$guide),
                       name = "event rank") +
    scale_size_continuous(range = c(1.5, 5), guide = "none") +
    labs(x = "Kernel \u03c3 (bp)",
         y = "Event position (bp, TSS-relative)",
         title = "Kernel \u03c3 sensitivity of deconvolved binding events",
         subtitle = "Flat lines = robust call | large drift = artefact of smoothing width") +
    theme_caspex() +
    theme(legend.position = "bottom")
}

# =============================================================================
# B.3  Event jackknife (leave-one-region-out)
# =============================================================================

#' Drop each region in turn, rerun the deconvolution, and count how often
#' each event position survives within `tol_bp` bp.
#'
#' Leave-one-out re-fits use the coverage-aware deconvolution, matching
#' the original call.
#' @param result (see function body).
#' @param tfs (see function body).
#' @param tol_bp (see function body).
#' @param weight_mode (see function body).
#' @param cov_floor (see function body).
#' @param edge_guard_frac (see function body).
#' @export
run_event_jackknife <- function(result, tfs = NULL, tol_bp = 100,
                                 weight_mode      = NULL,
                                 cov_floor        = NULL,
                                 edge_guard_frac  = NULL) {
  weight_mode      <- weight_mode      %||% result$weight_mode      %||% "mod_t"
  cov_floor        <- cov_floor        %||% result$cov_floor        %||% 0.05
  edge_guard_frac  <- edge_guard_frac  %||% result$edge_guard_frac  %||% 0.15
  if (is.null(tfs))
    tfs <- unique(result$binding_events$tf)
  if (length(tfs) == 0) return(invisible(NULL))

  rs <- names(result$pos_map)
  message("Jackknife across ", length(rs), " regions \u00d7 ",
          length(tfs), " TFs",
          "  [coverage-aware]",
          "...")

  rows <- list()
  for (tf in tfs) {
    obs <- result$binding_events[result$binding_events$tf == tf, ]
    hits <- if (!is.null(result$motif_results[[tf]]))
      result$motif_results[[tf]]$hits else integer(0)

    surv <- integer(nrow(obs))
    for (r in rs) {
      pos_sub <- result$pos_map[setdiff(names(result$pos_map), r)]
      ev <- predict_binding_events_coverage_aware(
              tf, result$long_data, pos_sub, hits,
              weight_mode = weight_mode, cov_floor = cov_floor,
              edge_guard_frac = edge_guard_frac)
      if (nrow(ev) == 0) next
      for (j in seq_len(nrow(obs))) {
        if (any(abs(ev$position - obs$position[j]) <= tol_bp))
          surv[j] <- surv[j] + 1L
      }
    }
    rows[[length(rows) + 1]] <- data.frame(
      tf = tf,
      position = obs$position,
      weight   = obs$weight,
      n_surv   = surv,
      n_drops  = length(rs),
      frac     = surv / length(rs),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(invisible(NULL))
  do.call(rbind, rows)
}

#' Plot results of run_event_jackknife()
#'
#' Visualises which TF events survive after dropping each region in turn.
#'
#' @return A ggplot.
#' @param jk_result (see function body).
#' @param top_n (see function body).
#' @export
plot_event_jackknife <- function(jk_result, top_n = 40) {
  if (is.null(jk_result) || nrow(jk_result) == 0)
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No jackknife data"))
  jk <- jk_result[order(jk_result$frac, jk_result$weight,
                         decreasing = TRUE), ]
  jk <- head(jk, top_n)
  jk$label <- sprintf("%s @ %+.0f bp", jk$tf, jk$position)
  jk$label <- factor(jk$label, levels = rev(jk$label))

  ggplot(jk, aes(x = frac, y = label, fill = frac)) +
    geom_col(alpha = 0.9) +
    geom_text(aes(label = sprintf("%d/%d", n_surv, n_drops)),
              hjust = -0.1, size = 2.6, color = COLS$neutral) +
    scale_x_continuous(limits = c(0, 1.15),
                       breaks = seq(0, 1, 0.25),
                       labels = scales::percent) +
    scale_fill_gradient(low = COLS$low, high = COLS$high, guide = "none") +
    labs(x = "Fraction of jackknife replicates in which the event survived",
         y = NULL,
         title = "Event robustness to leave-one-region-out",
         subtitle = paste0("Top-", top_n,
                           " events by survival fraction | tolerance: \u00b1100 bp")) +
    theme_caspex() +
    theme(panel.grid.major.y = element_blank())
}

# =============================================================================
# B.4  NNLS residual plot for a TF
# =============================================================================

#' Overlay observed signal, NNLS reconstruction, and residual for a TF.
#'
#' Three-panel diagnostic for a single TF: the observed CasPEX signal
#' \eqn{s(x)}, the NNLS reconstruction \eqn{X\beta} from JASPAR motif hits,
#' and the residual. A large residual = JASPAR motifs for this TF cannot
#' fully account for the observed signal at this locus (motif-orphan
#' enrichment, partner co-binding, or coverage-aware-only call).
#'
#' @param result Output of \code{\link{run_caspex}}.
#' @param tf_name TF symbol to plot.
#' @param kernel_sigma Gaussian labelling-radius sigma in bp (default 250).
#' @param weight_mode Region-weight mode. NULL inherits from
#'   \code{result$weight_mode}.
#' @return A ggplot.
#' @export
plot_nnls_residual <- function(result, tf_name, kernel_sigma = 250,
                                weight_mode = NULL) {
  weight_mode <- weight_mode %||% result$weight_mode %||% "mod_t"
  if (!tf_name %in% result$long_data$protein)
    stop(tf_name, " not in long_data")

  x_grid <- seq(-2500, 500, by = 5)
  sig <- build_caspex_signal(tf_name, result$long_data, result$pos_map,
                              x_grid, kernel_sigma, weight_mode)
  hits <- if (!is.null(result$motif_results[[tf_name]]))
    result$motif_results[[tf_name]]$hits else integer(0)
  hits <- hits[!is.na(hits) & hits >= min(x_grid) & hits <= max(x_grid)]

  reconstruction <- numeric(length(x_grid))
  if (length(hits) > 0 && requireNamespace("nnls", quietly = TRUE)) {
    X <- vapply(hits,
                function(m) exp(-0.5 * ((x_grid - m) / kernel_sigma)^2),
                numeric(length(x_grid)))
    if (length(hits) == 1) X <- matrix(X, ncol = 1)
    fit <- nnls::nnls(X, sig$y)
    reconstruction <- as.numeric(X %*% fit$x)
  }
  residual <- sig$y - reconstruction

  df <- rbind(
    data.frame(x = x_grid, y = sig$y,         lab = "Observed s(x)"),
    data.frame(x = x_grid, y = reconstruction, lab = "NNLS reconstruction X\u03b2"),
    data.frame(x = x_grid, y = residual,      lab = "Residual")
  )
  df$lab <- factor(df$lab,
                   levels = c("Observed s(x)",
                              "NNLS reconstruction X\u03b2",
                              "Residual"))
  ss_obs <- sum(sig$y^2)
  ss_res <- sum(residual^2)
  r2 <- if (ss_obs > 0) 1 - ss_res / ss_obs else NA_real_

  ggplot(df, aes(x = x, y = y)) +
    geom_area(aes(fill = lab), alpha = 0.4, color = NA) +
    geom_line(aes(color = lab), linewidth = 0.5) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = COLS$tss, linewidth = 0.6) +
    facet_wrap(~ lab, ncol = 1, scales = "free_y") +
    scale_color_manual(values = c("Observed s(x)"                = COLS$guide,
                                   "NNLS reconstruction X\u03b2" = COLS$low,
                                   "Residual"                    = COLS$high),
                       guide = "none") +
    scale_fill_manual(values = c("Observed s(x)"                = COLS$guide,
                                  "NNLS reconstruction X\u03b2" = COLS$low,
                                  "Residual"                    = COLS$high),
                      guide = "none") +
    scale_x_continuous(labels = scales::comma) +
    labs(x = "Position (bp, TSS-relative)", y = "Signal (a.u.)",
         title = paste0(tf_name, " \u2014 NNLS fit quality"),
         subtitle = sprintf("%d motif candidates | fit R\u00b2 = %.3f | \u03c3 = %d bp",
                            length(hits), r2, kernel_sigma)) +
    theme_caspex() +
    theme(strip.text = element_text(face = "bold"))
}

# =============================================================================
# C.1  Volcano plots per region
# =============================================================================

#' One volcano per region, TFs highlighted.
#' @param result (see function body).
#' @param pval_thresh (see function body).
#' @param lfc_thresh (see function body).
#' @param label_top (see function body).
#' @export
plot_volcano_per_region <- function(result, pval_thresh = 0.05,
                                     lfc_thresh = 1, label_top = 8) {
  ld <- result$long_data
  ld$logp <- -log10(pmax(ld$pval, 1e-16))
  ld$sig  <- ld$pval <= pval_thresh & abs(ld$lfc) >= lfc_thresh
  ld$kind <- ifelse(!ld$sig, "ns",
                    ifelse(ld$isTF, "TF (sig)", "non-TF (sig)"))
  ld$kind <- factor(ld$kind, levels = c("ns", "non-TF (sig)", "TF (sig)"))

  # Label top-N TFs per region by -log10(p)
  ld$label <- ""
  for (r in unique(ld$region)) {
    sub <- ld[ld$region == r & ld$isTF & ld$sig, ]
    sub <- sub[order(sub$logp, decreasing = TRUE), ]
    top <- head(sub$protein, label_top)
    ld$label[ld$region == r & ld$protein %in% top & ld$isTF & ld$sig] <-
      as.character(ld$protein[ld$region == r & ld$protein %in% top &
                                ld$isTF & ld$sig])
  }

  ggplot(ld, aes(x = lfc, y = logp)) +
    geom_point(aes(color = kind, size = kind, alpha = kind)) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
               linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_hline(yintercept = -log10(pval_thresh),
               linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_text(aes(label = label), vjust = -0.8, size = 2.6,
              color = COLS$neutral, fontface = "bold",
              check_overlap = TRUE) +
    facet_wrap(~ region, scales = "free") +
    scale_color_manual(values = c(ns = "grey80",
                                   `non-TF (sig)` = COLS$mid,
                                   `TF (sig)`     = COLS$high),
                       name = NULL) +
    scale_size_manual(values = c(ns = 0.6, `non-TF (sig)` = 1.1,
                                  `TF (sig)` = 1.6),
                      guide = "none") +
    scale_alpha_manual(values = c(ns = 0.3, `non-TF (sig)` = 0.75,
                                   `TF (sig)` = 0.95),
                      guide = "none") +
    labs(x = "logFC", y = "-log10(p)",
         title = "Per-region volcano plots",
         subtitle = sprintf("dashed lines: |logFC| \u2265 %g and p \u2264 %g",
                            lfc_thresh, pval_thresh)) +
    theme_caspex() +
    theme(legend.position = "bottom")
}

# =============================================================================
# C.2  Region correlation heatmap
# =============================================================================

#' Pearson correlation heatmap of region logFC profiles.
#' @param result (see function body).
#' @param method (see function body).
#' @export
plot_region_correlation <- function(result, method = "pearson") {
  ld <- result$long_data
  wide <- reshape(ld[, c("protein", "region", "lfc")],
                  idvar = "protein", timevar = "region",
                  direction = "wide")
  mat <- as.matrix(wide[, -1])
  colnames(mat) <- sub("^lfc\\.", "", colnames(mat))
  mat <- mat[complete.cases(mat), , drop = FALSE]
  cor_mat <- cor(mat, method = method)

  # To data.frame
  df <- as.data.frame(as.table(cor_mat))
  names(df) <- c("R1", "R2", "r")

  ggplot(df, aes(x = R1, y = R2, fill = r)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", r)),
              size = 3, color = "grey20") +
    scale_fill_gradient2(low = COLS$low, mid = "white", high = COLS$high,
                         midpoint = 0, limits = c(-1, 1),
                         name = method) +
    labs(x = NULL, y = NULL,
         title = "Inter-region logFC correlation",
         subtitle = paste0(nrow(mat),
                           " proteins with complete observations | ",
                           method)) +
    coord_fixed() +
    theme_caspex() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
}

# =============================================================================
# C.3  P-value distribution per region
# =============================================================================

#' Faceted p-value histogram — should be roughly uniform with a peak near 0.
#' @param result (see function body).
#' @param binwidth (see function body).
#' @export
plot_pval_histograms <- function(result, binwidth = 0.02) {
  ld <- result$long_data
  ggplot(ld, aes(x = pval)) +
    geom_histogram(binwidth = binwidth, fill = COLS$guide, color = "white",
                   alpha = 0.85) +
    geom_vline(xintercept = 0.05, linetype = "dashed",
               color = COLS$high, linewidth = 0.4) +
    facet_wrap(~ region) +
    labs(x = "p-value", y = "# proteins",
         title = "Per-region p-value distribution (QC)",
         subtitle = "Expected: uniform over (0.05, 1], with excess near 0") +
    theme_caspex()
}

# =============================================================================
# C.4  Motif-strength vs NNLS-weight
# =============================================================================

#' Scatter of PWM log-odds score against NNLS beta per motif hit.
#'
#' Faceted by TF. Reinforces the point that the bubble area on the main
#' binding-deconvolution plot is the NNLS coefficient \eqn{\beta}, not the
#' PWM match score: two motif hits with the same PWM strength can land at
#' very different \eqn{\beta} because NNLS is reflecting the observed
#' CasPEX signal, not the static motif quality.
#'
#' @param result Output of \code{\link{run_caspex}}.
#' @param tfs Optional character vector of TF symbols to render. NULL =
#'   every TF that has both called events and a motif scan.
#' @param kernel_sigma Gaussian labelling-radius sigma in bp (default 250).
#' @param weight_mode Region-weight mode. NULL inherits from
#'   \code{result$weight_mode}.
#' @return A ggplot.
#' @export
plot_motif_vs_nnls <- function(result, tfs = NULL,
                                kernel_sigma = 250,
                                weight_mode = NULL) {
  weight_mode <- weight_mode %||% result$weight_mode %||% "mod_t"
  if (is.null(tfs))
    tfs <- intersect(unique(result$binding_events$tf),
                      names(result$motif_results))
  if (length(tfs) == 0)
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No motif/event data"))

  x_grid <- seq(-2500, 500, by = 5)
  rows <- list()
  for (tf in tfs) {
    mr <- result$motif_results[[tf]]
    if (is.null(mr)) next
    hits <- mr$hits
    if (length(hits) == 0) next

    # Recompute PWM score per hit on the promoter sequence
    pwm <- mr$pwm$pwm; L <- mr$pwm$len
    if (is.null(pwm) || is.null(L)) next
    seq_chars <- strsplit(result$promoter_info$seq, "")[[1]]
    rev_map   <- c(A = "T", T = "A", G = "C", C = "G", N = "N")
    rev_chars <- rev(rev_map[seq_chars])
    tss_i <- result$promoter_info$tss_offset + 1
    # fwd scores indexed by 1-based window start
    fwd_scores <- score_pwm_positions(seq_chars, pwm)
    rev_scores <- score_pwm_positions(rev_chars, pwm)
    n <- length(seq_chars)

    fwd_pos_to_tssrel <- function(i) i - tss_i
    rev_pos_to_tssrel <- function(i) n - (i + L - 2) - tss_i

    # For each hit position (TSS-relative), best PWM score across strands
    pwm_score <- vapply(hits, function(h) {
      fi <- which(abs(fwd_pos_to_tssrel(seq_along(fwd_scores)) - h) <= 2)
      ri <- which(abs(rev_pos_to_tssrel(seq_along(rev_scores)) - h) <= 2)
      cands <- c(fwd_scores[fi], rev_scores[ri])
      if (length(cands) == 0) return(NA_real_)
      max(cands)
    }, numeric(1))

    # NNLS \u03b2 per hit
    if (!requireNamespace("nnls", quietly = TRUE)) next
    sig <- build_caspex_signal(tf, result$long_data, result$pos_map,
                                x_grid, kernel_sigma, weight_mode)
    valid <- hits >= min(x_grid) & hits <= max(x_grid)
    hits_v <- hits[valid]; pwm_score_v <- pwm_score[valid]
    if (length(hits_v) == 0) next
    X <- vapply(hits_v,
                function(m) exp(-0.5 * ((x_grid - m) / kernel_sigma)^2),
                numeric(length(x_grid)))
    if (length(hits_v) == 1) X <- matrix(X, ncol = 1)
    fit <- nnls::nnls(X, sig$y)

    rows[[length(rows) + 1]] <- data.frame(
      tf = tf, hit_pos = hits_v,
      pwm_score = pwm_score_v,
      beta = fit$x,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0)
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No motif/beta pairs"))
  df <- do.call(rbind, rows)
  df$survived <- df$beta > 0

  ggplot(df, aes(x = pwm_score, y = beta)) +
    geom_point(aes(color = survived), size = 2, alpha = 0.8) +
    facet_wrap(~ tf, scales = "free") +
    scale_color_manual(values = c(`TRUE` = COLS$high, `FALSE` = "grey70"),
                       labels = c(`TRUE` = "\u03b2 > 0 (kept)",
                                  `FALSE` = "\u03b2 = 0"),
                       name = NULL) +
    labs(x = "JASPAR PWM log-odds score",
         y = "NNLS coefficient \u03b2",
         title = "Motif match strength vs NNLS weight",
         subtitle = "Same PWM score, different \u03b2 \u2192 NNLS is reflecting the signal, not the PWM") +
    theme_caspex() +
    theme(legend.position = "bottom")
}

# =============================================================================
# A.5  TF motif co-occurrence matrix
# =============================================================================

#' Heat-map of TF-pair co-occurrence of predicted binding events.
#'
#' For every pair of TFs with called events, count how many events of one
#' TF fall within `tol_bp` of an event of the other TF. The resulting
#' symmetric matrix (proportion-scaled) is a candidate-cofactor screen:
#' pairs lighting up have overlapping binding footprints on this promoter.
#'
#' Mode-agnostic: works identically for default and coverage-aware events
#' because it only consumes the `tf` / `position` columns of
#' `result$binding_events`.
#' @param result (see function body).
#' @param tol_bp (see function body).
#' @param min_events (see function body).
#' @param top_peak_tfs (see function body).
#' @param top_tfs (see function body).
#' @export
plot_tf_cooccurrence <- function(result, tol_bp = 50, min_events = 1,
                                  top_peak_tfs = Inf,
                                  top_tfs = 50) {
  ev <- result$binding_events
  if (is.null(ev) || nrow(ev) == 0)
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No events"))
  # Roster construction:
  #   (a) EVERY TF with at least one motif-anchored event
  #       (motif_based = TRUE). No cap -- motif-anchored predictions
  #       carry the PWM-evidence weight the paper highlights (e.g.
  #       FOXP2-FOXP4 co-binding on the FOXP2 promoter), and
  #       capping them at top-N-by-weight previously dropped real
  #       calls whenever a few peak-driven outliers had higher
  #       summed weight.
  #   (b) Every peak-driven TF (motif_based = FALSE), ranked by
  #       summed event weight. Default `top_peak_tfs = Inf` keeps
  #       the full peak-driven roster; pass a finite integer to
  #       cap if the resulting heatmap gets too dense.
  #   (c) UNION with `result$detail_tfs` if present -- the Plot 10
  #       deck's TF roster. Guarantees roster-coherence: any TF
  #       with a per-TF detail page on Plot 10 also appears on A5.
  # `top_tfs` is retained as a legacy fallback for results that
  # pre-date the `motif_based` column (no motif_based -> roster
  # falls back to top-N by total weight).
  n_ev  <- table(ev$tf)
  with_events <- names(n_ev)[n_ev >= min_events]

  if (!is.null(ev$motif_based)) {
    motif_tfs <- unique(as.character(ev$tf[ev$motif_based %in% TRUE]))
    motif_tfs <- intersect(motif_tfs, with_events)
    ev_peak   <- ev[!(ev$motif_based %in% TRUE), , drop = FALSE]
    if (nrow(ev_peak) > 0L && top_peak_tfs > 0L) {
      peak_w <- tapply(ev_peak$weight, ev_peak$tf, sum)
      peak_w <- peak_w[!is.na(peak_w)]
      top_peak <- names(head(sort(peak_w, decreasing = TRUE), top_peak_tfs))
      top_peak <- intersect(top_peak, with_events)
    } else {
      top_peak <- character(0)
    }
    primary <- union(motif_tfs, top_peak)
  } else {
    # Legacy result with no motif_based column -- fall back to historical
    # top-N-by-summed-weight roster.
    tot_w <- tapply(ev$weight, ev$tf, sum)
    tot_w <- tot_w[with_events]
    primary <- names(head(sort(tot_w, decreasing = TRUE), top_tfs))
  }

  deck_tfs <- intersect(as.character(result$detail_tfs %||% character(0)),
                         with_events)
  keep <- union(primary, deck_tfs)

  ev   <- ev[ev$tf %in% keep, ]
  tfs  <- sort(unique(as.character(ev$tf)))
  n     <- length(tfs)
  if (n < 2)
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "Need \u2265 2 TFs with events"))

  M <- matrix(0L, n, n, dimnames = list(tfs, tfs))
  for (i in seq_len(n)) {
    pi_ <- ev$position[ev$tf == tfs[i]]
    for (j in seq_len(n)) {
      if (i == j) next
      pj_ <- ev$position[ev$tf == tfs[j]]
      # Count i-events whose nearest j-event is within tol_bp
      M[i, j] <- sum(vapply(pi_, function(p)
        any(abs(pj_ - p) <= tol_bp), logical(1)))
    }
  }
  # Symmetric proportion of i-events co-located with any j-event
  denom <- pmax(as.integer(n_ev[tfs]), 1L)
  frac  <- sweep(M, 1, denom, FUN = "/")
  # Symmetrize by max so the heat-map reads consistently in both directions
  frac_sym <- pmax(frac, t(frac))
  diag(frac_sym) <- NA

  # \u2500\u2500 Hi-C-style upright triangle layout \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
  # The symmetrised matrix has M[i,j] == M[j,i], so the upper and lower
  # halves carry identical information. Rather than rendering as a
  # half-filled square (one triangle empty), rotate the upper triangle
  # 45\u00b0 so the diagonal becomes a horizontal baseline and the triangle
  # stands upright with its apex pointing up \u2014 the convention used in
  # Hi-C contact-map browsers, which makes the "near-diagonal" (low
  # genomic-distance, here: low alphabetical-index distance) cells the
  # most prominent and lets the eye scan along the baseline to
  # individual TFs.
  #
  # Coordinate rotation: cell (i, j) with j >= i maps to
  #   (cx, cy) = ((i+j)/2, (j-i)/2)
  # so the diagonal (i == j) sits at y = 0, the row "one step off the
  # diagonal" (j == i+1) sits at y = 0.5, etc. Each unit-square cell
  # becomes a unit-side rhombus in the rotated frame, which we draw
  # with geom_polygon (4 vertices: top, right, bottom, left).
  df <- as.data.frame(as.table(frac_sym))
  names(df) <- c("tf_i", "tf_j", "frac")
  df$i_idx <- as.integer(factor(df$tf_i, levels = tfs))
  df$j_idx <- as.integer(factor(df$tf_j, levels = tfs))
  df <- df[df$j_idx >= df$i_idx, , drop = FALSE]
  df$cell_id <- seq_len(nrow(df))
  df$cx <- (df$i_idx + df$j_idx) / 2
  df$cy <- (df$j_idx - df$i_idx) / 2

  # Expand each cell into its 4 polygon vertices (rhombus). vapply'd
  # for speed on larger n; per-cell row groups via `cell_id` so
  # geom_polygon connects the right vertices.
  half <- 0.5
  poly_df <- df[rep(seq_len(nrow(df)), each = 4L), c("cell_id", "frac",
                                                      "cx", "cy"),
                 drop = FALSE]
  vx_off <- rep(c(0,  half, 0, -half), times = nrow(df))
  vy_off <- rep(c(half, 0, -half, 0), times = nrow(df))
  poly_df$x <- poly_df$cx + vx_off
  poly_df$y <- poly_df$cy + vy_off
  rownames(poly_df) <- NULL

  # \u2500\u2500 Diagonal-edge TF labels (replace bottom labels) \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
  # Each TF k corresponds to one cell on each of the three triangle
  # edges:
  #   bottom (diagonal):  (k, 0)            self-coloc (rendered black)
  #   LEFT edge:          ((1+k)/2,  (k-1)/2)  i=1 row, j=k column
  #   RIGHT edge:         ((k+n)/2,  (n-k)/2)  i=k row, j=n column
  # We label the diagonal edges (not the bottom \u2014 the user asked for
  # axis-style labels along the two upward edges of the triangle and
  # told us the bottom labels can go away if the edges look good).
  # Labels sit just outside each edge and extend PERPENDICULAR to the
  # edge (up-left from the left edge, up-right from the right edge).
  # Perpendicular orientation guarantees adjacent labels don't overlap
  # \u2014 successive labels are spaced ~0.7 units apart along the edge,
  # but their text bodies fan out perpendicular so each label has its
  # own sliver of space.
  # Per-TF source classification for label coloring:
  # "motif-anchored" if the TF has at least one event with motif_based=TRUE
  # in the filtered set, "peak-driven" otherwise. Lets the reader
  # distinguish JASPAR-evidence-backed TFs (e.g. FOXP2, FOXP4) from
  # motif-free peak-detection fallbacks at a glance on the heatmap.
  tf_is_motif <- if (!is.null(ev$motif_based)) {
    vapply(tfs, function(tf)
      any(ev$motif_based[ev$tf == tf] %in% TRUE), logical(1))
  } else {
    rep(TRUE, length(tfs))
  }
  tf_source <- ifelse(tf_is_motif, "motif-anchored", "peak-driven")

  outward_offset <- 0.45
  left_lbl <- data.frame(
    tf     = tfs,
    source = tf_source,
    x      = ((1 + seq_along(tfs)) / 2) - outward_offset * 0.707,
    y      = ((seq_along(tfs) - 1) / 2) + outward_offset * 0.707,
    angle  = -45,
    hjust  = 1,
    stringsAsFactors = FALSE)
  right_lbl <- data.frame(
    tf     = tfs,
    source = tf_source,
    x      = ((seq_along(tfs) + n) / 2) + outward_offset * 0.707,
    y      = ((n - seq_along(tfs)) / 2) + outward_offset * 0.707,
    angle  = 45,
    hjust  = 0,
    stringsAsFactors = FALSE)

  ggplot(poly_df, aes(x = x, y = y, group = cell_id, fill = frac)) +
    geom_polygon(color = "white", linewidth = 0.2) +
    # `na.value = "black"` paints the baseline (diagonal) row of
    # rhombuses \u2014 same-TF self-coloc is by definition 100% and isn't
    # a data point. Black gives the triangle a strong horizontal
    # foundation visible at the bottom of the plot.
    scale_fill_gradient(low = "white", high = COLS$high, na.value = "black",
                        labels = scales::percent, name = "co-loc",
                        limits = c(0, 1)) +
    # LEFT-edge labels: angle = -45\u00b0 + hjust = 1 means the text body
    # starts at the anchor and extends UP-LEFT (perpendicular to the
    # left edge, away from the triangle).
    geom_text(data = left_lbl,
              aes(x = x, y = y, label = tf, angle = angle,
                  hjust = hjust, color = source),
              inherit.aes = FALSE,
              size = 1.8, vjust = 0.5) +
    # RIGHT-edge labels: angle = +45\u00b0 + hjust = 0 makes text extend
    # UP-RIGHT (perpendicular to the right edge, away from the
    # triangle).
    geom_text(data = right_lbl,
              aes(x = x, y = y, label = tf, angle = angle,
                  hjust = hjust, color = source),
              inherit.aes = FALSE,
              size = 1.8, vjust = 0.5) +
    scale_color_manual(
      values = c(`motif-anchored` = "grey10",
                 `peak-driven`    = "#c66a2d"),
      name   = "TF source",
      drop   = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    # Axis-title-style annotations at the apex telling the reader
    # which diagonal direction encodes which matrix axis.
    annotate("text",
             x = (1 + n) / 2 - (n - 1) / 4,
             y = (n - 1) / 2 + 1.4,
             label = "TF_i \u2197", angle = 45,
             hjust = 0.5, vjust = 0,
             size = 3, color = "grey25", fontface = "italic") +
    annotate("text",
             x = (1 + n) / 2 + (n - 1) / 4,
             y = (n - 1) / 2 + 1.4,
             label = "\u2196 TF_j", angle = -45,
             hjust = 0.5, vjust = 0,
             size = 3, color = "grey25", fontface = "italic") +
    # x-axis ticks/labels suppressed (per-TF labels are now on the
    # diagonal edges, not the baseline).
    scale_x_continuous(breaks = NULL, labels = NULL,
                       expand = expansion(add = c(3.5, 3.5))) +
    scale_y_continuous(breaks = NULL, labels = NULL,
                       expand = expansion(add = c(0.1, 3.0))) +
    coord_fixed(clip = "off") +
    labs(x = NULL, y = NULL,
         title = "TF co-occurrence of binding events",
         subtitle = paste0("% of events within \u00b1", tol_bp,
                           " bp of another TF's event | roster = all motif-anchored TFs (n=",
                           sum(tf_is_motif), ") + ",
                           sum(!tf_is_motif), " peak-driven TFs",
                           if (length(deck_tfs) > 0)
                             paste0(" \u222a ", length(deck_tfs),
                                    " Plot 10 deck TFs") else "",
                           "  \u00b7  Hi-C-style triangle (matrix is symmetric); baseline in black = same-TF self-coloc"),
         caption = NULL) +
    theme_caspex() +
    theme(axis.text   = element_blank(),
          axis.ticks  = element_blank(),
          panel.grid  = element_blank(),
          plot.margin = margin(t = 5.5, r = 16, b = 5.5, l = 16,
                                unit = "pt"))
}

# =============================================================================
# A.6  Ranked event table (confidence-scored)
# =============================================================================

#' Rank every called binding event by a composite confidence score and
#' emit a barplot + CSV of the top-N events.
#'
#' Score:   conf = z_beta + z_surv
#'   z_beta = per-TF-normalized β   (how strong is this event for this TF)
#'   z_surv = jackknife survival fraction (B.3)           [default 1 if absent]
#'
#' The jackknife result is accepted as an optional second argument. If
#' provided, it will be merged by (tf, position) within ±tol_bp so events
#' inherit their survival fraction; otherwise survival defaults to 1 and
#' the ranking reduces to per-TF-normalized β.
#' @param result (see function body).
#' @param jk_result (see function body).
#' @param tol_bp (see function body).
#' @param top_n (see function body).
#' @export
rank_binding_events <- function(result, jk_result = NULL, tol_bp = 100,
                                 top_n = 50) {
  ev <- result$binding_events
  if (is.null(ev) || nrow(ev) == 0)
    return(invisible(NULL))
  out <- ev
  # Per-TF max-normalized \u03b2: places a TF's strongest event at 1
  tf_max <- tapply(out$weight, out$tf, max, na.rm = TRUE)
  out$beta_norm <- as.numeric(out$weight / pmax(tf_max[as.character(out$tf)], 1e-9))
  # Attach jackknife survival if we have it
  if (!is.null(jk_result) && nrow(jk_result) > 0) {
    out$surv <- vapply(seq_len(nrow(out)), function(i) {
      m <- jk_result$tf == out$tf[i] &
           abs(jk_result$position - out$position[i]) <= tol_bp
      if (!any(m)) return(NA_real_)
      mean(jk_result$frac[m], na.rm = TRUE)
    }, numeric(1))
  } else {
    out$surv <- NA_real_
  }
  # Final confidence. Treat missing survival as 1.0 so events without
  # jackknife support are not penalised relative to ones that also lack it.
  surv_fill <- ifelse(is.na(out$surv), 1, out$surv)
  out$confidence <- out$beta_norm + surv_fill
  out <- out[order(out$confidence, decreasing = TRUE), , drop = FALSE]
  top <- head(out, top_n)
  top$label <- sprintf("%s @ %+.0f", top$tf, top$position)
  top$label <- factor(top$label, levels = rev(top$label))

  p <- ggplot(top, aes(x = confidence, y = label, fill = beta_norm)) +
    geom_col(alpha = 0.9) +
    geom_text(aes(label = sprintf("\u03b2* %.2f | surv %s",
                                   beta_norm,
                                   ifelse(is.na(surv), "\u2014",
                                          sprintf("%.2f", surv)))),
              hjust = -0.05, size = 2.6, color = COLS$neutral) +
    scale_fill_gradient(low = COLS$guide, high = COLS$high, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(x = "Confidence = per-TF-normalized \u03b2 + jackknife survival",
         y = NULL,
         title = paste0("Top-", top_n, " ranked binding events"),
         subtitle = "Higher = stronger + more reproducible (coverage-aware)") +
    theme_caspex() +
    theme(panel.grid.major.y = element_blank())
  list(ranked = out, plot = p)
}

# =============================================================================
# D.1  Coverage-rescue audit scatter
# =============================================================================

#' Scatter of event \eqn{\beta} vs local coverage, coloured by gRNA distance.
#'
#' Each point is one called event from a coverage-aware run. Rescued calls
#' — the ones that depend on dividing by a small \eqn{C(x)} to survive —
#' land in the top-left (high \eqn{\beta}, low local_coverage) and are
#' coloured by how far they sit from the nearest gRNA cut site. This is
#' the single most direct sanity check on the coverage correction: it
#' makes "this call only exists because cov_floor clamped the
#' denominator" visible in one view instead of requiring a cross-reference
#' against \eqn{C(x)}.
#'
#' @param result Output of \code{\link{run_caspex}}. Must have
#'   \code{local_coverage} and \code{distance_to_nearest_grna} columns in
#'   \code{result$binding_events} (always present for coverage-aware runs).
#' @param top_label Number of top-\eqn{\beta} events to label. Default 15.
#' @return A ggplot.
#' @export
plot_coverage_rescue_scatter <- function(result, top_label = 15) {
  ev <- result$binding_events
  need <- c("local_coverage", "distance_to_nearest_grna")
  if (is.null(ev) || nrow(ev) == 0 || !all(need %in% names(ev)))
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No coverage-aware event table"))

  ev$rescued <- ev$local_coverage <= (result$cov_floor %||% 0.05) * 1.5
  # rank by \u03b2 descending to pick top labels
  ord <- order(ev$weight, decreasing = TRUE)
  ev$label <- ""
  ev$label[ord[seq_len(min(top_label, nrow(ev)))]] <-
    sprintf("%s @ %+.0f",
            ev$tf[ord[seq_len(min(top_label, nrow(ev)))]],
            ev$position[ord[seq_len(min(top_label, nrow(ev)))]])

  ggplot(ev, aes(x = local_coverage, y = weight)) +
    geom_vline(xintercept = (result$cov_floor %||% 0.05),
               linetype = "dashed", color = COLS$high, linewidth = 0.3) +
    annotate("text",
             x = (result$cov_floor %||% 0.05),
             y = Inf,
             label = " cov_floor",
             color = COLS$high, hjust = 0, vjust = 1.4, size = 3) +
    geom_point(aes(color = distance_to_nearest_grna,
                   size  = weight,
                   shape = rescued),
               alpha = 0.85, stroke = 0.3) +
    geom_text(aes(label = label), vjust = -0.9,
              size = 2.6, fontface = "bold", color = COLS$neutral,
              check_overlap = TRUE) +
    scale_color_gradient(low = COLS$guide, high = COLS$high,
                         name = "dist to gRNA (bp)") +
    scale_size_continuous(range = c(1.2, 6), guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16),
                       labels = c(`TRUE` = "near floor (rescued)",
                                  `FALSE` = "above floor"),
                       name = NULL) +
    labs(x = "local coverage C(event_pos)",
         y = "event weight \u03b2",
         title = "Coverage-rescue audit",
         subtitle = paste0("cov_floor = ", result$cov_floor %||% 0.05,
                           " | triangles sit near the floor and were ",
                           "amplified by the s/C correction")) +
    theme_caspex()
}

# =============================================================================
# D.2  cov_floor sensitivity sweep (coverage-aware mode only)
# =============================================================================

#' Re-score events at multiple cov_floor values and track event stability.
#'
#' Analog of B.2 (sigma sensitivity), but sweeping the coverage-floor
#' parameter that controls how aggressively distal / gap binders are
#' rescued. Lower floor = more distal rescues but more tail noise.
#' Robust calls drift little across the sweep; calls that appear only at
#' cov_floor ≤ 0.02 (say) should be treated as floor-sensitive.
#' @param result (see function body).
#' @param floors (see function body).
#' @param tfs (see function body).
#' @param weight_mode (see function body).
#' @param edge_guard_frac (see function body).
#' @export
run_covfloor_sensitivity <- function(result,
                                      floors = c(0.02, 0.05, 0.10, 0.20),
                                      tfs    = NULL,
                                      weight_mode = NULL,
                                      edge_guard_frac = NULL) {
  weight_mode     <- weight_mode     %||% result$weight_mode     %||% "mod_t"
  edge_guard_frac <- edge_guard_frac %||% result$edge_guard_frac %||% 0.15
  if (is.null(tfs)) {
    # Same selection as B2: every TF with motif-anchored binding events,
    # looked up against motif_results UNION motif_results_extra so the
    # augment-pool TFs (motif_scan_pool='spatial_all') are included.
    be <- result$binding_events
    motif_tfs_with_events <- if (!is.null(be$motif_based))
      unique(be$tf[be$motif_based %in% TRUE]) else unique(be$tf)
    all_motif_names <- union(names(result$motif_results),
                              names(result$motif_results_extra %||% list()))
    tfs <- intersect(motif_tfs_with_events, all_motif_names)
  }
  if (length(tfs) == 0) return(invisible(NULL))

  .get_hits <- function(tf) {
    h <- result$motif_results[[tf]]$hits
    if (is.null(h)) h <- result$motif_results_extra[[tf]]$hits
    if (is.null(h)) integer(0) else h
  }

  message("cov_floor sensitivity: ", length(floors), " floors \u00d7 ",
          length(tfs), " TFs...")
  rows <- list()
  for (tf in tfs) {
    hits <- .get_hits(tf)
    for (fl in floors) {
      ev <- predict_binding_events_coverage_aware(
        tf, result$long_data, result$pos_map, hits,
        weight_mode = weight_mode, cov_floor = fl,
        edge_guard_frac = edge_guard_frac)
      if (nrow(ev) == 0) next
      ev <- ev[order(ev$weight, decreasing = TRUE), , drop = FALSE]
      top <- head(ev, 3)
      rows[[length(rows) + 1]] <- data.frame(
        tf = tf, cov_floor = fl,
        rank = seq_len(nrow(top)),
        position = top$position,
        weight   = top$weight,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) return(invisible(NULL))
  do.call(rbind, rows)
}

#' Plot results of run_covfloor_sensitivity()
#'
#' Per-event amplitude across cov_floor values; highlights events whose
#' \eqn{\beta} is sensitive to the floor choice.
#'
#' @return A ggplot.
#' @param cf_result (see function body).
#' @export
plot_covfloor_sensitivity <- function(cf_result) {
  if (is.null(cf_result) || nrow(cf_result) == 0)
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No cov_floor-sensitivity data"))
  ggplot(cf_result, aes(x = cov_floor, y = position,
                         group = interaction(tf, rank),
                         color = factor(rank))) +
    geom_line(alpha = 0.5) +
    geom_point(aes(size = weight), alpha = 0.85) +
    facet_wrap(~ tf, ncol = 3, scales = "free_y") +
    scale_color_manual(values = c(`1` = COLS$high, `2` = COLS$mid,
                                   `3` = COLS$guide),
                       name = "event rank") +
    scale_size_continuous(range = c(1.5, 5), guide = "none") +
    scale_x_log10() +
    labs(x = "cov_floor (log scale)",
         y = "event position (bp, TSS-relative)",
         title = "cov_floor sensitivity of coverage-aware events",
         subtitle = "Flat lines = robust call | drift to the left = floor-sensitive rescue") +
    theme_caspex() +
    theme(legend.position = "bottom")
}

# =============================================================================
# D.3  s(x) / C(x) / \u03b2(x) per-TF stack (coverage-aware mode only)
# =============================================================================

#' Three-panel per-TF diagnostic: signal s(x), coverage C(x), ratio β(x).
#'
#' Analog of B.4 (NNLS residual) for coverage-aware mode. Stacking the
#' three panels with called events overlaid on the β panel shows exactly
#' what the s/C correction is doing at each event: where labeling
#' opportunity was scarce versus where β is simply tracking strong
#' enrichment.
#' @param result (see function body).
#' @param tf_name (see function body).
#' @param kernel_sigma (see function body).
#' @param weight_mode (see function body).
#' @export
plot_coverage_stack <- function(result, tf_name, kernel_sigma = NULL,
                                 weight_mode = NULL) {
  weight_mode  <- weight_mode  %||% result$weight_mode  %||% "mod_t"
  kernel_sigma <- kernel_sigma %||% result$kernel_sigma %||% 300

  if (!tf_name %in% result$long_data$protein)
    stop(tf_name, " not in long_data")

  # Use the run's actual window. Newer runs persist `upstream`/`downstream`
  # on the result; older results don't, so fall back to inferring the
  # window from the gRNA layout (pad \u00b12\u03c3 around the cut-site span so the
  # Gaussian tails fit). Last-resort fallback matches the legacy default
  # but is reached only when neither inference path yields anything.
  pos_r <- as.numeric(result$pos_map[!is.na(result$pos_map)])
  upstream   <- result$upstream
  downstream <- result$downstream
  if (is.null(upstream) || is.null(downstream)) {
    if (length(pos_r) > 0) {
      # Window has to span [TSS=0, every guide] plus 2\u03c3 padding on each
      # side so the Gaussian kernel tails fit. Wrap pos_r with 0 so the
      # window always contains the TSS itself even when all guides sit on
      # one side of it (e.g. Mackenzie FOXP2 \u2014 all guides are downstream).
      pad        <- 2L * kernel_sigma
      span_lo    <- min(c(0, pos_r))
      span_hi    <- max(c(0, pos_r))
      upstream   <- ceiling(-span_lo + pad)
      downstream <- ceiling( span_hi + pad)
    } else {
      upstream   <- 2500L
      downstream <- 500L
    }
  }
  x_grid <- seq(-upstream, downstream, by = 5)

  sig <- build_caspex_signal(tf_name, result$long_data, result$pos_map,
                              x_grid, kernel_sigma, weight_mode)
  # C(x) using the same helper the pipeline uses internally. compute_coverage
  # returns a list with $y; older fallback branch returned a bare numeric, so
  # normalise to a numeric vector here.
  cov_obj <- if (exists("compute_coverage")) {
    compute_coverage(result$pos_map, x_grid, kernel_sigma)
  } else {
    list(y = Reduce("+", lapply(result$pos_map,
      function(p) exp(-0.5 * ((x_grid - p) / kernel_sigma)^2))))
  }
  cov <- if (is.list(cov_obj)) cov_obj$y else cov_obj

  cov_floor       <- result$cov_floor       %||% 0.05
  edge_guard_frac <- result$edge_guard_frac %||% cov_floor
  floor_val       <- cov_floor * max(cov)
  beta_curve      <- sig$y / pmax(cov, floor_val)
  # Match the engine's support mask: \u03b2 is only trustworthy where C(x) is
  # comfortably above the clamp floor. Outside that region the s/C ratio
  # explodes against the floor and produces artifacts that have nothing to
  # do with binding signal. Zero those out for the plot.
  support_floor_val <- max(cov_floor, edge_guard_frac) * max(cov)
  beta_curve[cov <= support_floor_val] <- 0

  ev <- result$binding_events[result$binding_events$tf == tf_name, ,
                              drop = FALSE]

  df <- rbind(
    data.frame(x = x_grid, y = sig$y,       panel = "s(x) \u2014 signal"),
    data.frame(x = x_grid, y = cov,         panel = "C(x) \u2014 coverage"),
    data.frame(x = x_grid, y = beta_curve,  panel = "\u03b2(x) = s/C")
  )
  df$panel <- factor(df$panel,
                     levels = c("s(x) \u2014 signal",
                                "C(x) \u2014 coverage",
                                "\u03b2(x) = s/C"))
  ev_df <- if (nrow(ev) > 0)
    data.frame(x = ev$position, weight = ev$weight,
               panel = factor("\u03b2(x) = s/C",
                              levels = levels(df$panel)))
  else NULL

  p <- ggplot(df, aes(x = x, y = y)) +
    geom_area(aes(fill = panel), alpha = 0.35, color = NA) +
    geom_line(aes(color = panel), linewidth = 0.5) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = COLS$tss, linewidth = 0.6) +
    facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    scale_color_manual(values = c("s(x) \u2014 signal"     = COLS$guide,
                                   "C(x) \u2014 coverage"  = COLS$low,
                                   "\u03b2(x) = s/C"  = COLS$high),
                       guide = "none") +
    scale_fill_manual(values = c("s(x) \u2014 signal"      = COLS$guide,
                                  "C(x) \u2014 coverage"   = COLS$low,
                                  "\u03b2(x) = s/C"   = COLS$high),
                      guide = "none") +
    scale_x_continuous(labels = scales::comma) +
    labs(x = "position (bp, TSS-relative)", y = NULL,
         title = paste0(tf_name,
                        " \u2014 coverage-aware decomposition"),
         subtitle = sprintf("\u03c3 = %d bp | cov_floor = %g | %d event(s)",
                             kernel_sigma, result$cov_floor %||% 0.05,
                             nrow(ev))) +
    theme_caspex() +
    theme(strip.text = element_text(face = "bold"))
  if (!is.null(ev_df))
    p <- p + geom_point(data = ev_df,
                        aes(x = x, y = 0, size = weight),
                        inherit.aes = FALSE,
                        shape = 21, fill = COLS$high, color = "black",
                        stroke = 0.4, alpha = 0.9) +
      scale_size_continuous(range = c(2, 6), guide = "none")
  p
}

# =============================================================================
# F.  Convenience wrapper
# =============================================================================

#' Produce every extra plot as a separate PDF in `out_dir`.
#'
#' Auto-detects whether `result` came from a default-mode or coverage-aware
#'
#'   Steps A1-A6, B1-B3, C1-C3, D1-D3 (coverage-aware diagnostics).
#'
#' @param result        From run_caspex()
#' @param out_dir       Output directory (created if missing)
#' @param n_perm        Number of permutations for the null distribution
#' @param sigmas        Kernel widths for the sensitivity grid
#' @param cov_floors    cov_floor sweep used by D2 (coverage mode only)
#' @param one_pager_tfs TFs for which to emit a one-pager (default: top-10
#'                      composite TFs union every motif-scanned TF)
#' @param skip          Character vector of step names to skip. Any of:
#'                      "one_pager","family","event_density","scatter",
#'                      "cooccurrence","ranked_events",
#'                      "permutation","sigma","jackknife","residual",
#'                      "volcano","correlation","pvalhist","motif_vs_beta",
#'                      "cov_rescue","cov_floor_sweep","cov_stack"
#' @return An invisible list of all generated objects
#' @export
run_caspex_extras <- function(result,
                              out_dir        = "caspex_output/extras",
                              n_perm         = 500,
                              sigmas         = c(100, 200, 250, 300, 500),
                              cov_floors     = c(0.02, 0.05, 0.10, 0.20),
                              one_pager_tfs  = NULL,
                              skip           = character(0)) {
  if (missing(result) || is.null(result))
    stop("run_caspex_extras() needs the object returned by run_caspex(). ",
         "Call run_caspex_extras(result) \u2014 do not pass NULL.")
  for (f in c("long_data", "spatial_df", "pos_map"))
    if (is.null(result[[f]]))
      stop("result$", f, " is NULL. Did run_caspex() finish successfully?")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir_abs <- normalizePath(out_dir, mustWork = FALSE)
  mode_lbl <- paste0("coverage-aware (cov_floor=",
                     result$cov_floor %||% 0.05, ")")
  message("\n=== CasPEX extras: writing to ", out_dir_abs, " ===")
  message("    binding mode: ", mode_lbl)

  # `out` is an environment, not a list, so step expressions evaluated by
  # safely() can mutate it via `out$x <- val` (in-place reference semantics).
  # A plain list would require `<<-`, which doesn't work here: the step
  # expression is captured as a promise whose environment is run_caspex_extras's
  # eval frame, so `<<-` skips that frame and walks up to globalenv looking for
  # `out` \u2014 fails with "object 'out' not found". Environment fix sidesteps that.
  out         <- new.env(parent = emptyenv())
  step_status <- list()                   # name -> "ok" | "skipped" | "failed: ..."
  do_step     <- function(name) !(name %in% skip)

  # Wrap each step so one failure can't kill the rest. `expr` must write
  # its own file(s); on error we record the reason and continue.
  safely <- function(name, expr) {
    if (!do_step(name)) { step_status[[name]] <<- "skipped"; return(invisible(NULL)) }
    tryCatch({
      force(expr)
      step_status[[name]] <<- "ok"
    }, error = function(e) {
      step_status[[name]] <<- paste0("failed: ", conditionMessage(e))
      message("  !! step '", name, "' failed: ", conditionMessage(e))
      try(grDevices::dev.off(), silent = TRUE)  # close any half-open PDF
    })
  }

  # --- A.1 per-TF one-pagers (multi-page PDF)
  safely("one_pager", {
    if (is.null(one_pager_tfs)) {
      # Default: top-10 composite TFs \u222a every TF that was motif-scanned.
      # `motif_results` already contains the full union of the three
      # selection buckets (common + shared-focal + region-specific), so
      # this automatically produces a one-pager for every TF that shows
      # up on any deck in the main pipeline.
      top  <- head(as.character(result$spatial_df$protein), 10)
      scan <- if (!is.null(result$motif_results))
        names(result$motif_results) else character(0)
      one_pager_tfs <- unique(c(top, scan))
    }
    one_pager_tfs <- intersect(one_pager_tfs, result$long_data$protein)
    message("  [A1] one-pagers for ", length(one_pager_tfs), " TFs")
    pdf(file.path(out_dir, "A1_tf_one_pagers.pdf"), width = 11, height = 9)
    for (tf in one_pager_tfs) {
      p <- tryCatch(plot_tf_one_pager(result, tf),
                    error = function(e) {
                      message("    skipping ", tf, ": ", conditionMessage(e))
                      NULL
                    })
      if (!is.null(p)) print(p)
    }
    dev.off()
    out$one_pager_tfs <- one_pager_tfs
  })

  # --- A.2 TF-family enrichment (may take time; queries JASPAR)
  safely("family", {
    message("  [A2] TF-family enrichment (JASPAR queries)")
    p <- plot_tf_family_enrichment(result)
    ggsave(file.path(out_dir, "A2_tf_family.pdf"), p, width = 10, height = 5.5)
    out$family <- p
  })

  # --- A.3 event density
  safely("event_density", {
    message("  [A3] event density")
    p <- plot_event_density(result)
    ggsave(file.path(out_dir, "A3_event_density.pdf"), p, width = 10, height = 4)
    out$event_density <- p
  })

  # --- A.4 composite vs specificity
  safely("scatter", {
    message("  [A4] composite vs specificity")
    p <- plot_composite_vs_specificity(result)
    ggsave(file.path(out_dir, "A4_composite_vs_specificity.pdf"),
           p, width = 8.5, height = 6.5)
    out$scatter <- p
  })

  # --- A.5 TF-pair co-occurrence heatmap (mode-agnostic)
  safely("cooccurrence", {
    message("  [A5] TF-pair co-occurrence")
    p <- plot_tf_cooccurrence(result)
    ggsave(file.path(out_dir, "A5_tf_cooccurrence.pdf"),
           p, width = 9, height = 8)
    out$cooccurrence <- p
  })

  # --- B.1 permutation null
  safely("permutation", {
    message("  [B1] permutation null (B=", n_perm, ")")
    perm <- run_permutation_null(result, n_perm = n_perm)
    p <- plot_permutation_null(perm)
    ggsave(file.path(out_dir, "B1_permutation_null.pdf"), p, width = 9, height = 7)
    write.csv(perm$summary,
              file.path(out_dir, "B1_permutation_null_summary.csv"),
              row.names = FALSE)
    out$permutation <- perm
  })

  # --- B.2 sigma sensitivity
  safely("sigma", {
    message("  [B2] sigma sensitivity")
    sg <- run_sigma_sensitivity(result, sigmas = sigmas)
    if (!is.null(sg)) {
      out$sigma <- sg
      .save_sensitivity_paginated(
        sg, plot_sigma_sensitivity,
        out_path     = file.path(out_dir, "B2_sigma_sensitivity.pdf"),
        tfs_per_page = 12)
      write.csv(sg,
                file.path(out_dir, "B2_sigma_sensitivity.csv"),
                row.names = FALSE)
    }
  })

  # --- B.3 event jackknife
  safely("jackknife", {
    message("  [B3] event jackknife")
    jk <- run_event_jackknife(result)
    if (!is.null(jk)) {
      p <- plot_event_jackknife(jk)
      ggsave(file.path(out_dir, "B3_event_jackknife.pdf"),
             p, width = 8, height = 10)
      write.csv(jk, file.path(out_dir, "B3_event_jackknife.csv"),
                row.names = FALSE)
      out$jackknife <- jk
    }
  })

  # --- A.6 ranked events \u2014 consumes B.3 if present. Placed AFTER jackknife
  # so the confidence score can include the survival fraction; if jackknife
  # failed or was skipped, the score falls back to per-TF-normalized beta.
  safely("ranked_events", {
    message("  [A6] ranked event table")
    rk <- rank_binding_events(result, jk_result = out$jackknife)
    if (!is.null(rk)) {
      ggsave(file.path(out_dir, "A6_ranked_events.pdf"),
             rk$plot, width = 9, height = 10)
      write.csv(rk$ranked, file.path(out_dir, "A6_ranked_events.csv"),
                row.names = FALSE)
      out$ranked_events <- rk
    }
  })

  # --- B.4 NNLS residual \u2014 one page per TF with motif hits
  safely("residual", {
    tfs_r <- intersect(unique(result$binding_events$tf),
                        names(result$motif_results))
    if (length(tfs_r) > 0) {
      message("  [B4] NNLS residual for ", length(tfs_r), " TFs")
      pdf(file.path(out_dir, "B4_nnls_residual.pdf"), width = 9, height = 7)
      for (tf in tfs_r) {
        p <- tryCatch(plot_nnls_residual(result, tf),
                      error = function(e) NULL)
        if (!is.null(p)) print(p)
      }
      dev.off()
      out$residual_tfs <- tfs_r
    } else {
      message("  [B4] skipped: no TFs with both events and motif hits")
    }
  })

  # --- C.1 volcano per region
  safely("volcano", {
    message("  [C1] volcano per region")
    p <- plot_volcano_per_region(result)
    ggsave(file.path(out_dir, "C1_volcano_per_region.pdf"),
           p, width = 12, height = 8)
    out$volcano <- p
  })

  # --- C.2 region correlation
  safely("correlation", {
    message("  [C2] region correlation")
    p <- plot_region_correlation(result)
    ggsave(file.path(out_dir, "C2_region_correlation.pdf"),
           p, width = 7, height = 6)
    out$correlation <- p
  })

  # --- C.3 p-value histograms
  safely("pvalhist", {
    message("  [C3] p-value histograms")
    p <- plot_pval_histograms(result)
    ggsave(file.path(out_dir, "C3_pval_histograms.pdf"),
           p, width = 11, height = 7)
    out$pvalhist <- p
  })

  # --- C.4 motif score vs NNLS beta
  safely("motif_vs_beta", {
    message("  [C4] motif score vs NNLS beta")
    p <- plot_motif_vs_nnls(result)
    ggsave(file.path(out_dir, "C4_motif_vs_nnls_beta.pdf"),
           p, width = 12, height = 9)
    out$motif_vs_beta <- p
  })

  # --- D.1 coverage-rescue audit scatter
  safely("cov_rescue", {
    message("  [D1] coverage-rescue audit scatter")
    p <- plot_coverage_rescue_scatter(result)
    ggsave(file.path(out_dir, "D1_coverage_rescue.pdf"),
           p, width = 9, height = 7)
    out$cov_rescue <- p
  })

  # --- D.2 cov_floor sensitivity sweep
  safely("cov_floor_sweep", {
    message("  [D2] cov_floor sensitivity (floors: ",
            paste(cov_floors, collapse = ", "), ")")
    cf <- run_covfloor_sensitivity(result, floors = cov_floors)
    if (!is.null(cf)) {
      out$cov_floor_sweep <- cf
      .save_sensitivity_paginated(
        cf, plot_covfloor_sensitivity,
        out_path     = file.path(out_dir, "D2_covfloor_sensitivity.pdf"),
        tfs_per_page = 12)
      write.csv(cf, file.path(out_dir, "D2_covfloor_sensitivity.csv"),
                row.names = FALSE)
    }
  })

  # --- D.3 s(x) / C(x) / beta(x) per-TF stack (auto-skipped for default mode)
  safely("cov_stack", {
    tfs_s <- intersect(unique(result$binding_events$tf),
                        names(result$motif_results))
    if (length(tfs_s) > 0) {
      message("  [D3] s/C/beta stack for ", length(tfs_s), " TFs")
      pdf(file.path(out_dir, "D3_coverage_stack.pdf"), width = 9, height = 9)
      for (tf in tfs_s) {
        p <- tryCatch(plot_coverage_stack(result, tf),
                      error = function(e) {
                        message("    skipping ", tf, ": ",
                                conditionMessage(e))
                        NULL
                      })
        if (!is.null(p)) print(p)
      }
      dev.off()
      out$cov_stack_tfs <- tfs_s
    } else {
      message("  [D3] skipped: no TFs with both events and motif hits")
    }
  })

  # Final tally \u2014 what worked, what didn't, what was skipped
  message("\n--- Extras summary -----------------------------------------------")
  for (nm in names(step_status))
    message(sprintf("  %-14s : %s", nm, step_status[[nm]]))
  files_written <- list.files(out_dir, full.names = FALSE)
  message("  files written  : ", length(files_written))
  if (length(files_written) > 0)
    message("    ", paste(files_written, collapse = ", "))
  else
    message("    (none \u2014 check the errors above)")
  message("------------------------------------------------------------------\n")

  # Convert the env-backed `out` to a list before returning. Downstream code
  # (and existing call sites) expect list semantics, e.g. `extras$one_pager_tfs`.
  out_list                <- as.list(out)
  out_list$step_status    <- step_status
  out_list$files_written  <- files_written
  out_list$out_dir        <- out_dir_abs
  invisible(out_list)
}

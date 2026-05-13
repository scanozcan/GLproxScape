# =============================================================================
# caspex_epigenetic.R
#
# Epigenetic-factor binding-zone analysis. Parallel module to the TF binding-
# event pipeline in caspex_analysis.R, but tailored to chromatin-associated
# factors that lack JASPAR motifs:
#
#   * No motif scanning, no per-motif NNLS, no motif anchoring.
#   * Bubble-style point calls are replaced with horizontal "binding zones"
#     (contiguous bp ranges where β(x) > zone_frac · max(β)). This honestly
#     represents what the proximity-labelling kernel resolves at σ = 300 bp
#     for spreading marks (PRC2, BRD4, KDM domains, etc.) — a domain, not a
#     position. Compare to plot_binding_deconvolution()'s point-bubble lane.
#   * ChIP-Atlas validation is the primary lane (not supplementary), since
#     for chromatin factors public ChIP-seq is the directly comparable
#     readout.
#
# This module piggybacks on caspex_analysis.R's signal-building primitives
# (build_caspex_signal, compute_coverage, compute_region_weight, COLS,
# theme_caspex) and on caspex_chipatlas.R's run_chipatlas_scan(). Both must
# already be sourced into the global env before this file is loaded.
#
# Sourced explicitly from runner scripts (e.g. myers_2018_reanalysis/1-myers.R)
# alongside the other caspex_*.R modules. Not auto-loaded by caspex_analysis.R
# — keeps the TF deck path independent of the epigenetic deck path.
# =============================================================================

# `%||%` is defined once in R/utils-internal.R and visible to every R/
# file in the package namespace.

# =============================================================================
# SECTION 1: Per-factor zone detection
# =============================================================================

#' Predict binding zones for an epigenetic factor.
#'
#' Parallel to predict_binding_events_coverage_aware() but emits ZONES rather
#' than point events. For each contiguous run of β(x) above
#' `zone_frac` · max(β), report the bp interval, peak position inside the
#' zone, peak β, distance from peak to nearest gRNA, and the number of
#' detected regions whose logFC is positive within ±2σ of the zone (a loose
#' "supporting region count" diagnostic — not used as a filter, just for
#' reading the deck).
#'
#' Why zones, not bubbles. Chromatin-associated factors typically don't have
#' a sequence-specific anchor point; their binding "spreads" across a domain
#' (H3K27me3 over Polycomb domains, BRD4 over superenhancers, KDM5 over
#' H3K4me3 plateaus). Reporting a single bubble bp coordinate at the β-peak
#' overstates spatial precision. The zone bar honestly says "the proximity
#' signal is consistent with binding anywhere within this bp range."
#'
#' Why σ stays at 300 bp. The biotinylation radius of APEX2-bound dCas9 is a
#' physical property of the labelling reaction, not of the factor's
#' biological function. Increasing σ for spreading marks would be a
#' sleight-of-hand: pretending the *measurement* resolves a wider region
#' when in fact it doesn't. Instead we keep σ = 300 bp (matching the TF
#' path) and lower zone_frac to 0.3 (vs. 0.5 in the TF zone detector). That
#' preserves the kernel-honest signal but lets neighbouring near-peak β
#' values merge into broader domains, which is the correct way to
#' communicate spreading-mark biology under a fixed labelling radius.
#'
#' @return data.frame with columns: tf, zone_start, zone_end, zone_width,
#'   peak_position, peak_beta, distance_to_nearest_grna, n_regions_supporting.
#' @param tf_name (see function body).
#' @param long_data (see function body).
#' @param pos_map (see function body).
#' @param x_grid (see function body).
#' @param kernel_sigma (see function body).
#' @param zone_frac (see function body).
#' @param inner_zone_frac (see function body).
#' @param centroid_frac (see function body).
#' @param weight_mode (see function body).
#' @param cov_floor (see function body).
#' @param edge_guard_frac (see function body).
#' @param max_grna_distance (see function body).
#' @param edge_grna_weight_cap (see function body).
#' @noRd
predict_binding_zones_epigenetic <- function(
    tf_name, long_data, pos_map,
    x_grid          = seq(-2500, 500, by = 5),
    kernel_sigma    = 300,
    zone_frac       = 0.3,
    inner_zone_frac = 0.7,
    centroid_frac   = 0.7,
    weight_mode     = "lfc_signed",
    cov_floor       = 0.05,
    edge_guard_frac = 0.25,
    max_grna_distance    = NULL,
    edge_grna_weight_cap = NULL) {
  # `centroid_frac`: fraction of the zone's max regional logFC that a region
  # must reach to count as a centroid in addition to being a local maximum
  # among the regions inside the zone. The diamond(s) on the deck mark
  # *region-specific* centroids — i.e. the gRNA position(s) where the per-
  # region logFC is locally highest within the zone, not the kernel-summed
  # s(x) argmax (which biases toward multi-guide overlap zones even when
  # one region has the standout signal). For a monotonically-decreasing
  # zone (one strong region surrounded by weaker ones, like GATAD2A's R5)
  # this gives one centroid at the strongest region. For a two-peak
  # pattern within one zone (R5=high, R3=low, R2=high), both R5 and R2
  # qualify and get their own diamonds. centroid_frac = 0.7 keeps the
  # threshold reasonably tight: a secondary peak must be at least 70% of
  # the primary peak's logFC to be visualised separately, otherwise it's
  # treated as a flank of the primary.
  # `inner_zone_frac` defines a tighter "focal core" inside each outer zone.
  # The outer zone (zone_frac, default 0.3) captures the bp range where
  # β > 0.3 · max(β) — this is intentionally wide (~3.1σ ≈ 930 bp at σ=300
  # for a single-source signal) so it honestly represents the chromatin-
  # domain extent the labelling kernel can resolve. Inside each zone we
  # additionally walk outward from the zone's peak position until
  # β drops below inner_zone_frac · max(β within zone), and report that
  # contiguous run as the focal core (default 0.7 → ~1.7σ ≈ 510 bp single-
  # source baseline). Visually: light outer bar for broad domain, dark
  # inner bar for focal core. Avoids forcing the reader to choose between
  # "broad-domain biology" and "focal peak resolution" framings.
  empty <- data.frame(
    tf = character(), zone_start = numeric(), zone_end = numeric(),
    zone_width = numeric(), peak_position = numeric(),
    peak_beta = numeric(),
    # Inner-core schema. Single-value columns (core_start, core_end,
    # core_width) carry the FIRST core for back-compat with consumers
    # written against the single-core era; multi-core info lives in the
    # comma-separated string columns (core_starts, core_ends, core_widths)
    # plus n_cores. A zone with one core has 1-element strings; a zone
    # with multiple sub-peaks above the inner threshold has multi-element
    # strings, e.g. "-1900,-300" / "-1100,500".
    core_start = numeric(), core_end = numeric(), core_width = numeric(),
    core_starts = character(), core_ends = character(),
    core_widths = character(), n_cores = integer(),
    centroid_positions = character(),
    centroid_lfcs      = character(),
    n_centroids        = integer(),
    distance_to_nearest_grna = numeric(),
    n_regions_supporting = integer(),
    stringsAsFactors = FALSE)
  if (is.null(max_grna_distance)) max_grna_distance <- kernel_sigma

  # Reuse the TF-path signal building so we share the exact same
  # weight_mode plumbing, support_mask logic, and edge-guard semantics.
  sig    <- build_caspex_signal(tf_name, long_data, pos_map, x_grid,
                                 kernel_sigma, weight_mode)
  s_grid <- sig$y
  if (max(s_grid) <= 0) return(empty)
  cov_obj <- compute_coverage(pos_map, x_grid, kernel_sigma)
  c_grid  <- cov_obj$y
  if (max(c_grid) <= 0) return(empty)

  pos_r <- sort(as.numeric(pos_map[!is.na(pos_map)]))
  support_floor_val <- max(cov_floor, edge_guard_frac) * max(c_grid)
  support_mask      <- c_grid > support_floor_val

  # Per-region weights for this TF (used by both region-specific centroid
  # detection and the n_regions_supporting diagnostic). Computed once
  # here so both downstream blocks see the same filtered table — earlier
  # versions duplicated the subset in two places, which is fine but
  # error-prone if the filtering rule ever drifts.
  detected <- long_data[long_data$protein == tf_name, , drop = FALSE]
  detected$pos <- as.numeric(pos_map[detected$region])
  detected <- detected[!is.na(detected$pos) & !is.na(detected$lfc), ,
                       drop = FALSE]

  # ── Zone detection runs on s(x), NOT β = s(x)/C(x) ─────────────────────
  # This is the deliberate departure from the TF path. For the TF deck,
  # β = s/C is the right object to threshold: it answers "where is per-
  # guide enrichment density highest?" — which is what motif-anchored
  # scoring needs. But for broad chromatin factors with relatively
  # uniform per-region logFCs (e.g. KDM2A on hTERT with logFC≈0.2–0.4
  # across all 5 regions), β amplifies at lone-guide edges: at R1 alone
  # at +1100, β ≈ R1's logFC because numerator and denominator both
  # collapse to a single guide's contribution; in the multi-guide
  # interior at ~0 bp, β is the kernel-weighted average of all logFCs,
  # which is *lower* than any single one. The upshot is the β-peak
  # lands at the lone-guide edge even though the actual labelling
  # intensity peaks in the interior — so the zone bar and the s(x)
  # curve drawn in the upper panel disagree, which is misleading.
  #
  # s(x) directly answers "where is labelling intensity high?" — which
  # is what proximity labelling fundamentally measures, and which is
  # what the upper-panel curve already shows. Using it for zone detection
  # makes the bar visually consistent with the curve. The support mask
  # is preserved so events outside the kernel-resolved trust region are
  # still zeroed; edge artefacts past the gRNA tile are unaffected.
  signal <- s_grid
  signal[!support_mask] <- 0
  signal_max <- max(signal)
  if (signal_max <= 0) return(empty)

  threshold <- zone_frac * signal_max
  above     <- signal > threshold
  if (!any(above)) return(empty)

  rle_a   <- rle(above)
  ends    <- cumsum(rle_a$lengths)
  starts  <- c(1L, head(ends + 1L, -1L))
  zone_s  <- starts[rle_a$values]
  zone_e  <- ends[rle_a$values]

  zones <- data.frame(
    tf         = tf_name,
    zone_start = x_grid[zone_s],
    zone_end   = x_grid[zone_e],
    stringsAsFactors = FALSE)
  zones$zone_width <- zones$zone_end - zones$zone_start
  # `peak_beta` (column name retained for downstream-CSV continuity)
  # stores the kernel-summed peak — max(s(x)) inside the zone — and is
  # what drives the inner-core fill colour gradient. This is a property
  # of the smooth signal curve.
  zones$peak_beta  <- mapply(function(s, e) max(signal[s:e]),
                              zone_s, zone_e)

  # ── Region-specific centroids ──────────────────────────────────────────
  # The diamond(s) drawn on the deck mark *calculated* bp positions where
  # the per-region logFC pattern peaks. Each centroid is identified by
  # finding regions inside the zone whose logFC is a local maximum among
  # in-zone regions AND is at least centroid_frac × max(zone_lfcs). For
  # each such centroid REGION we then compute a CALCULATED bp coordinate
  # by taking the logFC-weighted average of that region and its
  # immediate positive-logFC neighbours in the zone. The diamond sits at
  # this calculated bp, NOT at the gRNA position — analogous in spirit
  # to the no-motif TF path emitting bubble positions at find_local_
  # maxima(s(x)) bp coordinates rather than snapping to gRNAs.
  #
  # This addresses the bias problems in both prior alternatives: it does
  # NOT use β = s/C (which puts the argmax at the outermost guide when
  # logFCs are uniform — KDM2A failure mode), and it does NOT use s(x)
  # argmax alone (which puts the argmax in the multi-guide kernel-
  # overlap zone even when one region clearly dominates — GATAD2A
  # failure mode). The weighted-average-of-neighbours is data-driven:
  # the centroid is pulled toward whichever side carries more positive
  # signal, equals the gRNA coordinate only when both immediate
  # neighbours are zero-or-negative, and produces interpretable
  # multi-centroid output for two-peak zones.
  #
  # Examples:
  #   GATAD2A (R5=0.85, R4=0.50, R3=0.42, R2=0.47, R1≈0): single
  #   centroid at R5 (only one passing 0.7 × 0.85 threshold), calculated
  #   position = (0.85·(−380) + 0.50·(−50)) / 1.35 = −258. Diamond at
  #   −258 (between R5 and R4, weighted toward R5).
  #
  #   Hypothetical two-peak (R5=0.85, R4=0.20, R3=0.10, R2=0.80):
  #   two centroids at R5 and R2. R5 calculated = (0.85·(−380) +
  #   0.20·(−50)) / 1.05 = −317. R2 calculated = (0.10·200 + 0.80·360)
  #   / 0.90 = +342. Two diamonds at −317 and +342.
  #
  # The `centroid_positions` / `centroid_lfcs` columns are comma-
  # separated strings so the deck plot can parse them and emit one
  # diamond per centroid; `peak_position` records the *primary*
  # (highest-logFC) centroid's calculated bp as a single value for
  # back-compat / scalar-column consumers.
  detected_pos <- detected[detected$lfc > 0, , drop = FALSE]
  detected_pos <- detected_pos[order(detected_pos$pos), , drop = FALSE]
  zones$peak_position      <- NA_real_
  zones$centroid_positions <- NA_character_
  zones$centroid_lfcs      <- NA_character_
  zones$n_centroids        <- 0L
  for (zi in seq_len(nrow(zones))) {
    zs <- zones$zone_start[zi]
    ze <- zones$zone_end[zi]
    in_zone <- detected_pos$pos >= zs & detected_pos$pos <= ze
    zr <- detected_pos[in_zone, , drop = FALSE]
    if (nrow(zr) == 0) {
      # No positive-logFC region inside the zone — fall back to the
      # s(x) argmax position so we still emit a meaningful marker.
      zones$peak_position[zi] <- {
        s_idx <- zone_s[zi]
        e_idx <- zone_e[zi]
        x_grid[s_idx + which.max(signal[s_idx:e_idx]) - 1L]
      }
      zones$centroid_positions[zi] <- as.character(zones$peak_position[zi])
      zones$centroid_lfcs[zi]      <- NA_character_
      zones$n_centroids[zi]        <- 1L
      next
    }
    zr <- zr[order(zr$pos), , drop = FALSE]
    nz <- nrow(zr)
    zone_max_lfc <- max(zr$lfc)
    # Local-maximum test: a region is a local max if its logFC is at
    # least as large as both immediate neighbours within the zone (with
    # endpoints counted as one-sided). For a single-region zone all
    # rows trivially qualify.
    is_lmax <- vapply(seq_len(nz), function(i) {
      left_ok  <- (i == 1L)  || zr$lfc[i] >= zr$lfc[i - 1L]
      right_ok <- (i == nz)  || zr$lfc[i] >= zr$lfc[i + 1L]
      left_ok && right_ok
    }, logical(1))
    is_centroid <- is_lmax & (zr$lfc >= centroid_frac * zone_max_lfc)
    if (!any(is_centroid)) {
      # Defensive: should never trigger because the zone-max is itself a
      # local max and is by construction >= centroid_frac × itself, but
      # guard against floating-point edge cases by promoting whichever
      # row has the highest logFC.
      is_centroid <- seq_len(nz) == which.max(zr$lfc)
    }
    # Calculated centroid bp position. For each centroid REGION, compute a
    # logFC-weighted average of (centroid + immediate left/right positive-
    # logFC neighbours in zone) so the diamond lands at a *calculated*
    # bp coordinate instead of snapping to the gRNA position. Pulled
    # toward whichever neighbour carries more signal; equals the gRNA
    # coordinate only when neighbours have zero positive contribution.
    # See the no-motif TF path (predict_binding_events_coverage_aware →
    # find_local_maxima on s(x)) — same intent (interpolated bp position
    # informed by neighbouring data), but here we use a per-region
    # weighted average instead of a curve local-max because the goal is
    # to centre on region-specific signal, not on kernel-summation peaks.
    centroid_idx <- which(is_centroid)
    c_pos <- vapply(centroid_idx, function(i) {
      nbhd_idx <- unique(c(if (i > 1L) i - 1L, i,
                           if (i < nz) i + 1L))
      if (length(nbhd_idx) == 0L) return(zr$pos[i])
      weighted.mean(zr$pos[nbhd_idx], zr$lfc[nbhd_idx])
    }, numeric(1))
    # Centroid LOGFC stays as the region's own logFC (not the weighted
    # average) — this is the value used for primary-centroid ranking and
    # for the centroid_lfcs CSV column, where the user should see the
    # actual measurement of the centroid region, not a derived statistic.
    c_lfc <- zr$lfc[centroid_idx]
    primary <- which.max(c_lfc)
    zones$peak_position[zi]      <- c_pos[primary]
    zones$centroid_positions[zi] <- paste(round(c_pos, 1), collapse = ",")
    zones$centroid_lfcs[zi]      <- paste(sprintf("%.3f", c_lfc),
                                           collapse = ",")
    zones$n_centroids[zi]        <- length(c_pos)
  }
  # Focal-core bounds within each zone. For each zone, find ALL contiguous
  # runs of s(x) > inner_zone_frac · zone_max — one inner core per
  # suprathreshold sub-peak. A zone with a single sharp peak gets one
  # core; a zone with two well-separated sub-peaks (separated by a valley
  # that dips below the inner threshold) gets two cores; etc. The
  # earlier single-core walker started at the zone's global max and
  # walked outward until the first valley, which missed any secondary
  # sub-peaks above threshold (the SSRP1 case, where R5/R6/R7 cluster
  # and R3/R2 cluster both pass 0.70·max but the valley between them
  # at ~−700 dips below). The core_starts / core_ends / core_widths
  # columns are comma-separated bp strings so the deck plot can parse
  # them and emit one geom_rect per core; core_start / core_end /
  # core_width still carry the FIRST core as scalars for back-compat.
  core_results <- lapply(seq_along(zone_s), function(zi) {
    s <- zone_s[zi]; e <- zone_e[zi]
    zb        <- signal[s:e]
    zone_max  <- max(zb)
    inner_thr <- inner_zone_frac * zone_max
    above_inner <- zb > inner_thr
    if (!any(above_inner)) {
      return(list(starts = NA_character_, ends = NA_character_,
                  widths = NA_character_, n = 0L,
                  first_start = NA_real_, first_end = NA_real_,
                  first_width = NA_real_))
    }
    # rle() over the boolean mask gives all contiguous TRUE runs; we
    # only keep runs where rle$values is TRUE.
    rle_a       <- rle(above_inner)
    ends_idx    <- cumsum(rle_a$lengths)
    starts_idx  <- c(1L, head(ends_idx + 1L, -1L))
    sub_s       <- starts_idx[rle_a$values]
    sub_e       <- ends_idx[rle_a$values]
    core_xs     <- x_grid[s + sub_s - 1L]
    core_xe     <- x_grid[s + sub_e - 1L]
    core_widths <- core_xe - core_xs
    list(
      starts      = paste(round(core_xs,     1), collapse = ","),
      ends        = paste(round(core_xe,     1), collapse = ","),
      widths      = paste(round(core_widths, 1), collapse = ","),
      n           = length(core_xs),
      first_start = core_xs[1],
      first_end   = core_xe[1],
      first_width = core_widths[1])
  })
  zones$core_starts <- vapply(core_results, function(x) x$starts, character(1))
  zones$core_ends   <- vapply(core_results, function(x) x$ends,   character(1))
  zones$core_widths <- vapply(core_results, function(x) x$widths, character(1))
  zones$n_cores     <- vapply(core_results, function(x) x$n,      integer(1))
  zones$core_start  <- vapply(core_results, function(x) x$first_start, numeric(1))
  zones$core_end    <- vapply(core_results, function(x) x$first_end,   numeric(1))
  zones$core_width  <- vapply(core_results, function(x) x$first_width, numeric(1))
  zones$distance_to_nearest_grna <- vapply(zones$peak_position,
    function(p) min(abs(p - pos_r)), numeric(1))

  # Supporting-regions count (diagnostic, not a filter). Counts detected
  # regions whose logFC is positive AND whose position falls within ±2σ
  # of the zone bounds. Tells the reader whether the zone is broadly
  # supported across the gRNA tile or driven by a single guide. Reuses
  # `detected` from the top of the function (single source of truth).
  zones$n_regions_supporting <- mapply(function(s, e) {
    if (nrow(detected) == 0) return(0L)
    sum(detected$pos >= (x_grid[s] - 2 * kernel_sigma) &
        detected$pos <= (x_grid[e] + 2 * kernel_sigma) &
        detected$lfc > 0)
  }, zone_s, zone_e)

  # Geometric edge cap on the peak (mirrors the TF path's max_grna_distance
  # filter — a zone whose peak β sits more than `max_grna_distance` away
  # from any gRNA is in the kernel-tail region where labelling intensity
  # is mathematically determined by one boundary guide, not by genuine
  # binding-zone signal).
  if (is.finite(max_grna_distance)) {
    zones <- zones[zones$distance_to_nearest_grna <= max_grna_distance, ,
                    drop = FALSE]
  }
  # Edge-gRNA weight-cap (mirrors the TF path). Drops zones whose peak β
  # is dominated by a single boundary guide's Gaussian tail.
  if (nrow(zones) > 0 && !is.null(edge_grna_weight_cap) &&
      is.finite(edge_grna_weight_cap) && length(pos_r) >= 2) {
    left_r  <- min(pos_r); right_r <- max(pos_r)
    d   <- outer(zones$peak_position, pos_r, FUN = function(p, r) p - r)
    w   <- exp(-0.5 * (d / kernel_sigma)^2)
    wsum <- rowSums(w); wsum[wsum <= 0] <- 1
    idx_l <- which(pos_r == left_r)[1]
    idx_r <- which(pos_r == right_r)[1]
    edge_frac <- pmax(w[, idx_l] / wsum, w[, idx_r] / wsum)
    zones <- zones[edge_frac <= edge_grna_weight_cap, , drop = FALSE]
  }

  if (nrow(zones) == 0) return(empty)
  zones[order(zones$peak_beta, decreasing = TRUE), , drop = FALSE]
}

# =============================================================================
# SECTION 2: Per-factor zone deck plot
# =============================================================================

#' Per-factor epigenetic zone deck page (Plot 10 analogue).
#'
#' Renders one page for a single epigenetic factor with three lanes:
#'   1. Upper panel: CasPEX signal curve s(x) + per-region logFC stems.
#'   2. Middle lane: zone bars (horizontal coloured rectangles spanning
#'      [zone_start, zone_end]) with a diamond at peak_position. Replaces
#'      the bubble + motif-tick lane from the TF deck.
#'   3. Lower lane: ChIP-Atlas SRX rows (when peaks are passed in).
#'
#' No motif tick lane (epigenetic factors typically have none); no MLE
#' position track (MLE is a TF-specific diagnostic that disambiguates the
#' "low logFC could be far binder" case for sequence-specific anchors).
#' @param tf_name (see function body).
#' @param long_data (see function body).
#' @param pos_map (see function body).
#' @param zones_df (see function body).
#' @param kernel_sigma (see function body).
#' @param upstream (see function body).
#' @param downstream (see function body).
#' @param weight_mode (see function body).
#' @param zone_frac (see function body).
#' @param inner_zone_frac (see function body).
#' @param cov_floor (see function body).
#' @param edge_guard_frac (see function body).
#' @param chipatlas_peaks (see function body).
#' @param peak_signal_range (see function body).
#' @export
plot_epigenetic_zone_deck <- function(
    tf_name, long_data, pos_map, zones_df,
    kernel_sigma    = 300,
    upstream        = 2500,
    downstream      = 500,
    weight_mode     = "lfc_signed",
    zone_frac       = 0.3,
    inner_zone_frac = 0.7,
    cov_floor       = 0.05,
    edge_guard_frac = 0.25,
    # `peak_signal_range`: optional 2-vector c(min, max) defining the
    # colourbar limits in the inner-core fill gradient. When NULL (default),
    # ggplot picks the per-plot range, so different factors render with
    # incomparable colour scales — convenient for solo plotting, misleading
    # for a deck where the reader is comparing across pages. When supplied,
    # the gradient is locked to this range across the whole deck (and could
    # be made cross-locus by the caller). Out-of-range values are squished
    # to the nearest endpoint via scales::squish, so a single high-amplitude
    # outlier doesn't blow the scale out for everyone else.
    peak_signal_range = NULL,
    chipatlas_peaks = NULL) {
  x_grid <- seq(-upstream, downstream, by = 5)
  sig    <- build_caspex_signal(tf_name, long_data, pos_map, x_grid,
                                 kernel_sigma, weight_mode)
  rd      <- sig$region_data
  sig_max <- max(c(sig$y, rd$lfc, 1), na.rm = TRUE)

  # Match plot_binding_deconvolution's window-clipping rationale: scope
  # the visible x-range to where the gRNA tile can actually constrain
  # signal (one σ past the outermost guides on each side).
  pos_r_detect <- sort(as.numeric(pos_map[!is.na(pos_map)]))
  if (length(pos_r_detect) >= 1) {
    left_cut  <- min(pos_r_detect) - kernel_sigma
    right_cut <- max(pos_r_detect) + kernel_sigma
  } else {
    left_cut  <- -upstream
    right_cut <-  downstream
  }
  sig_df <- data.frame(x = x_grid, y = sig$y)

  # ── Lower-panel lane y-positions ────────────────────────────────────────
  neg_floor   <- min(c(0, rd$lfc), na.rm = TRUE)
  track_gap   <- 0.06 * sig_max
  track_top   <- neg_floor - track_gap
  zone_y_top  <- track_top - 0.04 * sig_max
  zone_y_bot  <- track_top - 0.18 * sig_max
  zone_y_mid  <- (zone_y_top + zone_y_bot) / 2
  track_bot   <- zone_y_bot - 0.04 * sig_max

  # ChIP-Atlas band height scales with #SRX (mirrors plot_binding_deconvolution).
  ca_rows <- if (!is.null(chipatlas_peaks) && nrow(chipatlas_peaks) > 0)
    unique(chipatlas_peaks$srx) else character(0)
  n_ca_rows <- length(ca_rows)
  ca_band_h <- if (n_ca_rows == 0) 0 else min(0.50, 0.025 * n_ca_rows + 0.06)
  ca_gap    <- if (n_ca_rows == 0) 0 else 0.03 * sig_max
  ca_top    <- track_bot - ca_gap
  ca_bot    <- ca_top - ca_band_h * sig_max

  p <- ggplot() +
    # ── upper panel: signal + region logFCs ────────────────────────────
    geom_area(data = sig_df, aes(x = x, y = y),
              fill = COLS$guide, alpha = 0.25) +
    geom_line(data = sig_df, aes(x = x, y = y),
              color = COLS$neutral, linewidth = 0.4) +
    geom_segment(data = rd, aes(x = pos, xend = pos, y = 0, yend = lfc),
                 color = COLS$neutral, linewidth = 0.5, alpha = 0.7) +
    geom_point(data = rd, aes(x = pos, y = lfc),
               color = COLS$neutral, size = 3) +
    # Region labels: above for positive lfc, to the right side for
    # negative lfc — matches the TF deck's plot_binding_deconvolution
    # behaviour so the two decks stay visually consistent and so
    # negative-lfc labels don't drop into the zone-bar lane below.
    geom_text(data = transform(rd,
                                 .vj = ifelse(lfc >= 0, -0.9,  0.5),
                                 .hj = ifelse(lfc >= 0,  0.5, -0.3)),
              aes(x = pos, y = lfc, label = region,
                  vjust = .vj, hjust = .hj),
              size = 2.8, fontface = "bold", color = COLS$neutral) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = COLS$tss, linewidth = 0.7) +
    # ── lower-panel container strip ─────────────────────────────────────
    annotate("rect",
             xmin = -upstream, xmax = downstream,
             ymin = track_bot, ymax = track_top,
             fill = "grey97", color = NA)

  # ── Zone bars (two-tone: broad outer + focal inner core) ──────────────
  # Outer rect at [zone_start, zone_end] with low alpha + grey fill marks
  # the broad chromatin-domain extent (β > zone_frac · max). Inner rect at
  # [core_start, core_end] with full alpha + peak_β-coloured fill marks
  # the focal core (β > inner_zone_frac · max within zone). For factors
  # with truly broad uniform binding, outer ≈ inner; for factors with a
  # sharp focal peak inside a broader detectable zone, the dark inner bar
  # sits inside a light outer bar, communicating both pieces of
  # information honestly under the σ=300 bp resolution constraint.
  if (!is.null(zones_df) && nrow(zones_df) > 0) {
    z <- zones_df
    z$y_mid    <- zone_y_mid
    z$y_bottom <- zone_y_bot + 0.005 * sig_max
    z$y_top    <- zone_y_top - 0.005 * sig_max
    p <- p +
      annotate("rect",
               xmin = -upstream, xmax = downstream,
               ymin = zone_y_bot, ymax = zone_y_top,
               fill = "grey94", color = NA) +
      # Outer (broad domain) bar — light grey, no peak-β colouring (the
      # outer extent is a geometric statement about kernel-resolved support
      # range, not a magnitude statement).
      geom_rect(data = z,
                aes(xmin = zone_start, xmax = zone_end,
                    ymin = y_bottom,    ymax = y_top),
                fill = "grey78", color = "grey55",
                linewidth = 0.25, alpha = 0.55)
    # Inner (focal core) bar(s) — peak-β colour gradient, full alpha.
    # Each zone may have one or more cores depending on whether the
    # signal has multiple suprathreshold sub-peaks separated by
    # sub-threshold valleys. The new schema stores them as comma-
    # separated strings in `core_starts` / `core_ends`; we parse those
    # and expand to a long-format core_rects df so ggplot draws one
    # geom_rect per core. Older zones_df produced before this multi-
    # core change carry only `core_start` / `core_end` (single bp per
    # zone) — fall back to that when the new columns are absent.
    has_multi_core <- all(c("core_starts", "core_ends") %in% names(z))
    has_single_core <- all(c("core_start",  "core_end")  %in% names(z))
    core_rects <- if (has_multi_core) {
      lst <- lapply(seq_len(nrow(z)), function(i) {
        cs_str <- z$core_starts[i]
        ce_str <- z$core_ends[i]
        if (is.na(cs_str) || is.na(ce_str) || !nzchar(cs_str) || !nzchar(ce_str))
          return(NULL)
        cs_vec <- suppressWarnings(as.numeric(strsplit(cs_str, ",", fixed = TRUE)[[1]]))
        ce_vec <- suppressWarnings(as.numeric(strsplit(ce_str, ",", fixed = TRUE)[[1]]))
        if (length(cs_vec) == 0L || length(cs_vec) != length(ce_vec)) return(NULL)
        # Drop zero-width and NA entries
        keep <- is.finite(cs_vec) & is.finite(ce_vec) & (ce_vec > cs_vec)
        if (!any(keep)) return(NULL)
        data.frame(core_start = cs_vec[keep], core_end = ce_vec[keep],
                   y_bottom   = z$y_bottom[i], y_top    = z$y_top[i],
                   peak_beta  = z$peak_beta[i])
      })
      do.call(rbind, Filter(Negate(is.null), lst))
    } else if (has_single_core) {
      tmp <- z[!is.na(z$core_start) & !is.na(z$core_end) &
                z$core_end > z$core_start, , drop = FALSE]
      if (nrow(tmp) > 0)
        data.frame(core_start = tmp$core_start, core_end = tmp$core_end,
                   y_bottom   = tmp$y_bottom,   y_top    = tmp$y_top,
                   peak_beta  = tmp$peak_beta)
      else NULL
    } else NULL
    if (!is.null(core_rects) && nrow(core_rects) > 0) {
      # Build the gradient once, with optional global limits so the
      # colourbar is comparable across factor pages within a deck.
      # `scales::squish` handles values outside the range by pinning to
      # the nearest endpoint so a single outlier factor doesn't poison
      # the scale.
      fill_scale <- if (!is.null(peak_signal_range) &&
                        length(peak_signal_range) == 2L &&
                        all(is.finite(peak_signal_range))) {
        scale_fill_gradient(low = "#E6F0F5", high = COLS$high,
                            name = "peak signal",
                            limits = peak_signal_range,
                            oob    = scales::squish,
                            guide  = "colorbar")
      } else {
        scale_fill_gradient(low = "#E6F0F5", high = COLS$high,
                            name = "peak signal", guide = "colorbar")
      }
      p <- p +
        geom_rect(data = core_rects,
                  aes(xmin = core_start, xmax = core_end,
                      ymin = y_bottom,   ymax = y_top,
                      fill = peak_beta),
                  color = "black", linewidth = 0.3, alpha = 0.90) +
        fill_scale
    }
    # Centroid markers — one diamond per region-specific centroid inside
    # each zone. `centroid_positions` is a comma-separated string of bp
    # values (single value for single-centroid zones, multiple for
    # multi-peak zones). Parse and expand to a long-format data.frame so
    # ggplot draws one diamond per centroid. Falls back to peak_position
    # if the column is missing (legacy zones_df).
    centroid_df <- if ("centroid_positions" %in% names(z)) {
      lst <- lapply(seq_len(nrow(z)), function(i) {
        cp <- z$centroid_positions[i]
        if (is.na(cp) || !nzchar(cp)) return(NULL)
        pos <- suppressWarnings(as.numeric(strsplit(cp, ",",
                                                     fixed = TRUE)[[1]]))
        pos <- pos[is.finite(pos)]
        if (length(pos) == 0) return(NULL)
        data.frame(centroid_pos = pos, y_mid = z$y_mid[i])
      })
      do.call(rbind, Filter(Negate(is.null), lst))
    } else {
      data.frame(centroid_pos = z$peak_position, y_mid = z$y_mid)
    }
    if (!is.null(centroid_df) && nrow(centroid_df) > 0) {
      # Centroid markers — black, with thin white border for contrast
      # against either the red inner-core fill or the grey outer
      # background. Shape 23 (filled diamond with separate stroke)
      # gives us a fill + outline pair so the diamond is unambiguously
      # the topmost element regardless of what's underneath. Larger
      # size + thicker stroke than before so the marker reads at deck
      # zoom levels too.
      p <- p +
        geom_point(data = centroid_df,
                   aes(x = centroid_pos, y = y_mid),
                   shape = 23, size = 4.0,
                   fill = "black", color = "white", stroke = 0.6) +
        geom_text(data = centroid_df,
                  aes(x = centroid_pos, y = y_mid,
                      label = sprintf("%+.0f", centroid_pos)),
                  size = 2.4, fontface = "bold", vjust = 2.4,
                  color = "black")
    }
  }

  # ── ChIP-Atlas stacked sub-lane ───────────────────────────────────────
  if (n_ca_rows > 0) {
    row_h <- (ca_top - ca_bot) / n_ca_rows
    ord   <- order(suppressWarnings(as.integer(sub("SRX", "", ca_rows))),
                   decreasing = TRUE)
    srx_levels <- ca_rows[ord]
    ca <- chipatlas_peaks
    ca$srx <- factor(ca$srx, levels = srx_levels)
    ca$y_mid <- ca_top - (as.integer(ca$srx) - 0.5) * row_h
    ca$xs <- pmax(ca$start_rel, -upstream)
    ca$xe <- pmin(ca$end_rel,    downstream)
    ca$xs <- pmin(ca$xs, ca$xe - 5)
    p <- p +
      annotate("rect",
               xmin = -upstream, xmax = downstream,
               ymin = ca_bot, ymax = ca_top,
               fill = "grey98", color = NA) +
      geom_segment(data = ca,
                   aes(x = xs, xend = xe, y = y_mid, yend = y_mid),
                   color = COLS$tss, linewidth = 0.5, alpha = 0.65,
                   lineend = "butt") +
      annotate("text",
               x = downstream, y = (ca_top + ca_bot) / 2,
               label = {
                 ns <- attr(chipatlas_peaks, "n_srx_scanned")
                 sig <- isTRUE(attr(chipatlas_peaks, "is_special_interest"))
                 base <- if (is.null(ns) || is.na(ns))
                   sprintf("ChIP-Atlas  \u00b7  %d SRX", n_ca_rows)
                 else
                   sprintf("ChIP-Atlas  \u00b7  %d SRX / %d experiments",
                           n_ca_rows, ns)
                 if (sig) paste0(base, "*") else base
               },
               hjust = 1.02, vjust = -0.6, size = 2.3,
               color = "grey35", fontface = "italic")
  }

  y_label <- switch(weight_mode,
    z          = "CasPEX signal (z-weighted, signed z from p-value)",
    mod_t      = "CasPEX signal (mod-t weighted)",
    lfc_pos    = "CasPEX signal (positive logFC-weighted)",
    lfc_signed = "CasPEX signal (signed logFC-weighted)",
    lfc_x_negp = "CasPEX signal (logFC \u00d7 -log10 p-weighted)",
    paste0("CasPEX signal (", weight_mode, "-weighted)"))

  p +
    scale_y_continuous(
      expand = expansion(mult = c(0.08, 0.08)),
      breaks = function(lim) pretty(c(0, lim[2]))) +
    labs(x = "Position (bp, TSS-relative)",
         y = y_label,
         title = paste0(tf_name, " \u2014 epigenetic binding zones"),
         subtitle = sprintf(
           "%d zone(s) | kernel \u03c3 = %d bp  \u00b7  outer bar = \u03b2 > %.2f\u00b7max(\u03b2) (broad domain) | inner bar = \u03b2 > %.2f\u00b7max(\u03b2) within zone (focal core)",
           # Avoid `%||%` here — caspex_analysis.R's definition does
           # `is.na(a[[1]])` which returns a length-N logical when `a` is
           # a data.frame, tripping the `||` short-circuit. Use an
           # explicit null check.
           if (is.null(zones_df)) 0L else nrow(zones_df),
           kernel_sigma, zone_frac, inner_zone_frac)) +
    theme_caspex()
}

# =============================================================================
# SECTION 3: ChIP-Atlas overlap diagnostic (per-zone)
# =============================================================================

#' Compute ChIP-Atlas overlap statistics for each zone.
#'
#' For each zone, calculates:
#'   * peak_bp_in_zone — total bp of any ChIP-Atlas peak overlapping the zone
#'     (clipped to zone bounds, union across all SRXs)
#'   * peak_coverage_frac — peak_bp_in_zone / zone_width
#'   * n_srx_with_peak — number of distinct SRX experiments contributing
#'     at least one peak that overlaps the zone
#'   * any_overlap — boolean shortcut for downstream filtering
#'
#' Run AFTER fetch — purely a per-zone summary; does not alter zones_df.
#'
#' @param zones_df       output of predict_binding_zones_epigenetic()
#'                       (concatenated across all factors).
#' @param chipatlas_res  named list tf -> peaks data.frame (from
#'                       run_chipatlas_scan).
#' @return zones_df with the four overlap columns appended.
#' @noRd
compute_chipatlas_overlap <- function(zones_df, chipatlas_res) {
  zones_df$peak_bp_in_zone   <- 0L
  zones_df$peak_coverage_frac <- 0
  zones_df$n_srx_with_peak    <- 0L
  zones_df$any_overlap        <- FALSE
  if (is.null(chipatlas_res) || nrow(zones_df) == 0) return(zones_df)
  for (i in seq_len(nrow(zones_df))) {
    tf  <- zones_df$tf[i]
    pks <- chipatlas_res[[tf]]
    if (is.null(pks) || nrow(pks) == 0) next
    zs  <- zones_df$zone_start[i]
    ze  <- zones_df$zone_end[i]
    # Per-peak overlap with this zone (clipped to zone bounds)
    ov_starts <- pmax(pks$start_rel, zs)
    ov_ends   <- pmin(pks$end_rel,   ze)
    ov_lens   <- pmax(ov_ends - ov_starts, 0)
    keep      <- ov_lens > 0
    if (!any(keep)) next
    # Union of overlapping intervals (collapse overlapping peaks to bp)
    iv  <- data.frame(s = ov_starts[keep], e = ov_ends[keep])
    iv  <- iv[order(iv$s), , drop = FALSE]
    merged_s <- iv$s[1]; merged_e <- iv$e[1]; tot <- 0
    for (j in seq_len(nrow(iv))[-1]) {
      if (iv$s[j] <= merged_e) {
        merged_e <- max(merged_e, iv$e[j])
      } else {
        tot <- tot + (merged_e - merged_s)
        merged_s <- iv$s[j]; merged_e <- iv$e[j]
      }
    }
    tot <- tot + (merged_e - merged_s)
    zones_df$peak_bp_in_zone[i]    <- as.integer(tot)
    zones_df$peak_coverage_frac[i] <- tot / max(ze - zs, 1)
    zones_df$n_srx_with_peak[i]    <- length(unique(pks$srx[keep]))
    zones_df$any_overlap[i]        <- TRUE
  }
  zones_df
}

# =============================================================================
# SECTION 3b: Epigenetic complexes — loader + per-complex stacked plot
# =============================================================================

#' Load EpiGenes complexes from a CSV pair.
#'
#' Expected files (curated from EpiFactors-style annotations):
#'   * `epigenes_main_csv` — gene-level table with `HGNC_symbol` and
#'     `UniProt_ID` columns. Used to build a UniProt_ID → HGNC_symbol
#'     lookup so each complex's member proteins can be resolved to
#'     canonical HGNC names regardless of how the source file spells
#'     them in the protein column.
#'   * `complex_csv` — complex-level table with `Complex_name`,
#'     `Group_name`, `UniProt_ID` (comma-separated, may contain `?`,
#'     `+`, `|`, parens for "uncertain" / "one-of" annotations) and
#'     descriptive columns (`Function`, `Target`, etc.).
#'
#' Parses each complex's UniProt_ID list (stripping `?`, `+`, parens;
#' splitting `|` as alternatives), maps to HGNC, and optionally filters
#' to detected proteins. Returns one element per complex, sorted by
#' how many of its members were detected (descending).
#'
#' @param complex_csv      path to EpiGenes_complexes.csv.
#' @param main_csv         path to EpiGenes_main.csv (for UniProt → HGNC).
#' @param detected_proteins optional character vector of HGNC symbols
#'                         present in the proteomics data; used to
#'                         flag which members are detected at this
#'                         locus (other members still listed).
#' @return list of complex records:
#'   list(complex_id, complex_name, group_name, function_, target,
#'        members_total, members_hgnc, members_detected)
#' @noRd
load_epigenes_complexes <- function(complex_csv, main_csv,
                                     detected_proteins = NULL) {
  if (!file.exists(main_csv))
    stop("EpiGenes_main.csv not found: ", main_csv)
  if (!file.exists(complex_csv))
    stop("EpiGenes_complexes.csv not found: ", complex_csv)

  main_df <- read.csv(main_csv, stringsAsFactors = FALSE,
                       fileEncoding = "UTF-8", check.names = FALSE)
  if (!all(c("HGNC_symbol", "UniProt_ID") %in% names(main_df)))
    stop("EpiGenes_main.csv missing HGNC_symbol or UniProt_ID columns.")
  upid_to_hgnc <- setNames(toupper(main_df$HGNC_symbol),
                            toupper(main_df$UniProt_ID))
  upid_to_hgnc <- upid_to_hgnc[nzchar(names(upid_to_hgnc)) &
                                 nzchar(upid_to_hgnc)]

  cmplx_df <- read.csv(complex_csv, stringsAsFactors = FALSE,
                        fileEncoding = "UTF-8", check.names = FALSE)
  if (!all(c("Complex_name", "UniProt_ID") %in% names(cmplx_df)))
    stop("EpiGenes_complexes.csv missing Complex_name or UniProt_ID.")

  parse_upids <- function(s) {
    if (is.null(s) || is.na(s) || !nzchar(s)) return(character(0))
    # Strip uncertainty markers (?, +) and parentheses; treat | as
    # an alternative separator equivalent to comma.
    s2 <- gsub("[?()+]", "", s)
    s2 <- gsub("\\|", ",", s2)
    parts <- strsplit(s2, ",", fixed = TRUE)[[1]]
    parts <- toupper(trimws(parts))
    parts[nzchar(parts)]
  }

  detected_norm <- if (!is.null(detected_proteins))
    toupper(detected_proteins) else NULL

  out <- lapply(seq_len(nrow(cmplx_df)), function(i) {
    upids <- parse_upids(cmplx_df$UniProt_ID[i])
    hgnc_members <- unname(upid_to_hgnc[upids])
    hgnc_members <- hgnc_members[!is.na(hgnc_members) & nzchar(hgnc_members)]
    hgnc_members <- unique(hgnc_members)
    detected_members <- if (!is.null(detected_norm))
      intersect(hgnc_members, detected_norm) else hgnc_members
    list(
      complex_id       = cmplx_df$Id[i],
      complex_name     = cmplx_df$Complex_name[i],
      group_name       = cmplx_df$Group_name[i],
      function_        = cmplx_df$Function[i],
      target           = cmplx_df$Target[i],
      alternative_name = cmplx_df$Alternative_name[i],
      members_total    = length(hgnc_members),
      members_hgnc     = hgnc_members,
      members_detected = detected_members)
  })
  # Order: most-detected complexes first, ties broken by total members
  # ascending (smaller complexes more readable).
  ord <- order(-vapply(out, function(x) length(x$members_detected),
                        integer(1)),
                vapply(out, function(x) x$members_total, integer(1)))
  out[ord]
}

#' Per-complex stacked binding-zone page (analogue of histone-marks page).
#'
#' One row per member protein of an epigenetic complex, stacked
#' top-to-bottom in the same x-axis frame as the per-factor pages
#' (bp TSS-relative, gRNA arrowheads, TSS dashed line). Each row
#' shows that member's outer zone bars (light grey) and inner-core
#' bars (peak-signal coloured) at this locus — the same visual
#' vocabulary as the per-factor pages, but compressed into one row
#' per member so the reader can see at a glance whether members of
#' the complex co-localise (suggesting the complex acts as a unit
#' here) or scatter independently.
#'
#' Detected members render with their zones; not-detected members
#' (still listed in the complex's roster but absent from the
#' proteomics data) render with their label and an "(n/d)" annotation
#' so absence is visually distinct from absence-of-zones.
#' @param complex_info (see function body).
#' @param zones_df (see function body).
#' @param sig_df_per_member (see function body).
#' @param pos_map (see function body).
#' @param gene_info (see function body).
#' @param upstream (see function body).
#' @param downstream (see function body).
#' @param kernel_sigma (see function body).
#' @param peak_signal_range (see function body).
#' @export
plot_epigenetic_complex_locus <- function(
    complex_info, zones_df, sig_df_per_member, pos_map, gene_info,
    upstream         = 2500,
    downstream       = 500,
    kernel_sigma     = 300,
    peak_signal_range = NULL) {

  detected <- complex_info$members_detected
  not_detected <- setdiff(complex_info$members_hgnc, detected)
  # Order rows: detected first (alphabetical), then not-detected
  # (alphabetical) at the bottom so they don't visually crowd the
  # data rows.
  ordered <- c(sort(detected), sort(not_detected))
  if (!length(ordered)) return(NULL)
  n <- length(ordered)

  bar_h <- 0.7
  row_top <- 1.0
  row_bot <- row_top - n
  grna_y  <- row_top + 0.45

  # gRNA arrowhead row at panel top
  pos_r <- sort(as.numeric(pos_map[!is.na(pos_map)]))
  grna_df <- data.frame(x = pos_r, y = grna_y)

  p <- ggplot() +
    annotate("rect",
             xmin = -upstream, xmax = downstream,
             ymin = row_bot, ymax = row_top,
             fill = "grey97", color = NA) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = COLS$tss, linewidth = 0.7) +
    geom_point(data = grna_df, aes(x = x, y = y),
               shape = 25, size = 3, fill = COLS$guide,
               color = COLS$neutral, stroke = 0.4)

  # Build per-row outer + inner core rectangles by parsing the
  # zones_df rows for each member. zones_df was produced by
  # predict_binding_zones_epigenetic and may contain multiple zones
  # per TF; each zone may have multiple cores (core_starts /
  # core_ends comma-separated strings).
  outer_rects <- list(); core_rects <- list(); centroid_pts <- list()

  for (i in seq_along(ordered)) {
    tf <- ordered[i]
    y_mid <- row_top - i + 0.5
    if (!(tf %in% detected)) next
    z <- zones_df[zones_df$tf == tf, , drop = FALSE]
    if (nrow(z) == 0) next
    # Outer zone bars (one per zone row)
    outer_rects[[length(outer_rects) + 1]] <- data.frame(
      tf = tf, y_mid = y_mid,
      xmin = z$zone_start, xmax = z$zone_end,
      ymin = y_mid - bar_h / 2, ymax = y_mid + bar_h / 2)
    # Inner cores (parse comma-separated)
    for (zi in seq_len(nrow(z))) {
      cs_str <- z$core_starts[zi]; ce_str <- z$core_ends[zi]
      if (is.na(cs_str) || is.na(ce_str) || !nzchar(cs_str)) next
      cs <- suppressWarnings(as.numeric(strsplit(cs_str, ",", fixed = TRUE)[[1]]))
      ce <- suppressWarnings(as.numeric(strsplit(ce_str, ",", fixed = TRUE)[[1]]))
      keep <- is.finite(cs) & is.finite(ce) & ce > cs
      if (!any(keep)) next
      core_rects[[length(core_rects) + 1]] <- data.frame(
        tf = tf, y_mid = y_mid,
        xmin = cs[keep], xmax = ce[keep],
        ymin = y_mid - bar_h / 2 + 0.05,
        ymax = y_mid + bar_h / 2 - 0.05,
        peak_beta = z$peak_beta[zi])
    }
    # Centroid diamonds (parse comma-separated)
    for (zi in seq_len(nrow(z))) {
      cp_str <- z$centroid_positions[zi]
      if (is.na(cp_str) || !nzchar(cp_str)) next
      cps <- suppressWarnings(as.numeric(strsplit(cp_str, ",", fixed = TRUE)[[1]]))
      cps <- cps[is.finite(cps)]
      if (!length(cps)) next
      centroid_pts[[length(centroid_pts) + 1]] <- data.frame(
        tf = tf, x = cps, y_mid = y_mid)
    }
  }
  outer_df    <- if (length(outer_rects))    do.call(rbind, outer_rects)    else NULL
  core_df     <- if (length(core_rects))     do.call(rbind, core_rects)     else NULL
  centroid_df <- if (length(centroid_pts))   do.call(rbind, centroid_pts)   else NULL

  if (!is.null(outer_df) && nrow(outer_df) > 0) {
    p <- p + geom_rect(data = outer_df,
                       aes(xmin = xmin, xmax = xmax,
                           ymin = ymin, ymax = ymax),
                       fill = "grey78", color = "grey55",
                       linewidth = 0.25, alpha = 0.55)
  }
  if (!is.null(core_df) && nrow(core_df) > 0) {
    fill_scale <- if (!is.null(peak_signal_range) &&
                      length(peak_signal_range) == 2L &&
                      all(is.finite(peak_signal_range))) {
      scale_fill_gradient(low = "#E6F0F5", high = COLS$high,
                          name = "peak signal",
                          limits = peak_signal_range,
                          oob    = scales::squish,
                          guide  = "colorbar")
    } else {
      scale_fill_gradient(low = "#E6F0F5", high = COLS$high,
                          name = "peak signal", guide = "colorbar")
    }
    p <- p + geom_rect(data = core_df,
                       aes(xmin = xmin, xmax = xmax,
                           ymin = ymin, ymax = ymax,
                           fill = peak_beta),
                       color = "black", linewidth = 0.25, alpha = 0.90) +
      fill_scale
  }
  if (!is.null(centroid_df) && nrow(centroid_df) > 0) {
    # Black diamond with white outline — same convention as the
    # per-factor zone deck. Drawn last so it always sits on top of
    # the inner-core fills and remains visible at any zoom.
    p <- p + geom_point(data = centroid_df,
                        aes(x = x, y = y_mid),
                        shape = 23, size = 3.0,
                        fill = "black", color = "white", stroke = 0.5)
  }

  # Right-edge member labels — detected members in bold, not-detected
  # in faded grey + "(n/d)" annotation.
  label_df <- data.frame(
    tf = ordered,
    y_mid = row_top - seq_along(ordered) + 0.5,
    label = ifelse(ordered %in% detected,
                   ordered,
                   paste0(ordered, "  (n/d)")),
    is_det = ordered %in% detected,
    stringsAsFactors = FALSE)
  p <- p +
    geom_text(data = label_df,
              aes(x = downstream, y = y_mid, label = label,
                  color = is_det, fontface = ifelse(is_det, "bold", "plain")),
              hjust = -0.05, size = 2.8) +
    scale_color_manual(values = c(`TRUE` = "grey20", `FALSE` = "grey55"),
                        guide = "none")

  # Title / subtitle — complex info
  title_str <- paste0(
    gene_info$name %||% "(locus)", "  \u2014  ",
    complex_info$group_name, " \u00b7 ", complex_info$complex_name,
    "  (", length(detected), "/", length(ordered),
    " members detected)")
  subtitle_parts <- c()
  if (!is.null(complex_info$function_) && nzchar(complex_info$function_) &&
      complex_info$function_ != "#")
    subtitle_parts <- c(subtitle_parts, complex_info$function_)
  if (!is.null(complex_info$target) && nzchar(complex_info$target) &&
      complex_info$target != "#")
    subtitle_parts <- c(subtitle_parts, paste0("target: ", complex_info$target))
  if (!is.null(complex_info$alternative_name) &&
      nzchar(complex_info$alternative_name) &&
      complex_info$alternative_name != "#")
    subtitle_parts <- c(subtitle_parts,
                         paste0("aka ", complex_info$alternative_name))
  subtitle_str <- if (length(subtitle_parts) > 0)
    paste(subtitle_parts, collapse = "  \u00b7  ") else ""

  p +
    scale_y_continuous(
      limits = c(row_bot - 0.3, grna_y + 0.4),
      expand = c(0, 0)) +
    scale_x_continuous(
      expand = expansion(mult = c(0.05, 0.28))) +
    coord_cartesian(clip = "off") +
    labs(x = "Position (bp, TSS-relative)", y = NULL,
         title = title_str, subtitle = subtitle_str) +
    theme_caspex() +
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin  = margin(t = 5.5, r = 16, b = 5.5, l = 5.5,
                                 unit = "pt"))
}

# =============================================================================
# SECTION 4: Locus-level histone-marks page
# =============================================================================

#' Compute union-of-peaks intervals across SRXs for a single mark/condition.
#'
#' Given a peaks data.frame (any number of rows from any number of SRXs),
#' returns a data.frame of merged non-overlapping intervals in TSS-relative
#' bp coordinates. "Union" = "any SRX with a peak at this position → mark
#' this position as covered". Intervals are in [start_rel, end_rel] format
#' so the plot lane can render one geom_segment per interval.
#' @param peaks_df (see function body).
#' @param x_lo (see function body).
#' @param x_hi (see function body).
#' @noRd
.histone_union_intervals <- function(peaks_df, x_lo, x_hi) {
  if (is.null(peaks_df) || nrow(peaks_df) == 0)
    return(data.frame(start_rel = numeric(), end_rel = numeric()))
  iv <- data.frame(start_rel = pmax(peaks_df$start_rel, x_lo),
                    end_rel   = pmin(peaks_df$end_rel,   x_hi))
  iv <- iv[iv$end_rel > iv$start_rel, , drop = FALSE]
  if (nrow(iv) == 0)
    return(data.frame(start_rel = numeric(), end_rel = numeric()))
  iv <- iv[order(iv$start_rel), , drop = FALSE]
  out_s <- numeric(0); out_e <- numeric(0)
  cur_s <- iv$start_rel[1]; cur_e <- iv$end_rel[1]
  for (j in seq_len(nrow(iv))[-1]) {
    if (iv$start_rel[j] <= cur_e) {
      cur_e <- max(cur_e, iv$end_rel[j])
    } else {
      out_s <- c(out_s, cur_s); out_e <- c(out_e, cur_e)
      cur_s <- iv$start_rel[j]; cur_e <- iv$end_rel[j]
    }
  }
  out_s <- c(out_s, cur_s); out_e <- c(out_e, cur_e)
  data.frame(start_rel = out_s, end_rel = out_e)
}

#' Locus-level histone-marks summary page (one PDF page).
#'
#' Two stacked sections of six rows each. Top section: cell-type-matched
#' (e.g. HEK293T) consensus peaks per mark. Bottom section: all-cell-types
#' aggregated consensus peaks per mark. Same x-axis as the per-factor pages
#' (bp TSS-relative, gRNA positions and TSS dashed line for context).
#'
#' Active marks (warm fill) on top of each section, then repressive marks
#' (cool fill). Empty rows still render with their label so absence-of-data
#' is visually distinct from absence-of-signal — readers see "this mark
#' wasn't measured in HEK293T at this locus" rather than the row vanishing.
#'
#' @param histone_data    list returned by fetch_histone_peaks_for_locus().
#' @param pos_map         passed-through gRNA position vector (for the gRNA
#'                        ticks at the top).
#' @param gene_info       passed-through gene info (for plot title).
#' @param upstream        bp upstream of TSS for the x-axis.
#' @param downstream      bp downstream of TSS for the x-axis.
#' @export
plot_histone_marks_locus <- function(
    histone_data, pos_map, gene_info,
    upstream   = 2500,
    downstream = 500) {
  marks    <- histone_data$marks
  ct       <- histone_data$cell_type %||% "(unspecified)"
  # Active first, repressive last — fixed order regardless of input order
  # so every deck is visually comparable. Anything passed in but not
  # recognised falls into "other" and renders below the canonical six in
  # neutral grey.
  active_marks     <- c("H3K4me3", "H3K27ac", "H3K4me1", "H3K36me3")
  repressive_marks <- c("H3K27me3", "H3K9me3")
  ordered <- c(intersect(active_marks,     marks),
               intersect(repressive_marks, marks),
               setdiff(marks, c(active_marks, repressive_marks)))

  # Render rows top-to-bottom: matched section first (12 rows total: 6 +
  # divider + 6, with row index decreasing so the visual top-to-bottom
  # order matches the requested layout).
  n_marks <- length(ordered)
  matched_y_top <- 1.00
  matched_y_bot <- matched_y_top - n_marks
  divider_y     <- matched_y_bot - 0.5
  all_y_top     <- divider_y - 0.5
  all_y_bot     <- all_y_top - n_marks

  fill_for <- function(mark) {
    if (mark %in% active_marks)     return("#E63946")  # warm (COLS$high)
    if (mark %in% repressive_marks) return("#2A9D8F")  # cool (COLS$low)
    return("grey60")
  }

  # Build one row of intervals for each (section, mark). A section-mark
  # combination with no peaks gets an empty intervals data.frame; the
  # row label still draws from `row_meta` even when the bar is empty.
  rows_intervals <- list()
  rows_meta      <- data.frame(
    section = character(), mark = character(),
    y_mid   = numeric(),   n_srx = integer(),
    fill    = character(), stringsAsFactors = FALSE)

  x_lo <- -upstream; x_hi <- downstream
  for (i in seq_along(ordered)) {
    mark <- ordered[i]
    y_mid_matched <- matched_y_top - i + 0.5
    y_mid_all     <- all_y_top     - i + 0.5
    iv_m <- .histone_union_intervals(histone_data$matched[[mark]], x_lo, x_hi)
    iv_a <- .histone_union_intervals(histone_data$all[[mark]],     x_lo, x_hi)
    if (nrow(iv_m) > 0) {
      iv_m$y_mid <- y_mid_matched
      iv_m$mark  <- mark
      iv_m$section <- "matched"
      rows_intervals[[length(rows_intervals) + 1]] <- iv_m
    }
    if (nrow(iv_a) > 0) {
      iv_a$y_mid <- y_mid_all
      iv_a$mark  <- mark
      iv_a$section <- "all"
      rows_intervals[[length(rows_intervals) + 1]] <- iv_a
    }
    rows_meta <- rbind(rows_meta, data.frame(
      section = c("matched", "all"),
      mark    = c(mark, mark),
      y_mid   = c(y_mid_matched, y_mid_all),
      n_srx   = c(histone_data$n_srx_matched[[mark]] %||% 0L,
                  histone_data$n_srx_all[[mark]]     %||% 0L),
      fill    = c(fill_for(mark), fill_for(mark)),
      stringsAsFactors = FALSE))
  }
  intervals_df <- if (length(rows_intervals) > 0)
    do.call(rbind, rows_intervals) else
    data.frame(start_rel = numeric(), end_rel = numeric(),
               y_mid = numeric(), mark = character(), section = character())

  # gRNA tick row at top of plot (above the matched section).
  pos_r <- sort(as.numeric(pos_map[!is.na(pos_map)]))
  grna_y <- matched_y_top + 0.4
  grna_df <- data.frame(x = pos_r, y = grna_y)

  bar_h <- 0.7  # bar height per row (in y-units; row spacing is 1)

  p <- ggplot() +
    # ── Section background bands ────────────────────────────────────────
    annotate("rect",
             xmin = -upstream, xmax = downstream,
             ymin = matched_y_bot, ymax = matched_y_top,
             fill = "grey97", color = NA) +
    annotate("rect",
             xmin = -upstream, xmax = downstream,
             ymin = all_y_bot, ymax = all_y_top,
             fill = "grey97", color = NA) +
    # ── Section labels (left-edge italics) ──────────────────────────────
    annotate("text",
             x = -upstream, y = (matched_y_top + matched_y_bot) / 2,
             label = sprintf("Cell-type-matched (%s)", ct),
             hjust = 1.05, size = 3, fontface = "italic", color = "grey25") +
    annotate("text",
             x = -upstream, y = (all_y_top + all_y_bot) / 2,
             label = "All cell types aggregated",
             hjust = 1.05, size = 3, fontface = "italic", color = "grey25") +
    # ── TSS line + gRNAs ────────────────────────────────────────────────
    geom_vline(xintercept = 0, linetype = "dashed",
               color = COLS$tss, linewidth = 0.7) +
    geom_point(data = grna_df, aes(x = x, y = y),
               shape = 25, size = 3, fill = COLS$guide,
               color = COLS$neutral, stroke = 0.4)

  # Per-row peak intervals
  if (nrow(intervals_df) > 0) {
    intervals_df$fill <- fill_for_vec(intervals_df$mark, active_marks, repressive_marks)
    p <- p +
      geom_rect(data = intervals_df,
                aes(xmin = start_rel, xmax = end_rel,
                    ymin = y_mid - bar_h / 2,
                    ymax = y_mid + bar_h / 2,
                    fill = fill),
                color = "grey20", linewidth = 0.15, alpha = 0.85) +
      scale_fill_identity()
  }

  # Mark labels at right edge (always render — empty rows still show
  # their label so absence-of-data is visible).
  p <- p +
    geom_text(data = rows_meta,
              aes(x = downstream, y = y_mid,
                  label = sprintf("%s  (n=%d)", mark, n_srx)),
              hjust = -0.05, size = 2.8, fontface = "bold",
              color = "grey20")

  # gRNA region labels (faint, above the gRNA tick row)
  if (length(pos_r) > 0) {
    valid <- pos_map[!is.na(pos_map)]
    grna_lbl_df <- data.frame(x = as.numeric(valid), label = names(valid))
    p <- p +
      geom_text(data = grna_lbl_df,
                aes(x = x, y = grna_y + 0.15, label = label),
                size = 2.4, fontface = "bold", color = COLS$neutral)
  }

  p +
    scale_y_continuous(
      limits = c(all_y_bot - 0.3, matched_y_top + 0.7),
      expand = c(0, 0)) +
    scale_x_continuous(
      # Generous horizontal expansion on BOTH sides so the section labels
      # ("Cell-type-matched (HEK293T)" on the left, ~26 chars) and the
      # row labels with 4-5 digit SRX counts on the right
      # ("H3K27ac (n=11827)") aren't clipped at the panel edges. With
      # the default 0.10/0.10 expansion the bottom-row labels lost their
      # closing paren and the left section header lost "Cell-".
      expand = expansion(mult = c(0.28, 0.28))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = "Position (bp, TSS-relative)",
      y = NULL,
      title = paste0(
        gene_info$name %||% "(locus)",
        " \u2014 histone-mark landscape"),
      subtitle = sprintf(
        "Top: %s-matched ChIP-Atlas peaks per mark, union across SRXs  |  Bottom: all cell types aggregated, union across newest-N SRXs",
        ct)) +
    theme_caspex() +
    theme(panel.grid = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          # Give the plot device a bit of horizontal slack so the off-
          # panel section/row labels don't run into the figure edge.
          plot.margin  = margin(t = 5.5, r = 16, b = 5.5, l = 16,
                                 unit = "pt"))
}

# Helper: vectorised version of fill_for() used inside ggplot data — same
# colour palette but operates on a character vector. Defined at top level
# so plot_histone_marks_locus() can reference it inside the geom layer.
#' Recycle a fill vector to a target length, padding with grey.
#' @param marks (see function body).
#' @param active_marks (see function body).
#' @param repressive_marks (see function body).
#' @noRd
fill_for_vec <- function(marks, active_marks, repressive_marks) {
  ifelse(marks %in% active_marks,     "#E63946",
  ifelse(marks %in% repressive_marks, "#2A9D8F", "grey60"))
}

# =============================================================================
# SECTION 5: Top-level orchestrator — run_caspex_epigenetic()
# =============================================================================

#' Run the epigenetic-factor binding-zone analysis on a CasPEX result.
#'
#' Designed to be called AFTER `run_caspex()` on the same locus, taking the
#' returned result list as input. Reuses `result$long_data`, `result$pos_map`,
#' `result$gene_info`, `result$promoter_info`, `result$kernel_sigma`,
#' `result$weight_mode` — no inputs are reloaded from disk. Mirrors the
#' `run_caspex_extras()` ergonomic pattern.
#'
#' What this function does, in order:
#'   1. Filter result$long_data to the supplied epigenetic factor list.
#'      Reports the detection-overlap rate (factors in list ∩ proteins
#'      detected in long_data) — the only proteins we can model.
#'   2. Run the same spatial model (run_spatial_model) used by run_caspex(),
#'      but with tfs_only = FALSE so the TFDatabase filter doesn't drop
#'      epigenetic factors that happen to be missing the isTF flag.
#'   3. Per surviving factor, call predict_binding_zones_epigenetic() to
#'      get its zone bars.
#'   4. Optionally fetch ChIP-Atlas peaks for each factor in the deck
#'      (default ON; same wrapper as the TF path).
#'   5. Compute per-zone ChIP-Atlas overlap statistics
#'      (compute_chipatlas_overlap).
#'   6. Generate a per-factor PDF deck and write zones / spatial CSVs.
#'
#'
#' @param result Output of \code{\link{run_caspex}} on the same locus.
#'   Reused: \code{long_data}, \code{pos_map}, \code{gene_info},
#'   \code{promoter_info}, \code{kernel_sigma}, \code{weight_mode}.
#' @param epigenetic_factors Character vector of HGNC symbols for the
#'   factors to model (e.g. read from \code{EpiGenes_main.csv}). The
#'   intersection with \code{result$long_data$protein} is what actually
#'   gets analysed; missing factors are reported but skipped.
#' @param out_dir Output directory; created if absent.
#'
#'
#' @param zone_frac Outer-zone threshold for zone detection:
#'   \eqn{\beta(x) > zone\_frac \cdot max(\beta)} defines the broad
#'   chromatin-domain extent. Default 0.3 (lower than the TF path's
#'   0.5 so near-peak beta values merge into a single broad domain).
#' @param inner_zone_frac Inner-core threshold within each outer zone:
#'   contiguous runs above \eqn{inner\_zone\_frac \cdot zone\_peak\_beta}
#'   become darker inner bars marking the most concentrated binding.
#'   Default 0.7.
#' @param centroid_frac Threshold on the per-zone regional logFC that a
#'   candidate centroid region must clear (relative to the zone's max
#'   regional logFC). Default 0.7.
#'
#'
#' @param kernel_sigma Labelling-radius sigma in bp. NULL = inherit from
#'   \code{result$kernel_sigma} (typically 300).
#' @param weight_mode Region-weight mode. NULL = inherit from
#'   \code{result$weight_mode} (typically "z").
#' @param cov_floor Coverage-denominator floor (see \code{\link{run_caspex}}).
#'   Default 0.05.
#' @param edge_guard_frac In-support beta mask (see \code{\link{run_caspex}}).
#'   Default 0.25.
#' @param edge_grna_weight_cap Optional boundary-gRNA weight-share cap
#'   (see \code{\link{run_caspex}}). NULL disables (default).
#' @param max_grna_distance Hard geometric cap on event-to-nearest-guide
#'   distance (see \code{\link{run_caspex}}). NULL (default).
#' @param upstream,downstream bp window around TSS. Defaults 2500 / 500.
#'   Should match the TF run's window so coordinates line up.
#'
#'
#' @param min_n_regions Minimum regions a factor must be detected in to
#'   enter the spatial model. Default 2.
#' @param subtract_tf_overlap If TRUE, drop factors that were also in
#'   \code{result$motif_results} (i.e. motif-scanned in the TF deck) so
#'   the epigenetic deck doesn't duplicate them. Default FALSE - render
#'   both decks for dual-class proteins (GATA, FOX, KLF, KRAB-ZNFs) so
#'   the reader can compare bubble vs. zone framings.
#' @param detail_top_n Max number of per-factor detail pages (default 50;
#'   ranked by max peak beta).
#' @param peak_signal_range Optional 2-vector \code{c(min, max)} for the
#'   deck-wide colour-bar limits on the inner-core fill. NULL (default)
#'   auto-computes from the union of all rendered factors' peak betas
#'   so every page shares one colour scale. Override when stitching
#'   multiple decks (e.g. shared range across hTERT + MYC).
#'
#'
#' @param chipatlas Fetch and overlay ChIP-Atlas peaks per factor
#'   (default TRUE). Cache shared with the TF run.
#' @param chipatlas_threshold One of \code{"05"} (Q<1e-5, default),
#'   \code{"10"}, or \code{"20"}.
#' @param chipatlas_max_experiments Cap on SRX experiments per factor
#'   (default 100).
#' @param special_interest_gene Optional character vector of factor
#'   symbols whose ChIP-Atlas scan bypasses the per-factor cap (same
#'   semantics as in \code{\link{run_caspex}}).
#' @param special_interest_cap Optional integer cap on the
#'   special-interest SRX count. NULL = scan all SRX (default).
#' @param chipatlas_quiet Suppress per-SRX download messages (default TRUE).
#'
#'
#' @param histone_marks_pdf If TRUE (default), generate a single-page
#'   \code{histone_marks.pdf} showing the chromatin-state landscape at
#'   this locus (active vs. repressive ChIP-Atlas peaks).
#' @param histone_marks Character vector of mark names to render.
#'   Default \code{c("H3K4me3", "H3K27ac", "H3K4me1", "H3K36me3",
#'   "H3K27me3", "H3K9me3")}.
#' @param histone_cell_type Cell-type substring (case- and punctuation-
#'   insensitive) used to filter the top section to matched SRXs.
#'   Default \code{"HEK293T"}.
#' @param histone_max_experiments_matched Cap on cell-type-matched SRX
#'   per mark (default 50).
#' @param histone_max_experiments_all Cap on all-cell-type SRX per
#'   mark in the bottom section (default 50).
#'
#'
#' @param epigenetic_complexes_csv Optional path to
#'   \code{EpiGenes_complexes.csv}. When supplied together with
#'   \code{epigenes_main_csv}, generates an \code{epigenetic_complexes.pdf}
#'   with one page per complex showing all members' zone bars
#'   side-by-side.
#' @param epigenes_main_csv Optional path to \code{EpiGenes_main.csv}
#'   (required if \code{epigenetic_complexes_csv} is set).
#' @param complex_min_detected Minimum detected members for a complex
#'   to get a page (default 2).
#' @param complex_max_pages Cap on number of complex pages (default 50,
#'   ranked by detected-member count).
#'
#'
#' @param save_plots Write PDFs to \code{out_dir} (default TRUE).
#' @param plot_width,plot_height PDF dimensions in inches. Defaults 12 x 8.
#'
#' @return Invisibly, a list with: \code{spatial_df_epi}, \code{zones_df},
#'   \code{chipatlas_peaks}, \code{factors_present}, \code{factors_missing},
#'   \code{kernel_sigma}, \code{zone_frac}, \code{weight_mode}.
#' @export
run_caspex_epigenetic <- function(
    result,
    epigenetic_factors = NULL,
    zone_frac        = 0.3,
    inner_zone_frac  = 0.7,
    centroid_frac    = 0.7,
    kernel_sigma     = NULL,
    weight_mode      = NULL,
    cov_floor        = 0.05,
    edge_guard_frac  = 0.25,
    edge_grna_weight_cap = NULL,
    max_grna_distance    = NULL,
    chipatlas               = TRUE,
    chipatlas_threshold     = "05",
    chipatlas_max_experiments = 100,
    # `special_interest_gene`: same semantics as in run_caspex(). Listed
    # factors bypass the chipatlas_max_experiments cap. Cache is shared with
    # the TF run, so factors already downloaded there are free here.
    # `special_interest_cap`: NULL = scan ALL SRX for special-interest
    # factors. Set to e.g. 250 to take a top-N slice instead.
    special_interest_gene   = NULL,
    special_interest_cap    = NULL,
    chipatlas_quiet         = TRUE,
    upstream         = NULL,
    downstream       = NULL,
    detail_top_n     = 50,
    min_n_regions    = 2,
    subtract_tf_overlap = FALSE,
    # `peak_signal_range`: 2-vector c(min, max) controlling the deck-wide
    # colour-bar limits on the inner-core fill. NULL (default) auto-computes
    # the range from the union of all rendered factors' peak_beta values,
    # so every page in this deck shares one consistent colour scale.
    # Override when stitching multiple decks together (e.g. one shared
    # range across hTERT + MYC for cross-locus comparison) by passing the
    # union of both deck's ranges explicitly.
    peak_signal_range   = NULL,
    # ── Locus-level histone-marks summary page ───────────────────────────
    # When TRUE, fetches ChIP-Atlas peaks for the canonical six histone
    # marks (or whatever `histone_marks` lists) and saves a single-page
    # PDF (`histone_marks.pdf`) showing the chromatin-state landscape at
    # this locus. Top section = peaks from cell-type-matched SRXs only;
    # bottom section = peaks aggregated across all cell-type SRXs (cap
    # `max_experiments_all`). Page is locus-level (factor-agnostic) so the
    # reader sees the chromatin-state context once per locus instead of
    # repeated on every per-factor page.
    histone_marks_pdf       = TRUE,
    histone_marks           = c("H3K4me3", "H3K27ac", "H3K4me1",
                                 "H3K36me3", "H3K27me3", "H3K9me3"),
    histone_cell_type       = "HEK293T",
    histone_max_experiments_matched = 50,
    histone_max_experiments_all     = 50,
    # ── Epigenetic complexes (EpiGenes-derived) deck ─────────────────────
    # When `epigenetic_complexes_csv` and `epigenes_main_csv` are both
    # supplied, an `epigenetic_complexes.pdf` is generated with one page
    # per known epigenetic complex (BAF, NuRF, PRC2, etc.). Each page
    # stacks one row per member protein (detected ones with their zone
    # bars, not-detected ones as labelled empty rows) so the reader sees
    # at a glance whether the members of a complex co-localise at this
    # locus. Filtered to complexes with >= `complex_min_detected`
    # members detected in the proteomics; capped at `complex_max_pages`
    # pages, ranked by detected-member count.
    epigenetic_complexes_csv = NULL,
    epigenes_main_csv        = NULL,
    complex_min_detected     = 2,
    complex_max_pages        = 50,
    out_dir,
    plot_width       = 12,
    plot_height      = 8,
    save_plots       = TRUE) {

  if (missing(out_dir) || is.null(out_dir))
    stop("run_caspex_epigenetic(): `out_dir` is required.")

  # Inherit defaults from the result object so the analysis is consistent
  # with the upstream TF run unless the caller explicitly overrides.
  # Critical for `upstream` / `downstream`: the zone-detection x_grid is
  # built as seq(-upstream, downstream, by = 5), and gRNAs landing outside
  # that range produce empty signal on every factor.  Mackenzie FOXP2 sits
  # at +980..+1194 — well past the legacy +500 default — so without this
  # inheritance every factor returns 0 zones.
  if (is.null(kernel_sigma))
    kernel_sigma <- (result$kernel_sigma %||% 300)
  if (is.null(weight_mode))
    weight_mode  <- (result$weight_mode  %||% "z")
  if (is.null(upstream))
    upstream     <- (result$upstream     %||% 2500)
  if (is.null(downstream))
    downstream   <- (result$downstream   %||% 500)
  if (is.null(max_grna_distance))
    max_grna_distance <- kernel_sigma

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # ── Resolve epigenetic-factor universe + complexes paths from bundled
  # databases when not user-supplied.  Keeps the simple
  # `run_caspex_epigenetic(result, out_dir = ...)` call self-contained;
  # users who want to pin a custom list still pass `epigenetic_factors = ...`
  # (and/or the two complex csv paths) and override the defaults.
  if (is.null(epigenetic_factors)) {
    epi_path <- system.file("extdata/databases/EpiGenes_main.csv",
                             package = "GLproxScape")
    if (nzchar(epi_path) && file.exists(epi_path)) {
      epi_df  <- utils::read.csv(epi_path, stringsAsFactors = FALSE)
      sym_col <- if ("HGNC_symbol" %in% names(epi_df))
        "HGNC_symbol" else names(epi_df)[1]
      epigenetic_factors <- toupper(trimws(epi_df[[sym_col]]))
      epigenetic_factors <- epigenetic_factors[
        nzchar(epigenetic_factors) & !is.na(epigenetic_factors)]
      message("Using bundled epigenetic-factor list (",
              length(epigenetic_factors), " entries) from ",
              basename(epi_path))
    }
  }
  if (is.null(epigenes_main_csv)) {
    p <- system.file("extdata/databases/EpiGenes_main.csv",
                      package = "GLproxScape")
    if (nzchar(p) && file.exists(p)) epigenes_main_csv <- p
  }
  if (is.null(epigenetic_complexes_csv)) {
    p <- system.file("extdata/databases/EpiGenes_complexes.csv",
                      package = "GLproxScape")
    if (nzchar(p) && file.exists(p)) epigenetic_complexes_csv <- p
  }

  message("=== CasPEX Epigenetic Zone Predictor ===")
  message("Gene: ", result$gene_info$name,
          "  |  factors requested: ", length(epigenetic_factors))

  # ── Detection-overlap diagnostic ───────────────────────────────────────
  detected_proteins <- unique(result$long_data$protein)
  factors_present   <- intersect(epigenetic_factors, detected_proteins)
  factors_missing   <- setdiff(epigenetic_factors, detected_proteins)
  message(sprintf(
    "  detected in long_data: %d / %d  (%.1f%%)",
    length(factors_present), length(epigenetic_factors),
    100 * length(factors_present) / max(length(epigenetic_factors), 1)))
  if (length(factors_missing) > 0) {
    msg_n <- min(10, length(factors_missing))
    message("  not detected in this run (first ", msg_n, " of ",
            length(factors_missing), "): ",
            paste(factors_missing[seq_len(msg_n)], collapse = ", "),
            if (length(factors_missing) > msg_n) ", ..." else "")
  }

  # ── Optional dual-class subtraction ────────────────────────────────────
  if (isTRUE(subtract_tf_overlap) && !is.null(result$motif_results)) {
    tf_motif_set <- names(result$motif_results)
    pre <- length(factors_present)
    factors_present <- setdiff(factors_present, tf_motif_set)
    if (pre != length(factors_present)) {
      message("  subtract_tf_overlap=TRUE: dropped ", pre - length(factors_present),
              " factor(s) already on TF deck.")
    }
  }

  if (length(factors_present) == 0) {
    message("  no epigenetic factors with usable data; aborting.")
    return(invisible(NULL))
  }

  # ── Spatial model on the epigenetic-factor subset ──────────────────────
  ld_epi <- result$long_data[result$long_data$protein %in% factors_present, ,
                              drop = FALSE]
  # Force isTF = TRUE locally so run_spatial_model's tfs_only filter does
  # not silently drop these proteins. We're explicitly opting in by name.
  ld_epi$isTF <- TRUE
  spatial_df_epi <- run_spatial_model(
    ld_epi, result$pos_map,
    tfs_only    = FALSE,
    pval_thresh = 0.05,
    min_regions = min_n_regions,
    min_lfc     = 0,
    weight_mode = weight_mode)
  message("  spatial-modelled epigenetic factors: ", nrow(spatial_df_epi))
  if (nrow(spatial_df_epi) == 0) {
    message("  no factors survived spatial filter; aborting.")
    return(invisible(NULL))
  }

  # ── Per-factor zone detection ──────────────────────────────────────────
  message(sprintf(
    "\nDetecting binding zones (kernel \u03c3=%d bp, zone_frac=%.2f)...",
    kernel_sigma, zone_frac))
  x_grid <- seq(-upstream, downstream, by = 5)
  zones_list <- lapply(spatial_df_epi$protein, function(tf) {
    z <- predict_binding_zones_epigenetic(
      tf, ld_epi, result$pos_map,
      x_grid          = x_grid,
      kernel_sigma    = kernel_sigma,
      zone_frac       = zone_frac,
      inner_zone_frac = inner_zone_frac,
      centroid_frac   = centroid_frac,
      weight_mode     = weight_mode,
      cov_floor       = cov_floor,
      edge_guard_frac = edge_guard_frac,
      max_grna_distance    = max_grna_distance,
      edge_grna_weight_cap = edge_grna_weight_cap)
    if (nrow(z) > 0) message("  ", tf, ": ", nrow(z), " zone(s)")
    z
  })
  zones_df <- do.call(rbind, zones_list)
  if (is.null(zones_df)) zones_df <- data.frame()
  message("  Total zones: ", nrow(zones_df),
          if (nrow(zones_df) > 0)
            paste0("  (across ", length(unique(zones_df$tf)), " factor(s))")
          else "")

  factors_with_zones <- if (nrow(zones_df) > 0) unique(zones_df$tf)
                                              else character(0)

  # Top-N selection for the deck (rank by max peak β across each factor's
  # zones — same spirit as detail_top_n on the TF path).
  detail_factors <- factors_with_zones
  if (is.finite(detail_top_n) && length(detail_factors) > detail_top_n) {
    top_betas <- sapply(detail_factors, function(f)
      max(zones_df$peak_beta[zones_df$tf == f]))
    detail_factors <- detail_factors[order(top_betas, decreasing = TRUE)][
      seq_len(detail_top_n)]
  }

  # ── ChIP-Atlas fetch ───────────────────────────────────────────────────
  chipatlas_res <- NULL
  if (isTRUE(chipatlas) && length(detail_factors) > 0 &&
      exists("run_chipatlas_scan", mode = "function")) {
    chipatlas_res <- tryCatch(
      run_chipatlas_scan(
        detail_factors,
        gene_info       = result$gene_info,
        promoter_info   = result$promoter_info,
        upstream        = upstream,
        downstream      = downstream,
        threshold       = chipatlas_threshold,
        max_experiments = chipatlas_max_experiments,
        special_interest_gene = special_interest_gene,
        special_interest_cap  = special_interest_cap,
        quiet           = chipatlas_quiet),
      error = function(e) {
        message("  ChIP-Atlas fetch failed: ", conditionMessage(e))
        NULL
      })
  }

  # ── Per-zone overlap diagnostic ────────────────────────────────────────
  if (nrow(zones_df) > 0) {
    zones_df <- compute_chipatlas_overlap(zones_df, chipatlas_res)
  }

  # ── Per-factor deck PDF ────────────────────────────────────────────────
  if (isTRUE(save_plots) && length(detail_factors) > 0) {
    deck_path <- file.path(out_dir, "epigenetic_zone_deck.pdf")
    # Compute deck-wide peak-signal range so every page shares the same
    # colourbar. Auto-compute from the rendered factors' peak_beta column
    # unless the caller passed one in. Float-pad the lower bound to 0 if
    # the auto range starts above 0 — the colourbar reads more naturally
    # when the low end of the gradient corresponds to "no signal" rather
    # than "lowest detected signal across this run".
    deck_range <- peak_signal_range
    if (is.null(deck_range) && nrow(zones_df) > 0) {
      pk <- zones_df$peak_beta[zones_df$tf %in% detail_factors]
      pk <- pk[is.finite(pk)]
      if (length(pk) > 0) {
        deck_range <- c(0, max(pk, na.rm = TRUE))
      }
    }
    if (!is.null(deck_range))
      message(sprintf("  Deck colour scale: peak signal in [%.3f, %.3f]",
                      deck_range[1], deck_range[2]))
    deck_plots <- lapply(detail_factors, function(tf) {
      plot_epigenetic_zone_deck(
        tf,
        long_data         = ld_epi,
        pos_map           = result$pos_map,
        zones_df          = zones_df[zones_df$tf == tf, , drop = FALSE],
        kernel_sigma      = kernel_sigma,
        upstream          = upstream,
        downstream        = downstream,
        weight_mode       = weight_mode,
        zone_frac         = zone_frac,
        inner_zone_frac   = inner_zone_frac,
        cov_floor         = cov_floor,
        edge_guard_frac   = edge_guard_frac,
        peak_signal_range = deck_range,
        chipatlas_peaks   = if (!is.null(chipatlas_res))
                              chipatlas_res[[tf]] else NULL)
    })
    # `.safe_pdf` (defined module-level in caspex_analysis.R) wraps
    # the multi-page write in a tryCatch so a stale-file-lock failure
    # ("dev.off() : write failed") doesn't abort the whole run.
    if (exists(".safe_pdf", mode = "function")) {
      .safe_pdf(deck_path,
                width  = plot_width,
                height = max(5, plot_height * 0.6),
                plots  = deck_plots)
    } else {
      pdf(deck_path,
          width  = plot_width,
          height = max(5, plot_height * 0.6))
      for (pl in deck_plots) print(pl)
      invisible(dev.off())
    }
    message("  Zone deck saved: ", deck_path)
  }

  # ── Locus-level histone-marks page ────────────────────────────────────
  histone_data <- NULL
  if (isTRUE(histone_marks_pdf) && length(histone_marks) > 0 &&
      exists("fetch_histone_peaks_for_locus", mode = "function")) {
    message(sprintf(
      "\nFetching histone marks (%d marks, cell type = %s)...",
      length(histone_marks), histone_cell_type %||% "(no filter)"))
    histone_data <- tryCatch(
      fetch_histone_peaks_for_locus(
        marks         = histone_marks,
        gene_info     = result$gene_info,
        promoter_info = result$promoter_info,
        upstream      = upstream,
        downstream    = downstream,
        cell_type     = histone_cell_type,
        max_experiments_matched = histone_max_experiments_matched,
        max_experiments_all     = histone_max_experiments_all,
        threshold     = chipatlas_threshold,
        quiet         = chipatlas_quiet),
      error = function(e) {
        message("  Histone fetch failed: ", conditionMessage(e))
        NULL
      })
    if (isTRUE(save_plots) && !is.null(histone_data)) {
      hist_path <- file.path(out_dir, "histone_marks.pdf")
      hist_plot <- plot_histone_marks_locus(
        histone_data = histone_data,
        pos_map      = result$pos_map,
        gene_info    = result$gene_info,
        upstream     = upstream,
        downstream   = downstream)
      if (exists(".safe_pdf", mode = "function")) {
        .safe_pdf(hist_path,
                  width  = plot_width,
                  height = max(6, plot_height * 0.8),
                  plots  = list(hist_plot))
      } else {
        pdf(hist_path,
            width  = plot_width,
            height = max(6, plot_height * 0.8))
        print(hist_plot)
        invisible(dev.off())
      }
      message("  Histone-marks page saved: ", hist_path)
    }
  }

  # ── Epigenetic-complexes deck ─────────────────────────────────────────
  complexes_pdf <- NULL
  if (!is.null(epigenetic_complexes_csv) &&
      !is.null(epigenes_main_csv) &&
      file.exists(epigenetic_complexes_csv) &&
      file.exists(epigenes_main_csv) &&
      exists("load_epigenes_complexes", mode = "function")) {

    detected_proteins <- unique(result$long_data$protein)
    message(sprintf(
      "\nLoading epigenetic complexes from %s ...",
      basename(epigenetic_complexes_csv)))
    complexes <- tryCatch(
      load_epigenes_complexes(
        complex_csv = epigenetic_complexes_csv,
        main_csv    = epigenes_main_csv,
        detected_proteins = detected_proteins),
      error = function(e) {
        message("  Complex loader failed: ", conditionMessage(e))
        NULL
      })

    if (!is.null(complexes) && length(complexes) > 0) {
      n_total <- length(complexes)
      keep <- vapply(complexes,
        function(c) length(c$members_detected) >= complex_min_detected,
        logical(1))
      complexes <- complexes[keep]
      message(sprintf(
        "  %d / %d complexes have >= %d detected member(s)",
        length(complexes), n_total, complex_min_detected))
      if (length(complexes) > complex_max_pages)
        complexes <- complexes[seq_len(complex_max_pages)]

      if (length(complexes) > 0 && isTRUE(save_plots)) {
        # Reuse the deck-wide peak-signal range (same colour scale as
        # the per-factor pages) so a "peak signal" colour means the
        # same thing whether you're looking at a per-factor page or a
        # complex page.
        complex_pages <- lapply(complexes, function(ci) {
          plot_epigenetic_complex_locus(
            complex_info     = ci,
            zones_df         = zones_df[zones_df$tf %in% ci$members_detected,
                                         , drop = FALSE],
            sig_df_per_member = NULL,   # reserved for future enrichment overlay
            pos_map          = result$pos_map,
            gene_info        = result$gene_info,
            upstream         = upstream,
            downstream       = downstream,
            kernel_sigma     = kernel_sigma,
            peak_signal_range = if (exists("deck_range",
                                             inherits = FALSE)) deck_range
                                 else NULL)
        })
        # Drop NULLs (defensive)
        complex_pages <- Filter(Negate(is.null), complex_pages)
        if (length(complex_pages) > 0) {
          complexes_pdf <- file.path(out_dir, "epigenetic_complexes.pdf")
          if (exists(".safe_pdf", mode = "function")) {
            .safe_pdf(complexes_pdf,
                      width  = plot_width,
                      height = max(6, plot_height * 0.7),
                      plots  = complex_pages)
          } else {
            pdf(complexes_pdf,
                width  = plot_width,
                height = max(6, plot_height * 0.7))
            for (pl in complex_pages) print(pl)
            invisible(dev.off())
          }
          message("  Complexes deck saved: ", complexes_pdf,
                  " (", length(complex_pages), " complex(es))")

          # CSV summary of the deck (one row per complex, listing
          # detected vs. not-detected members, for downstream
          # filtering / inspection).
          summary_df <- do.call(rbind, lapply(complexes, function(ci) {
            data.frame(
              complex_name     = ci$complex_name,
              group_name       = ci$group_name,
              n_members_total  = ci$members_total,
              n_members_detected = length(ci$members_detected),
              members_detected = paste(ci$members_detected, collapse = ";"),
              members_not_detected = paste(
                setdiff(ci$members_hgnc, ci$members_detected),
                collapse = ";"),
              function_        = ci$function_,
              target           = ci$target,
              stringsAsFactors = FALSE)
          }))
          write.csv(summary_df,
                    file.path(out_dir, paste0(result$gene_info$name,
                              "_epigenetic_complexes.csv")),
                    row.names = FALSE)
        }
      }
    }
  }

  # ── CSV outputs ────────────────────────────────────────────────────────
  gene_name <- result$gene_info$name
  if (nrow(zones_df) > 0) {
    zones_df$gene            <- gene_name
    zones_df$kernel_sigma    <- kernel_sigma
    zones_df$zone_frac       <- zone_frac
    zones_df$inner_zone_frac <- inner_zone_frac
    zones_df$centroid_frac   <- centroid_frac
    zones_path <- file.path(out_dir,
      paste0(gene_name, "_epigenetic_zones.csv"))
    write.csv(zones_df[order(zones_df$tf, -zones_df$peak_beta), ],
              zones_path, row.names = FALSE)
    message("  Zones CSV: ", zones_path)
  }
  spatial_path <- file.path(out_dir,
    paste0(gene_name, "_epigenetic_spatial.csv"))
  write.csv(spatial_df_epi, spatial_path, row.names = FALSE)
  message("  Spatial CSV: ", spatial_path)

  # ── Run summary ────────────────────────────────────────────────────────
  message("\n--- Epigenetic run summary --------------------------------------")
  message(" Factors requested  : ", length(epigenetic_factors))
  message(" Detected in run    : ", length(factors_present))
  message(" Spatial-modelled   : ", nrow(spatial_df_epi))
  message(" Factors with zones : ", length(factors_with_zones))
  message(" Total zones        : ", nrow(zones_df))
  if (!is.null(chipatlas_res)) {
    n_with_peaks <- sum(vapply(chipatlas_res,
      function(x) !is.null(x) && nrow(x) > 0, logical(1)))
    message(" ChIP-Atlas         : ", n_with_peaks, "/",
            length(chipatlas_res), " factor(s) had peaks in window")
  }
  if (!is.null(histone_data)) {
    n_marks_matched <- sum(vapply(histone_data$matched,
      function(x) !is.null(x) && nrow(x) > 0, logical(1)))
    n_marks_all <- sum(vapply(histone_data$all,
      function(x) !is.null(x) && nrow(x) > 0, logical(1)))
    message(" Histone marks      : matched=", n_marks_matched,
            "/", length(histone_data$marks),
            "  all=",     n_marks_all,
            "/", length(histone_data$marks),
            " (cell_type=", histone_data$cell_type %||% "NULL", ")")
  }
  message("------------------------------------------------------------------")
  message("Done.")

  invisible(list(
    spatial_df_epi   = spatial_df_epi,
    zones_df         = zones_df,
    chipatlas_peaks  = chipatlas_res,
    histone_data     = histone_data,
    factors_present  = factors_present,
    factors_missing  = factors_missing,
    factors_with_zones = factors_with_zones,
    kernel_sigma     = kernel_sigma,
    zone_frac        = zone_frac,
    inner_zone_frac  = inner_zone_frac,
    centroid_frac    = centroid_frac,
    peak_signal_range = if (exists("deck_range",
                                    inherits = FALSE)) deck_range else NULL,
    weight_mode      = weight_mode))
}

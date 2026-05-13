# =============================================================================
# Engine — gene lookup, gRNA matching, signal modelling, motif scan, and
# the deconvolution pipeline (run_caspex). Plotting helpers and the
# top-level user-facing entry points live here too. Internal utilities
# (`%||%`, `COLS`, `.safe_pdf`) are in R/utils-internal.R; ChIP-Atlas
# integration is in R/caspex_chipatlas.R; chromatin-factor zone path is
# in R/caspex_epigenetic.R; diagnostic plots in R/caspex_extras.R.
# =============================================================================

# Sticky flag historically used by the script-loading model to track
# whether caspex_chipatlas.R had been sourced. In package context every
# R/ file is loaded by the package machinery, so this is always TRUE.
# Kept as a constant binding (not @export) so existing internal guards
# `isTRUE(.caspex_chipatlas_loaded)` continue to evaluate correctly
# without code churn.
.caspex_chipatlas_loaded <- TRUE

# ── Palette & theme ───────────────────────────────────────────────────────────
# COLS lives in R/utils-internal.R.

#' GLproxScape default ggplot2 theme
#'
#' Minimal theme tuned for the package's plot decks (bold titles,
#' grey-50 subtitles, no minor grid lines, light-grey major grid).
#'
#' @param base_size base font size in points (default 11).
#' @return A ggplot2 theme object.
#' @noRd
theme_caspex <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1),
      plot.subtitle    = element_text(color = "grey50", size = base_size - 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92"),
      legend.position  = "right",
      strip.text       = element_text(face = "bold")
    )
}

# =============================================================================
# SECTION 1: Data loading
# =============================================================================

#' Load CasPEX enrichment files for all regions
#'
#' @param data_files   Named character vector: names = region IDs (e.g. "R1"),
#'                     values = file paths.
#' @param tf_universe  Optional character vector of HGNC symbols (e.g. read
#'                     from TFLibrary.txt). When supplied, the `isTF`
#'                     column is set by membership in this vector. When
#'                     NULL (default), falls back to a `TFDatabase` /
#'                     `is_tf` column in the input file if present, else
#'                     `isTF = FALSE` for all rows. Case-insensitive.
#' @param epi_universe Optional character vector of HGNC symbols (e.g.
#'                     read from EpiGenes_main.csv). When supplied, the
#'                     `isEpi` column is set by membership. NULL → FALSE.
#'                     Case-insensitive.
#' @return Long-format data frame with columns: protein, region, lfc, pval,
#'                     t_stat, isTF, isEpi.
#' @noRd
load_region_data <- function(data_files,
                              tf_universe  = NULL,
                              epi_universe = NULL) {
  tf_set  <- if (!is.null(tf_universe))  toupper(tf_universe[nzchar(tf_universe)])  else character(0)
  epi_set <- if (!is.null(epi_universe)) toupper(epi_universe[nzchar(epi_universe)]) else character(0)
  out <- vector("list", length(data_files))
  for (i in seq_along(data_files)) {
    region <- names(data_files)[i]
    path   <- data_files[[i]]
    if (!file.exists(path)) stop("File not found: ", path)
    df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)

    lfc_col  <- grep("^logFC$|logfc|log2fc|log_fc", names(df), ignore.case = TRUE, value = TRUE)[1]
    pval_col <- grep("P\\.Value|p_value|pvalue|^p$", names(df), ignore.case = TRUE, value = TRUE)[1]
    name_col <- grep("^name$|gene|protein|symbol", names(df), ignore.case = TRUE, value = TRUE)[1]
    tf_col   <- grep("TFDatabase|tf_db|is_tf", names(df), ignore.case = TRUE, value = TRUE)[1]
    t_col    <- grep("^t$|t[._]stat|moderated[._]t", names(df),
                     ignore.case = TRUE, value = TRUE)[1]

    if (any(is.na(c(lfc_col, pval_col, name_col)))) {
      stop("Could not identify required columns in ", path,
           "\nFound: ", paste(names(df), collapse=", "))
    }

    proteins_upper <- toupper(df[[name_col]])

    # Resolve isTF: tf_universe wins if supplied; else fall back to a
    # TFDatabase-style column in the file if present; else FALSE.
    isTF_vec <- if (length(tf_set) > 0) {
      proteins_upper %in% tf_set
    } else if (!is.na(tf_col)) {
      tolower(df[[tf_col]]) == "exist"
    } else {
      rep(FALSE, nrow(df))
    }

    isEpi_vec <- if (length(epi_set) > 0) {
      proteins_upper %in% epi_set
    } else {
      rep(FALSE, nrow(df))
    }

    out[[i]] <- data.frame(
      protein = df[[name_col]],
      region  = region,
      lfc     = as.numeric(df[[lfc_col]]),
      pval    = as.numeric(df[[pval_col]]),
      t_stat  = if (!is.na(t_col)) as.numeric(df[[t_col]]) else NA_real_,
      isTF    = isTF_vec,
      isEpi   = isEpi_vec,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

#' Load the CasPEX inputs manifest from an inputs/ folder
#'
#' Parses a tab-separated manifest file (default `inputs/grnas.tsv`) listing
#' per-region gRNA protospacer and logFC file, and returns the `grnas` +
#' `data_files` named vectors expected by `run_caspex()`. This replaces the
#' old hand-written vectors in the driver script.
#'
#' The manifest must have a header row with columns: `region`, `sequence`,
#' `data_file`. Lines beginning with `#` are treated as comments; blank
#' lines are ignored. `data_file` is resolved relative to the inputs
#' directory — so the manifest and the Region*.txt files live side by side.
#'
#' A `sequence` cell with value NA / "NA" / empty is kept as NA_character_
#' so that region still contributes its logFC data but no gRNA cut-site
#' tick is drawn.
#'
#' @param inputs_dir Path to the inputs directory (default `"inputs"`).
#' @param manifest   Manifest filename within `inputs_dir`
#'                   (default `"grnas.tsv"`).
#' @return List with:
#'   * `grnas`: named character vector (names = region IDs, values =
#'     protospacers). NA for regions without a gRNA sequence.
#'   * `data_files`: named character vector (names = region IDs, values =
#'     absolute paths to per-region logFC files).
#'   * `inputs_dir`: resolved path to the inputs directory.
#'   * `manifest_path`: resolved path to the manifest file.
#'
#' @examples
#' # Use the bundled FOXP2 example dataset
#' inputs_dir <- system.file("extdata/examples/foxp2_mackenzie",
#'                           package = "GLproxScape")
#' inputs <- load_caspex_inputs(inputs_dir)
#' str(inputs)
#'
#' \dontrun{
#' # Then feed into the pipeline:
#' run_caspex(gene       = "FOXP2",
#'            transcript = "ENST00000901759",
#'            grnas      = inputs$grnas,
#'            data_files = inputs$data_files,
#'            upstream   = 200,
#'            downstream = 2000,
#'            out_dir    = tempfile("foxp2_run_"))
#' }
#' @export
load_caspex_inputs <- function(inputs_dir = "inputs",
                                manifest   = "grnas.tsv") {
  if (!dir.exists(inputs_dir))
    stop("Inputs directory not found: ", inputs_dir,
         "\nExpected a folder containing `", manifest, "` and the per-region ",
         "logFC files (e.g. Region1.txt, Region2.txt, ...).")

  manifest_path <- file.path(inputs_dir, manifest)
  if (!file.exists(manifest_path))
    stop("Manifest not found: ", manifest_path,
         "\nExpected a tab-separated file with header `region\\tsequence\\tdata_file`.")

  # Read, stripping comment lines (#...) and blank rows. read.delim handles
  # the tab-split; comment.char = "#" drops commented lines including the
  # provenance header block.
  df <- read.delim(manifest_path,
                   header           = TRUE,
                   sep              = "\t",
                   comment.char     = "#",
                   stringsAsFactors = FALSE,
                   strip.white      = TRUE,
                   blank.lines.skip = TRUE)

  required <- c("region", "sequence", "data_file")
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0)
    stop("Manifest ", manifest_path, " is missing required column(s): ",
         paste(missing, collapse = ", "),
         "\nFound: ", paste(names(df), collapse = ", "))

  # Drop any fully-blank rows that survived read.delim
  df <- df[nzchar(trimws(df$region)), , drop = FALSE]
  if (nrow(df) == 0)
    stop("Manifest ", manifest_path, " has no data rows after stripping comments.")

  # Validate region uniqueness
  dup <- df$region[duplicated(df$region)]
  if (length(dup) > 0)
    stop("Duplicate region IDs in manifest: ", paste(unique(dup), collapse = ", "))

  # Normalise sequence: "", "NA", "na" → NA_character_
  seq_raw <- trimws(df$sequence)
  seq_raw[seq_raw == "" | toupper(seq_raw) == "NA"] <- NA_character_

  # Resolve data_file paths relative to inputs_dir and verify they exist
  df$data_file <- trimws(df$data_file)
  if (any(!nzchar(df$data_file)))
    stop("Manifest ", manifest_path,
         " has row(s) with empty data_file. Every region needs a logFC file.")
  data_paths <- file.path(inputs_dir, df$data_file)
  not_found  <- data_paths[!file.exists(data_paths)]
  if (length(not_found) > 0)
    stop("data_file(s) referenced by manifest not found under ", inputs_dir, ":\n  ",
         paste(not_found, collapse = "\n  "))

  grnas      <- setNames(seq_raw, df$region)
  data_files <- setNames(data_paths, df$region)

  message("Loaded CasPEX inputs manifest: ", manifest_path)
  message("  Regions: ", nrow(df),
          " | with gRNA: ", sum(!is.na(grnas)),
          " | without gRNA: ", sum(is.na(grnas)))

  list(grnas         = grnas,
       data_files    = data_files,
       inputs_dir    = normalizePath(inputs_dir, mustWork = TRUE),
       manifest_path = normalizePath(manifest_path, mustWork = TRUE))
}

# =============================================================================
# SECTION 2: Ensembl gene lookup + promoter sequence
# =============================================================================

#' Query Ensembl REST API
#' @param path (see function body).
#' @param params (see function body).
#' @param accept (see function body).
#' @noRd
ensembl_get <- function(path, params = list(), accept = "application/json") {
  url <- modify_url("https://rest.ensembl.org", path = path, query = params)
  # IMPORTANT: do NOT hardcode Content-Type to application/json on every
  # request. Content-Type describes a request BODY (which a GET doesn't
  # have), but Ensembl's REST router was treating the hardcoded JSON
  # content-type as a hint to serve JSON even when Accept said text/plain
  # — wrapping /sequence/region responses in a JSON envelope of the form
  #   {"MOLECULE":"DNA","QUERY":"...","seq":"ATCG...",...}
  # which then got blindly fed into gRNA matching as if it were raw
  # sequence, producing a uniform Δbp drift in pos_map between runs as
  # the JSON serializer reordered fields. Send Accept only.
  res <- GET(url, add_headers(Accept = accept), timeout(20))
  if (status_code(res) == 200) {
    if (accept == "text/plain") return(trimws(content(res, "text", encoding = "UTF-8")))
    return(content(res, "parsed", simplifyVector = FALSE))
  }
  warning("Ensembl API returned ", status_code(res), " for: ", path)
  NULL
}

#' Look up gene coordinates from Ensembl
#'
#' By default the TSS is anchored on the Ensembl canonical transcript
#' (`transcript = "canonical"`). This is deterministic across REST calls and
#' matches the "canonical promoter" frame used in most TF / ChIP studies.
#' Pass `transcript = "ENST..."` to pin to a specific transcript ID, or
#' `transcript = NA` to fall back to the legacy gene-level outermost
#' transcript boundary (Ensembl's gene$start / gene$end). The legacy
#' boundary is the union over all annotated transcripts and can drift
#' between Ensembl releases / load-balanced cluster nodes — which produces
#' uniform TSS-frame shifts (e.g. all gRNA positions translate by the same
#' Δ bp), so it is no longer the default.
#'
#' @param gene_name  HGNC symbol (e.g. "MYC")
#' @param species    Ensembl species string (default "homo_sapiens")
#' @param transcript Transcript anchor for the TSS:
#'                   * "canonical" (default): use the Ensembl canonical transcript
#'                   * "ENSTxxxxxxxxx": pin to the named transcript
#'                   * NA: legacy gene-level start/end (subject to drift; not recommended)
#' @return List: name, chr, strand, tss, start, end, species, transcript_id
#' @export
lookup_gene <- function(gene_name, species = "homo_sapiens",
                         transcript = "canonical") {
  # NA is the explicit escape hatch into legacy gene-level coords.
  # Treat it the same as the old `NULL` default (skips the transcript
  # resolution branch below).
  if (length(transcript) == 1 && is.na(transcript)) transcript <- NULL
  message("  Looking up gene: ", gene_name,
          if (!is.null(transcript))
            paste0(" (transcript=", transcript, ")") else "")
  js <- ensembl_get(paste0("/lookup/symbol/", species, "/", gene_name),
                    params = list(expand = 1, `content-type` = "application/json"))
  if (is.null(js)) stop("Gene '", gene_name, "' not found in Ensembl (", species, ")")

  start         <- js$start
  end           <- js$end
  transcript_id <- NA_character_

  if (!is.null(transcript)) {
    tx_list <- js$Transcript
    if (is.null(tx_list) || length(tx_list) == 0)
      stop("No Transcript field in Ensembl lookup for ", gene_name,
           " \u2014 cannot resolve transcript='", transcript, "'.")
    chosen <- NULL
    if (identical(transcript, "canonical")) {
      for (tx in tx_list) {
        if (isTRUE(tx$is_canonical == 1)) { chosen <- tx; break }
      }
      if (is.null(chosen))
        stop("No canonical transcript flagged for ", gene_name,
             " (", length(tx_list), " transcripts returned).")
    } else {
      for (tx in tx_list) {
        if (identical(tx$id, transcript)) { chosen <- tx; break }
      }
      if (is.null(chosen))
        stop("Transcript '", transcript, "' not found among ",
             length(tx_list), " transcripts for ", gene_name, ". ",
             "Available IDs: ",
             paste(vapply(tx_list, function(tx) tx$id, character(1)),
                   collapse = ", "))
    }
    message("    using transcript ", chosen$id,
            " (", chosen$start, "-", chosen$end, ")")
    start         <- chosen$start
    end           <- chosen$end
    transcript_id <- chosen$id
  }

  # Attach species so downstream consumers (fetch_promoter_seq) can hit the
  # correct Ensembl assembly without being passed the parameter separately.
  list(
    name          = js$display_name,
    chr           = js$seq_region_name,
    strand        = js$strand,
    tss           = if (js$strand == 1) start else end,
    start         = start,
    end           = end,
    species       = species,
    transcript_id = transcript_id
  )
}

#' Fetch promoter sequence from Ensembl (strand-aware, TSS-centred)
#'
#' @param gene_info  Output of lookup_gene()
#' @param upstream   bp upstream of TSS to include
#' @param downstream bp downstream of TSS to include
#' @param species    Ensembl species token for the `/sequence/region/<species>`
#'                   endpoint. Defaults to `gene_info$species` (set by
#'                   lookup_gene()), falling back to "homo_sapiens" for
#'                   backward compatibility with any hand-built gene_info
#'                   lists that pre-date the species field. Previously this
#'                   was hardcoded to "human", which silently routed every
#'                   non-human gene to the human genome.
#' @return List: seq (character), tss_offset (int), seq_start, seq_end, strand
#' @export
fetch_promoter_seq <- function(gene_info, upstream = 2500, downstream = 500,
                                species = NULL) {
  # Resolve species: explicit arg wins, else pull from gene_info (set by
  # lookup_gene()), else default to human for legacy call sites. Kept as an
  # `if`-chain rather than `%||%` because caspex_analysis.R doesn't declare
  # that operator and we shouldn't force an R ≥ 4.4 dependency.
  if (is.null(species)) species <- gene_info$species
  if (is.null(species) || !nzchar(species)) species <- "homo_sapiens"
  message("  Fetching promoter sequence (", upstream, " up / ", downstream,
          " down, species=", species, ")...")
  tss    <- gene_info$tss
  strand <- gene_info$strand
  chr    <- gene_info$chr

  if (strand == 1) {
    seq_start <- max(1, tss - upstream)
    seq_end   <- tss + downstream
  } else {
    seq_start <- max(1, tss - downstream)
    seq_end   <- tss + upstream
  }

  seq <- ensembl_get(
    paste0("/sequence/region/", species, "/", chr, ":",
           seq_start, "..", seq_end, ":", strand),
    accept = "text/plain"
  )
  if (is.null(seq) || nchar(seq) < 10)
    stop("Failed to fetch sequence for ", gene_info$name,
         " (species=", species, ")")

  message("  Got ", nchar(seq), " bp")
  list(seq        = toupper(seq),
       tss_offset = upstream,        # TSS is at position upstream + 1 (1-based)
       seq_start  = seq_start,
       seq_end    = seq_end,
       strand     = strand)
}

# =============================================================================
# SECTION 3: gRNA → cut site matching
# =============================================================================

#' Reverse complement of a DNA string
#' @param s (see function body).
#' @noRd
rc <- function(s) {
  map <- c(A="T",T="A",G="C",C="G",N="N")
  chars <- strsplit(toupper(s), "")[[1]]
  paste(rev(map[chars]), collapse = "")
}

#' Strip PAM from protospacer if present
#'
#' Handles common cases:
#'   3' NGG PAM  (SpCas9)
#'   5' CCN  PAM (reverse complement form)
#' @param g (see function body).
#' @noRd
strip_pam <- function(g) {
  g <- toupper(trimws(g))
  g <- sub("[ACGT]GG$", "", g)   # 3' NGG
  g <- sub("^CC[ACGT]", "", g)   # 5' CCN
  g
}

#' Find all occurrences of a pattern in a string (1-based positions)
#' @param pattern (see function body).
#' @param text (see function body).
#' @noRd
str_find_all <- function(pattern, text) {
  hits <- c()
  start <- 1
  repeat {
    m <- regexpr(pattern, substr(text, start, nchar(text)), fixed = TRUE)
    if (m == -1) break
    hits  <- c(hits, start + m - 1)
    start <- start + m   # advance past current hit
  }
  hits
}

#' Match a single gRNA protospacer to the promoter and return TSS-relative cut site(s)
#'
#' Cas9 cuts between nt17 and nt18 of the 20nt protospacer (3 bp upstream of PAM).
#' For a + strand hit at position p: cut = p + 16
#' For a - strand hit at position p: cut = p + len(guide) - 4
#'
#' @param grna_seq       Raw protospacer string (PAM optional)
#' @param promoter_info  Output of fetch_promoter_seq()
#' @return List: positions (TSS-relative bp vector), strands, n_hits, status message
#' @noRd
match_grna <- function(grna_seq, promoter_info) {
  g <- strip_pam(grna_seq)
  if (nchar(g) < 15) return(list(positions = NA, strands = NA, n_hits = 0,
                                  msg = "gRNA too short after PAM removal"))
  prom  <- promoter_info$seq
  g_rc  <- rc(g)
  tss_i <- promoter_info$tss_offset + 1   # 1-based index of TSS in seq

  fwd_hits <- str_find_all(g,    prom)
  rev_hits <- str_find_all(g_rc, prom)

  positions <- c()
  strands   <- c()

  for (p in fwd_hits) {
    cut <- p + 16                  # 0-based between nt17 and nt18
    positions <- c(positions, cut - tss_i)
    strands   <- c(strands, "+")
  }
  for (p in rev_hits) {
    cut <- p + nchar(g) - 4        # symmetric cut on rev strand
    positions <- c(positions, cut - tss_i)
    strands   <- c(strands, "-")
  }

  n <- length(positions)
  if (n == 0) return(list(positions = NA, strands = NA, n_hits = 0,
                           msg = paste0("No match for: ", g)))
  msg <- if (n == 1) sprintf("%+d bp (%s strand)", positions[1], strands[1]) else
    sprintf("%d hits; using first at %+d bp", n, positions[1])
  list(positions = positions, strands = strands, n_hits = n, msg = msg)
}

#' Match all gRNAs and return a named position vector
#'
#' @param grnas Named character vector: names = region IDs, values = protospacers
#' @param promoter_info Output of fetch_promoter_seq()
#' @return Named numeric vector of TSS-relative cut sites (NA for unmatched)
#' @noRd
match_all_grnas <- function(grnas, promoter_info) {
  message("\nMatching gRNA sequences:")
  pos <- setNames(rep(NA_real_, length(grnas)), names(grnas))
  for (r in names(grnas)) {
    g <- grnas[[r]]
    if (is.na(g) || nchar(trimws(g)) == 0) { message("  ", r, ": skipped (empty)"); next }
    result <- match_grna(g, promoter_info)
    message("  ", r, ": ", result$msg)
    if (result$n_hits > 0) pos[[r]] <- result$positions[1]
  }
  pos
}

# =============================================================================
# SECTION 4: Spatial triangulation model
# =============================================================================

#' Signed z-score derived from a two-sided p-value, signed by the effect direction.
#'
#' Inverts p = 2 * (1 - \Phi(|z|))  ->  |z| = \Phi^{-1}(1 - p/2).
#'
#' Numerical note: we evaluate the inverse via the *upper tail* form
#'   qnorm(p/2, lower.tail = FALSE)
#' rather than the algebraically-equivalent qnorm(1 - p/2).  For p ≲ 1e-16 the
#' subtraction `1 - p/2` rounds to exactly 1.0 in double precision, and
#' qnorm(1) = Inf — which is why earlier versions of this function had to
#' clip pval at 1e-16, collapsing every "very significant" region (p<=1e-16)
#' onto the same |z| ≈ 8.13.  The upper-tail form avoids the catastrophic
#' cancellation entirely, so we can push the floor down to the smallest
#' representable double and let z range up to ~37.5.  This matters for
#' synthesized p-value pipelines that hard-cap at 1e-300 (e.g. the Myers
#' 2018 reanalysis loader) — without this fix, dozens of strongly enriched
#' proteins per region all received identical z weights and the spatial
#' deconvolution lost their relative ranking.
#'
#' `pval_floor` exists *only* to guard against p == 0 (which would give Inf).
#'
#' Rationale: a z-score is the "signal per standard error" — the statistically
#' principled weight for inverse-variance-weighted spatial means.  It compresses
#' extreme p-values (p = 1e-10 -> z ~ 6.4, not 10), so one spectacular outlier
#' cannot dominate a centroid computed from a handful of regions.
#' @param lfc (see function body).
#' @param pval (see function body).
#' @param pval_floor (see function body).
#' @noRd
signed_z_from_p <- function(lfc, pval, pval_floor = .Machine$double.xmin) {
  pval <- pmax(pmin(pval, 1), pval_floor)
  z    <- qnorm(pval / 2, lower.tail = FALSE)   # always non-negative; stable in upper tail
  s    <- sign(lfc)
  s[is.na(s) | s == 0] <- 1                     # treat NA/zero-effect as positive
  s * z
}

#' Compute per-region weights given a data frame with columns lfc, pval, and
#' optionally t_stat (moderated t from limma::topTable).
#'
#' Modes:
#'   "z"          — signed z from p-value.  **Default and recommended.**
#'                  This is the approach used in the original ATP7B run; it
#'                  keeps the weight scale consistent across datasets that
#'                  may or may not carry a `t_stat` column. See note below.
#'   "mod_t"      — moderated t if t_stat column is present & non-NA,
#'                  else signed z-score from p-value. Kept for backward
#'                  compatibility / forensic comparison; not recommended
#'                  for production runs because the moderated-t scale is
#'                  much wider than z (large |t| values inflate the
#'                  smoothed-signal magnitude in build_caspex_signal()
#'                  without changing the rank-ordering of binding events).
#'   "lfc_x_negp" — legacy logFC x -log10(p).  Kept for backward compat.
#'   "lfc_pos"    — pmax(lfc, 0).  Effect size only, no significance.
#'   "lfc_signed" — raw lfc (allows negative).
#' @param df (see function body).
#' @param mode (see function body).
#' @noRd
compute_region_weight <- function(df, mode = c("z", "mod_t", "lfc_x_negp",
                                                "lfc_pos", "lfc_signed")) {
  mode <- match.arg(mode)
  has_t <- "t_stat" %in% names(df) && any(!is.na(df$t_stat))
  switch(mode,
    mod_t      = if (has_t) as.numeric(df$t_stat)
                 else        signed_z_from_p(df$lfc, df$pval),
    z          = signed_z_from_p(df$lfc, df$pval),
    lfc_x_negp = pmax(df$lfc, 0) *
                   (-log10(pmax(df$pval, 1e-300) + 1e-10)),
    lfc_pos    = pmax(df$lfc, 0),
    lfc_signed = as.numeric(df$lfc)
  )
}

#' Compute TSS-relative binding centroid for a single protein
#'
#' Weights per region are chosen via `weight_mode` (see compute_region_weight).
#' For centroid estimation only the *positive* part of the weight is used
#' (enrichment only).  This preserves the original interpretation ("where on
#' the promoter is this TF enriched?") while swapping in a statistically
#' sounder quantity.
#'
#' Centroid  = Σ(w⁺ × pos) / Σ(w⁺)
#' Spread    = weighted SD of positions around the centroid
#' Composite = mean(w⁺) × log1p(n_regions)
#'
#' @param protein_name  Character
#' @param long_data     Long-format data frame from load_region_data()
#' @param pos_map       Named numeric vector from match_all_grnas()
#' @param pval_thresh   Significance cutoff
#' @param weight_mode   Passed to compute_region_weight() (default "z")
#' @return List with centroid, spread, composite, n_regions, etc.; or NULL
#' @param min_lfc (see function body).
#' @noRd
compute_spatial <- function(protein_name, long_data, pos_map,
                             pval_thresh = 0.05,
                             min_lfc     = 0,
                             weight_mode = "z") {
  df <- long_data[long_data$protein == protein_name, ]
  df <- df[df$region %in% names(pos_map), ]
  df$pos <- pos_map[df$region]
  df <- df[!is.na(df$pos) & !is.na(df$lfc), ]
  if (nrow(df) == 0) return(NULL)

  df$lfc <- pmax(df$lfc, min_lfc)
  w_raw  <- compute_region_weight(df, mode = weight_mode)
  df$w   <- pmax(w_raw, 0)              # enrichment only
  # Floor to keep spread finite if every weight is ~0
  total_w <- sum(df$w)
  if (total_w < 1e-8) return(NULL)

  centroid <- sum(df$pos * df$w) / total_w
  spread   <- sqrt(sum(df$w * (df$pos - centroid)^2) / total_w)
  spread   <- max(spread, 30)

  list(
    protein   = protein_name,
    isTF      = any(df$isTF),
    centroid  = centroid,
    spread    = spread,
    composite = mean(df$w) * log1p(nrow(df)),
    n_regions = nrow(df),
    mean_lfc  = mean(df$lfc),
    min_pval  = min(df$pval, na.rm = TRUE),
    sig_any   = any(df$pval <= pval_thresh, na.rm = TRUE),
    regions   = df
  )
}

#' Run spatial model on all proteins (or TFs only)
#'
#' @param long_data   From load_region_data()
#' @param pos_map     From match_all_grnas()
#' @param tfs_only    Logical: restrict to TFDatabase == "exist"
#' @param pval_thresh p-value threshold
#' @param min_regions Only keep proteins detected in >= this many regions
#' @param min_lfc     Floor applied to logFC before weighting (default 0 =
#'                    disabled, equivalent to legacy behaviour). Forwarded to
#'                    compute_spatial() so callers can clip mildly negative
#'                    logFC values without bypassing the wrapper.
#' @param weight_mode Region-weight mode; see compute_region_weight()
#' @return Data frame of spatial predictions, sorted by composite score
#' @noRd
run_spatial_model <- function(long_data, pos_map,
                               tfs_only    = TRUE,
                               pval_thresh = 0.05,
                               min_regions = 2,
                               min_lfc     = 0,
                               weight_mode = "z") {
  message("\nRunning spatial model...")
  if (tfs_only) long_data <- long_data[long_data$isTF, ]
  proteins <- unique(long_data$protein)
  has_t    <- "t_stat" %in% names(long_data) && any(!is.na(long_data$t_stat))
  message("  Proteins to model: ", length(proteins))
  message("  Weight mode     : ", weight_mode,
          if (weight_mode == "mod_t")
            paste0("  (using ", if (has_t) "moderated t from input"
                                else       "signed z from p-value", ")")
          else "")

  results <- lapply(proteins, function(prot) {
    compute_spatial(prot, long_data, pos_map, pval_thresh,
                    min_lfc     = min_lfc,
                    weight_mode = weight_mode)
  })
  results <- Filter(Negate(is.null), results)
  results <- Filter(function(x) x$n_regions >= min_regions, results)

  df <- do.call(rbind, lapply(results, function(x) {
    data.frame(protein   = x$protein,
               isTF      = x$isTF,
               centroid  = round(x$centroid, 1),
               spread    = round(x$spread, 1),
               composite = round(x$composite, 4),
               n_regions = x$n_regions,
               mean_lfc  = round(x$mean_lfc, 4),
               min_pval  = x$min_pval,
               sig_any   = x$sig_any,
               stringsAsFactors = FALSE)
  }))
  df[order(df$composite, decreasing = TRUE), ]
}

# =============================================================================
# SECTION 5: JASPAR motif fetch + PWM scan
# =============================================================================

# In-session PWM cache. Keeps (tf_name, host) -> PWM object across repeated
# fetches in the same R session. Critical for deterministic A/B comparisons
# (e.g., `1-try-compare.R`) where two `run_caspex()` calls would otherwise
# each hit the JASPAR API independently. For genes with multiple matrix
# profiles (HOX family is the canonical case), two independent API calls
# can return DIFFERENT matrices — different PWM → different score
# distribution → different hit counts. Caching pins the choice to whatever
# the first call resolved, so the two compare runs share identical motifs.
# `.caspex_pwm_cache` lives in the package namespace (defined in
# R/utils-internal.R, parent = emptyenv()). Within one loaded package
# session it persists across all run_caspex() calls, so multi-matrix TFs
# (HOX family etc.) resolve to a stable matrix pick on every call.

#' Fetch position weight matrix from JASPAR for a TF name
#'
#' JASPAR REST API uses `tf_name` (case-sensitive, exact) or `search`
#' (substring) as query parameters. `version` is not a valid filter — the
#' latest-version profile is selected by ordering via `-version`.
#' We try `tf_name` first, then fall back to `search` if empty.
#'
#' Results are cached in `.caspex_pwm_cache` keyed by (TF, host) so
#' repeated calls within one R session return the identical PWM. This
#' keeps A/B runs (default vs coverage-aware) reproducible.
#'
#' @param tf_name  HGNC symbol
#' @param host     JASPAR API host (default current canonical mirror)
#' @return List with id, name, pwm (4×L matrix, rows=ACGT), length; or NULL
#' @export
fetch_jaspar_pwm <- function(tf_name, host = "https://jaspar.elixir.no") {
  # Cache lookup — key by uppercased TF + host to avoid trivial mismatches.
  cache_key <- paste0(toupper(tf_name), "@", host)
  if (exists(cache_key, envir = .caspex_pwm_cache, inherits = FALSE))
    return(get(cache_key, envir = .caspex_pwm_cache, inherits = FALSE))
  api_search <- function(params) {
    # Paginate to cover the whole CORE vertebrates set if needed
    out <- list()
    page <- 1
    repeat {
      res <- GET(paste0(host, "/api/v1/matrix/"),
                 query = c(params, list(format    = "json",
                                         page_size = 100,
                                         page      = page,
                                         order     = "-version")),
                 timeout(20))
      if (status_code(res) != 200) break
      js <- content(res, "parsed", simplifyVector = FALSE)
      if (is.null(js$results) || length(js$results) == 0) break
      out <- c(out, js$results)
      if (is.null(js$`next`)) break
      page <- page + 1
      if (page > 30) break    # safety cap (3000 matrices)
    }
    out
  }

  # Accept matrix names that match the TF symbol, ignoring JASPAR
  # variant suffixes like "(var.2)" and case.
  clean_name <- function(n) toupper(sub("\\s*\\(.*$", "", n))
  target     <- toupper(tf_name)

  pick_match <- function(results) {
    if (length(results) == 0) return(NULL)
    hit_names <- vapply(results, function(r) clean_name(r$name), character(1))
    keep      <- which(hit_names == target)
    if (length(keep) == 0) return(NULL)     # <-- strict: no match -> NULL
    # Pick highest version among matching matrices; break ties deterministically
    # by matrix_id (lexicographic, higher wins). Without the tiebreak the pick
    # depends on the JASPAR API's result ordering, which is server-dependent
    # and can differ across processes — making motif hit counts drift across
    # fresh R sessions for TFs that have multiple matrices at the same max
    # version (e.g. HOX family: MA0898.1, MA1502.1, MA1557.1 all at v=1).
    versions <- vapply(results[keep],
                       function(r) suppressWarnings(as.numeric(r$version)),
                       numeric(1))
    versions[is.na(versions)] <- 0
    mat_ids  <- vapply(results[keep],
                       function(r) as.character(r$matrix_id), character(1))
    # order() is deterministic — sort by (-version, -matrix_id) and take first
    ord <- order(-versions, -xtfrm(mat_ids))
    results[[keep[ord[1]]]]
  }

  # 1) Exact tf_name match in CORE/vertebrates
  chosen <- pick_match(api_search(list(tf_name    = tf_name,
                                       collection = "CORE",
                                       tax_group  = "vertebrates")))
  # 2) Broaden: drop tax_group (some vertebrate TFs are tagged elsewhere)
  if (is.null(chosen))
    chosen <- pick_match(api_search(list(tf_name    = tf_name,
                                          collection = "CORE")))
  # 3) Fall back to substring search (name, synonyms, ID)
  if (is.null(chosen))
    chosen <- pick_match(api_search(list(search     = tf_name,
                                          collection = "CORE")))

  if (is.null(chosen)) {
    assign(cache_key, NULL, envir = .caspex_pwm_cache)
    return(NULL)
  }
  mat_id <- chosen$matrix_id

  # Step 2: fetch the actual PFM
  res2 <- GET(paste0(host, "/api/v1/matrix/", mat_id, "/"),
              query = list(format = "json"), timeout(15))
  if (status_code(res2) != 200) {
    # Don't cache transient HTTP failures — let the caller retry next time.
    return(NULL)
  }
  mat  <- content(res2, "parsed", simplifyVector = FALSE)
  pfm  <- mat$pfm  # list with keys A, C, G, T

  # Convert PFM → log-odds PWM  (pseudocount 0.25, bg = 0.25 each)
  bases <- c("A", "C", "G", "T")
  m     <- sapply(bases, function(b) unlist(pfm[[b]]))   # L × 4
  m     <- t(m)                                           # 4 × L
  m     <- m + 0.25
  m     <- sweep(m, 2, colSums(m), "/")                  # normalise cols
  pwm   <- log2(m / 0.25)                                 # log-odds vs uniform bg

  pwm_obj <- list(id = mat_id, name = tf_name, pwm = pwm, len = ncol(pwm))
  assign(cache_key, pwm_obj, envir = .caspex_pwm_cache)
  pwm_obj
}

#' Score a DNA sequence with a PWM
#'
#' Vectorized sliding-window scorer. Replaces the legacy O(ns·L) double
#' for-loop with an O(L) rowSums over a ns×L matrix of per-position PWM
#' scores. Non-ACGT bases (IUPAC ambiguity codes, Ns, soft-masked repeats)
#' contribute 0 at that column position, matching the legacy `if (!is.na(b))`
#' skip behaviour exactly. Runtime on 3 kb promoters × 12-mer PWMs drops
#' from ~100 ms per call to ~1 ms, which compounds across run_motif_scan
#' (two strands × dozens of TFs).
#' @param seq_chars (see function body).
#' @param pwm (see function body).
#' @noRd
score_pwm_positions <- function(seq_chars, pwm) {
  L  <- ncol(pwm)
  n  <- length(seq_chars)
  ns <- n - L + 1
  if (ns <= 0) return(numeric(0))
  base_idx <- c(A = 1L, C = 2L, G = 3L, T = 4L)

  # Numeric 1..4 index per base; NA for anything outside ACGT. Downstream
  # we replace NAs with a column-wise 0 contribution via a dummy row 5.
  idx <- base_idx[seq_chars]

  # Build a ns × L matrix where row i is positions [i, i+L-1] of `idx`.
  # `embed()` returns rows in reverse column order (j = L, L-1, ..., 1), so
  # we flip left-to-right to align column j with PWM column j.
  win <- embed(idx, L)[, L:1, drop = FALSE]

  # Augment PWM with a dummy 5th row of zeros so NA-index lookups become
  # "contribute 0 here" — same behaviour as the legacy `if (!is.na(b))`.
  pwm_aug <- rbind(pwm, rep(0, L))
  win[is.na(win)] <- 5L

  # For each window column j, `pwm_aug[win[, j], j]` is a length-ns vector of
  # contributions. Build the ns × L score matrix by column, then rowSums.
  score_mat <- vapply(seq_len(L),
                      function(j) pwm_aug[win[, j], j],
                      numeric(ns))
  if (!is.matrix(score_mat)) score_mat <- matrix(score_mat, ncol = L)
  rowSums(score_mat)
}

#' Scan both strands of a promoter sequence for PWM hits
#'
#' @param promoter_info  Output of fetch_promoter_seq()
#' @param pwm_obj        Output of fetch_jaspar_pwm()
#' @param threshold_frac Fraction of max possible score to use as cutoff
#' @return Numeric vector of TSS-relative hit positions, with attributes
#'         `scores` (parallel raw log-odds score per hit), `max_score`,
#'         `min_score`, `score_frac` (per-hit fraction of max log-odds, in
#'         [threshold_frac, 1]). Existing callers that only consume the
#'         positions vector are unaffected; downstream code that wants
#'         score-weighted NNLS reads `attr(hits, "score_frac")`.
#' @noRd
scan_sequence <- function(promoter_info, pwm_obj, threshold_frac = 0.80) {
  if (is.null(pwm_obj)) return(integer(0))
  pwm    <- pwm_obj$pwm
  L      <- pwm_obj$len
  seq    <- promoter_info$seq
  tss_i  <- promoter_info$tss_offset + 1

  max_score <- sum(apply(pwm, 2, max))
  min_score <- sum(apply(pwm, 2, min))
  threshold <- min_score + threshold_frac * (max_score - min_score)

  fwd_chars <- strsplit(seq, "")[[1]]
  rev_map   <- c(A="T", T="A", G="C", C="G", N="N")
  rev_chars <- rev(rev_map[fwd_chars])

  fwd_scores <- score_pwm_positions(fwd_chars, pwm)
  rev_scores <- score_pwm_positions(rev_chars, pwm)

  fwd_idx <- which(fwd_scores >= threshold)
  rev_idx <- which(rev_scores >= threshold)
  # Reverse strand: position in fwd coords = nchar(seq) - (rev_pos + L - 2)
  n        <- nchar(seq)
  fwd_hits <- fwd_idx
  rev_hits <- n - (rev_idx + L - 2)

  # Score per hit (forward scores at fwd_idx; reverse scores at rev_idx)
  fwd_score_vals <- fwd_scores[fwd_idx]
  rev_score_vals <- rev_scores[rev_idx]

  raw_pos    <- c(fwd_hits, rev_hits)
  raw_scores <- c(fwd_score_vals, rev_score_vals)
  ord        <- order(raw_pos)
  raw_pos    <- raw_pos[ord]
  raw_scores <- raw_scores[ord]
  # Deduplicate identical positions (forward and reverse can coincide for
  # palindromic motifs); keep the higher-scoring orientation.
  if (anyDuplicated(raw_pos)) {
    keep_idx <- !duplicated(raw_pos) |
      (duplicated(raw_pos) &
         raw_scores > c(-Inf, raw_scores[-length(raw_scores)]))
    # Simpler: per-position, keep max score
    agg <- tapply(raw_scores, raw_pos, max)
    raw_pos    <- as.integer(names(agg))
    raw_scores <- as.numeric(agg)
  }

  positions  <- raw_pos - tss_i  # TSS-relative
  score_frac <- (raw_scores - min_score) / (max_score - min_score)

  attr(positions, "scores")     <- raw_scores
  attr(positions, "max_score")  <- max_score
  attr(positions, "min_score")  <- min_score
  attr(positions, "score_frac") <- score_frac
  positions
}

#' Fetch PWMs and scan sequence for a vector of TF names
#'
#' @param tf_names       Character vector of TF symbols
#' @param promoter_info  Output of fetch_promoter_seq()
#' @param threshold_frac Score threshold
#' @return Named list: each element has $hits (TSS-relative positions),
#'         $score_frac (per-hit fraction of max log-odds, parallel to
#'         $hits), $pwm, and $n_hits.
#' @noRd
run_motif_scan <- function(tf_names, promoter_info, threshold_frac = 0.80) {
  message("\nRunning JASPAR motif scan for ", length(tf_names), " TFs...")
  out <- list()
  for (tf in tf_names) {
    message("  ", tf, " ...", appendLF = FALSE)
    pwm <- tryCatch(fetch_jaspar_pwm(tf),
                    error = function(e) { message(" error: ", conditionMessage(e)); NULL })
    if (is.null(pwm)) { message(" no motif found"); next }
    hits <- scan_sequence(promoter_info, pwm, threshold_frac)
    score_frac <- attr(hits, "score_frac") %||% rep(NA_real_, length(hits))
    # Strip attributes so `hits` consumed elsewhere stays a plain integer
    # vector (subset-safe). The scores travel alongside in $score_frac.
    plain_hits <- as.integer(hits)
    message(" [", pwm$id, "] ", length(plain_hits), " hit(s)")
    out[[tf]] <- list(pwm        = pwm,
                      hits       = plain_hits,
                      score_frac = score_frac,
                      n_hits     = length(plain_hits))
  }
  out
}

# =============================================================================
# SECTION 5.5: Motif-constrained binding deconvolution
#
# The per-TF centroid from SECTION 4 is a single weighted midpoint. That
# collapses two distinct binding events (e.g., R1=3, R2=0, R3=2) into one
# phantom peak at R2. This section builds a continuous CasPEX signal along
# the promoter and decomposes it onto JASPAR motif positions via NNLS so
# multiple bindings — and bindings that fall *between* two adjacent
# enriched regions — can be separated and localised.
# =============================================================================

#' Continuous CasPEX signal for a TF along x_grid.
#'
#' Each gRNA cut site contributes a Gaussian labeling kernel weighted by the
#' region's logFC (positive only, by default). The combined curve is what a
#' proximity-labeling experiment would produce if binding happened everywhere
#' proportional to the contribution map.
#'
#' @param tf_name       Protein symbol
#' @param long_data     From load_region_data()
#' @param pos_map       From match_all_grnas()
#' @param x_grid        bp grid (TSS-relative)
#' @param kernel_sigma  Labeling kernel σ in bp (default 300; APEX-ish radius)
#' @param weight_mode   "z" (default — signed z from p, ATP7B-style), "mod_t",
#'                      "lfc_pos", "lfc_signed", "lfc_x_negp" (= "weighted",
#'                      legacy).
#'                      Only the positive part is retained for signal building
#'                      — negative weights would subtract, which doesn't model
#'                      a proximity-labeling intensity map.
#' @return list(x, y, region_data, kernel_sigma)
#' @noRd
build_caspex_signal <- function(tf_name, long_data, pos_map,
                                x_grid       = seq(-2500, 500, by = 5),
                                kernel_sigma = 300,
                                weight_mode  = c("z", "mod_t", "lfc_pos",
                                                  "lfc_signed", "lfc_x_negp",
                                                  "weighted")) {
  weight_mode <- match.arg(weight_mode)
  # `weighted` is a legacy alias for `lfc_x_negp`
  if (weight_mode == "weighted") weight_mode <- "lfc_x_negp"

  # --- All regions from pos_map (with valid positions) — used for lollipops.
  # Every region the guide layout covers should appear on the plot, even if
  # the protein was undetected in that region or reported lfc = NA. Undetected
  # / NA entries are plotted at lfc = 0 (dot on baseline, no bar) so the
  # viewer can see "measured but no enrichment" vs. a silent drop.
  valid_regions <- names(pos_map)[!is.na(pos_map)]
  all_df <- data.frame(region  = valid_regions,
                       pos     = as.numeric(pos_map[valid_regions]),
                       lfc     = 0,
                       protein = tf_name,
                       stringsAsFactors = FALSE)

  # Detected rows for this protein (raw, pre-NA filter)
  detected <- long_data[long_data$protein == tf_name &
                          long_data$region %in% valid_regions, ]
  if (nrow(detected) > 0) {
    # If the same (protein, region) appears multiple times, keep the first
    # non-NA lfc; otherwise fall back to the first row.
    keep_idx <- vapply(valid_regions, function(r) {
      hits <- which(detected$region == r)
      if (length(hits) == 0) return(NA_integer_)
      nn <- hits[!is.na(detected$lfc[hits])]
      if (length(nn) > 0) nn[1] else hits[1]
    }, integer(1))
    has_row <- !is.na(keep_idx)
    if (any(has_row)) {
      lfc_vec <- detected$lfc[keep_idx[has_row]]
      lfc_vec[is.na(lfc_vec)] <- 0
      all_df$lfc[match(valid_regions[has_row], all_df$region)] <- lfc_vec
    }
  }

  # --- Signal construction uses only detected rows with non-NA lfc/pos.
  df_sig <- detected
  if (nrow(df_sig) > 0) {
    df_sig$pos <- as.numeric(pos_map[df_sig$region])
    df_sig <- df_sig[!is.na(df_sig$pos) & !is.na(df_sig$lfc), ]
  }

  y <- numeric(length(x_grid))
  if (nrow(df_sig) > 0) {
    w_raw <- compute_region_weight(df_sig, mode = weight_mode)
    # Signal is an intensity map; clip to non-negative.
    w     <- pmax(w_raw, 0)
    df_sig$w <- w
    for (i in seq_len(nrow(df_sig))) {
      y <- y + w[i] * exp(-0.5 * ((x_grid - df_sig$pos[i]) / kernel_sigma)^2)
    }
  }

  list(x = x_grid, y = y, region_data = all_df, kernel_sigma = kernel_sigma)
}

#' CasPEX labeling-opportunity (coverage) map.
#'
#' Returns C(x) = Σ_r G(x − pos_r; σ) — the UNWEIGHTED Gaussian coverage
#' produced by the gRNA cut-site layout alone. Physically: the cumulative
#' biotin-labeling "opportunity" at each genomic position, regardless of
#' what proteins happen to be bound there. Depleted-weight gRNAs still
#' contribute to C(x) because labeling is a physical property of the
#' cut site, not a property of proteomic enrichment.
#'
#' Used by `predict_binding_events_coverage_aware()` to invert the
#' labeling decay: if s(x) is the enrichment-intensity map and C(x) is
#' the labeling-opportunity map, then
#'
#'     β(x) = s(x) / max(C(x), cov_floor · max(C))
#'
#' is the "per-unit-labeling enrichment" — an occupancy estimate that
#' corrects for the geometric attenuation of biotinylation with distance.
#' In the overlap region between two gRNAs, both s and C are boosted, so
#' β recovers the underlying weight; in gRNA-poor tails, both go to zero
#' and the `cov_floor` term caps the amplification so noise can't be
#' inflated indefinitely.
#'
#' Also useful as a standalone diagnostic overlay on the deconvolution
#' plot, to show the viewer where labeling is strong vs. where a call
#' was coverage-rescued.
#'
#' @param pos_map       Named bp positions (from match_all_grnas)
#' @param x_grid        bp grid (TSS-relative)
#' @param kernel_sigma  Labeling kernel σ (default 300)
#' @return list(x, y, pos, kernel_sigma) — `y` is C(x)
#' @noRd
compute_coverage <- function(pos_map,
                             x_grid       = seq(-2500, 500, by = 5),
                             kernel_sigma = 300) {
  pos <- as.numeric(pos_map[!is.na(pos_map)])
  y   <- numeric(length(x_grid))
  for (p in pos) y <- y + exp(-0.5 * ((x_grid - p) / kernel_sigma)^2)
  list(x = x_grid, y = y, pos = pos, kernel_sigma = kernel_sigma)
}

#' Simple local-maxima peak detection with minimum separation.
#' Greedy selection in descending height; preserves valleys.
#' @param x (see function body).
#' @param y (see function body).
#' @param min_height (see function body).
#' @param min_dist (see function body).
#' @noRd
find_local_maxima <- function(x, y, min_height = 0, min_dist = 150) {
  n <- length(y)
  if (n < 3) return(integer(0))
  peaks <- which(y[2:(n - 1)] > y[1:(n - 2)] &
                   y[2:(n - 1)] > y[3:n]) + 1
  peaks <- peaks[y[peaks] >= min_height]
  if (length(peaks) <= 1) return(peaks)
  ord <- peaks[order(y[peaks], decreasing = TRUE)]
  kept <- integer(0)
  for (p in ord) {
    if (length(kept) == 0 || all(abs(x[p] - x[kept]) >= min_dist))
      kept <- c(kept, p)
  }
  sort(kept)
}

#' Predict binding events with COVERAGE-NORMALIZED per-motif scoring.
#'
#' Motivation ---------------------------------------------------------------
#' CasPEX biotinylation intensity falls off with distance from each gRNA
#' cut site (APEX-like Gaussian decay, σ ~ 300 bp). So the per-region
#' enrichment `w_r` we observe is the product of TRUE OCCUPANCY at nearby
#' loci × the LABELING EFFICIENCY of the nearest gRNA(s) at that distance.
#' The smoothed-signal path (`predict_binding_events`) inherits that decay:
#' a real binder sitting halfway between two gRNAs, or in a gRNA-poor gap,
#' generates a `s(x)` that is geometrically smaller than one sitting under
#' a gRNA, even if both are equally occupied. NNLS on s(x) then discards
#' the midpoint / gap events because they look like weak shoulders.
#'
#' Biologically, that's wrong. We want to invert the labeling decay, not
#' absorb it. The key quantity is
#'
#'     β(x) = s(x) / max(C(x), floor)
#'
#' where
#'     s(x) = Σ_r w_r⁺ · G(x − pos_r; σ)         (enrichment intensity)
#'     C(x) = Σ_r       G(x − pos_r; σ)           (labeling opportunity)
#'
#' In words: "per-unit-labeling enrichment at x". At x directly under a
#' gRNA, C ≈ 1 (with small tail contributions from others), so β ≈ w_r.
#' At a midpoint between two equally-enriched gRNAs, both s and C roughly
#' double in the overlap region, and β recovers the individual region's
#' weight — i.e. a motif there is called with the same strength as one
#' under a gRNA. In a gRNA-poor tail both s and C approach zero; the
#' `floor = cov_floor · max(C)` term caps the amplification so noise
#' can't be inflated indefinitely.
#'
#' Worked example
#' --------------
#' R1 at −500, R2 at −1000, σ = 300 bp, w_1 = w_2 = 1.
#' Motifs at −450, −700 (midpoint), −1050:
#'   m = −450 : s = 0.98·1 + 0.089·1 = 1.07, C = 1.07, β = 1.00
#'   m = −700 : s = 0.61·1 + 0.61·1 = 1.22, C = 1.22, β = 1.00
#'   m = −1050: s = 0.089·1 + 0.98·1 = 1.07, C = 1.07, β = 1.00
#' All three are called with equal weight — which is what you asked for.
#'
#' Contrast: the default smoothed-s(x) NNLS will typically put weight on
#' −450 and −1050 and clip the −700 motif, because a 2-gRNA signal can
#' be reconstructed by 2 basis columns and NNLS picks the sparse fit.
#'
#' Guardrails --------------------------------------------------------------
#'   * `cov_floor` (default 0.05) is a *relative* floor: the denominator
#'     is `max(C(x), cov_floor · max(C))`. Effectively limits inverse-
#'     coverage amplification to ~1/cov_floor × the max-coverage value.
#'   * `min_weight_frac` defines the above-threshold zones on β(x):
#'     zones are contiguous intervals where β(x) > frac · max(β(x)).
#'   * `merge_dist` collapses nearby survivors (amplitude-weighted centroid).
#'   * Motif-less TFs fall back to peak detection on the RAW signal s(x)
#'     (same as the default path), with the in-support mask applied so any
#'     signal outside the guide layout is zeroed before peak-picking. The
#'     earlier β(x)-based fallback tended to mark window edges as peaks
#'     because pmax(C, floor_val) clamps near the tiled-region boundaries
#'     and amplifies residual s there.
#'
#' Zone-based detection (motif branch)
#' -----------------------------------
#' A pure point-β-at-motif rule (score each motif's β, keep the ones above
#' `frac · max(β_at_motifs)`) is brittle on slopes: a motif sitting on
#' the shoulder of a strong peak can fall below threshold even when β
#' there is clearly elevated above the floor. Instead we compute β(x) on
#' the GRID, find contiguous above-threshold zones, and emit EVERY motif
#' that falls inside any such zone. A zone with no motif still emits one
#' bubble at its peak so real signal without a JASPAR hit is not lost.
#' This matches the detection-style question "is this motif inside a
#' plausible binding zone?" rather than "does this motif sit at a β peak?".
#'
#' Interpretation
#' --------------
#' Unlike the NNLS paths (smoothed or per-region), each motif here is
#' evaluated INDEPENDENTLY — there is no sparsity pressure. Every motif
#' inside a plausible binding zone is surfaced, subject only to the zone
#' threshold (`min_weight_frac`) and the coverage floor (`cov_floor`).
#'
#' @return data.frame(tf, position, weight, motif_based, n_motifs_merged,
#'                    distance_to_nearest_grna, local_coverage)

# Note: an earlier compute_mle_positions() function lived here, providing
# a generative-model MLE position track for Plot 10. It was removed
# because the diagnostic added visual noise without changing the
# bubble-placement story in any actionable way — the bubble at β-peak
# already conveys the relevant per-event information, and the MLE band
# duplicated information that's already available from the s(x) curve
# and the per-region logFC stems. Removed: the function definition, the
# `mle_position` parameter on every call site
# (predict_binding_events_coverage_aware → predict_all_binding_events →
# plot_binding_deconvolution → run_caspex), the has_mle layout block in
# Plot 10, the schema columns added by add_diag, and the result-list
# bookkeeping. The Wild bootstrap "Position support" infrastructure
# remains in place under `position_stability` for back-compatibility but
# is not rendered.

#' Predict binding events with coverage-aware deconvolution
#'
#' Coverage-normalised counterpart to predict_binding_events(): operates
#' on \eqn{\beta(x) = s(x)/C(x)} inside a zone-based detector that emits
#' one bubble per qualifying motif inside a contiguous above-threshold
#' region of \eqn{\beta}. See run_caspex() for the full parameter wiring.
#'
#' @param tf_name (see function body).
#' @param long_data (see function body).
#' @param pos_map (see function body).
#' @param motif_hits (see function body).
#' @param x_grid (see function body).
#' @param kernel_sigma (see function body).
#' @param min_weight_frac (see function body).
#' @param min_peak_dist (see function body).
#' @param merge_dist (see function body).
#' @param weight_mode (see function body).
#' @param cov_floor (see function body).
#' @param motif_scores (see function body).
#' @param motif_score_weight (see function body).
#' @param edge_guard_frac (see function body).
#' @param zone_peak_frac (see function body).
#' @param max_events_per_tf (see function body).
#' @param merge_position (see function body).
#' @param max_grna_distance (see function body).
#' @param edge_grna_weight_cap (see function body).
#' @param position_stability (see function body).
#' @param n_bootstrap (see function body).
#' @noRd
predict_binding_events_coverage_aware <- function(
    tf_name, long_data, pos_map, motif_hits,
    # Optional per-motif PWM score-fraction in [0,1] (parallel to motif_hits).
    # Used only when motif_score_weight != "none".
    motif_scores    = NULL,
    # `motif_score_weight`: how to incorporate PWM score into the per-motif
    # weight emitted from each above-threshold zone.
    #   "none"   - default, legacy behaviour: weight = β(motif_pos).
    #   "linear" - weight = β(motif_pos) × score_frac.
    #   "log"    - weight = β(motif_pos) × 2^(score_frac − 1)
    #              (compresses the difference; weak motifs get 0.5×, strong
    #               motifs get 1.0×).
    motif_score_weight = c("none", "linear", "log"),
    x_grid          = seq(-2500, 500, by = 5),
    kernel_sigma    = 300,
    min_weight_frac = 0.15,
    min_peak_dist   = 150,
    merge_dist      = 100,
    weight_mode     = "z",
    cov_floor       = 0.05,
    # `edge_guard_frac`: fraction-of-max-coverage threshold that defines the
    # in-support region. β is trusted only where c_grid > edge_guard_frac ×
    # max(c_grid); elsewhere it is zeroed before zone detection and peak
    # picking. MUST be ≥ cov_floor — the clamp only prevents β from going
    # to infinity, while this guard prevents β from being evaluated in the
    # low-denominator transition band just inside the clamp (where a real
    # peak-finder / zone-finder can still see the ramp as a spurious peak
    # or zone). Default 0.25 (= 5× default cov_floor); raised from 0.15
    # because σ=300 extends the single-gRNA tail far enough that 0.15
    # leaked edge events ~1.4σ west of the westernmost guide. At 0.25 the
    # trust boundary sits ~1σ outside an isolated guide, matching the
    # geometric assumption that β is uninformative in the single-gRNA
    # tail (numerator and denominator are the same Gaussian and cancel).
    # Raise further (0.30) for sparser gRNA layouts; lower (0.15–0.20)
    # for densely tiled regions where multi-guide overlap lifts max(C).
    edge_guard_frac = 0.25,
    # ---- readability guards (prevent HOXB6-style plateau floods) ----
    # `zone_peak_frac`: within each above-threshold zone, only keep motifs
    # whose local β is >= zone_peak_frac × max(β) inside that zone. This
    # filters shoulders of broad plateaus while leaving motifs on genuinely
    # elevated slopes (e.g. GATA6 R7-R6 ramp, β > 0.8 × zone peak) alone.
    # Set to 0 to disable and emit every zone motif as before.
    zone_peak_frac  = 0.50,
    # `max_events_per_tf`: hard cap on bubbles per TF, top-N by weight.
    # Prevents plots from becoming unreadable on promiscuous PWMs with
    # broadly elevated signal (HOXB6 had 273 candidates pre-cap). Set to
    # Inf to disable. Applied AFTER merging so merged centroids are ranked
    # on their combined amplitude.
    max_events_per_tf = 30,
    # `merge_position`: how to report the position for a merged cluster of
    # motif candidates. "argmax" (default, new behaviour) snaps the bubble
    # to the strongest motif in the cluster — so every motif-based bubble
    # sits exactly on a real motif tick, matching what the no-motif
    # fallback does at an s(x) peak. "centroid" keeps the legacy
    # amplitude-weighted mean position, which can land between motifs when
    # 2-3 of them chain within merge_dist. Motif-less-zone bubbles (single-
    # candidate clusters emitted at a zone's β peak) are unaffected — they
    # carry their peak position regardless of this setting.
    merge_position = c("argmax", "centroid"),
    # `max_grna_distance`: hard geometric cap in bp on how far a called
    # event can sit from the nearest gRNA. Events whose
    # `distance_to_nearest_grna` exceeds this value are dropped AFTER the
    # merge step. Motivation: in the single-gRNA tail (positions west of
    # the westernmost guide or east of the easternmost), β = s/C is
    # mathematically flat because numerator and denominator share the
    # same Gaussian shape and cancel — so the zone detector can emit a
    # swarm of bubbles all carrying the same β, even at positions where
    # the labeling model provides no geometric information to distinguish
    # them. This cap is the belt-and-suspenders companion to
    # `edge_guard_frac`: it scales automatically with σ, which the
    # relative-coverage mask does not. Default NULL resolves to
    # `kernel_sigma` at runtime (so "events within one σ of a guide").
    # Set `Inf` to disable. Set to a smaller multiple of σ (e.g. 0.75·σ)
    # to be stricter about tail leakage.
    max_grna_distance = NULL,
    # `edge_grna_weight_cap`: fraction-of-Gaussian-weight-sum above which an
    # outermost gRNA (westernmost OR easternmost) is considered to be
    # "dominating" an event's local signal. Events where either boundary
    # guide contributes more than this fraction of the summed per-gRNA
    # Gaussian weight at the event position are dropped.
    #
    # Motivates this: `max_grna_distance` gates events by distance to the
    # NEAREST guide, which protects against tail leakage past the outermost
    # gRNA, but does nothing for interior events that still get much of
    # their β from the westernmost / easternmost guide's kernel tail. In
    # ATP7B at σ=300, R7(-1945) contributes ~37% of the summed Gaussian
    # weight at position -1620 even though R6(-1507) is only 113 bp away —
    # producing a huge "edge-inflated" bubble at a position that is fully
    # inside the trusted support mask.
    #
    # 0.30 is a reasonable starting point — it keeps events that are
    # solidly backed by interior guides while rejecting those where a
    # single boundary guide accounts for more than a third of the local
    # per-gRNA weight sum.  NULL (default) disables the filter entirely,
    # preserving legacy behaviour.
    edge_grna_weight_cap = NULL,
    # `position_stability` (coverage-aware only): per-event diagnostic that
    # measures how robust a called event's position is to noise perturbation
    # of the per-region observed weights.
    #
    # NOT a confidence interval on the TF's true binding location. A high-β
    # event can still be a distant binder (high enrichment despite distance);
    # this metric captures the position-spread the detector itself would
    # produce if the same dataset had landed slightly differently under
    # heteroskedastic noise. A wide spread means the regional support pattern
    # would have placed the bubble somewhere else under modest perturbation;
    # a narrow spread means the position is determined by the data, not by
    # noise.
    #
    #   "none"           — default; no bootstrap; output unchanged.
    #   "wild_bootstrap" — Rademacher Wild bootstrap on the residuals of a
    #                      forward-model NNLS fit at the called event
    #                      positions. Preserves the gRNA layout and design
    #                      matrix (which is correct for our small N=5–10
    #                      regions) while resampling the per-region weight
    #                      noise. For each bootstrap draw, perturbed weights
    #                      are re-injected into long_data via the lfc field
    #                      (with t_stat = NA, weight_mode = "lfc_signed") and
    #                      the inner pipeline (zone detection + clustering +
    #                      filters) is re-run. Each original event is then
    #                      matched to its nearest bootstrap event within 2σ,
    #                      and the matched-position distribution is
    #                      summarised as pos_stab_low/high (2.5/97.5%
    #                      quantiles), pos_stab_median, and
    #                      bootstrap_survival_frac (= matched draws /
    #                      n_bootstrap).
    position_stability = c("none", "wild_bootstrap"),
    # Number of Wild bootstrap draws when `position_stability` is enabled.
    # 200 is a sensible default for 95% CI-style summaries; bump to 500–1000
    # for tighter tails or stable-binders for a paper figure.
    n_bootstrap        = 200L
) {
  merge_position     <- match.arg(merge_position)
  position_stability <- match.arg(position_stability)
  motif_score_weight <- match.arg(motif_score_weight)
  # Resolve NULL → kernel_sigma. Doing this in the body (not as a default
  # expression) keeps the dependence explicit and survives callers that
  # pass `max_grna_distance = NULL` intentionally.
  if (is.null(max_grna_distance)) max_grna_distance <- kernel_sigma
  empty <- data.frame(
    tf = character(), position = numeric(), weight = numeric(),
    motif_based = logical(), n_motifs_merged = integer(),
    distance_to_nearest_grna = numeric(),
    local_coverage = numeric(),
    stringsAsFactors = FALSE
  )

  # ── Wild bootstrap position stability (Myers benchmarking) ──────────────
  # Defined as an inner closure so BOTH return paths in this function — the
  # motif-less fallback and the motif-based zone branch — get the same column
  # schema for `out`. Without this, no-motif TFs returned `ev` early without
  # `pos_stab_*`, while motif-anchored TFs got the columns, and rbind across
  # TFs in predict_all_binding_events failed with "numbers of columns of
  # arguments do not match".
  #
  # See `position_stability` param docs above for the framing — these
  # columns describe how robust each bubble's POSITION is to noise
  # perturbation, not a confidence interval on the TF's true binding
  # location. The procedure preserves the gRNA layout (design matrix) and
  # only resamples the per-region weight residual via Rademacher signs:
  #
  #   1. Fit a_k ≥ 0 by NNLS so that  Σ_k a_k · G(pos_r - p_k; σ)  ≈  w_obs[r],
  #      where w_obs[r] is the per-region weight that build_caspex_signal()
  #      used (i.e. compute_region_weight() with the active weight_mode,
  #      clipped at 0). p_k are the called event positions in `out`.
  #   2. Compute residuals  ε[r] = w_obs[r] − ŵ[r].
  #   3. For b = 1…n_bootstrap, draw signs  s_b[r] ∈ {−1, +1}  i.i.d., set
  #      perturbed weights  w_b[r] = max(ŵ[r] + s_b[r] · ε[r], 0),  inject
  #      them into a long_data clone via the lfc column with t_stat = NA,
  #      and recursively call this same function with weight_mode =
  #      "lfc_signed" (so build_caspex_signal recovers exactly w_b) and
  #      position_stability = "none" (recursion guard).
  #   4. For each original event, record the position of the nearest
  #      bootstrap event within 2σ (no match → no contribution).
  #   5. Append pos_stab_{low,high,median} (2.5 / 97.5 / 50% quantiles of
  #      matched positions) and bootstrap_survival_frac (= matched draws /
  #      n_bootstrap).
  do_pos_stability <- function(out) {
    # When stability is off / disabled / nothing to bootstrap, return as-is
    # (no extra columns). predict_all_binding_events filters nrow(ev) > 0
    # before rbind, so empty results never reach the column-mismatch path.
    if (position_stability != "wild_bootstrap" || !isTRUE(n_bootstrap > 0) ||
        nrow(out) == 0) {
      return(out)
    }
    # Initialise columns up front so callers see a consistent schema even
    # when the bootstrap can't actually run (e.g. fewer than 2 detected
    # regions for this TF, or NNLS fails). This is the rbind-across-TFs
    # guarantee.
    out$pos_stab_low            <- NA_real_
    out$pos_stab_high           <- NA_real_
    out$pos_stab_median         <- NA_real_
    out$bootstrap_survival_frac <- NA_real_
    detected <- long_data[long_data$protein == tf_name, , drop = FALSE]
    detected$pos <- as.numeric(pos_map[detected$region])
    detected <- detected[!is.na(detected$pos) & !is.na(detected$lfc), ,
                         drop = FALSE]
    if (nrow(detected) < 2) return(out)
    w_obs   <- pmax(compute_region_weight(detected, mode = weight_mode), 0)
    pos_obs <- detected$pos
    # Forward design matrix at called positions: rows = regions, cols = events
    A <- exp(-0.5 * outer(pos_obs, out$position,
                          FUN = function(p, q) (p - q) / kernel_sigma)^2)
    fit_a <- tryCatch(nnls::nnls(A, w_obs)$x,
                      error = function(e) NULL)
    if (is.null(fit_a) || length(fit_a) != nrow(out)) return(out)
    w_hat <- as.numeric(A %*% fit_a)
    eps   <- w_obs - w_hat
    n_r   <- length(w_obs)
    n_ev  <- nrow(out)
    boot_pos <- vector("list", n_ev)
    for (k in seq_len(n_ev)) boot_pos[[k]] <- numeric(0)
    n_match <- integer(n_ev)
    match_tol <- 2 * kernel_sigma
    for (b in seq_len(n_bootstrap)) {
      sgn <- sample(c(-1, 1), n_r, replace = TRUE)
      w_b <- pmax(w_hat + sgn * eps, 0)
      ld_b <- detected
      ld_b$lfc <- w_b
      if ("t_stat"  %in% names(ld_b)) ld_b$t_stat  <- NA_real_
      if ("pval"    %in% names(ld_b)) ld_b$pval    <- NA_real_
      if ("p_value" %in% names(ld_b)) ld_b$p_value <- NA_real_
      # Stitch back to the long_data schema so the recursive call's subset
      # on `protein == tf_name` finds these rows.
      ld_b$protein <- tf_name
      # Recursive call with the same params except (a) weight_mode forced
      # to "lfc_signed" so build_caspex_signal reproduces w_b exactly, and
      # (b) position_stability = "none" so we don't bootstrap-of-bootstraps.
      ev_b <- tryCatch(
        predict_binding_events_coverage_aware(
          tf_name, ld_b, pos_map, motif_hits,
          motif_scores       = motif_scores,
          motif_score_weight = motif_score_weight,
          x_grid          = x_grid,
          kernel_sigma    = kernel_sigma,
          min_weight_frac = min_weight_frac,
          min_peak_dist   = min_peak_dist,
          merge_dist      = merge_dist,
          weight_mode     = "lfc_signed",
          cov_floor       = cov_floor,
          edge_guard_frac = edge_guard_frac,
          zone_peak_frac  = zone_peak_frac,
          max_events_per_tf = max_events_per_tf,
          merge_position  = merge_position,
          max_grna_distance = max_grna_distance,
          edge_grna_weight_cap = edge_grna_weight_cap,
          position_stability = "none",
          n_bootstrap     = 0L
        ),
        error = function(e) NULL)
      if (is.null(ev_b) || nrow(ev_b) == 0) next
      # Match each original event to its nearest bootstrap event within 2σ.
      # Allow many-to-one matching (we only need a position estimate per
      # original bubble; double-counting a shared bootstrap event is fine
      # because it's a per-event diagnostic).
      for (k in seq_len(n_ev)) {
        d  <- abs(ev_b$position - out$position[k])
        j  <- which.min(d)
        if (length(j) && d[j] <= match_tol) {
          boot_pos[[k]] <- c(boot_pos[[k]], ev_b$position[j])
          n_match[k]    <- n_match[k] + 1L
        }
      }
    }
    out$pos_stab_low    <- vapply(boot_pos, function(v)
      if (length(v) >= 3) unname(quantile(v, 0.025)) else NA_real_,
      numeric(1))
    out$pos_stab_high   <- vapply(boot_pos, function(v)
      if (length(v) >= 3) unname(quantile(v, 0.975)) else NA_real_,
      numeric(1))
    out$pos_stab_median <- vapply(boot_pos, function(v)
      if (length(v) >= 1) median(v) else NA_real_,
      numeric(1))
    out$bootstrap_survival_frac <- n_match / n_bootstrap
    out
  }

  # Diagnostic augmentation wrapper. Both return branches (motif-based
  # and motif-less fallback) route through this so the column schema is
  # consistent for the predict_all_binding_events rbind. Currently only
  # the Wild-bootstrap "Position support" diagnostic lives here (see
  # `do_pos_stability`); kept as a wrapper so future per-event
  # diagnostics can be threaded in without re-touching both return paths.
  add_diag <- function(events) {
    do_pos_stability(events)
  }

  # Enrichment-intensity map s(x). Negative weights clipped to zero (same
  # contract as the default path); we are modelling labeling INTENSITY.
  sig   <- build_caspex_signal(tf_name, long_data, pos_map, x_grid,
                                kernel_sigma, weight_mode)
  s_grid <- sig$y
  if (max(s_grid) <= 0) return(empty)

  # Labeling-opportunity map C(x). UNWEIGHTED: every gRNA contributes
  # regardless of its proteomic enrichment direction, because labeling is
  # a physical property of the gRNA cut site — depleted gRNAs still emit
  # biotin, they just don't capture proteins preferentially.
  cov_obj <- compute_coverage(pos_map, x_grid, kernel_sigma)
  c_grid  <- cov_obj$y
  if (max(c_grid) <= 0) return(empty)

  # gRNA positions (for diagnostics: distance_to_nearest_grna)
  pos_r <- sort(as.numeric(pos_map[!is.na(pos_map)]))

  # Relative floor on the denominator — caps amplification. Setting
  # floor_val as a fraction of max(C) makes `cov_floor` dimensionless and
  # independent of how many gRNAs you tiled.
  floor_val <- cov_floor * max(c_grid)

  # Filter motifs to the requested window so downstream behaviour matches
  # the default path. Keep motif_scores parallel.
  in_range <- !is.na(motif_hits) &
              motif_hits >= min(x_grid) &
              motif_hits <= max(x_grid)
  if (!is.null(motif_scores) && length(motif_scores) == length(in_range)) {
    motif_scores <- motif_scores[in_range]
  } else {
    motif_scores <- NULL  # length mismatch -> ignore
  }
  motif_hits <- motif_hits[in_range]
  # Build a fast lookup from position -> score_frac (used inside the
  # zone loop below to weight emitted bubbles when motif_score_weight
  # != "none").
  score_lookup <- if (!is.null(motif_scores))
    setNames(pmax(0, pmin(1, motif_scores)), as.character(motif_hits))
  else NULL

  # In-support mask: β is trusted only where RAW C(x) exceeds a multiple of
  # the clamp floor. The clamp prevents β from going to infinity, but it
  # does NOT prevent β from being evaluated in the low-denominator band
  # just inside the clamp — where a tiny residual s(x) divided by a tiny
  # c_grid(x) produces an inflated β that peak-finders and zone-finders
  # happily mark as a hit at the west edge. Setting the support floor to
  # `edge_guard_frac × max(C)` (default 0.25, i.e. 5× the default 0.05
  # cov_floor) pushes the trust boundary well off the clamp transition so
  # the ramp lives in the masked-out zone and can't form a zone on its own.
  # `edge_guard_frac` is capped at no less than cov_floor (the clamp floor)
  # so reducing one without the other never re-opens the old failure mode.
  support_floor_val <- max(cov_floor, edge_guard_frac) * max(c_grid)
  support_mask      <- c_grid > support_floor_val

  # ── Motif-less fallback ────────────────────────────────────────────────
  # Peak detection on the RAW signal s(x) — identical in spirit to the
  # default path. Previously this branch peak-picked on β(x) = s/max(C,floor),
  # which tends to spike near the window edges: C(x) trails off fast at the
  # tiled-region boundaries, pmax(C, floor_val) clamps to floor_val, and any
  # tiny residual s there gets divided by a tiny denominator and amplified
  # into a "peak" at the grid edge. The support_mask helped but did not
  # fully eliminate these artefacts — a real peak-finder still sees the
  # upward slope of β ramping into the clamped zone. Peak-picking on s(x)
  # (which is already smoothed by the Gaussian kernel) gives the same
  # robust no-motif behaviour as the default path; we still zero out s
  # outside the in-support region so anything beyond the guide layout is
  # off the table by construction.
  if (length(motif_hits) == 0) {
    s_masked <- s_grid
    s_masked[!support_mask] <- 0
    if (max(s_masked) <= 0) return(empty)
    peaks  <- find_local_maxima(x_grid, s_masked,
                                 min_height = min_weight_frac * max(s_masked),
                                 min_dist   = min_peak_dist)
    if (length(peaks) == 0) return(empty)
    ev <- data.frame(
      tf              = tf_name,
      position        = x_grid[peaks],
      weight          = s_masked[peaks],
      motif_based     = FALSE,
      n_motifs_merged = 1L,
      stringsAsFactors = FALSE
    )
    ev$distance_to_nearest_grna <- vapply(ev$position,
                                          function(p) min(abs(p - pos_r)),
                                          numeric(1))
    ev$local_coverage <- approx(x_grid, c_grid, xout = ev$position,
                                 rule = 2)$y
    # Geometric safety net (see `max_grna_distance` param docs). In the
    # no-motif fallback, `support_mask` already zeroes s(x) outside the
    # edge_guard_frac band, so the hard distance cap is mostly redundant
    # here — but we still apply it for consistency with the motif branch,
    # and to keep behaviour deterministic if a caller disables
    # edge_guard_frac while leaving max_grna_distance set.
    if (is.finite(max_grna_distance)) {
      ev <- ev[ev$distance_to_nearest_grna <= max_grna_distance, ,
               drop = FALSE]
    }
    if (nrow(ev) == 0) return(empty)
    ev <- ev[order(ev$weight, decreasing = TRUE), , drop = FALSE]
    if (is.finite(max_events_per_tf) && nrow(ev) > max_events_per_tf) {
      ev <- ev[seq_len(max_events_per_tf), , drop = FALSE]
    }
    # Schema-consistency: route through the diagnostics closure (MLE +
    # bootstrap) even on the no-motif branch so the column set matches
    # the motif-based branch (NA-filled when each diagnostic is off,
    # populated when on). Without this the rbind-across-TFs in
    # predict_all_binding_events fails when a run mixes motif-anchored
    # and no-motif TFs.
    return(add_diag(ev))
  }

  # ── Motif-based: zone-based detection on β(x) ─────────────────────────
  # Rationale: scoring motifs one-by-one against max(β_at_motifs) is brittle
  # on slopes — a motif sitting on the shoulder of a strong peak can fall
  # below threshold even when the local β is clearly elevated above floor.
  # Biological question is "is this motif inside a plausible binding zone?",
  # not "does this motif sit at a β peak?". So: find contiguous above-
  # threshold intervals of β(x) on the grid, and emit EVERY motif inside
  # them. Zones with no motif still emit a single peak bubble so true
  # signal never goes uncalled.
  y_corr <- s_grid / pmax(c_grid, floor_val)
  # Suppress β outside the in-support region — any β value produced by the
  # floor clamp in near-zero-coverage tails is not a real binding signal.
  y_corr[!support_mask] <- 0
  y_max  <- max(y_corr)
  if (y_max <= 0) return(empty)
  threshold <- min_weight_frac * y_max

  # Contiguous above-threshold runs via rle() on the boolean mask.
  # support_mask is already baked into y_corr (zeroed out), so `above`
  # automatically stays inside the in-support region.
  above <- y_corr > threshold
  if (!any(above)) return(empty)
  rle_a  <- rle(above)
  ends   <- cumsum(rle_a$lengths)
  starts <- c(1L, head(ends + 1L, -1L))
  zone_s <- starts[rle_a$values]
  zone_e <- ends[rle_a$values]

  # Collect per-motif (and motif-less-zone) candidates across all zones.
  pos_cand <- numeric(0)
  w_cand   <- numeric(0)
  c_cand   <- numeric(0)
  mb_cand  <- logical(0)     # motif_based flag per candidate
  for (zi in seq_along(zone_s)) {
    a <- x_grid[zone_s[zi]]
    b <- x_grid[zone_e[zi]]
    in_zone <- motif_hits[motif_hits >= a & motif_hits <= b]
    if (length(in_zone) > 0) {
      # Emit motifs inside the zone whose local β is at least
      # `zone_peak_frac` × (zone's own peak β). Motifs on broad plateaus
      # with near-peak β all pass (GATA6 R7-R6 ramp: β ≈ 0.85 · peak,
      # passes); motifs on the low flanks of a sharp peak get filtered
      # (HOXB6 plateau shoulders at β ≈ 0.3 · peak, drop).
      s_k <- approx(x_grid, s_grid, xout = in_zone, rule = 2)$y
      c_k <- approx(x_grid, c_grid, xout = in_zone, rule = 2)$y
      b_k <- s_k / pmax(c_k, floor_val)
      # Zone-local peak: use the grid max inside the zone, not just the
      # motif-position max, so we don't bias the reference when motifs
      # happen to sit away from the true β crest.
      idxs <- zone_s[zi]:zone_e[zi]
      zone_peak_b <- max(y_corr[idxs])
      keep_zone <- which(b_k >= zone_peak_frac * zone_peak_b)
      if (length(keep_zone) > 0) {
        kept_pos <- in_zone[keep_zone]
        kept_w   <- b_k[keep_zone]
        # Optional PWM-score weighting: scale the emitted weight by a
        # function of the per-motif PWM score-fraction. Default "none"
        # leaves weights as β(motif_pos), preserving legacy behaviour.
        if (motif_score_weight != "none" && !is.null(score_lookup)) {
          sf <- score_lookup[as.character(kept_pos)]
          sf[is.na(sf)] <- 1     # missing-score motif -> no penalty
          factor <- switch(motif_score_weight,
                           linear = sf,
                           log    = 2^(sf - 1))
          kept_w <- kept_w * factor
        }
        pos_cand <- c(pos_cand, kept_pos)
        w_cand   <- c(w_cand,   kept_w)
        c_cand   <- c(c_cand,   c_k[keep_zone])
        mb_cand  <- c(mb_cand,  rep(TRUE, length(keep_zone)))
      } else {
        # All motifs in this zone were below the zone-local threshold —
        # fall back to a single peak bubble so the zone isn't orphaned.
        pk <- idxs[which.max(y_corr[idxs])]
        pos_cand <- c(pos_cand, x_grid[pk])
        w_cand   <- c(w_cand,   y_corr[pk])
        c_cand   <- c(c_cand,   c_grid[pk])
        mb_cand  <- c(mb_cand,  FALSE)
      }
    } else {
      # No motif inside this zone — emit one bubble at the zone's peak so
      # a real β peak without a JASPAR hit still gets surfaced.
      idxs <- zone_s[zi]:zone_e[zi]
      pk   <- idxs[which.max(y_corr[idxs])]
      pos_cand <- c(pos_cand, x_grid[pk])
      w_cand   <- c(w_cand,   y_corr[pk])
      c_cand   <- c(c_cand,   c_grid[pk])
      mb_cand  <- c(mb_cand,  FALSE)
    }
  }

  if (length(pos_cand) == 0) return(empty)

  # Merge closely-spaced candidates (amplitude-weighted centroid), same
  # merge rule as the default path so event counts stay comparable. A
  # merged cluster is motif_based if ANY of its members was motif-based.
  # NOTE: `used` is indexed by position in `pos_cand`, and `ord` is a
  # permutation of pos_cand indices. We MUST skip on `used[ord[j]]`, not
  # on `used[j]` — otherwise a candidate absorbed into an earlier cluster
  # can block an unrelated higher-weight candidate whose rank happens to
  # coincide with the absorbed index (the R7–R6 miss was exactly this).
  ord  <- order(w_cand, decreasing = TRUE)
  used <- rep(FALSE, length(pos_cand))
  events <- list()
  for (j in seq_along(ord)) {
    idx <- ord[j]
    if (used[idx]) next
    cluster_mask <- !used & abs(pos_cand - pos_cand[idx]) <= merge_dist
    cluster_idx  <- which(cluster_mask)
    if (length(cluster_idx) == 0) next   # defensive: should never trigger
    used[cluster_idx] <- TRUE
    w_sum  <- sum(w_cand[cluster_idx])
    # Position rule depends on merge_position.
    # - "argmax" : position of the highest-weighted candidate in the
    #   cluster. For motif_based clusters, that's the strongest motif's
    #   exact bp coordinate, so the bubble sits on a motif tick. For a
    #   zone-fallback single-candidate cluster (no motif survived), the
    #   top is the zone β-peak and argmax reduces to that peak.
    # - "centroid" : legacy amplitude-weighted mean position. Can land
    #   between 2-3 motifs that chained within merge_dist.
    top_in_cluster <- cluster_idx[which.max(w_cand[cluster_idx])]
    p_avg  <- if (merge_position == "argmax") {
      pos_cand[top_in_cluster]
    } else {
      sum(pos_cand[cluster_idx] * w_cand[cluster_idx]) / w_sum
    }
    d_nn   <- min(abs(p_avg - pos_r))
    # Interpolated C(x) at the merged centroid — more faithful than the
    # average of cluster members, which would be biased by the weights.
    c_loc  <- approx(x_grid, c_grid, xout = p_avg, rule = 2)$y
    events[[length(events) + 1]] <- data.frame(
      tf              = tf_name,
      position        = round(p_avg, 1),
      weight          = w_sum,
      motif_based     = any(mb_cand[cluster_idx]),
      n_motifs_merged = sum(mb_cand[cluster_idx]),
      distance_to_nearest_grna = round(d_nn, 1),
      local_coverage           = round(c_loc, 4),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, events)
  # Geometric safety net: drop events whose nearest gRNA is farther than
  # `max_grna_distance`. In the single-gRNA tail beyond the westernmost /
  # easternmost guide, β = s/C is flat (same Gaussian cancels top and
  # bottom) and the zone detector can spray bubbles all the way to the
  # edge_guard_frac mask. This cap clips those at a σ-scaled distance
  # from the guide layout — complementary to edge_guard_frac, which is
  # relative to max(C) and doesn't automatically rescale when kernel_sigma
  # changes.
  if (is.finite(max_grna_distance)) {
    out <- out[out$distance_to_nearest_grna <= max_grna_distance, ,
               drop = FALSE]
  }
  # Edge-gRNA dominance filter. For each candidate event position p, compute
  # the per-gRNA Gaussian weight w_i(p) = exp(-0.5 · ((p - r_i)/σ)²) and the
  # fractional contribution of each gRNA: f_i(p) = w_i(p) / Σ_j w_j(p). If
  # EITHER the westernmost (min r_i) or easternmost (max r_i) guide has
  # f_i(p) > edge_grna_weight_cap, drop the event — it's being inflated by
  # a boundary guide's kernel tail even though it may sit inside the
  # support mask and within max_grna_distance of some interior guide.
  # Applies only in coverage-aware paths that actually defined such a cap;
  # NULL keeps legacy behaviour.
  if (nrow(out) > 0 && !is.null(edge_grna_weight_cap) &&
      is.finite(edge_grna_weight_cap) && length(pos_r) >= 2) {
    left_r  <- min(pos_r)
    right_r <- max(pos_r)
    # Vectorised weight matrix: rows = events, cols = gRNAs
    d   <- outer(out$position, pos_r, FUN = function(p, r) p - r)
    w   <- exp(-0.5 * (d / kernel_sigma)^2)
    wsum <- rowSums(w)
    wsum[wsum <= 0] <- 1   # defensive: avoid divide-by-zero
    idx_left  <- which(pos_r == left_r)[1]
    idx_right <- which(pos_r == right_r)[1]
    frac_left  <- w[, idx_left]  / wsum
    frac_right <- w[, idx_right] / wsum
    edge_frac  <- pmax(frac_left, frac_right)
    keep <- edge_frac <= edge_grna_weight_cap
    if (!all(keep)) {
      message(sprintf(
        "    edge_grna_weight_cap(%.2f): dropped %d/%d %s event(s)",
        edge_grna_weight_cap, sum(!keep), length(keep), tf_name))
    }
    out <- out[keep, , drop = FALSE]
  }
  if (nrow(out) == 0) return(empty)
  out <- out[order(out$weight, decreasing = TRUE), , drop = FALSE]
  # Top-N cap (readability backstop). Applied AFTER merging so merged
  # centroids compete on their combined amplitude, not on per-motif β.
  if (is.finite(max_events_per_tf) && nrow(out) > max_events_per_tf) {
    out <- out[seq_len(max_events_per_tf), , drop = FALSE]
  }

  # Schema-consistent diagnostics pass — same closure the no-motif
  # branch uses. When position_stability is off this is a no-op; when
  # on it appends pos_stab_low/high/median + bootstrap_survival_frac.
  add_diag(out)
}

#' Batch: predict binding events for a TF list, using motif_results when present.
#'

#' @param cov_floor         Relative
#'   floor on the denominator: amplification is capped at ~ 1/cov_floor ×
#'   the max-coverage value. Default 0.05 (~20× amplification cap).
#' @param tfs (see function body).
#' @param long_data (see function body).
#' @param pos_map (see function body).
#' @param motif_results (see function body).
#' @param x_grid (see function body).
#' @param kernel_sigma (see function body).
#' @param min_weight_frac (see function body).
#' @param min_peak_dist (see function body).
#' @param merge_dist (see function body).
#' @param weight_mode (see function body).
#' @param edge_guard_frac (see function body).
#' @param zone_peak_frac (see function body).
#' @param max_events_per_tf (see function body).
#' @param max_grna_distance (see function body).
#' @param edge_grna_weight_cap (see function body).
#' @param n_bootstrap (see function body).
#' @param merge_position (see function body).
#' @param position_stability (see function body).
#' @param motif_score_weight (see function body).
#' @noRd
predict_all_binding_events <- function(tfs, long_data, pos_map, motif_results,
                                        x_grid          = seq(-2500, 500, by = 5),
                                        kernel_sigma    = 300,
                                        min_weight_frac = 0.15,
                                        min_peak_dist   = 150,
                                        merge_dist      = 100,
                                        weight_mode     = "z",
                                                                            cov_floor        = 0.05,
                                        edge_guard_frac  = 0.25,
                                        zone_peak_frac    = 0.50,
                                        max_events_per_tf = 30,
                                        # Coverage-aware only. See
                                        # predict_binding_events_coverage_aware()
                                        # for the full description.
                                        merge_position    = c("argmax", "centroid"),
                                        max_grna_distance = NULL,
                                        edge_grna_weight_cap = NULL,
                                        # Coverage-aware only. Wild bootstrap
                                        # position-stability diagnostic. See
                                        # predict_binding_events_coverage_aware()
                                        # for the full description / framing.
                                        position_stability = c("none", "wild_bootstrap"),
                                        n_bootstrap        = 200L,
                                        # PWM-score weighting of bubble size.
                                        # See run_caspex() / predict_binding_events*()
                                        # docstrings for details.
                                        motif_score_weight = c("none", "linear", "log")) {
  position_stability <- match.arg(position_stability)
  motif_score_weight <- match.arg(motif_score_weight)
  mode_tag <- "coverage-normalized per-motif (s/C)"
  message("\nPredicting binding events (", mode_tag, ", \u03c3=",
          kernel_sigma, " bp)...")
  # PWM-score weighting of per-motif basis columns / per-zone \u03b2 bubbles.
  # Default "none" is silent (legacy behaviour). When the user opts into
  # "linear" or "log" we surface it once, here, so the active scoring scheme
  # is visible alongside the binding-path mode it modifies.
  if (motif_score_weight != "none") {
    weight_desc <- switch(
      motif_score_weight,
      linear = paste0("linear (column / bubble weight = score_frac, ",
                      "i.e. raw fraction of max log-odds in [thresh, 1])"),
      log    = paste0("log (column / bubble weight = 1 + log2(1 + score_frac), ",
                      "compresses score range)"))
    message("  PWM-score weighting: ", motif_score_weight,
            " \u2014 ", weight_desc)
  }
  out <- list()
  for (tf in tfs) {
    hits <- if (!is.null(motif_results) && tf %in% names(motif_results))
      motif_results[[tf]]$hits else integer(0)
    scores <- if (!is.null(motif_results) && tf %in% names(motif_results))
      motif_results[[tf]]$score_frac else NULL
    ev <- predict_binding_events_coverage_aware(
        tf, long_data, pos_map, hits,
        motif_scores       = scores,
        motif_score_weight = motif_score_weight,
        kernel_sigma    = kernel_sigma,
        min_weight_frac = min_weight_frac,
        min_peak_dist   = min_peak_dist,
        merge_dist      = merge_dist,
        weight_mode     = weight_mode,
        cov_floor       = cov_floor,
        edge_guard_frac = edge_guard_frac,
        zone_peak_frac    = zone_peak_frac,
        max_events_per_tf = max_events_per_tf,
        merge_position    = merge_position,
        max_grna_distance = max_grna_distance,
        edge_grna_weight_cap = edge_grna_weight_cap,
        position_stability = position_stability,
        n_bootstrap        = n_bootstrap,
        x_grid          = x_grid
      )
    if (nrow(ev) > 0) {
      tag <- if (length(hits) > 0)
        paste0(" [motif-anchored, ", length(hits), " candidates]")
      else " [no motif; peak detection]"
      message("  ", tf, ": ", nrow(ev), " event(s)", tag)
      out[[tf]] <- ev
    }
  }
  if (length(out) == 0) {
    base_cols <- data.frame(tf = character(), position = numeric(),
                            weight = numeric(), motif_based = logical(),
                            n_motifs_merged = integer(),
                            stringsAsFactors = FALSE)
    base_cols$distance_to_nearest_grna <- numeric(0)
    base_cols$local_coverage <- numeric(0)
    return(base_cols)
  }
  # rbind across TFs — coverage-aware returns extra cols, smoothed doesn't.
  # Unify schema so the CSV always has the same columns regardless of which
  # branch was taken for each TF (motif-less fallback on coverage mode still
  # carries the extra cols because we tagged them in the fallback).
  do.call(rbind, out)
}

#' Per-TF deconvolution detail plot: signal + region logFCs + predicted events.
#'
#' Layout (top -> bottom):
#'   y > 0  : CasPEX signal ribbon (histogram-style), per-region logFC lollipops
#'   y = 0  : baseline
#'   y < 0  : "called peak track" — motif ticks (upper sub-lane), event circles
#'            sized by NNLS weight (lower sub-lane). Events sit *under* the
#'            histogram directly below their predicted position, with a thin
#'            dotted connector up to the signal peak they explain.
#' @param tf_name (see function body).
#' @param long_data (see function body).
#' @param pos_map (see function body).
#' @param motif_hits (see function body).
#' @param kernel_sigma (see function body).
#' @param min_weight_frac (see function body).
#' @param upstream (see function body).
#' @param downstream (see function body).
#' @param weight_mode (see function body).

#' @param cov_floor (see function body).
#' @param edge_guard_frac (see function body).
#' @param zone_peak_frac (see function body).
#' @param max_events_per_tf (see function body).
#' @param merge_position (see function body).
#' @param max_grna_distance (see function body).
#' @param edge_grna_weight_cap (see function body).
#' @param n_bootstrap (see function body).
#' @param motif_score_weight (see function body).
#' @param position_stability (see function body).
#' @param chipatlas_peaks (see function body).
#' @param motif_scores (see function body).
#' @export
plot_binding_deconvolution <- function(tf_name, long_data, pos_map, motif_hits,
                                        kernel_sigma    = 300,
                                        min_weight_frac = 0.15,
                                        upstream        = 2500,
                                        downstream      = 500,
                                        weight_mode     = "z",
                                        cov_floor       = 0.05,
                                        edge_guard_frac = 0.25,
                                        zone_peak_frac    = 0.50,
                                        max_events_per_tf = 30,
                                        merge_position    = c("argmax", "centroid"),
                                        max_grna_distance = NULL,
                                        edge_grna_weight_cap = NULL,
                                        # Coverage-aware only. Wild bootstrap
                                        # position-stability diagnostic. See
                                        # predict_binding_events_coverage_aware()
                                        # for the framing — when enabled, a
                                        # "Position support" sub-track renders
                                        # below the event bubbles, plotting
                                        # the central 95% bootstrap range and
                                        # median for each event. Off by default
                                        # because bootstrap is N×slower.
                                        position_stability = c("none", "wild_bootstrap"),
                                        n_bootstrap        = 200L,
                                        # ChIP-Atlas overlay: data.frame as returned
                                        # by fetch_chipatlas_peaks() for THIS TF
                                        # (cols srx, cell_type, start_rel, end_rel,
                                        # ...). NULL skips the sub-lane.
                                        chipatlas_peaks   = NULL,
                                        # PWM-score weighting of bubble size
                                        # (mirror of run_caspex / predict_*'s
                                        # motif_score_weight + motif_scores).
                                        # Without these, the deck plot would
                                        # silently rerun the deconvolution in
                                        # motif_score_weight="none" mode and
                                        # disagree with result$binding_events
                                        # whenever the parent run used
                                        # "linear" or "log".
                                        motif_scores      = NULL,
                                        motif_score_weight = c("none",
                                                                "linear",
                                                                "log")) {
  position_stability <- match.arg(position_stability)
  motif_score_weight <- match.arg(motif_score_weight)
  x_grid <- seq(-upstream, downstream, by = 5)
  sig    <- build_caspex_signal(tf_name, long_data, pos_map, x_grid,
                                 kernel_sigma, weight_mode)
  events <- predict_binding_events_coverage_aware(
      tf_name, long_data, pos_map, motif_hits,
      motif_scores       = motif_scores,
      motif_score_weight = motif_score_weight,
      kernel_sigma    = kernel_sigma,
      min_weight_frac = min_weight_frac,
      weight_mode     = weight_mode,
      cov_floor       = cov_floor,
      edge_guard_frac = edge_guard_frac,
      zone_peak_frac    = zone_peak_frac,
      max_events_per_tf = max_events_per_tf,
      merge_position    = merge_position,
      max_grna_distance = max_grna_distance,
      edge_grna_weight_cap = edge_grna_weight_cap,
      position_stability = position_stability,
      n_bootstrap        = n_bootstrap,
      x_grid          = x_grid)
  rd      <- sig$region_data
  sig_max <- max(c(sig$y, rd$lfc, 1), na.rm = TRUE)
  # Clip motif ticks to the gRNA-supported region (one kernel σ past the
  # outermost guides on each side). Past that envelope, CasPEX has no
  # signal to pair with sequence-level motifs, so showing ticks there is
  # visual noise that invites misreading. Without this clip, the ggplot
  # panel auto-expands whenever the ChIP-Atlas band is drawn at the full
  # [-upstream, downstream] width, which makes far-upstream motif clusters
  # (e.g. the NR2F1 cluster past -2000 on ATP7B where R7 ≈ -1900 is the
  # westernmost guide) suddenly visible even though no events are — or
  # ever could be — called in that region. The other mini-browser decks
  # (07/09/12) intentionally show the full motif window; this detail plot
  # scopes to what the detector can actually resolve.
  pos_r_detect <- sort(as.numeric(pos_map[!is.na(pos_map)]))
  if (length(pos_r_detect) >= 1) {
    left_cut  <- min(pos_r_detect) - kernel_sigma
    right_cut <- max(pos_r_detect) + kernel_sigma
  } else {
    left_cut  <- -upstream
    right_cut <-  downstream
  }
  motif_hits_in <- motif_hits[!is.na(motif_hits) &
                                motif_hits >= max(-upstream, left_cut) &
                                motif_hits <= min(downstream, right_cut)]
  sig_df <- data.frame(x = x_grid, y = sig$y)

  # Called-peak track y-positions. Float the whole strip BELOW the most-
  # negative logFC lollipop so depletion bars are never buried by motif
  # ticks or event circles. When all lfc are >= 0 this collapses to the
  # previous layout (strip at y < 0).
  neg_floor   <- min(c(0, rd$lfc), na.rm = TRUE)
  track_gap   <- 0.06 * sig_max          # breathing room below lowest lollipop
  track_top   <- neg_floor - track_gap
  motif_y_top <- track_top - 0.04 * sig_max
  motif_y_bot <- track_top - 0.10 * sig_max
  event_y     <- track_top - 0.20 * sig_max
  track_bot   <- track_top - 0.28 * sig_max

  # ChIP-Atlas sub-lane sits BELOW the event bubbles, stacked one row per
  # SRX experiment. Height scales with the number of experiments but is
  # capped so the deck stays readable. If no peaks were passed, these
  # y-coords collapse to track_bot and the lane draws nothing.
  ca_rows <- if (!is.null(chipatlas_peaks) && nrow(chipatlas_peaks) > 0) {
    unique(chipatlas_peaks$srx)
  } else character(0)
  n_ca_rows <- length(ca_rows)
  # Total height for the ChIP-Atlas band (as a fraction of sig_max).
  ca_band_h <- if (n_ca_rows == 0) 0 else min(0.50, 0.025 * n_ca_rows + 0.06)
  ca_gap    <- if (n_ca_rows == 0) 0 else 0.03 * sig_max
  ca_top    <- track_bot - ca_gap
  ca_bot    <- ca_top - ca_band_h * sig_max

  p <- ggplot() +
    # --- upper panel: signal + region logFCs --------------------------------
    geom_area(data = sig_df, aes(x = x, y = y),
              fill = COLS$guide, alpha = 0.25) +
    geom_line(data = sig_df, aes(x = x, y = y),
              color = COLS$neutral, linewidth = 0.4) +
    geom_segment(data = rd, aes(x = pos, xend = pos, y = 0, yend = lfc),
                 color = COLS$neutral, linewidth = 0.5, alpha = 0.7) +
    geom_point(data = rd, aes(x = pos, y = lfc),
               color = COLS$neutral, size = 3) +
    # Region labels (R1, R2, ...). Positive-logFC stems get the label
    # ABOVE the head (vjust = -0.9). Negative-logFC stems get the label
    # to the RIGHT side (vjust = 0.5 = vertically centred on the head;
    # hjust = -0.3 = nudged right of the head) — this keeps the label
    # out of the called-peak strip below the baseline, which previously
    # got cluttered when several depletion lollipops stacked their
    # labels into the same band.
    geom_text(data = transform(rd,
                                 .vj = ifelse(lfc >= 0, -0.9,  0.5),
                                 .hj = ifelse(lfc >= 0,  0.5, -0.3)),
              aes(x = pos, y = lfc, label = region,
                  vjust = .vj, hjust = .hj),
              size = 2.8, fontface = "bold",
              color = COLS$neutral) +
    # --- baseline + TSS -----------------------------------------------------
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = COLS$tss, linewidth = 0.7) +
    # --- lower panel container (subtle shaded strip) ------------------------
    # Strip floats below the most-negative lfc lollipop (track_top); this
    # keeps depletion bars visible above the called-peak lane.
    annotate("rect",
             xmin = -upstream, xmax = downstream,
             ymin = track_bot, ymax = track_top,
             fill = "grey97", color = NA)

  # JASPAR motif ticks — upper sub-lane of the called-peak track
  if (length(motif_hits_in) > 0) {
    p <- p + geom_segment(data = data.frame(x = motif_hits_in),
                           aes(x = x, xend = x,
                               y = motif_y_bot, yend = motif_y_top),
                           color = COLS$mid, linewidth = 0.4, alpha = 0.7)
  }

  # Predicted binding events — lower sub-lane. Position the circle on
  # the event_y baseline, with a dotted connector up to the signal value
  # at that x so the viewer sees "this peak → this call".
  if (nrow(events) > 0) {
    ev <- events
    ev$y_base   <- event_y
    # signal height at the event position (for the dotted connector)
    ev$y_signal <- approx(x_grid, sig$y, xout = ev$position, rule = 2)$y
    p <- p +
      geom_segment(data = ev,
                   aes(x = position, xend = position,
                       y = y_base, yend = y_signal),
                   color = COLS$neutral, linewidth = 0.3,
                   linetype = "dotted", alpha = 0.5) +
      geom_point(data = ev,
                 aes(x = position, y = y_base,
                     size = weight, fill = motif_based),
                 shape = 21, color = "black", stroke = 0.5, alpha = 0.95) +
      # Bubble bp-coordinate label, centred INSIDE the bubble (vjust =
      # 0.5, hjust = 0.5). Black + bold so it stays legible on top of
      # the size-and-fill-coded circle regardless of motif_based status.
      # Bare integer (no "bp" suffix) — keeps neighbouring labels short
      # when events are closely spaced; x-axis title already declares
      # the unit.
      geom_text(data = ev,
                aes(x = position, y = y_base,
                    label = sprintf("%+.0f", position)),
                size = 2.6, fontface = "bold",
                vjust = 0.5, hjust = 0.5,
                color = "black")
  }

  # --- ChIP-Atlas stacked-experiment sub-lane ---------------------------------
  # Each SRX occupies a thin row in the (ca_bot, ca_top) band. Peak intervals
  # are drawn as short horizontal bars clipped to the promoter window. A thin
  # separator strip above and a row label on the right-hand margin help read
  # the stack when it's crowded.
  if (n_ca_rows > 0) {
    row_h <- (ca_top - ca_bot) / n_ca_rows
    # Newest-first by SRX number (matches the cap ordering upstream).
    ord   <- order(suppressWarnings(as.integer(sub("SRX", "", ca_rows))),
                   decreasing = TRUE)
    srx_levels <- ca_rows[ord]
    ca <- chipatlas_peaks
    ca$srx <- factor(ca$srx, levels = srx_levels)
    ca$y_mid <- ca_top - (as.integer(ca$srx) - 0.5) * row_h
    # Clip to the plotted window so a peak that spills past doesn't drag the
    # lane off-axis. Also widen 1-bp intervals to at least ~10 bp so they
    # render visibly.
    ca$xs <- pmax(ca$start_rel, -upstream)
    ca$xe <- pmin(ca$end_rel,    downstream)
    ca$xs <- pmin(ca$xs, ca$xe - 5)   # ensure non-zero width
    p <- p +
      # Background band so the sub-lane is visible even when peaks are sparse.
      annotate("rect",
               xmin = -upstream, xmax = downstream,
               ymin = ca_bot,    ymax = ca_top,
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

  # Y-axis label tracks the active weight_mode rather than hard-coding
  # "logFC-weighted". Previously the label always read "logFC-weighted" even
  # when the curve was actually z-weighted (the default since the ATP7B-style
  # default switch) or t-weighted, which made magnitude reads of the curve
  # confusing — e.g., a z-weighted Myers curve peaking at ~100 looked like
  # logFC=100 even though the per-region lfc stems on the same plot show 1-3.
  y_label <- switch(weight_mode,
    z          = "CasPEX signal (z-weighted, signed z from p-value)",
    mod_t      = "CasPEX signal (mod-t weighted)",
    lfc_pos    = "CasPEX signal (positive logFC-weighted)",
    lfc_signed = "CasPEX signal (signed logFC-weighted)",
    lfc_x_negp = "CasPEX signal (logFC \u00d7 -log10 p-weighted)",
    paste0("CasPEX signal (", weight_mode, "-weighted)")
  )

  p +
    scale_fill_manual(values = c(`TRUE` = COLS$high, `FALSE` = COLS$low),
                      labels = c(`TRUE` = "motif-anchored",
                                 `FALSE` = "no motif (peak)"),
                      name = "Event type", drop = FALSE) +
    scale_size_continuous(range = c(3, 9), guide = "none") +
    scale_y_continuous(
      expand = expansion(mult = c(0.08, 0.08)),
      # Hide negative ticks/labels — the "called peak" lane is a legend-like
      # strip, not a numeric reading
      breaks = function(lim) pretty(c(0, lim[2]))
    ) +
    labs(x = "Position (bp, TSS-relative)",
         y = y_label,
         title = paste0(tf_name, " \u2014 binding event deconvolution",
                        " (coverage-aware)"),
         subtitle = sprintf(
           "%d event(s) | kernel \u03c3 = %d bp | %d motif candidate(s)  \u00b7  called-peak track below baseline",
           nrow(events), kernel_sigma, length(motif_hits_in))) +
    theme_caspex()
}

# =============================================================================
# SECTION 6: Plots
# =============================================================================

#' Plot 1 — gRNA positions ruler
#' @param pos_map (see function body).
#' @param gene_info (see function body).
#' @param upstream (see function body).
#' @param downstream (see function body).
#' @export
plot_grna_positions <- function(pos_map, gene_info, upstream, downstream) {
  valid <- pos_map[!is.na(pos_map)]
  df    <- data.frame(region = names(valid), pos = as.numeric(valid))

  ggplot(df, aes(x = pos, y = 1)) +
    annotate("rect",
             xmin = -upstream, xmax = downstream,
             ymin = 0.80, ymax = 1.20,
             fill = "#EEF4FB", color = "#5B8DB8", linewidth = 0.4) +
    geom_segment(aes(xend = pos, y = 0.72, yend = 1.28),
                 color = COLS$guide, linewidth = 1.4) +
    geom_point(size = 5, color = COLS$guide) +
    geom_text(aes(label = paste0(region, "\n", pos, " bp")),
              vjust = -0.5, size = 3, fontface = "bold",
              color = COLS$guide, lineheight = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = COLS$tss, linewidth = 1) +
    annotate("text", x = 0, y = 0.60, label = "TSS",
             color = COLS$tss, size = 3.8, fontface = "bold") +
    scale_y_continuous(limits = c(0.45, 1.55)) +
    labs(x = "Position (bp, relative to TSS)", y = NULL,
         title = paste0("gRNA cut sites \u2014 ", gene_info$name),
         subtitle = paste0("chr", gene_info$chr, "  |  ",
                           ifelse(gene_info$strand == 1, "+", "-"), " strand  |  TSS = ",
                           format(gene_info$tss, big.mark = ","))) +
    theme_caspex() +
    theme(axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank())
}

#' Plot 2 — Spatial binding track (histogram-style)
#'
#' Each TF shown as a Gaussian ridge centred on the predicted binding centroid,
#' height proportional to composite score, width proportional to binding spread.
#' Bars stacked in a genomic coordinate space.
#' @param spatial_df (see function body).
#' @param pos_map (see function body).
#' @param top_n (see function body).
#' @param upstream (see function body).
#' @param downstream (see function body).
#' @export
plot_spatial_track <- function(spatial_df, pos_map,
                                top_n = 25, upstream = 2500, downstream = 500) {
  df  <- head(spatial_df[order(spatial_df$composite, decreasing = TRUE), ], top_n)
  pos <- as.numeric(pos_map[!is.na(pos_map)])
  x   <- seq(-upstream, downstream, length.out = 800)

  # Build ribbon data for each TF
  ribbons <- lapply(seq_len(nrow(df)), function(i) {
    row     <- df[i, ]
    sigma   <- max(row$spread, 50)
    density <- exp(-0.5 * ((x - row$centroid) / sigma)^2)
    density <- density * row$composite   # height = composite score
    data.frame(x        = x,
               y_lo     = i - 0.45,
               y_hi     = i - 0.45 + density / max(density) * 0.88,
               protein  = row$protein,
               centroid = row$centroid,
               composite = row$composite,
               rank     = i)
  })
  rib_df <- do.call(rbind, ribbons)

  # Centroid points
  cent_df <- data.frame(
    rank    = seq_len(nrow(df)),
    centroid = df$centroid,
    composite = df$composite,
    protein = df$protein,
    sig     = df$sig_any
  )

  # TF labels on y-axis
  y_labels <- setNames(seq_len(nrow(df)), df$protein)

  ggplot() +
    # sgRNA guide lines
    geom_vline(xintercept = pos,
               linetype = "dotted", color = COLS$guide, linewidth = 0.5, alpha = 0.7) +
    # TSS line
    geom_vline(xintercept = 0,
               linetype = "dashed", color = COLS$tss, linewidth = 0.8) +
    # Ribbons (Gaussian binding density per TF)
    geom_ribbon(data = rib_df,
                aes(x = x, ymin = y_lo, ymax = y_hi,
                    fill = composite, group = protein),
                alpha = 0.65) +
    # Baseline per TF
    geom_hline(data = cent_df,
               aes(yintercept = rank - 0.45),
               color = "grey80", linewidth = 0.2) +
    # Centroid dot
    geom_point(data = cent_df,
               aes(x = centroid, y = rank - 0.45 + 0.45,
                   color = composite, shape = sig),
               size = 2.5) +
    scale_fill_gradient2(low  = COLS$low,
                         mid  = COLS$mid,
                         high = COLS$high,
                         midpoint = median(df$composite),
                         name = "Composite\nscore") +
    scale_color_gradient2(low  = COLS$low,
                          mid  = COLS$mid,
                          high = COLS$high,
                          midpoint = median(df$composite),
                          guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 16),
                       labels = c(`TRUE` = paste0("p\u2264sig"), `FALSE` = "ns"),
                       name = "Significance") +
    scale_y_continuous(breaks = seq_len(nrow(df)),
                       labels = df$protein,
                       expand = c(0.01, 0.01)) +
    scale_x_continuous(limits = c(-upstream, downstream), expand = c(0, 0)) +
    annotate("text", x = 5, y = nrow(df) + 0.3,
             label = "TSS", color = COLS$tss, size = 3, hjust = 0) +
    labs(x = "Position (bp, relative to TSS)", y = NULL,
         title = "Predicted TF binding zones",
         subtitle = paste0("Gaussian ridge centred on weighted spatial centroid  |  ",
                           "Height = composite enrichment score  |  ",
                           "Width = binding spread")) +
    theme_caspex(base_size = 10) +
    theme(panel.grid.major.x = element_line(color = "grey90"),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 8, face = "bold"))
}

#' Plot 3 — Per-TF enrichment heatmap across regions
#' @param long_data (see function body).
#' @param spatial_df (see function body).
#' @param pos_map (see function body).
#' @param top_n (see function body).
#' @export
plot_heatmap <- function(long_data, spatial_df, pos_map, top_n = 30) {
  tfs   <- head(spatial_df$protein, top_n)
  df    <- long_data[long_data$protein %in% tfs &
                       long_data$region %in% names(pos_map), ]

  # Add position labels to region axis
  region_labels <- sapply(names(pos_map), function(r) {
    if (!is.na(pos_map[r])) paste0(r, "\n(", pos_map[r], "bp)") else r
  })

  df$protein <- factor(df$protein, levels = rev(tfs))
  df$region  <- factor(df$region,
                        levels = names(pos_map),
                        labels = region_labels[names(pos_map)])

  ggplot(df, aes(x = region, y = protein, fill = lfc)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = ifelse(pval <= 0.05, sprintf("%.1f", lfc), "")),
              size = 2.6, color = "white", fontface = "bold") +
    scale_fill_gradient2(low = "#3A86FF", mid = "grey97", high = "#E63946",
                         midpoint = 0, na.value = "grey92",
                         name = "logFC") +
    labs(x = NULL, y = NULL,
         title = "Enrichment heatmap: top TFs \u00d7 regions",
         subtitle = "Bold values = p \u2264 0.05") +
    theme_caspex(base_size = 10) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 9))
}

#' Plot 4 — GC content track for the promoter
#' @param promoter_info (see function body).
#' @param window (see function body).
#' @param step (see function body).
#' @noRd
plot_gc_track <- function(promoter_info, window = 50, step = 10) {
  seq <- promoter_info$seq
  n   <- nchar(seq)
  tss <- promoter_info$tss_offset + 1
  pos <- seq(1, n - window, by = step)

  gc <- sapply(pos, function(i) {
    s <- strsplit(substr(seq, i, i + window - 1), "")[[1]]
    mean(s %in% c("G", "C"))
  })

  df <- data.frame(pos = pos - tss, gc = gc)
  overall_gc <- mean(strsplit(seq, "")[[1]] %in% c("G", "C")) * 100

  ggplot(df, aes(x = pos, y = gc)) +
    geom_area(fill = COLS$guide, alpha = 0.4) +
    geom_line(color = COLS$neutral, linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = COLS$tss, linewidth = 0.7) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                       limits = c(0, 1)) +
    labs(x = "Position (bp, TSS-relative)", y = "GC content",
         title = sprintf("Promoter GC content  |  %d bp  |  Overall: %.1f%%  |  Window: %d bp",
                         n, overall_gc, window)) +
    theme_caspex()
}

#' Plot 5 — Binding prediction track with JASPAR motif overlay
#'
#' One horizontal lane per TF. Inside each lane:
#'   - Shaded Gaussian ribbon = CasPEX spatial enrichment prediction
#'   - Vertical ticks below = JASPAR motif hit positions
#'   - Color convergence between ribbon peak and motif cluster = high-confidence site
#' @param spatial_df (see function body).
#' @param motif_results (see function body).
#' @param pos_map (see function body).
#' @param promoter_info (see function body).
#' @param upstream (see function body).
#' @param downstream (see function body).
#' @export
plot_motif_track <- function(spatial_df, motif_results, pos_map,
                              promoter_info, upstream = 2500, downstream = 500) {
  tfs <- intersect(names(motif_results), spatial_df$protein)
  if (length(tfs) == 0) {
    message("No overlap between motif scan TFs and spatial model. Skipping motif track.")
    return(NULL)
  }

  x      <- seq(-upstream, downstream, length.out = 1000)
  pos    <- as.numeric(pos_map[!is.na(pos_map)])
  n_tfs  <- length(tfs)
  colors <- c("#E63946","#2A9D8F","#F4A261","#457B9D","#6D6875",
              "#C77DFF","#E9C46A","#264653","#B5838D","#A8DADC")

  p <- ggplot() +
    # Promoter span background
    annotate("rect",
             xmin = -upstream, xmax = downstream,
             ymin = 0.3, ymax = n_tfs + 0.7,
             fill = "grey98", color = NA) +
    # sgRNA positions
    geom_vline(xintercept = pos,
               linetype = "dotted", color = COLS$guide,
               linewidth = 0.4, alpha = 0.6) +
    # TSS
    geom_vline(xintercept = 0,
               linetype = "dashed", color = COLS$tss, linewidth = 0.9)

  # Lane separators
  for (i in seq_len(n_tfs)) {
    p <- p + annotate("rect",
                      xmin = -upstream, xmax = downstream,
                      ymin = i - 0.5, ymax = i + 0.5,
                      fill = if (i %% 2 == 0) "grey96" else "white",
                      color = NA, alpha = 0.6)
  }

  # Per-lane layout (genome-browser style — ribbon on top, motif ticks below):
  #   [i-0.50 ------------------------ i+0.50]   lane extent
  #   [i-0.42 ... i-0.22]                        motif sub-lane (ticks BELOW)
  #   [i-0.15 ... i+0.45]                        enrichment ribbon (ABOVE ticks)
  ribbon_base <- -0.15
  ribbon_h    <-  0.60
  motif_y0    <- -0.42
  motif_y1    <- -0.22

  for (i in seq_along(tfs)) {
    tf  <- tfs[i]
    col <- colors[((i - 1) %% length(colors)) + 1]
    sp  <- spatial_df[spatial_df$protein == tf, ]

    # Enrichment ribbon + centroid marker (confined to THIS TF's lane only)
    if (nrow(sp) > 0) {
      sigma   <- max(sp$spread[1], 60)
      density <- exp(-0.5 * ((x - sp$centroid[1]) / sigma)^2)
      rib_df  <- data.frame(x = x,
                             y_lo = i + ribbon_base,
                             y_hi = i + ribbon_base + density * ribbon_h)
      # Centroid tick spans only the ribbon height (not the motif sub-lane)
      cent_df <- data.frame(xc = sp$centroid[1],
                             y0 = i + ribbon_base,
                             y1 = i + ribbon_base + ribbon_h)
      p <- p +
        geom_ribbon(data = rib_df,
                    aes(x = x, ymin = y_lo, ymax = y_hi),
                    fill = col, alpha = 0.22) +
        geom_segment(data = cent_df,
                     aes(x = xc, xend = xc, y = y0, yend = y1),
                     color = col, linewidth = 0.7, alpha = 0.7,
                     linetype = "solid")
    }

    # JASPAR motif ticks — sub-lane BELOW the ribbon
    hits <- motif_results[[tf]]$hits
    if (length(hits) > 0) {
      hits_in <- hits[hits >= -upstream & hits <= downstream]
      if (length(hits_in) > 0) {
        # Faint sub-lane strip to frame the motif ticks
        p <- p + annotate("rect",
                          xmin = -upstream, xmax = downstream,
                          ymin = i + motif_y0, ymax = i + motif_y1,
                          fill = col, alpha = 0.06, color = NA)
        # Bake lane y into the data so `i` isn't late-bound by aes()
        tick_df <- data.frame(x  = hits_in,
                              y0 = i + motif_y0,
                              y1 = i + motif_y1,
                              yc = i + (motif_y0 + motif_y1) / 2)
        p <- p +
          geom_segment(data = tick_df,
                       aes(x = x, xend = x, y = y0, yend = y1),
                       color = col, linewidth = 0.7, alpha = 0.9) +
          geom_point(data = tick_df,
                     aes(x = x, y = yc),
                     color = col, size = 1.1, alpha = 0.9)
      }
    }

    # TF name label (left margin)
    p <- p + annotate("text",
                      x = -upstream - 30, y = i,
                      label = tf, hjust = 1, vjust = 0.5,
                      size = 3.2, color = col, fontface = "bold")
  }

  # sgRNA region labels along top
  valid_pos <- pos_map[!is.na(pos_map)]
  for (r in names(valid_pos)) {
    p <- p + annotate("text",
                      x = valid_pos[[r]], y = n_tfs + 0.55,
                      label = r, size = 2.6, color = COLS$guide,
                      angle = 90, vjust = 0.5, hjust = 0)
  }

  p +
    annotate("text", x = 10, y = n_tfs + 0.55,
             label = "TSS", color = COLS$tss, size = 3, hjust = 0) +
    scale_x_continuous(limits    = c(-upstream - 200, downstream),
                       expand    = c(0, 0),
                       labels    = comma) +
    scale_y_continuous(limits    = c(0, n_tfs + 1),
                       expand    = c(0, 0)) +
    labs(x = "Position (bp, relative to TSS)", y = NULL,
         title = "TF binding prediction: CasPEX enrichment + JASPAR motif overlay",
         subtitle = paste0("Upper ribbon = spatial enrichment (Gaussian)  |  ",
                           "Lower sub-lane ticks = JASPAR motif hits  |  ",
                           "Vertical line = predicted binding centroid")) +
    theme_caspex(base_size = 11) +
    theme(axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          plot.margin = margin(5, 10, 5, upstream / 10 + 30))
}

#' Generic genome track for a subset of TFs, with optional JASPAR overlay
#'
#' One lane per TF, stacked vertically. Each lane shows:
#'   * Gaussian enrichment ribbon centred on the spatial centroid
#'   * Lane-local centroid tick
#'   * Optional motif hit ticks inside the lane (only if motif_results provided)
#' Global reference lines (TSS dashed; sgRNA dotted) span the full plot.
#'
#' @param tfs            Character vector of TF symbols to plot (in order).
#' @param spatial_df     From run_spatial_model()
#' @param pos_map        From match_all_grnas()
#' @param promoter_info  From fetch_promoter_seq()
#' @param motif_results  From run_motif_scan(), or NULL to hide motif ticks
#' @param lane_labels    Optional named char vec (names = TF symbol, values =
#'                       label to show on the left). Defaults to TF symbol.
#' @param upstream       bp upstream of TSS
#' @param downstream     bp downstream of TSS
#' @param title,subtitle Plot text (subtitle auto-generated if NULL)
#' @return ggplot, or NULL if no TFs to plot
#' @param chipatlas_peaks (see function body).
#' @export
plot_tf_track <- function(tfs, spatial_df, pos_map, promoter_info,
                          motif_results = NULL,
                          lane_labels   = NULL,
                          upstream = 2500, downstream = 500,
                          title = "Predicted TF binding zones",
                          subtitle = NULL,
                          # ChIP-Atlas overlay: named list (TF -> data.frame of
                          # windowed peaks, cols start_rel/end_rel/srx/...), or
                          # NULL to skip the sub-lane.
                          chipatlas_peaks = NULL) {
  tfs <- as.character(tfs)
  tfs <- tfs[tfs %in% as.character(spatial_df$protein)]
  tfs <- unique(tfs)
  if (length(tfs) == 0) return(NULL)

  if (is.null(lane_labels) || length(lane_labels) == 0) {
    lane_labels <- setNames(tfs, tfs)
  } else {
    lbl <- lane_labels[tfs]
    lbl[is.na(lbl) | lbl == ""] <- tfs[is.na(lbl) | lbl == ""]
    lane_labels <- lbl
  }

  x      <- seq(-upstream, downstream, length.out = 1000)
  pos    <- as.numeric(pos_map[!is.na(pos_map)])
  n_tfs  <- length(tfs)
  colors <- c("#E63946","#2A9D8F","#F4A261","#457B9D","#6D6875",
              "#C77DFF","#E9C46A","#264653","#B5838D","#A8DADC")

  if (is.null(subtitle)) {
    subtitle <- if (is.null(motif_results)) {
      "Spatial enrichment only  |  shaded ribbon = Gaussian binding zone"
    } else {
      paste0("Upper ribbon = spatial enrichment (Gaussian)  |  ",
             "Lower sub-lane ticks = JASPAR motif hits  |  ",
             "Vertical line = binding centroid",
             if (!is.null(chipatlas_peaks)) "  |  dark strip = ChIP-Atlas union peaks" else "")
    }
  }

  p <- ggplot() +
    annotate("rect",
             xmin = -upstream, xmax = downstream,
             ymin = 0.3, ymax = n_tfs + 0.7,
             fill = "grey98", color = NA) +
    geom_vline(xintercept = pos,
               linetype = "dotted", color = COLS$guide,
               linewidth = 0.4, alpha = 0.45) +
    geom_vline(xintercept = 0,
               linetype = "dashed", color = COLS$tss, linewidth = 0.9)

  for (i in seq_len(n_tfs)) {
    p <- p + annotate("rect",
                      xmin = -upstream, xmax = downstream,
                      ymin = i - 0.5, ymax = i + 0.5,
                      fill = if (i %% 2 == 0) "grey96" else "white",
                      color = NA, alpha = 0.6)
  }

  # Per-lane layout (genome-browser style):
  #   [i-0.50 ------------------------ i+0.50]   lane extent
  #   [i-0.48 ... i-0.44]                        ChIP-Atlas union-peak strip
  #   [i-0.42 ... i-0.22]                        motif sub-lane (ticks BELOW)
  #   [i-0.15 ... i+0.45]                        enrichment ribbon (ABOVE ticks)
  ribbon_base <- -0.15   # ribbon baseline, relative to lane center (i)
  ribbon_h    <-  0.60   # max ribbon height
  motif_y0    <- -0.42   # motif tick bottom
  motif_y1    <- -0.22   # motif tick top
  chip_y0     <- -0.48   # ChIP-Atlas strip bottom
  chip_y1     <- -0.44   # ChIP-Atlas strip top

  # Interval union helper — collapses overlapping peak intervals to a minimal
  # set of non-overlapping ones. Input is a data.frame with start_rel/end_rel
  # columns; output is a 2-column data.frame (xs, xe).
  .union_intervals <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(data.frame(xs = numeric(0), xe = numeric(0)))
    iv <- data.frame(xs = pmin(df$start_rel, df$end_rel),
                     xe = pmax(df$start_rel, df$end_rel))
    iv <- iv[order(iv$xs), , drop = FALSE]
    merged_xs <- iv$xs[1]; merged_xe <- iv$xe[1]
    out_xs <- numeric(0); out_xe <- numeric(0)
    if (nrow(iv) >= 2) {
      for (k in 2:nrow(iv)) {
        if (iv$xs[k] <= merged_xe) {
          merged_xe <- max(merged_xe, iv$xe[k])
        } else {
          out_xs <- c(out_xs, merged_xs); out_xe <- c(out_xe, merged_xe)
          merged_xs <- iv$xs[k]; merged_xe <- iv$xe[k]
        }
      }
    }
    out_xs <- c(out_xs, merged_xs); out_xe <- c(out_xe, merged_xe)
    data.frame(xs = out_xs, xe = out_xe)
  }

  for (i in seq_along(tfs)) {
    tf  <- tfs[i]
    col <- colors[((i - 1) %% length(colors)) + 1]
    sp  <- spatial_df[spatial_df$protein == tf, ]

    if (nrow(sp) > 0) {
      sigma   <- max(sp$spread[1], 60)
      density <- exp(-0.5 * ((x - sp$centroid[1]) / sigma)^2)
      rib_df  <- data.frame(x = x,
                             y_lo = i + ribbon_base,
                             y_hi = i + ribbon_base + density * ribbon_h)
      # Centroid tick spans only the ribbon height (not into the motif lane)
      cent_df <- data.frame(xc = sp$centroid[1],
                             y0 = i + ribbon_base,
                             y1 = i + ribbon_base + ribbon_h)
      p <- p +
        geom_ribbon(data = rib_df,
                    aes(x = x, ymin = y_lo, ymax = y_hi),
                    fill = col, alpha = 0.22) +
        geom_segment(data = cent_df,
                     aes(x = xc, xend = xc, y = y0, yend = y1),
                     color = col, linewidth = 0.7, alpha = 0.7)
    }

    if (!is.null(motif_results) && tf %in% names(motif_results)) {
      hits <- motif_results[[tf]]$hits
      hits <- hits[!is.na(hits) & hits >= -upstream & hits <= downstream]
      if (length(hits) > 0) {
        # Faint sub-lane strip to frame the motif ticks
        p <- p + annotate("rect",
                          xmin = -upstream, xmax = downstream,
                          ymin = i + motif_y0, ymax = i + motif_y1,
                          fill = col, alpha = 0.06, color = NA)
        # Bake lane position into the data so ggplot doesn't late-bind `i`
        tick_df <- data.frame(x  = hits,
                              y0 = i + motif_y0,
                              y1 = i + motif_y1,
                              yc = i + (motif_y0 + motif_y1) / 2)
        p <- p +
          geom_segment(data = tick_df,
                       aes(x = x, xend = x, y = y0, yend = y1),
                       color = col, linewidth = 0.7, alpha = 0.9) +
          geom_point(data = tick_df,
                     aes(x = x, y = yc),
                     color = col, size = 1.1, alpha = 0.9)
      }
    }

    # ChIP-Atlas union-peak strip (compact — lanes are space-constrained here;
    # the full stacked-SRX view lives in Plot 10). Renders only when we were
    # handed a per-TF peak table; unions overlapping intervals across all
    # SRXs so the strip shows "where public ChIP-seq ever called a peak for
    # this TF in the promoter window".
    if (!is.null(chipatlas_peaks) && tf %in% names(chipatlas_peaks) &&
        !is.null(chipatlas_peaks[[tf]]) && nrow(chipatlas_peaks[[tf]]) > 0) {
      ivs <- .union_intervals(chipatlas_peaks[[tf]])
      # Clip to window
      ivs$xs <- pmax(ivs$xs, -upstream)
      ivs$xe <- pmin(ivs$xe,  downstream)
      ivs <- ivs[ivs$xe > ivs$xs, , drop = FALSE]
      if (nrow(ivs) > 0) {
        ivs$ymin <- i + chip_y0; ivs$ymax <- i + chip_y1
        p <- p + geom_rect(data = ivs,
                           aes(xmin = xs, xmax = xe,
                               ymin = ymin, ymax = ymax),
                           fill = "grey20", color = NA, alpha = 0.75)
      }
    }

    p <- p + annotate("text",
                      x = -upstream - 30, y = i,
                      label = lane_labels[tf], hjust = 1, vjust = 0.5,
                      size = 3.2, color = col, fontface = "bold")
  }

  valid_pos <- pos_map[!is.na(pos_map)]
  for (r in names(valid_pos)) {
    p <- p + annotate("text",
                      x = valid_pos[[r]], y = n_tfs + 0.55,
                      label = r, size = 2.6, color = COLS$guide,
                      angle = 90, vjust = 0.5, hjust = 0)
  }

  p +
    annotate("text", x = 10, y = n_tfs + 0.55,
             label = "TSS", color = COLS$tss, size = 3, hjust = 0) +
    scale_x_continuous(limits    = c(-upstream - 200, downstream),
                       expand    = c(0, 0),
                       labels    = comma) +
    scale_y_continuous(limits    = c(0, n_tfs + 1),
                       expand    = c(0, 0)) +
    labs(x = "Position (bp, relative to TSS)", y = NULL,
         title = title, subtitle = subtitle) +
    theme_caspex(base_size = 11) +
    theme(axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          plot.margin = margin(5, 10, 5, upstream / 10 + 30))
}

#' Plot 6 — Single TF inspector: enrichment bars by region position
#' @param tf_name (see function body).
#' @param long_data (see function body).
#' @param pos_map (see function body).
#' @param pval_thresh (see function body).
#' @export
plot_tf_profile <- function(tf_name, long_data, pos_map, pval_thresh = 0.05) {
  df <- long_data[long_data$protein == tf_name &
                    long_data$region %in% names(pos_map), ]
  df$pos <- as.numeric(pos_map[df$region])
  df <- df[order(df$pos), ]
  df$sig <- df$pval <= pval_thresh

  # Bar width = 80% of minimum spacing between guides
  positions <- sort(unique(df$pos[!is.na(df$pos)]))
  bar_w <- if (length(positions) > 1) min(diff(positions)) * 0.7 else 120

  ggplot(df, aes(x = pos, y = lfc, fill = sig)) +
    geom_col(width = bar_w, alpha = 0.85) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey40") +
    geom_text(aes(label = region,
                  y     = ifelse(lfc >= 0, lfc + 0.04, lfc - 0.04),
                  vjust = ifelse(lfc >= 0, 0, 1)),
              size = 3.2, fontface = "bold") +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = COLS$tss, linewidth = 0.7) +
    scale_fill_manual(values = c(`TRUE`  = COLS$high,
                                 `FALSE` = COLS$guide),
                      labels = c(`TRUE`  = paste0("p \u2264 ", pval_thresh),
                                 `FALSE` = paste0("p > ", pval_thresh)),
                      name = NULL) +
    labs(x     = "Position (bp, TSS-relative)",
         y     = "logFC",
         title = paste0(tf_name, " \u2014 enrichment by region")) +
    theme_caspex() +
    theme(legend.position = "bottom")
}

# =============================================================================
# SECTION 7: Main wrapper
# =============================================================================

#' Select TFs to scan for motifs: top-N by composite + top-N per-region specific
#'
#' "Common" TFs  : highest composite enrichment across all regions.
#' "Specific" TFs: for each region, the TFs whose logFC in that region most
#'                 exceeds their mean logFC across the other regions (region-
#'                 focal signal), restricted to pval <= pval_thresh and
#'                 positive logFC. Up to `top_specific` per region.
#'                 TFs already in `common_tfs` are removed from each region's
#'                 specific list so the two decks don't duplicate each other
#'                 (e.g. a strong global hit like GATA6 stays in common only).
#'
#' @param long_data     From load_region_data()
#' @param spatial_df    From run_spatial_model()  (sorted by composite desc)
#' @param pos_map       From match_all_grnas()
#' @param top_common    Number of global top-composite TFs (default 10)
#' @param top_specific  Number of region-specific TFs per region (default 10).
#'                      Was 5 historically; bumped to 10 so single-region
#'                      focal hits that sit a little lower in the per-region
#'                      specificity ranking still get surfaced on decks
#'                      08/09 instead of being cut off at the top five.
#' @param pval_thresh   P-value cutoff for region-specific selection
#' @return Named character vector of TFs (deduplicated), with an attr
#'         "source" listing common/region membership per TF.
#' @param top_shared (see function body).
#' @noRd
select_motif_tfs <- function(long_data, spatial_df, pos_map,
                             top_common   = 10,
                             top_shared   = 10,
                             top_specific = 10,
                             pval_thresh  = 0.05) {
  # ── Three-deck partition ────────────────────────────────────────────────
  # Every TF in the motif-scan set lands in exactly one of:
  #   1. common        — top-N by composite score (broadly strong);
  #                      rendered in decks 06/07.
  #   2. shared-focal  — significant (p ≤ thr & lfc > 0) in ≥2 regions, not
  #                      in common. Each TF is labelled with its full set
  #                      of significant regions (e.g. "R1+R5");
  #                      rendered in decks 11/12.
  #   3. region-spec   — significant in exactly 1 region, not in common;
  #                      rendered in decks 08/09, one page per region.
  # The three buckets are pairwise disjoint by construction and enforced
  # by post-condition assertions at the end of this function.
  # ────────────────────────────────────────────────────────────────────────

  # (1) Common: top N by composite
  common_tfs <- head(as.character(spatial_df$protein), top_common)

  # Restrict long_data to modelled TFs and regions with a matched gRNA
  modelled <- as.character(spatial_df$protein)
  ld <- long_data[long_data$protein %in% modelled &
                    long_data$region %in% names(pos_map), ]

  # Build the full candidate table: one row per (significant TF, region).
  # spec(p, R) = lfc_p,R − mean_{R'≠R}(lfc_p,R') captures focality.
  regions <- sort(unique(ld$region))
  cand_rows <- list()
  for (r in regions) {
    in_r  <- ld[ld$region == r, ]
    other <- ld[ld$region != r, ]
    if (nrow(in_r) == 0) next

    mean_other <- tapply(other$lfc, other$protein, mean, na.rm = TRUE)
    in_r$mean_other <- ifelse(in_r$protein %in% names(mean_other),
                              mean_other[in_r$protein], 0)
    in_r$spec <- in_r$lfc - in_r$mean_other

    in_r <- in_r[!is.na(in_r$pval) & in_r$pval <= pval_thresh &
                   !is.na(in_r$lfc) & in_r$lfc > 0, ]
    if (nrow(in_r) == 0) next
    in_r <- in_r[!(in_r$protein %in% common_tfs), ]
    if (nrow(in_r) == 0) next

    cand_rows[[r]] <- data.frame(tf     = as.character(in_r$protein),
                                  region = r,
                                  spec   = as.numeric(in_r$spec),
                                  stringsAsFactors = FALSE)
  }
  cand <- if (length(cand_rows) > 0) do.call(rbind, cand_rows)
          else data.frame(tf = character(), region = character(),
                          spec = numeric(), stringsAsFactors = FALSE)

  # Classify each TF by the number of regions in which it is significant.
  # This partitions the non-common candidate set into:
  #   shared_candidates   (n_sig_regions ≥ 2)
  #   specific_candidates (n_sig_regions == 1)
  tf_counts <- if (nrow(cand) > 0) table(cand$tf) else integer(0)
  shared_tfs_all    <- names(tf_counts)[tf_counts >= 2]
  specific_tfs_all  <- names(tf_counts)[tf_counts == 1]

  # (2) Shared-focal: one row per TF, ranked by total specificity across
  # its significant regions. We keep the top N.
  if (length(shared_tfs_all) > 0) {
    sh <- cand[cand$tf %in% shared_tfs_all, , drop = FALSE]
    tot <- aggregate(spec ~ tf, data = sh, FUN = sum)
    reg <- aggregate(region ~ tf, data = sh,
                     FUN = function(x) paste(sort(unique(x)), collapse = "+"))
    shared_summary <- merge(tot, reg, by = "tf")
    shared_summary <- shared_summary[order(shared_summary$spec,
                                            decreasing = TRUE), ]
    shared_summary <- head(shared_summary, top_shared)
    shared_tfs <- as.character(shared_summary$tf)
    shared_tag <- setNames(as.character(shared_summary$region), shared_tfs)
  } else {
    shared_summary <- data.frame(tf = character(0), spec = numeric(0),
                                  region = character(0),
                                  stringsAsFactors = FALSE)
    shared_tfs <- character(0)
    shared_tag <- setNames(character(0), character(0))
  }

  # (3) Region-specific: single-region TFs only (so no TF ever appears on
  # two pages of the region-specific deck). Top-N per region by spec.
  specific_list <- list()
  if (length(specific_tfs_all) > 0) {
    sp <- cand[cand$tf %in% specific_tfs_all, , drop = FALSE]
    for (r in regions) {
      sub <- sp[sp$region == r, , drop = FALSE]
      if (nrow(sub) == 0) next
      sub <- sub[order(sub$spec, decreasing = TRUE), ]
      specific_list[[r]] <- head(as.character(sub$tf), top_specific)
    }
  }
  specific_tfs <- unique(unlist(specific_list, use.names = FALSE))

  # Union (preserve this order so colouring is stable)
  all_tfs <- unique(c(common_tfs, shared_tfs, specific_tfs))

  # Attach provenance for logging / plot labels
  src <- lapply(all_tfs, function(tf) {
    tag <- c()
    if (tf %in% common_tfs) tag <- c(tag, "common")
    if (tf %in% shared_tfs) tag <- c(tag, paste0("shared:", shared_tag[tf]))
    for (r in names(specific_list))
      if (tf %in% specific_list[[r]]) tag <- c(tag, paste0("R:", r))
    paste(tag, collapse = ",")
  })
  attr(all_tfs, "source")   <- setNames(unlist(src), all_tfs)
  attr(all_tfs, "common")   <- common_tfs
  attr(all_tfs, "shared")   <- shared_tfs           # character vector
  attr(all_tfs, "shared_tag") <- shared_tag         # tf -> "R1+R5"
  attr(all_tfs, "specific") <- specific_list        # region -> TFs

  # Per-TF region tag for lane labels. Shared TFs carry their full region
  # set (e.g. "R1+R5"); region-specific TFs carry their single region.
  region_tag <- setNames(character(length(all_tfs)), all_tfs)
  for (tf in all_tfs) {
    if (tf %in% shared_tfs) {
      region_tag[tf] <- shared_tag[tf]
    } else {
      rs <- names(specific_list)[vapply(specific_list,
                                        function(v) tf %in% v, logical(1))]
      region_tag[tf] <- if (length(rs) > 0) paste(rs, collapse = ",") else ""
    }
  }
  attr(all_tfs, "region_tag") <- region_tag

  # ── Post-conditions: three disjoint buckets, single-region specific ───
  .check_disjoint <- function(a, b, a_name, b_name) {
    x <- intersect(a, b)
    if (length(x) > 0)
      warning("select_motif_tfs(): ", length(x),
              " TF(s) appear in BOTH ", a_name, " and ", b_name, ": ",
              paste(x, collapse = ", "),
              "\n  This should be impossible. Please file a bug.",
              call. = FALSE)
  }
  .check_disjoint(common_tfs, shared_tfs,   "common", "shared")
  .check_disjoint(common_tfs, specific_tfs, "common", "region-specific")
  .check_disjoint(shared_tfs, specific_tfs, "shared", "region-specific")

  .flat <- unlist(specific_list, use.names = FALSE)
  .dups <- unique(.flat[duplicated(.flat)])
  if (length(.dups) > 0) {
    warning("select_motif_tfs(): ", length(.dups),
            " TF(s) appear in MULTIPLE region-specific lists: ",
            paste(.dups, collapse = ", "),
            "\n  This should be impossible. Please file a bug.",
            call. = FALSE)
  }

  message("  Motif TF selection: ", length(common_tfs), " common + ",
          length(shared_tfs), " shared-focal + ",
          length(specific_tfs), " region-specific (",
          length(all_tfs), " unique total)")
  # Human-readable breakdown so the three decks' disjointness is auditable
  # from the console.
  message("    Common (deck 06/07): ",
          paste(common_tfs, collapse = ", "))
  for (r in names(specific_list)) {
    message("    Region-specific (deck 08/09) ", r, ": ",
            paste(specific_list[[r]], collapse = ", "))
  }
  if (length(shared_tfs) > 0) {
    message("    Shared-focal (deck 11/12):")
    for (tf in shared_tfs)
      message("      ", tf, "  [", shared_tag[tf], "]")
  } else {
    message("    Shared-focal (deck 11/12): (none)")
  }

  all_tfs
}

#' Run the full GLproxScape analysis pipeline
#'
#' End-to-end driver for spatial deconvolution of dCas9-APEX2 proximity
#' proteomics. Takes a manifest of guide RNAs and per-region enrichment
#' tables, anchors them to an Ensembl transcript's TSS, builds the
#' Gaussian labelling kernel + coverage map, runs the spatial model and
#' coverage-aware deconvolution against JASPAR PWMs, optionally overlays
#' independent ChIP-Atlas peaks, and writes the result deck (motif-track
#' PDF, deconvolution detail pages, CSVs) to \code{out_dir}.
#'
#'
#' @param gene HGNC gene symbol whose promoter window is analysed.
#' @param grnas Named character vector of protospacer sequences (17-23 bp,
#'   PAM optional). Names match the region IDs used by \code{data_files}.
#'   Set an element to \code{NA} for regions without a sequence to match
#'   (the region still contributes its logFC data but no cut-site tick).
#' @param data_files Named character vector of per-region proteomics
#'   table paths. Names align with \code{grnas}. Each file must carry
#'   a \code{logFC} column and a p-value column.
#' @param out_dir Output directory; created if absent. The deconvolution
#'   PDF, motif-track PDF, gRNA-positions PDF, and predictions CSVs all
#'   land here.
#'
#'
#' @param tf_universe Optional character vector of HGNC TF symbols (e.g.
#'   read from \code{TFLibrary.txt}). When supplied, sets the \code{isTF}
#'   flag on the loaded long_data and gates the spatial model when
#'   \code{tfs_only = TRUE}. NULL falls back to a \code{TFDatabase}-style
#'   column in the input files if present, else FALSE everywhere.
#' @param epi_universe Optional character vector of epigenetic-factor
#'   HGNC symbols (e.g. read from \code{EpiGenes_main.csv}). Sets
#'   \code{isEpi} on loaded rows; used downstream by
#'   \code{\link{run_caspex_epigenetic}} but does not affect the TF deck.
#' @param tfs_only If TRUE (default) the spatial model is restricted to
#'   rows with \code{isTF = TRUE}; FALSE includes every protein.
#'
#'
#' @param species Ensembl species token (default \code{"homo_sapiens"}).
#' @param transcript TSS-anchor selection. \code{"canonical"} (default)
#'   uses the Ensembl canonical transcript - deterministic across REST
#'   calls. \code{"ENST..."} pins to a specific transcript ID (e.g.
#'   Myers MYC P2: \code{"ENST00000377970"}). \code{NA} falls back to
#'   the legacy gene-level union boundary (subject to Ensembl annotation
#'   drift; not recommended).
#' @param upstream,downstream bp window around the TSS to fetch and
#'   model. Defaults 2500 / 500 - widen \code{downstream} when guides
#'   tile downstream of the canonical TSS (Mackenzie FOXP2 uses 200 / 2000).
#'
#'
#' @param pval_thresh Per-region p-value floor for a TF to count as
#'   significant in that region (default 0.05).
#' @param min_regions Minimum number of regions in which a TF must pass
#'   the p-value gate (default 2).
#' @param min_lfc Optional logFC floor applied before weighting in the
#'   spatial model. Default 0 (no floor). Set positive to ignore mildly
#'   negative logFC entirely rather than relying on \code{compute_region_weight}'s
#'   positive-weight clip.
#' @param top_n Number of top TFs to render on the spatial track plot
#'   (default 25).
#' @param motif_tfs Optional character vector of TF symbols to run the
#'   JASPAR motif scan on. NULL = the engine's own selection
#'   (\code{select_motif_tfs}: common + shared-focal + region-specific
#'   buckets, total ~44 TFs). Override if you want a custom set.
#' @param n_common,n_shared,n_specific Per-bucket caps for
#'   \code{select_motif_tfs}: top-N common across regions, top-N shared-
#'   focal (significant in >=2 regions), top-N per region-specific list.
#'   Defaults 20 / 20 / 20.
#'
#'
#' @param motif_thresh JASPAR PWM log-odds threshold as a fraction of the
#'   matrix max score (default 0.80). 0.75 is a common relaxation for
#'   long high-IC matrices; lower values surface more peripheral hits.
#' @param motif_scan_pool One of \code{"selected"} (default) or
#'   \code{"spatial_all"}. \code{"selected"} scans only the 44-TF
#'   select_motif_tfs union, matching the legacy pipeline. \code{"spatial_all"}
#'   additionally scans every TF in the spatial model that didn't make
#'   that cut; the extra results live in \code{result$motif_results_extra}
#'   and are consumed by the deconvolution deck when
#'   \code{deconv_min_motif_hits > 0}.
#' @param motif_score_weight How the PWM log-odds score per motif feeds
#'   into bubble amplitude: \code{"none"} (default; binary above-threshold
#'   filter only), \code{"linear"} (amplitude scaled by score_frac in
#'   [0,1]), \code{"log"} (amplitude scaled by 2^(score_frac - 1) -
#'   compresses the dynamic range). See the supplemental methods for the
#'   exact NNLS / zone-path math.
#'
#'
#' @param kernel_sigma Gaussian labelling kernel width in bp (default 300).
#'   Roughly matches APEX2's biotinylation radius along linear DNA.
#' @param min_weight_frac Minimum fraction of the local peak amplitude
#'   below which events are pruned (default 0.15).
#' @param min_peak_dist Minimum bp separation between local maxima in
#'   the no-motif fallback peak detector (default 150).
#' @param merge_dist Closely-spaced motif hits within this many bp are
#'   merged into a single amplitude-weighted cluster (default 100 ~ sigma/3).

#' @param cov_floor Relative floor on the labelling-coverage denominator:
#'   \eqn{\beta(x) = s(x) / \max(C(x), cov\_floor \cdot \max(C))}.
#'   Effective amplification is capped near \code{1 / cov_floor}.
#'   Default 0.05 (~20x cap). Lower amplifies distal signals more
#'   aggressively and also amplifies noise.
#' @param edge_guard_frac Fraction-of-max-coverage floor that defines
#'   the in-support region for beta. Set well above \code{cov_floor} so
#'   beta is not evaluated in the low-denominator transition band at the
#'   tile edges. Default 0.25 (5x default cov_floor).
#' @param zone_peak_frac Per-zone beta floor: motifs whose local beta is
#'   below \code{zone_peak_frac * zone_peak_beta} are dropped. 0 disables.
#'   Default 0.50.
#' @param max_events_per_tf Top-N cap on events per TF after merging
#'   (default 30; Inf disables).
#' @param merge_position How the position of a merged cluster is reported.
#'   \code{"argmax"} (default) snaps to the strongest motif's bp coord.
#'   \code{"centroid"} uses the amplitude-weighted mean - can land
#'   between motifs.
#' @param max_grna_distance Hard geometric cap in bp on how far a called
#'   event can sit from the nearest gRNA. NULL (default) resolves to
#'   \code{kernel_sigma}. Inf disables. Auto-scales with the kernel
#'   width; complements \code{edge_guard_frac}.
#' @param edge_grna_weight_cap Optional fraction in (0,1]; drops events
#'   where either boundary gRNA contributes more than this fraction of
#'   the local Gaussian weight sum. Suppresses "edge-bleed" events
#'   inflated by a single boundary guide's kernel tail. NULL (default)
#'   disables; 0.5 is a typical setting when needed.
#'
#'
#' @param position_stability \code{"none"} (default) or
#'   \code{"wild_bootstrap"}. When enabled, runs a Rademacher Wild
#'   bootstrap on the residuals of a forward-model NNLS fit at the called
#'   event positions and adds four columns to \code{binding_events}
#'   (pos_stab_low/high/median, bootstrap_survival_frac). NxSlower; only
#'   currently used as a benchmarking diagnostic.
#' @param n_bootstrap Wild bootstrap draw count when
#'   \code{position_stability = "wild_bootstrap"} (default 200L).
#'
#'
#' @param weight_mode Region-weight transform applied to the per-region
#'   enrichment when building the smoothed signal s(x). \code{"z"}
#'   (default) = signed z-score derived from p-value. Alternatives:
#'   \code{"mod_t"} (moderated t from a limma-style \code{t} column),
#'   \code{"lfc_pos"}, \code{"lfc_signed"}, \code{"lfc_x_negp"}.
#' @param signal_weight Backward-compat alias for \code{weight_mode}.
#'   If non-NULL, overrides \code{weight_mode} for signal building.
#'
#'
#' @param chipatlas FALSE (default; off so runs stay network-free) or
#'   TRUE to fetch and render public ChIP-seq peaks for every
#'   motif-scanned TF.
#' @param chipatlas_threshold One of \code{"05"} (Q<1e-5, default),
#'   \code{"10"}, or \code{"20"}.
#' @param chipatlas_max_experiments Cap on SRX experiments per TF
#'   (default 100, newest first).
#' @param special_interest_gene Optional character vector of TF symbols
#'   that bypass the per-TF SRX cap and have ALL (or
#'   \code{special_interest_cap}) of their available SRX scanned.
#'   Useful for biologically focal regulators where rare cell-type-
#'   specific peaks live in older SRX entries the cap excludes.
#' @param special_interest_cap Optional integer cap on the special-
#'   interest SRX count per TF. NULL (default) = scan all SRX.
#' @param chipatlas_quiet Suppress per-SRX download messages (default TRUE).
#'
#'
#' @param detail_top_n How many top TFs to render on the per-TF
#'   deconvolution detail PDF (default 100).
#' @param deconv_min_motif_hits Minimum number of JASPAR hits in the
#'   promoter window required for a TF to appear on the detail PDF.
#'   Default 0 (no filter); set > 0 to focus the deck on TFs with motif-
#'   grounded calls and exclude motif-less peak-detection fallbacks.
#' @param deconv_min_max_lfc Optional logFC floor on top of the motif-hit
#'   filter: a TF is kept iff its max per-region logFC across matched
#'   regions is >= this value. Default 0 (no filter). Useful when input
#'   p-values are unreliable; the two filters compose AND-style.
#'
#'
#' @param save_plots Write the PDF deck to \code{out_dir} (default TRUE).
#'   Set FALSE for a fast headless run that only computes the result list.
#' @param plot_width,plot_height PDF dimensions in inches. Defaults 10 x 8.
#'
#' @return Invisibly, a list with: spatial_df, long_data, pos_map,
#'   gene_info, promoter_info, motif_results, motif_results_extra,
#'   binding_events, plots, plus run metadata (weight_mode, cov_floor,
#'   edge_guard_frac, kernel_sigma, upstream, downstream,
#'   chipatlas_peaks, chipatlas_threshold).
#' @export
run_caspex <- function(
    gene,
    grnas,
    data_files,
    # Optional gene-class universes used to flag rows in the loaded
    # long_data:
    #   tf_universe  — e.g. read from TFLibrary.txt; sets `isTF`
    #   epi_universe — e.g. read from EpiGenes_main.csv; sets `isEpi`
    # When NULL, falls back to a `TFDatabase`-style column in the input
    # files if present, else FALSE for everything.
    tf_universe   = NULL,
    epi_universe  = NULL,
    out_dir       = "caspex_output",
    species       = "homo_sapiens",
    # Transcript anchor for the TSS reference frame. As of the canonical-
    # default change, runs are deterministic across sessions: the Ensembl
    # canonical transcript is a stable per-release anchor, whereas the
    # legacy gene-level start/end is the union over all annotated
    # transcripts and can drift between REST calls (uniform TSS-frame
    # shifts on the order of tens of bp).
    #   "canonical" (default): Ensembl canonical transcript — deterministic
    #                          and matches the promoter frame used in most
    #                          TF / ChIP studies.
    #   "ENST..."            : pin to a specific Ensembl transcript ID
    #                          (e.g. Myers 2018 MYC uses ENST00000377970
    #                          for the P2 promoter — overrides "canonical").
    #   NA                   : explicit legacy gene-level start/end. Subject
    #                          to Ensembl annotation drift; only use this
    #                          if you specifically need the union boundary.
    transcript    = "canonical",
    upstream      = 2500,
    downstream    = 500,
    pval_thresh   = 0.05,
    min_regions   = 2,
    # Floor applied to logFC before weighting in the spatial model. Default 0
    # = legacy behaviour (no floor). Set > 0 if you want the spatial centroid
    # to ignore mildly negative logFC values entirely rather than letting the
    # positive-weight clip in compute_region_weight handle them.
    min_lfc       = 0,
    top_n         = 25,
    motif_tfs     = NULL,
    n_common      = 20,
    n_shared      = 20,
    n_specific    = 20,
    motif_thresh  = 0.80,
    # `motif_score_weight`: how the JASPAR PWM score per motif feeds into
    # bubble size in the deconvolution. The score still acts as the binary
    # above-threshold filter regardless of this setting; this parameter
    # controls whether it ALSO weights the per-motif amplitude.
    #   "none"   - default; legacy behaviour. Above-threshold motifs are
    #              treated equally (unit-amplitude basis on the smoothed-
    #              s(x) NNLS path; raw β on the coverage-aware zone path).
    #   "linear" - per-motif amplitude is scaled by score_frac in [0,1],
    #              i.e. higher PWM score gets more amplitude per unit β
    #              (NNLS path) or per unit β-readout (coverage-aware).
    #   "log"    - amplitude scaled by 2^(score_frac - 1), monotonic but
    #              compresses the difference (0.5x at threshold, 1.0x at
    #              consensus).
    motif_score_weight = c("none", "linear", "log"),
    tfs_only      = TRUE,
    # Binding deconvolution parameters
    kernel_sigma     = 300,
    min_weight_frac  = 0.15,
    min_peak_dist    = 150,
    merge_dist       = 100,
    cov_floor        = 0.05,
    # `edge_guard_frac` (coverage-aware only): fraction-of-max-coverage floor
    # that defines the in-support region for β. Set well above `cov_floor`
    # so β is not evaluated in the low-denominator transition band at the
    # tiled-region edges (which would otherwise spawn spurious west-edge
    # peaks / zones). Default 0.25 = 5× default cov_floor (raised from
    # 0.15 when kernel_sigma went 250 → 300; wider kernel extends the
    # single-gRNA tail and re-opened the edge artifact at 0.15).
    edge_guard_frac  = 0.25,
    # Readability guards. `zone_peak_frac` = per-zone β floor
    # (motifs below frac × zone_peak β are dropped; 0 disables).
    # `max_events_per_tf` = top-N cap applied after merging (Inf disables).
    zone_peak_frac    = 0.50,
    max_events_per_tf = 30,
    # `merge_position` (coverage-aware only): "argmax" (default) reports
    # each merged cluster at the position of its strongest motif so
    # bubbles snap to a real motif tick. "centroid" keeps the legacy
    # amplitude-weighted mean, which can land between motifs.
    merge_position    = c("argmax", "centroid"),
    # `max_grna_distance` (coverage-aware only): hard geometric cap in bp
    # on how far a called event can sit from the nearest gRNA. NULL
    # (default) resolves to `kernel_sigma` at runtime — i.e. events must
    # be within one labeling σ of some guide. Inf disables. Complements
    # `edge_guard_frac`: the relative-coverage mask does not auto-rescale
    # with kernel_sigma, this one does.
    max_grna_distance = NULL,
    # `edge_grna_weight_cap` (coverage-aware only): fraction in (0, 1], or
    # NULL to disable (default). Drops events where either boundary gRNA
    # (westernmost `min(pos_r)` or easternmost `max(pos_r)`) contributes
    # more than this fraction of the local Gaussian weight sum
    # Σ_i exp(-0.5·((p - r_i)/σ)²). Suppresses "edge-bleed" events in the
    # tail beyond the outermost guides where a single boundary guide can
    # inflate the NNLS β at interior positions far from its own coordinate
    # (seen on ATP7B at p=-1620 where R7 @ -1945 contributes ~37%).
    # Complements `max_grna_distance` (which is a pure geometric cap on
    # position→nearest-guide): this one caps Gaussian weight-share.
    edge_grna_weight_cap = NULL,
    # `position_stability` (coverage-aware only): per-event diagnostic that
    # measures how robust each bubble's POSITION is to noise perturbation
    # of the per-region weights — NOT a confidence interval on the TF's
    # true binding location. See `predict_binding_events_coverage_aware()`
    # for the full framing. "wild_bootstrap" runs a Rademacher Wild
    # bootstrap on the residuals of a forward-model NNLS fit at the called
    # event positions; "none" (default) skips the diagnostic. When enabled,
    # binding_events gains four columns (pos_stab_low / pos_stab_high /
    # pos_stab_median / bootstrap_survival_frac) and the per-TF Plot 10
    # detail panel renders a "Position support" sub-track between the event
    # bubbles and the ChIP-Atlas band.
    #
    # Off by default because the bootstrap is N×slower and is currently
    # used as a benchmarking diagnostic on the Myers 2018 reanalysis only.
    position_stability   = c("none", "wild_bootstrap"),
    # Number of Wild bootstrap draws when `position_stability` is enabled.
    n_bootstrap          = 200L,
    # Region-weight mode — applied to BOTH spatial model and signal building.
    # "z" (default) = signed z-score derived from the p-value, ALWAYS, even
    # when the input files carry a `t` column. This matches the ATP7B-era
    # CasPEX behaviour (where no t-stat was available) and keeps the
    # smoothed-signal magnitude on a bounded scale (|z| ≲ 5) so build-
    # _caspex_signal()'s curve does not get visually dominated by a few
    # high-leverage replicates with small SE (moderated |t| can reach 15+).
    # Alternatives: "mod_t" (opt-in moderated-t from limma; falls back to
    # signed z when no `t` column is present), "lfc_pos", "lfc_signed",
    # "lfc_x_negp" (legacy).
    weight_mode      = "z",
    # Backward-compat alias; if set, overrides weight_mode for signal building
    signal_weight    = NULL,
    # ChIP-Atlas overlay (public ChIP-seq peaks as a supplementary track):
    #   chipatlas               - FALSE (default, off so runs stay network-free)
    #                             / TRUE to fetch + render peaks for every
    #                             motif-scanned TF.
    #   chipatlas_threshold     - "05" (Q<1e-5, default), "10", or "20".
    #   chipatlas_max_experiments - cap on SRXs per TF (default 100; newest first).
    #   chipatlas_quiet         - suppress per-SRX download messages.
    #   special_interest_gene   - character vector of TF symbols (e.g.
    #                             c("EP300","CTCF")) that should bypass the
    #                             chipatlas_max_experiments cap and have ALL
    #                             (or `special_interest_cap`) of their
    #                             available SRX scanned. Other TFs continue
    #                             to use the standard cap. Useful for
    #                             biologically focal regulators where a rare
    #                             cell-type-specific peak might live in the
    #                             older SRX entries the cap excludes.
    #   special_interest_cap    - optional integer cap on the special-interest
    #                             SRX count per TF. NULL (default) = scan ALL
    #                             SRX (e.g. all ~2200 CTCF SRXs). Set e.g.
    #                             250 to take a top-N slice that covers most
    #                             cell-type diversity without a multi-thousand
    #                             BED download.
    chipatlas               = FALSE,
    chipatlas_threshold     = "05",
    chipatlas_max_experiments = 100,
    special_interest_gene   = NULL,
    special_interest_cap    = NULL,
    chipatlas_quiet         = TRUE,
    detail_top_n     = 100,
    # `deconv_min_motif_hits` filters which TFs land on Plot 10
    # (binding_deconvolution.pdf). Default 0 = keep legacy behaviour of
    # rendering one page per motif-scanned TF (the common + region-specific
    # union from select_motif_tfs). Setting > 0 keeps only TFs whose JASPAR
    # PWM scan in the promoter window returned at least that many hits, which
    # focuses the deck on TFs whose binding case is grounded in motif evidence
    # rather than enrichment alone. Useful on datasets like Myers 2018 where
    # the predicted TF list is dominated by long high-IC ZNFs that often
    # return a handful (or zero) hits even at motif_thresh = 0.75 — those
    # pages are mostly empty and crowd out the TFs with real motif support.
    # TFs without any JASPAR matrix (absent from result$motif_results) are
    # treated as having 0 hits and are dropped whenever this threshold > 0.
    deconv_min_motif_hits = 0,
    # `deconv_min_max_lfc` adds an optional logFC floor on top of the motif-hit
    # filter for Plot 10. A TF is kept iff its maximum per-region logFC across
    # the matched regions (those in pos_map) is >= this value. Default 0 = no
    # filter (all TFs surviving the motif-hit step are plotted). Useful when
    # the upstream p-values are unreliable or saturated (e.g. Myers 2018, where
    # P.Value is synthesized from a robust z-score and floored at 1e-300, so
    # filtering on p loses its discriminating power) — a per-zone logFC cut
    # gives an interpretable proximity-enrichment criterion that is independent
    # of the synthesized p-value pipeline. Composes AND-style with
    # deconv_min_motif_hits: both must pass.
    deconv_min_max_lfc    = 0,
    # `motif_scan_pool` controls which TFs the JASPAR motif scan covers.
    #   "selected"    (default) — scan only the 44-TF select_motif_tfs union
    #                 (common + shared-focal + region-specific). Matches the
    #                 legacy pipeline; decks 06/07/08/09/11/12 stay narrative-
    #                 driven by spatial enrichment and Plot 10 can only filter
    #                 *within* that pre-narrowed set.
    #   "spatial_all" — *additionally* scan every TF in the spatial model
    #                 that did not make the select_motif_tfs cut (the
    #                 ~spatial_df$protein \ motif_tfs complement, ≈140 TFs
    #                 on a typical Myers run). The extra scan results live
    #                 in result$motif_results_extra and are *only* consumed
    #                 by Plot 10 when deconv_min_motif_hits > 0 — every other
    #                 deck still uses the unchanged motif_res so their layout
    #                 and TF roster do not shift. This lets Plot 10 surface
    #                 TFs with strong motif evidence at the locus that the
    #                 enrichment-driven select_motif_tfs missed (e.g. low-rank
    #                 TFs with a perfect JASPAR consensus in the window).
    # Cost: one extra fetch_jaspar_pwm() + scan_sequence() per non-selected
    # spatial TF; PWMs are cached on disk so subsequent runs are cheap.
    motif_scan_pool       = c("selected", "spatial_all"),
    save_plots    = TRUE,
    plot_width    = 10,
    plot_height   = 8
) {
  if (is.null(signal_weight)) signal_weight <- weight_mode
  position_stability <- match.arg(position_stability)
  motif_score_weight <- match.arg(motif_score_weight)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # ── Resolve TF / Epi universes from bundled databases when not user-supplied
  # Keeps the simple `run_caspex(gene, grnas, data_files)` call self-contained:
  # users who want to pin a custom universe still pass `tf_universe = ...` /
  # `epi_universe = ...` and override the defaults.
  if (is.null(tf_universe)) {
    tf_path <- system.file("extdata/databases/TFLibrary.txt",
                            package = "GLproxScape")
    if (nzchar(tf_path) && file.exists(tf_path)) {
      tf_universe <- toupper(trimws(readLines(tf_path)))
      tf_universe <- tf_universe[nzchar(tf_universe)]
      message("Using bundled TF universe (", length(tf_universe),
              " entries) from ", basename(tf_path))
    }
  }
  if (is.null(epi_universe)) {
    epi_path <- system.file("extdata/databases/EpiGenes_main.csv",
                             package = "GLproxScape")
    if (nzchar(epi_path) && file.exists(epi_path)) {
      epi_df <- utils::read.csv(epi_path, stringsAsFactors = FALSE)
      sym_col <- if ("HGNC_symbol" %in% names(epi_df))
        "HGNC_symbol" else names(epi_df)[1]
      epi_universe <- toupper(trimws(epi_df[[sym_col]]))
      epi_universe <- epi_universe[nzchar(epi_universe) & !is.na(epi_universe)]
      message("Using bundled Epi universe (", length(epi_universe),
              " entries) from ", basename(epi_path))
    }
  }

  message("=== CasPEX Binding Zone Predictor ===")
  message("Gene: ", gene, "  |  Regions: ", length(data_files))

  # ── 1. Load enrichment data ───────────────────────────────────────────────
  message("\nLoading enrichment data...")
  long_data <- load_region_data(data_files,
                                 tf_universe  = tf_universe,
                                 epi_universe = epi_universe)
  n_with_t  <- sum(!is.na(long_data$t_stat))
  has_t     <- n_with_t > 0
  # Resolve the human-readable label for the active weight source. With the
  # default `weight_mode = "z"` we always derive a signed z from the p-value,
  # regardless of whether the input carries a `t_stat` column — this is the
  # ATP7B-era approach and keeps the smoothed-signal magnitude bounded
  # (typical |z| ≲ 5 vs. moderated |t| reaching 15+ for tight replicates,
  # which inflates build_caspex_signal()'s curve without changing event
  # ranks). "mod_t" remains opt-in.
  weight_src <- if (weight_mode == "z") {
    if (has_t)
      "signed z-score (derived from p-value; t_stat column present but ignored \u2014 z-mode default)"
    else
      "signed z-score (derived from p-value; no `t` column in input)"
  } else if (weight_mode == "mod_t") {
    if (has_t)
      "moderated t (from input `t` column)"
    else
      "signed z-score (derived from p-value; no `t` column in input)"
  } else {
    weight_mode
  }
  message("  ", nrow(long_data), " rows loaded | ",
          sum(long_data$isTF), " TF rows | ",
          sum(long_data$isEpi %||% FALSE), " Epi rows | ",
          length(unique(long_data$protein)), " unique proteins | ",
          "t_stat populated: ", n_with_t, "/", nrow(long_data))
  message("  Weight source   : ", weight_src)

  # ── 2. Gene + promoter sequence ───────────────────────────────────────────
  gene_info     <- lookup_gene(gene, species, transcript = transcript)
  # Pass species explicitly so the sequence-fetch URL always matches the
  # lookup species, even if someone later calls fetch_promoter_seq() with a
  # hand-built gene_info that's missing the `species` field.
  promoter_info <- fetch_promoter_seq(gene_info, upstream, downstream,
                                       species = species)

  # ── 3. Match gRNAs ────────────────────────────────────────────────────────
  pos_map <- match_all_grnas(grnas, promoter_info)

  n_matched <- sum(!is.na(pos_map))
  if (n_matched < 2)
    stop("At least 2 gRNAs must match the promoter sequence. ",
         "Got ", n_matched, ". Check sequences and gene name.")

  # ── 4. Spatial model ─────────────────────────────────────────────────────
  spatial_df <- run_spatial_model(long_data, pos_map,
                                   tfs_only    = tfs_only,
                                   pval_thresh = pval_thresh,
                                   min_regions = min_regions,
                                   min_lfc     = min_lfc,
                                   weight_mode = weight_mode)
  message("  ", nrow(spatial_df), " TFs in spatial model")

  # ── 5. JASPAR motif scan ──────────────────────────────────────────────────
  if (is.null(motif_tfs))
    motif_tfs <- select_motif_tfs(long_data, spatial_df, pos_map,
                                  top_common   = n_common,
                                  top_shared   = n_shared,
                                  top_specific = n_specific,
                                  pval_thresh  = pval_thresh)

  motif_res <- run_motif_scan(motif_tfs, promoter_info, motif_thresh)

  # ── 5b. Optional augment-pool motif scan ──────────────────────────────────
  # When motif_scan_pool == "spatial_all", scan JASPAR for every TF in the
  # spatial model that did not make it into select_motif_tfs's 44-TF cut.
  # Results are kept in a SEPARATE list (`motif_res_extra`) so the existing
  # decks (06/07, 08/09, 11/12) — which size lanes/PDF heights from
  # length(motif_res) — are unaffected. Plot 10 will merge the two when
  # building its TF list (see deconv_min_motif_hits handling below).
  motif_scan_pool <- match.arg(motif_scan_pool)
  motif_res_extra <- list()
  if (motif_scan_pool == "spatial_all") {
    extra_tfs <- setdiff(as.character(spatial_df$protein),
                         as.character(motif_tfs))
    if (length(extra_tfs) > 0) {
      message("\n  Augment-pool motif scan: ", length(extra_tfs),
              " additional TFs from spatial_df (motif_scan_pool='spatial_all')")
      motif_res_extra <- run_motif_scan(extra_tfs, promoter_info, motif_thresh)
    }
  }

  # ── 5a. ChIP-Atlas peak overlay (optional) ────────────────────────────────
  # Pulls public ChIP-seq peaks for every motif-scanned TF from ChIP-Atlas,
  # windowed to the same promoter region. Used as a supplementary validation
  # track in plots 06/07, 08/09, 11/12 (union peak row per TF) and 10
  # (stacked per-experiment ticks below the bubble strip).
  chipatlas_res <- NULL
  if (isTRUE(chipatlas)) {
    if (!isTRUE(.caspex_chipatlas_loaded)) {
      warning("chipatlas=TRUE but caspex_chipatlas.R could not be sourced; ",
              "skipping ChIP-Atlas overlay.")
    } else {
      chipatlas_res <- run_chipatlas_scan(
        tfs              = motif_tfs,
        gene_info        = gene_info,
        promoter_info    = promoter_info,
        upstream         = upstream,
        downstream       = downstream,
        threshold        = chipatlas_threshold,
        max_experiments  = chipatlas_max_experiments,
        special_interest_gene = special_interest_gene,
        special_interest_cap  = special_interest_cap,
        quiet            = chipatlas_quiet
      )
    }
  }

  # ── 5b. Motif-constrained binding deconvolution ──────────────────────────
  x_grid_bind <- seq(-upstream, downstream, by = 5)
  binding_events <- predict_all_binding_events(
    tfs             = motif_tfs,
    long_data       = long_data,
    pos_map         = pos_map,
    motif_results   = motif_res,
    x_grid          = x_grid_bind,
    kernel_sigma    = kernel_sigma,
    min_weight_frac = min_weight_frac,
    min_peak_dist   = min_peak_dist,
    merge_dist      = merge_dist,
    weight_mode     = signal_weight,
    cov_floor        = cov_floor,
    edge_guard_frac  = edge_guard_frac,
    zone_peak_frac    = zone_peak_frac,
    max_events_per_tf = max_events_per_tf,
    merge_position    = merge_position,
    max_grna_distance = max_grna_distance,
    edge_grna_weight_cap = edge_grna_weight_cap,
    position_stability = position_stability,
    n_bootstrap        = n_bootstrap,
    motif_score_weight = motif_score_weight
  )
  message("  Total predicted events: ", nrow(binding_events),
          "  (across ", length(unique(binding_events$tf)), " TFs)")

  # ── 6. Plots ──────────────────────────────────────────────────────────────
  message("\nGenerating plots...")

  p_grna    <- plot_grna_positions(pos_map, gene_info, upstream, downstream)
  p_heat    <- plot_heatmap(long_data, spatial_df, pos_map, top_n = min(30, nrow(spatial_df)))
  p_gc      <- plot_gc_track(promoter_info)

  # ── Deck-provenance lookup (used for tagging binding_events $deck below) ──
  common_tfs   <- attr(motif_tfs, "common")
  shared_tfs   <- attr(motif_tfs, "shared")
  shared_tag   <- attr(motif_tfs, "shared_tag")
  specific_map <- attr(motif_tfs, "specific")
  region_tag   <- attr(motif_tfs, "region_tag")

  # If user supplied motif_tfs manually, those attrs won't exist; derive them
  if (is.null(common_tfs)) {
    common_tfs   <- head(as.character(spatial_df$protein), n_common)
    shared_tfs   <- character(0)
    shared_tag   <- setNames(character(0), character(0))
    specific_map <- list()
    region_tag   <- setNames(rep("", length(motif_tfs)), motif_tfs)
  }
  if (is.null(shared_tfs)) {
    shared_tfs <- character(0)
    shared_tag <- setNames(character(0), character(0))
  }
  specific_tfs <- setdiff(unique(unlist(specific_map, use.names = FALSE)),
                          character(0))

  # Initialise the binding-deconvolution deck TF roster up front so it's always
  # available for the result list (and downstream consumers like A5)
  # even when save_plots = FALSE. The save_plots block below
  # overwrites this with the actual filtered + ordered list when the
  # deck rendering runs.
  detail_tfs <- character(0)

  if (save_plots) {
    # Plots from earlier development (02_spatial_track, 05_motif_overlay,
    # 06/07_track_common*, 08/09_track_region_specific*, 11/12_track_shared*,
    # and the 00_summary composite) have been retired — they were diagnostic
    # artefacts that didn't pull their weight in the final figure set. Only
    # the load-bearing outputs are written: gRNA layout, per-region heatmap,
    # GC content, the binding-deconvolution deck, and the epigenetic deck
    # (in run_caspex_epigenetic).
    ggsave(file.path(out_dir, "grna_positions.pdf"),
           p_grna, width = plot_width, height = 4)
    ggsave(file.path(out_dir, "heatmap.pdf"),
           p_heat, width = 10, height = max(5, min(30, nrow(spatial_df)) * 0.28))
    ggsave(file.path(out_dir, "gc_content.pdf"),
           p_gc, width = plot_width, height = 3)

    # ── Binding deconvolution detail plots ─
    # Apply to *all* plotted TFs (the full motif_tfs set = common + region-specific),
    # not just the top-N by event weight. TFs with zero events are skipped so the
    # deck stays meaningful. Order: by total predicted event weight desc, then
    # alphabetical for TFs without any events.
    if (length(motif_tfs) > 0) {
      # Merged motif lookup for Plot 10 only. With motif_scan_pool='selected'
      # this is just `motif_res` and the augment branch is empty. With
      # motif_scan_pool='spatial_all' it carries every spatial-modeled TF.
      # Other decks read from `motif_res` directly and are not affected.
      motif_res_for_deck <- c(motif_res, motif_res_extra)

      # Base pool: the 44-TF select_motif_tfs union (legacy behaviour).
      all_plot_tfs <- as.character(motif_tfs)

      # Augment branch — only meaningful when both (a) the user opted into
      # the wider scan via motif_scan_pool='spatial_all' AND (b) a motif-hit
      # threshold was requested. With threshold = 0 the augmented pool would
      # add 100+ pages (every spatial TF that happens to have a PWM, even
      # 1-hit weak matches), defeating the purpose; require threshold > 0
      # for the augmentation to take effect.
      if (motif_scan_pool == "spatial_all" &&
          isTRUE(deconv_min_motif_hits > 0) &&
          length(motif_res_extra) > 0) {
        extra_n_hits <- vapply(motif_res_extra, function(m) {
          n <- m$n_hits
          if (is.null(n)) length(m$hits) else as.integer(n)
        }, integer(1))
        extra_pass <- names(motif_res_extra)[extra_n_hits >=
                                              deconv_min_motif_hits]
        new_extras <- setdiff(extra_pass, all_plot_tfs)
        if (length(new_extras) > 0) {
          message("  Augment-pool: adding ", length(new_extras),
                  " TF(s) with >= ", deconv_min_motif_hits,
                  " JASPAR hits that were not in motif_tfs (",
                  paste(head(new_extras, 8), collapse = ", "),
                  if (length(new_extras) > 8)
                    paste0(", ... +", length(new_extras) - 8, " more") else "",
                  ")")
        }
        all_plot_tfs <- unique(c(all_plot_tfs, extra_pass))
      }

      # Optional motif-hit filter — keep TFs whose JASPAR PWM scan returned at
      # least `deconv_min_motif_hits` hits in the promoter window. Default 0
      # is a no-op.
      #
      # IMPORTANT semantics: a TF that fetch_jaspar_pwm() couldn't resolve
      # (no JASPAR matrix in either the user-given alias or the species/
      # paralogue lookup) is *exempt* from this filter and kept. Such TFs
      # are common among CasPEX hits on Myers 2018 MYC — chromatin-associated
      # proteins (NONO, MTA1, SSRP1, GATAD2A/B, etc.) that don't have direct
      # sequence-specific DBDs and therefore aren't in JASPAR at all. We can
      # only fairly penalise a TF for "weak motif evidence" if it has a PWM
      # in the first place; without one, it falls through to the peak-
      # detection branch (already wired in build_caspex_signal()) and Plot 10
      # still shows its spatial deconvolution + ChIP-Atlas peak track. This
      # asymmetry only affects the original motif_tfs pool — the augment-pool
      # branch above can only add TFs that came through run_motif_scan
      # successfully (i.e., already have a PWM with >= threshold hits).
      #
      # Useful for datasets dominated by long high-IC ZNFs (e.g. Myers 2018
      # hTERT) where some predicted TFs have a PWM but happen to score
      # weakly at the locus, while others are PWM-less chromatin readers
      # that we still want on the deck.
      if (isTRUE(deconv_min_motif_hits > 0)) {
        has_pwm <- vapply(all_plot_tfs, function(tf) {
          tf %in% names(motif_res_for_deck)
        }, logical(1))
        hit_count <- vapply(all_plot_tfs, function(tf) {
          if (tf %in% names(motif_res_for_deck)) {
            n <- motif_res_for_deck[[tf]]$n_hits
            if (is.null(n)) length(motif_res_for_deck[[tf]]$hits)
            else as.integer(n)
          } else {
            NA_integer_
          }
        }, integer(1))
        # Keep: (a) TFs with no PWM in JASPAR (exempt — peak-detection mode);
        #       (b) TFs with a PWM and >= threshold hits.
        # Drop: TFs with a PWM but < threshold hits.
        keep <- is.na(hit_count) | hit_count >= deconv_min_motif_hits
        n_no_pwm    <- sum(!has_pwm)
        n_pwm_pass  <- sum(has_pwm & hit_count >= deconv_min_motif_hits,
                           na.rm = TRUE)
        n_pwm_fail  <- sum(has_pwm & hit_count <  deconv_min_motif_hits,
                           na.rm = TRUE)
        message("  Deconvolution motif-hit filter: keeping ", sum(keep),
                " / ", length(all_plot_tfs), " TFs (",
                n_pwm_pass, " with PWM >= ", deconv_min_motif_hits,
                " hits, ", n_no_pwm, " with no JASPAR PWM (exempt), ",
                n_pwm_fail, " with PWM dropped)")
        all_plot_tfs <- all_plot_tfs[keep]
      }

      # Optional per-zone logFC floor — keep only TFs whose maximum per-region
      # logFC across the matched regions (those in pos_map) is >= the threshold.
      # Composes AND-style with the motif-hit filter above. Useful when input
      # p-values are unreliable (saturated, synthesized, etc.) and a raw
      # proximity-enrichment criterion is the cleaner gate. TFs with no rows
      # in long_data for the matched regions (or all-NA lfc) are treated as
      # max_lfc = -Inf and dropped whenever this threshold > 0.
      if (isTRUE(deconv_min_max_lfc > 0) && length(all_plot_tfs) > 0) {
        ld_match <- long_data[long_data$region %in% names(pos_map) &
                                long_data$protein %in% all_plot_tfs &
                                !is.na(long_data$lfc), ]
        max_lfc <- if (nrow(ld_match) > 0)
          tapply(ld_match$lfc, ld_match$protein, max, na.rm = TRUE)
        else setNames(numeric(0), character(0))
        max_lfc_vec <- vapply(all_plot_tfs, function(tf) {
          v <- max_lfc[tf]
          if (is.null(v) || is.na(v)) -Inf else as.numeric(v)
        }, numeric(1))
        keep_lfc <- max_lfc_vec >= deconv_min_max_lfc
        message("  Deconvolution logFC filter: keeping ", sum(keep_lfc),
                " / ", length(all_plot_tfs), " TFs with max-zone logFC >= ",
                deconv_min_max_lfc, " (dropping ", sum(!keep_lfc), ")")
        all_plot_tfs <- all_plot_tfs[keep_lfc]
      }

      # Order is computed in two pools so augment-mode TFs don't get an
      # alphabetical-tail page slot they don't deserve:
      #
      #   Pool 1 — original motif_tfs (the select_motif_tfs cut, ~44 on
      #     Myers/hTERT). Sorted by total predicted-event weight desc using
      #     `binding_events`; motif_tfs without any events fall to the
      #     alphabetical tail of pool 1 (legacy behaviour).
      #
      #   Pool 2 — augment-pool TFs (those added via
      #     motif_scan_pool='spatial_all' that were NOT in motif_tfs).
      #     Sorted among themselves by JASPAR motif-hit count desc, with
      #     name asc as tiebreaker. Rationale: predict_all_binding_events
      #     only ran on motif_tfs, so augment-pool TFs aren't in
      #     binding_events and a true summed-β order across the whole deck
      #     would require on-the-fly event computation. Ordering by motif
      #     evidence is the cheap, honest tiebreaker for this second pool —
      #     "more JASPAR hits in the window" is a fair priority signal when
      #     we can't compare event totals directly.
      #
      # The deck therefore reads: best motif_tfs by event weight → motif_tfs
      # with no events (alpha) → augment-pool TFs by motif-hit count.
      if (length(all_plot_tfs) == 0) {
        detail_tfs <- character(0)
      } else {
        base_pool    <- intersect(as.character(motif_tfs), all_plot_tfs)
        augment_pool <- setdiff(all_plot_tfs, as.character(motif_tfs))

        # ---- Pool 1: legacy ordering by event weight (motif_tfs only) ----
        if (length(base_pool) > 0 && nrow(binding_events) > 0) {
          totals <- aggregate(weight ~ tf, data = binding_events, FUN = sum)
          totals <- totals[order(totals$weight, decreasing = TRUE), ]
          with_events  <- as.character(totals$tf)
          base_with    <- intersect(with_events, base_pool)
          base_without <- setdiff(base_pool, base_with)
          base_ordered <- c(base_with, sort(base_without))
        } else {
          base_ordered <- sort(base_pool)
        }

        # ---- Pool 2: augment-pool by motif-hit count desc, then alpha ----
        if (length(augment_pool) > 0) {
          aug_hits <- vapply(augment_pool, function(tf) {
            if (tf %in% names(motif_res_for_deck)) {
              n <- motif_res_for_deck[[tf]]$n_hits
              if (is.null(n)) length(motif_res_for_deck[[tf]]$hits)
              else as.integer(n)
            } else {
              0L
            }
          }, integer(1))
          # order() is stable: hits desc (negate) then name asc
          aug_ordered <- augment_pool[order(-aug_hits, augment_pool)]
        } else {
          aug_ordered <- character(0)
        }

        detail_tfs <- c(base_ordered, aug_ordered)
      }

      if (length(detail_tfs) == 0) {
        message("  Deconvolution plots: 0 TF(s) survived the filters ",
                "(deconv_min_motif_hits = ", deconv_min_motif_hits,
                ", deconv_min_max_lfc = ", deconv_min_max_lfc,
                ") \u2014 skipping binding_deconvolution.pdf")
      } else {
        # Augment-pool ChIP-Atlas top-up — survivors that weren't in the
        # original 44 motif_tfs haven't had their ChIP-Atlas peaks fetched
        # yet (run_chipatlas_scan above only ran on motif_tfs). Without
        # this top-up, augment-pool pages would render in Plot 10 missing
        # the ChIP-Atlas lane while their motif_tfs counterparts still
        # have it — silently dropping public-peak validation evidence.
        # Done AFTER all filters have collapsed all_plot_tfs -> detail_tfs
        # so we only spend ChIP-Atlas API/cache I/O on the TFs that will
        # actually be drawn. No-op when chipatlas=FALSE or the chipatlas
        # module wasn't loaded; also no-op when augment-pool is disabled
        # (motif_scan_pool='selected', augment_survivors is empty).
        if (isTRUE(chipatlas) && isTRUE(.caspex_chipatlas_loaded)) {
          augment_survivors <- setdiff(as.character(detail_tfs),
                                        as.character(motif_tfs))
          if (length(augment_survivors) > 0) {
            message("  Augment-pool ChIP-Atlas scan for ",
                    length(augment_survivors), " new TF(s) ",
                    "(deck survivors not in the original motif_tfs cut)")
            extra_chipatlas <- tryCatch(
              run_chipatlas_scan(
                tfs              = augment_survivors,
                gene_info        = gene_info,
                promoter_info    = promoter_info,
                upstream         = upstream,
                downstream       = downstream,
                threshold        = chipatlas_threshold,
                max_experiments  = chipatlas_max_experiments,
                special_interest_gene = special_interest_gene,
                special_interest_cap  = special_interest_cap,
                quiet            = chipatlas_quiet
              ),
              error = function(e) {
                message("    ChIP-Atlas augment-scan failed: ",
                        conditionMessage(e),
                        " \u2014 augment-pool TFs will render without peaks")
                list()
              })
            if (is.null(chipatlas_res)) chipatlas_res <- list()
            # Merge: any name overlap (shouldn't happen by construction,
            # since augment_survivors is disjoint from motif_tfs) is
            # resolved in favour of the new scan.
            chipatlas_res[names(extra_chipatlas)] <- extra_chipatlas
          }
        }

        detail_plots <- lapply(detail_tfs, function(tf) {
          # motif_res_for_deck = motif_res ∪ motif_res_extra (see filter block
          # above). With motif_scan_pool='selected' it equals motif_res, so
          # this is identical to the legacy path; with 'spatial_all' it lets
          # newly-augmented TFs render with their actual motif hits instead of
          # falling through to the empty integer(0) branch.
          hits <- if (tf %in% names(motif_res_for_deck))
            motif_res_for_deck[[tf]]$hits else integer(0)
          # Per-motif PWM scores parallel to `hits` — must be threaded
          # through so the deck plot's internal predict_*() rerun applies
          # the same motif_score_weight as the parent run. Without this,
          # the deck always rendered the "none" version even for runs
          # that used "linear" or "log", silently disagreeing with
          # binding_events.csv.
          scores <- if (tf %in% names(motif_res_for_deck))
            motif_res_for_deck[[tf]]$score_frac else NULL
          plot_binding_deconvolution(tf, long_data, pos_map, hits,
                                      kernel_sigma    = kernel_sigma,
                                      min_weight_frac = min_weight_frac,
                                      upstream        = upstream,
                                      downstream      = downstream,
                                      weight_mode     = signal_weight,
                                      cov_floor       = cov_floor,
                                      edge_guard_frac = edge_guard_frac,
                                      zone_peak_frac    = zone_peak_frac,
                                      max_events_per_tf = max_events_per_tf,
                                      merge_position    = merge_position,
                                      max_grna_distance = max_grna_distance,
                                      edge_grna_weight_cap = edge_grna_weight_cap,
                                      position_stability = position_stability,
                                      n_bootstrap        = n_bootstrap,
                                      motif_scores       = scores,
                                      motif_score_weight = motif_score_weight,
                                      chipatlas_peaks   = if (!is.null(chipatlas_res))
                                                            chipatlas_res[[tf]] else NULL)
        })

        # Multi-page PDF: 1 TF per page (cleaner than a giant patchwork grid now
        # that N can be 20+). Each page is one `plot_binding_deconvolution()` plot.
        # `.safe_pdf` handles the recurring "write failed" crash when the file
        # is locked by an open PDF viewer.
        .safe_pdf(file.path(out_dir, "binding_deconvolution.pdf"),
                  width  = plot_width,
                  height = max(5, plot_height * 0.6),
                  plots  = detail_plots)

        message("  Deconvolution plots: ", length(detail_tfs), " TF(s) (",
                length(intersect(as.character(unique(binding_events$tf)), detail_tfs)),
                " with predicted events)")
      }
    }

    message("\nPlots saved to: ", normalizePath(out_dir))
  }

  # ── 7. Write results tables ───────────────────────────────────────────────
  csv_path <- file.path(out_dir, paste0(gene, "_spatial_predictions.csv"))
  write.csv(spatial_df, csv_path, row.names = FALSE)

  # Binding events CSV (motif-constrained NNLS + peak-fallback predictions)
  if (nrow(binding_events) > 0) {
    binding_events$gene          <- gene
    binding_events$kernel_sigma  <- kernel_sigma
    binding_events$min_weight_frac <- min_weight_frac

    # Provenance: which selection deck did this TF come from?
    #   "common"               — top-composite across all regions (deck 06/07)
    #   "shared:R1+R5"         — significant in >=2 regions   (deck 11/12)
    #   "region-specific:R3"   — region-focal in R3           (deck 08/09)
    # By construction (see select_motif_tfs), the three buckets are disjoint,
    # so `deck` has exactly one "provenance type" per TF.
    .tag_deck <- function(tf) {
      if (tf %in% common_tfs) return("common")
      if (tf %in% shared_tfs) return(paste0("shared:", shared_tag[[tf]]))
      rs <- names(specific_map)[vapply(specific_map,
                                       function(v) tf %in% v, logical(1))]
      if (length(rs) > 0) return(paste0("region-specific:",
                                         paste(rs, collapse = ",")))
      NA_character_
    }
    binding_events$deck <- vapply(as.character(binding_events$tf),
                                   .tag_deck, character(1))

    events_path <- file.path(out_dir, paste0(gene, "_binding_events.csv"))
    write.csv(binding_events[order(binding_events$tf,
                                    -binding_events$weight), ],
              events_path, row.names = FALSE)
    message("Binding events saved to: ", events_path)
  }

  # Augment with motif hit counts
  spatial_df$motif_hits <- sapply(spatial_df$protein, function(tf) {
    if (tf %in% names(motif_res)) motif_res[[tf]]$n_hits else NA_integer_
  })

  message("Results saved to: ", csv_path)
  # Final reminder of the weight source — visible even if the earlier message
  # scrolled out of view.
  message("\n--- Run summary --------------------------------------------------")
  message(" Weight mode      : ", weight_mode)
  message(" Weight source    : ", weight_src)
  message(" Kernel sigma (bp): ", kernel_sigma)
  message(" Binding path     : ",
          paste0("coverage-normalized per-motif \u03b2=s/C (cov_floor=",
                 cov_floor, ")"))
  message(" PWM-score weight : ", motif_score_weight,
          if (motif_score_weight == "none")
            "  (column / bubble weight = 1; binary motif filter only)"
          else if (motif_score_weight == "linear")
            "  (column / bubble weight = score_frac)"
          else
            "  (column / bubble weight = 1 + log2(1 + score_frac))")
  message(" Spatial TFs      : ", nrow(spatial_df))
  message(" Binding events   : ", nrow(binding_events),
          " (", length(unique(binding_events$tf)), " TFs)")
  if (isTRUE(chipatlas) && !is.null(chipatlas_res)) {
    n_tf_hits <- sum(vapply(chipatlas_res,
                            function(x) !is.null(x) && nrow(x) > 0,
                            logical(1)))
    n_peaks   <- sum(vapply(chipatlas_res,
                            function(x) if (is.null(x)) 0L else nrow(x),
                            integer(1)))
    message(" ChIP-Atlas       : ", n_peaks, " peaks across ",
            n_tf_hits, "/", length(chipatlas_res),
            " TFs (threshold=", chipatlas_threshold, ")")
  }
  message("------------------------------------------------------------------")
  message("Done.")

  invisible(list(
    spatial_df     = spatial_df,
    long_data      = long_data,
    pos_map        = pos_map,
    gene_info      = gene_info,
    promoter_info  = promoter_info,
    motif_results  = motif_res,
    # Augment-pool motif scan results (named list, same shape as motif_results)
    # for the spatial_df TFs that were not in motif_tfs. Empty list when
    # motif_scan_pool='selected' (the default). Consumed only by Plot 10's
    # deconv_min_motif_hits filter; downstream code that iterates over
    # motif_results stays unchanged.
    motif_results_extra = motif_res_extra,
    motif_scan_pool     = motif_scan_pool,
    # The binding-deconvolution deck (binding_deconvolution.pdf) TF roster — every
    # TF that landed on a per-TF detail page after applying
    # `deconv_min_motif_hits` and `deconv_min_max_lfc`. Exposed on the
    # result so downstream consumers (notably the A5 TF co-occurrence
    # extras plot) can union-merge with their own selection criterion
    # and guarantee that any TF visible on Plot 10 is also visible on
    # A5. Empty character vector when save_plots = FALSE.
    detail_tfs     = detail_tfs,
    binding_events = binding_events,
    weight_mode    = weight_mode,
    weight_source  = weight_src,
    # Mode bookkeeping — lets caspex_extras.R auto-detect which analytical
    # path produced this result and rerun its jackknife / σ sweep / one-pagers
    # under the matching mode. Without these fields the extras would silently
    # run in default-mode NNLS even when `result` came from a coverage-aware
    # run, and would disagree with the events in `binding_events`.
    cov_floor        = cov_floor,
    edge_guard_frac  = edge_guard_frac,
    kernel_sigma     = kernel_sigma,
    # Window that the run actually used (TSS-relative bp). Stored so
    # downstream diagnostics — notably plot_coverage_stack / D3 — can
    # rebuild the correct x_grid instead of falling back to its old
    # hardcoded seq(-2500, 500, 5) which sits in the wrong place for any
    # promoter where guides land outside that range (e.g. Mackenzie FOXP2
    # has guides at +973..+1188, well past the legacy +500 boundary).
    upstream         = upstream,
    downstream       = downstream,
    # Wild bootstrap diagnostic — bookkeeping so callers / extras can detect
    # whether the position_stability columns were generated for this run.
    position_stability = position_stability,
    n_bootstrap        = n_bootstrap,
    binding_path     = paste0("coverage-normalized per-motif \u03b2=s/C (cov_floor=",
                              cov_floor, ", edge_guard_frac=",
                              edge_guard_frac, ")"),
    # ChIP-Atlas overlay (named list tf -> data.frame; NULL if chipatlas=FALSE).
    # Kept on the result so caspex_extras or interactive inspection can reuse
    # the windowed peak tables without re-hitting the network.
    chipatlas_peaks       = chipatlas_res,
    chipatlas_threshold   = if (isTRUE(chipatlas)) chipatlas_threshold else NULL,
    plots = list(grna = p_grna,
                 heat = p_heat,
                 gc   = p_gc)
  ))
}

# =============================================================================
# SECTION 8: Convenience helpers
# =============================================================================

#' Inspect a single TF from a completed run
#'
#' @param result   Output of run_caspex()
#' @param tf_name  TF symbol to inspect
#' @export
inspect_tf <- function(result, tf_name) {
  p <- plot_tf_profile(tf_name,
                        result$long_data,
                        result$pos_map)
  sp <- result$spatial_df[result$spatial_df$protein == tf_name, ]
  if (nrow(sp) == 0) {
    message(tf_name, " not found in spatial model")
  } else {
    cat("\n\u2500\u2500 Spatial estimate for", tf_name, "\u2500\u2500\n")
    cat("  Centroid  :", sprintf("%+.0f bp\n", sp$centroid))
    cat("  Spread    :", sprintf("%.0f bp\n",  sp$spread))
    cat("  Composite :", sprintf("%.4f\n",     sp$composite))
    cat("  Regions   :", sp$n_regions, "\n")
    cat("  Mean logFC:", sprintf("%.3f\n",     sp$mean_lfc))
    cat("  Min p-val :", format(sp$min_pval, digits = 3), "\n")
  }
  if (tf_name %in% names(result$motif_results)) {
    cat("  Motif hits:", result$motif_results[[tf_name]]$n_hits, "\n")
  }
  print(p)
  invisible(p)
}

#' Re-run motif scan with a different TF set or threshold
#' @param result (see function body).
#' @param tf_names (see function body).
#' @param threshold_frac (see function body).
#' @export
rescan_motifs <- function(result, tf_names, threshold_frac = 0.80) {
  run_motif_scan(tf_names, result$promoter_info, threshold_frac)
}

#' Audit which Ensembl transcript best matches your sgRNA layout
#'
#' Pre-flight diagnostic for any new GLproxScape analysis. For the
#' requested gene, queries Ensembl for ALL annotated transcripts and
#' checks which sgRNAs from a manifest fall inside each transcript's
#' promoter window. Reports per-transcript: number of sgRNAs matched,
#' their TSS-relative positions, canonical flag, biotype, and transcript
#' span; finishes with a recommended \code{transcript = "ENST..."}
#' anchor for \code{\link{run_caspex}}.
#'
#' Use BEFORE picking the \code{transcript} argument for
#' \code{\link{run_caspex}}, since alt-promoters of the same gene can
#' sit hundreds of kb apart (e.g. FOXP2's HEK293 "active TSS1" is
#' ~328 kb upstream of the canonical FOXP2-201 TSS, so a run anchored
#' on the canonical transcript would not find Mackenzie's sgRNAs at all).
#'
#' @param gene HGNC symbol (e.g. \code{"FOXP2"}, \code{"ABCB1"}).
#' @param manifest_path Path to the folder containing the sgRNA
#'   manifest. Resolved relative to the current working directory if
#'   not absolute. Default \code{"inputs"}.
#' @param manifest Manifest filename inside \code{manifest_path}.
#'   Default \code{"grnas.tsv"}.
#' @param species Ensembl species token. Default
#'   \code{"homo_sapiens"}.
#' @param upstream,downstream bp window around each candidate TSS in
#'   which to look for sgRNA matches. Deliberately wider than what
#'   \code{run_caspex()} uses, so transcripts with TSSs slightly
#'   offset from the guide tile still register matches.
#'   Defaults 5000 / 1500.
#' @param biotype_filter Character vector of Ensembl biotypes to keep.
#'   Default \code{"protein_coding"}. Set to NULL or \code{character(0)}
#'   to scan every biotype (includes \code{nonsense_mediated_decay},
#'   \code{retained_intron}, \code{processed_transcript}, etc.) —
#'   useful when your TF binds an alt-promoter Ensembl tags as
#'   non-coding.
#' @param verbose Print the per-transcript trace + recommendation
#'   table to the console (default TRUE). Set FALSE to suppress.
#' @return Invisibly, a data.frame with one row per transcript:
#'   columns \code{transcript}, \code{biotype}, \code{is_canonical},
#'   \code{tss}, \code{n_matched}, plus one column per sgRNA region
#'   holding its TSS-relative bp position (NA = no match). Sorted by
#'   \code{n_matched} descending.
#' @examples
#' \dontrun{
#' # Use the bundled FOXP2 example data
#' inputs_dir <- system.file("extdata/examples/foxp2_mackenzie",
#'                           package = "GLproxScape")
#' df <- check_transcripts(gene = "FOXP2", manifest_path = inputs_dir)
#' head(df, 5)
#' }
#' @export
check_transcripts <- function(gene,
                              manifest_path  = "inputs",
                              manifest       = "grnas.tsv",
                              species        = "homo_sapiens",
                              upstream       = 5000,
                              downstream     = 1500,
                              biotype_filter = "protein_coding",
                              verbose        = TRUE) {
  inputs <- load_caspex_inputs(manifest_path, manifest = manifest)
  grnas  <- inputs$grnas
  if (verbose) {
    cat("=== sgRNA manifest ===\n")
    for (i in seq_along(grnas))
      cat(sprintf("  %s : %s\n", names(grnas)[i], grnas[[i]]))
    cat(sprintf("\n=== Looking up %s in Ensembl (%s) ===\n", gene, species))
  }
  js <- ensembl_get(paste0("/lookup/symbol/", species, "/", gene),
                    params = list(expand = 1,
                                  `content-type` = "application/json"))
  if (is.null(js))
    stop("Gene '", gene, "' not found in Ensembl (", species, ")")
  tx_list <- js$Transcript
  if (is.null(tx_list) || length(tx_list) == 0)
    stop("No transcripts returned for ", gene)

  if (verbose)
    cat(sprintf("Gene span: chr%s:%d-%d, strand=%s, %d transcripts\n",
                js$seq_region_name, js$start, js$end,
                if (js$strand == 1) "+" else "-", length(tx_list)))

  if (length(biotype_filter) > 0) {
    bt   <- vapply(tx_list, function(tx) tx$biotype %||% "", character(1))
    keep <- bt %in% biotype_filter
    if (verbose)
      cat(sprintf("Filtering by biotype %s: %d / %d transcripts kept\n",
                  paste(biotype_filter, collapse = ","),
                  sum(keep), length(tx_list)))
    tx_list <- tx_list[keep]
  }

  make_gene_info <- function(tx, js)
    list(name = js$display_name, chr = js$seq_region_name,
         strand = tx$strand %||% js$strand,
         tss = if ((tx$strand %||% js$strand) == 1) tx$start else tx$end,
         start = tx$start, end = tx$end, species = species,
         transcript_id = tx$id)

  if (verbose)
    cat(sprintf("\n=== Testing %d transcript(s) with window upstream=%d / downstream=%d bp ===\n\n",
                length(tx_list), upstream, downstream))

  results <- vector("list", length(tx_list))
  for (i in seq_along(tx_list)) {
    tx <- tx_list[[i]]
    gi <- make_gene_info(tx, js)
    is_canon <- isTRUE(tx$is_canonical == 1)
    bio <- tx$biotype %||% NA_character_
    if (verbose)
      cat(sprintf("[%d/%d] %s%s  biotype=%s  span=%d-%d (%d bp)  TSS=%d\n",
                  i, length(tx_list), tx$id,
                  if (is_canon) "  [canonical]" else "",
                  bio, tx$start, tx$end, tx$end - tx$start + 1, gi$tss))
    prom <- tryCatch(suppressMessages(fetch_promoter_seq(gi, upstream,
                                                          downstream,
                                                          species = species)),
                     error = function(e) NULL)
    if (is.null(prom)) {
      results[[i]] <- list(tx = tx$id, biotype = bio, is_canon = is_canon,
                           tss = gi$tss, n_matched = NA_integer_,
                           positions = setNames(rep(NA_real_, length(grnas)),
                                                 names(grnas)))
      next
    }
    positions <- setNames(rep(NA_real_, length(grnas)), names(grnas))
    strands   <- setNames(rep(NA_character_, length(grnas)), names(grnas))
    for (r in names(grnas)) {
      res <- match_grna(grnas[[r]], prom)
      if (res$n_hits > 0) {
        positions[[r]] <- res$positions[1]
        strands[[r]]   <- res$strands[1]
      }
    }
    n_matched <- sum(!is.na(positions))
    if (verbose) {
      for (r in names(grnas)) {
        if (!is.na(positions[[r]]))
          cat(sprintf("    %-3s : %+d bp (%s strand)\n",
                      r, as.integer(positions[[r]]), strands[[r]]))
        else
          cat(sprintf("    %-3s : no match\n", r))
      }
      cat(sprintf("    -> %d / %d sgRNAs matched\n\n",
                  n_matched, length(grnas)))
    }
    results[[i]] <- list(tx = tx$id, biotype = bio, is_canon = is_canon,
                         tss = gi$tss, n_matched = n_matched,
                         positions = positions)
  }
  ord <- order(vapply(results, function(r) r$n_matched %||% -1L, integer(1)),
               decreasing = TRUE, na.last = TRUE)
  if (verbose) {
    cat("=== Summary (sorted by sgRNAs matched, DESC) ===\n")
    hdr <- sprintf("%-18s  %-22s  %-10s  %10s  %s",
                   "transcript", "biotype", "canonical", "TSS",
                   paste(sprintf("%6s", names(grnas)), collapse = "  "))
    cat(hdr, "\n"); cat(strrep("-", nchar(hdr)), "\n")
    for (i in ord) {
      r <- results[[i]]
      pos_str <- vapply(names(grnas), function(rg) {
        p <- r$positions[[rg]]
        if (is.na(p)) "    -" else sprintf("%+6d", as.integer(p))
      }, character(1))
      cat(sprintf("%-18s  %-22s  %-10s  %10d  %s\n",
                  r$tx, r$biotype %||% "-",
                  if (isTRUE(r$is_canon)) "yes" else "no",
                  r$tss, paste(pos_str, collapse = "  ")))
    }
    best <- which.max(vapply(results, function(r) r$n_matched %||% -1L,
                              integer(1)))
    cat("\n=== Recommendation ===\n")
    cat(sprintf("Best transcript: %s  (%d / %d sgRNAs matched)\n",
                results[[best]]$tx, results[[best]]$n_matched, length(grnas)))
    cat(sprintf("Use in run_caspex():  transcript = \"%s\"\n",
                results[[best]]$tx))
  }
  rows <- lapply(results, function(r) {
    base <- data.frame(transcript = r$tx,
                       biotype = r$biotype %||% NA_character_,
                       is_canonical = isTRUE(r$is_canon),
                       tss = r$tss,
                       n_matched = r$n_matched %||% NA_integer_,
                       stringsAsFactors = FALSE)
    for (rg in names(r$positions)) base[[rg]] <- r$positions[[rg]]
    base
  })
  invisible(do.call(rbind, rows[ord]))
}

# (Previously a script-style "CasPEX loaded; usage:" message printed at the
# end of source(). In package context, top-level message() calls are
# package-load side effects that R CMD check flags. Removed; usage lives
# in ?GLproxScape and the README. To re-add a startup banner, define a
# `.onAttach <- function(libname, pkgname) packageStartupMessage(...)`
# in this file — that's the package-correct hook.)

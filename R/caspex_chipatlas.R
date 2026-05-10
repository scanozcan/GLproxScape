# =============================================================================
# caspex_chipatlas.R — ChIP-Atlas peak backend for CasPEX
# -----------------------------------------------------------------------------
# Pulls public ChIP-seq peak calls from ChIP-Atlas (chip-atlas.dbcls.jp) and
# exposes them as a supplementary track beneath the predicted-event bubbles in
# plot 10 (binding deconvolution), plus a compact union-peak row beneath the
# JASPAR motif ticks in the mini-browser decks (06/07, 08/09, 11/12).
#
# Strategy:
#   1. Download experimentList.tab (one-time, ~50 MB, cached) -> maps SRX -> TF
#      symbol (antigen) and cell type. ChIP-Atlas publishes this at
#      https://chip-atlas.dbcls.jp/data/metadata/experimentList.tab.
#   2. For each requested TF, look up its SRX IDs in hg38 TFs class, take the
#      top-N most recent (by SRX number DESC), and fetch per-SRX BED files
#      from .../eachData/bed{threshold}/SRX*.{threshold}.bed. Each is small
#      (<1 MB typically), cached to R_user_dir("caspex","cache")/chipatlas/.
#   3. Filter every SRX BED to the promoter window and convert genomic
#      coordinates to TSS-relative bp.
#
# The "assembled" per-antigen path (Oth.ALL.05.{TF}.AllCell.bed) was rejected
# because files for well-studied TFs are 100-300 MB each; the per-SRX path
# downloads only what we show.
# =============================================================================

# `%||%` is defined once in R/utils-internal.R and visible to every R/
# file in the package namespace.

#' ChIP-Atlas cache directory.
#'
#' @return Path to the cache root (created on demand).
#' @noRd
.chipatlas_cache_dir <- function() {
  base <- tools::R_user_dir("caspex", which = "cache")
  dir  <- file.path(base, "chipatlas")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  dir
}

# ---- experimentList.tab ------------------------------------------------------

#' URL of the ChIP-Atlas experimentList.tab metadata.
#' @noRd
.chipatlas_experiment_list_url <- function() {
  "https://chip-atlas.dbcls.jp/data/metadata/experimentList.tab"
}

#' Download (if missing) and cache ChIP-Atlas experimentList.tab
#'
#' @param force logical; re-download even if cached.
#' @param quiet passed to download.file().
#' @return path to the cached TSV.
#' @export
download_chipatlas_experiment_list <- function(force = FALSE, quiet = FALSE) {
  fpath <- file.path(.chipatlas_cache_dir(), "experimentList.tab")
  # Sanity check existing cache. ChIP-Atlas experimentList.tab is ~300+ MB and
  # each line starts with an experiment accession — SRX (SRA), DRX (DDBJ), or
  # ERX (ENA) — followed by a tab. A previous failed download could leave an
  # HTML stub in the cache; we reject anything that doesn't match the
  # expected pattern.
  looks_like_exp_list <- function(first_line) {
    length(first_line) >= 1 &&
      grepl("^[SDE]RX[0-9]+\t", first_line[[1]])
  }
  if (file.exists(fpath) && !force) {
    ok_size <- file.size(fpath) > 50e6    # file is typically 300+ MB
    first   <- tryCatch(readLines(fpath, n = 1, warn = FALSE),
                        error = function(e) "")
    ok_head <- looks_like_exp_list(first)
    if (ok_size && ok_head) return(fpath)
    if (!quiet) message("  Cached experimentList.tab looks bad (",
                        format(file.size(fpath), big.mark = ","), " bytes, ",
                        "first line: ", substr(first[[1]], 1, 40), "...); ",
                        "re-downloading.")
    file.remove(fpath)
  }
  url <- .chipatlas_experiment_list_url()
  if (!quiet) message("  Downloading ChIP-Atlas experimentList (~300 MB, one-time)...")
  # ChIP-Atlas is served over HTTPS; mode="wb" for cross-platform safety.
  utils::download.file(url, fpath, mode = "wb", quiet = quiet)
  if (!file.exists(fpath) || file.size(fpath) < 50e6) {
    sz <- if (file.exists(fpath)) file.size(fpath) else 0
    if (file.exists(fpath)) file.remove(fpath)
    stop("experimentList.tab download produced a truncated file (",
         sz, " bytes); try again or check network.")
  }
  first <- readLines(fpath, n = 1, warn = FALSE)
  if (!looks_like_exp_list(first)) {
    file.remove(fpath)
    stop("experimentList.tab first line doesn't look like an ",
         "SRX/DRX/ERX record; server returned something unexpected. ",
         "First line was: ", substr(first[[1]], 1, 80))
  }
  fpath
}

# Parse into an in-memory data.frame. ChIP-Atlas layout (as of 2024+):
#   col1 = SRX ID         col2 = Genome (hg38 etc)
#   col3 = Antigen class  col4 = Antigen (TF symbol)
#   col5 = Cell type class col6 = Cell type
#   col7 = Cell type desc  col8 = Processing log (messy)
#   col9+ = per-threshold stats, GSE, title, attributes
# We only need cols 1-6.
#
# We use readLines + strsplit instead of read.table because experimentList.tab
# has variable column counts per row (some rows have 14 fields, some 20+) and
# read.table with col.names/flush/fill has a history of silently truncating
# rows or misaligning columns when the file has mixed widths. The line-based
# parse is bulletproof: we always take exactly the first 6 fields.
.chipatlas_experiment_list <- local({
  cached <- NULL
  function(force_reload = FALSE, verbose = TRUE) {
    if (!is.null(cached) && !force_reload) return(cached)
    fpath <- download_chipatlas_experiment_list(force = FALSE, quiet = !verbose)
    if (verbose) message("  Parsing experimentList.tab (",
                          format(file.size(fpath), big.mark = ","), " bytes)...")
    lines <- readLines(fpath, warn = FALSE, encoding = "UTF-8")
    if (verbose) message("    Read ", format(length(lines), big.mark = ","),
                         " lines.")
    # Split each line on tabs, take first 6 fields, pad if fewer.
    split6 <- function(x) {
      v <- strsplit(x, "\t", fixed = TRUE)[[1]]
      length(v) <- 6            # NA-pads short rows; truncates long rows
      v
    }
    mat <- do.call(rbind, lapply(lines, split6))
    df  <- data.frame(srx             = mat[, 1],
                      genome          = mat[, 2],
                      antigen_class   = mat[, 3],
                      antigen         = mat[, 4],
                      cell_type_class = mat[, 5],
                      cell_type       = mat[, 6],
                      stringsAsFactors = FALSE)
    if (verbose) {
      message("    Parsed ", format(nrow(df), big.mark = ","), " rows.")
      # Breakdown before filter so misparses are obvious.
      g_tbl <- sort(table(df$genome), decreasing = TRUE)
      message("    Top genomes: ",
              paste(sprintf("%s=%s", names(g_tbl)[seq_len(min(5, length(g_tbl)))],
                            format(as.integer(g_tbl[seq_len(min(5, length(g_tbl)))]),
                                   big.mark = ",")),
                    collapse = ", "))
      c_tbl <- sort(table(df$antigen_class), decreasing = TRUE)
      message("    Top antigen_class values: ",
              paste(sprintf("%s=%s", names(c_tbl)[seq_len(min(5, length(c_tbl)))],
                            format(as.integer(c_tbl[seq_len(min(5, length(c_tbl)))]),
                                   big.mark = ",")),
                    collapse = ", "))
    }
    # Keep hg38 TFs and Histones; drop DNase / ATAC / Bisulfite / etc. Try
    # several spellings because ChIP-Atlas has been known to shift labels
    # between releases. Histone rows are needed by the locus-level histone-
    # marks page (fetch_histone_peaks_for_locus); TF code paths query by
    # antigen NAME (e.g. "MYC") so adding histone rows doesn't change any
    # existing TF behaviour.
    keep_class_patterns <- c("TFs and others", "TFs_and_others",
                              "TF", "TFs", "Transcription factor",
                              "Histone", "Histones", "Histone modification")
    keep <- df$genome == "hg38" & df$antigen_class %in% keep_class_patterns
    if (verbose) message("    hg38 TF + Histone rows after filter: ",
                         format(sum(keep), big.mark = ","))
    if (sum(keep) == 0) {
      warning("ChIP-Atlas: 0 rows matched hg38 + (TF | Histone) filter. ",
              "Check the antigen_class breakdown above and update ",
              "keep_class_patterns in caspex_chipatlas.R.")
    }
    df <- df[keep, , drop = FALSE]
    cached <<- df
    df
  }
})

#' Experiment IDs for a given TF (antigen symbol) on hg38
#'
#' @param tf HGNC symbol; matched case-insensitively against `antigen`.
#' @return character vector of experiment IDs (SRX/DRX/ERX), ordered
#'   newest-first by numeric suffix within each prefix.
#' @noRd
chipatlas_srx_for_tf <- function(tf) {
  el <- .chipatlas_experiment_list()
  hit <- which(toupper(el$antigen) == toupper(tf))
  if (!length(hit)) return(character(0))
  ids <- el$srx[hit]
  # ChIP-Atlas stores experiments from SRA (SRX), DDBJ (DRX), and ENA (ERX).
  # Strip any 3-letter prefix before ordering so DRX1234567 sorts by 1234567.
  # Accession numbers within each prefix are monotonic in submission date.
  ord <- order(suppressWarnings(as.integer(sub("^[SDE]RX", "", ids))),
               decreasing = TRUE)
  ids[ord]
}

# ---- per-SRX BED fetch -------------------------------------------------------

#' URL of a per-SRX BED file at a given threshold.
#' @param srx (see function body).
#' @param genome (see function body).
#' @param threshold (see function body).
#' @noRd
.chipatlas_srx_bed_url <- function(srx, genome = "hg38", threshold = "05") {
  sprintf("https://chip-atlas.dbcls.jp/data/%s/eachData/bed%s/%s.%s.bed",
          genome, threshold, srx, threshold)
}

#' Download (if missing) a single SRX peak BED
#'
#' Returns the cached file path, or NULL on failure (404, truncated, HTML
#' error stub, etc.). `download.file()` returns 0 on success; we also
#' sniff the first line of the downloaded file because some ChIP-Atlas
#' 404 paths return a 200 status with an HTML body.
#' @param srx (see function body).
#' @param genome (see function body).
#' @param threshold (see function body).
#' @param force (see function body).
#' @param quiet (see function body).
#' @param max_retries (see function body).
#' @param timeout_sec (see function body).
#' @noRd
download_chipatlas_srx_bed <- function(srx, genome = "hg38", threshold = "05",
                                       force = FALSE, quiet = TRUE,
                                       max_retries = 3,
                                       timeout_sec = 600) {
  cache <- .chipatlas_cache_dir()
  sub   <- file.path(cache, genome, paste0("bed", threshold))
  dir.create(sub, showWarnings = FALSE, recursive = TRUE)
  fpath <- file.path(sub, sprintf("%s.%s.bed", srx, threshold))
  if (file.exists(fpath) && !force) return(fpath)
  url <- .chipatlas_srx_bed_url(srx, genome, threshold)

  # Per-SRX BEDs for well-studied TFs (EP300, CTCF, RAD21, ...) are routinely
  # 1-3 MB. R's default download.file() timeout is 60 s, which silently
  # truncates these on slower or congested networks - download.file then
  # returns non-zero status, the partial file is deleted, and the caller
  # sees NULL with no warning. We bump the session timeout for the duration
  # of this call and add up to `max_retries` retries with 2s/4s backoff
  # before giving up.
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(timeout_sec, old_timeout %||% 0))

  for (attempt in seq_len(max_retries)) {
    status <- tryCatch(
      suppressWarnings(
        utils::download.file(url, fpath, mode = "wb", quiet = quiet)),
      error = function(e) -1L
    )
    ok_status <- identical(as.integer(status), 0L) && file.exists(fpath)

    # Sniff first line: valid BED rows start with "chr" (or a digit for some
    # bacterial assemblies). HTML 404 stubs start with "<" or whitespace.
    ok_content <- FALSE
    if (ok_status) {
      first <- tryCatch(readLines(fpath, n = 1, warn = FALSE),
                        error = function(e) character(0))
      ok_head <- length(first) >= 1 &&
                 (grepl("^(chr|track|browser)", first[[1]]) ||
                  grepl("^[0-9]", first[[1]]))
      ok_content <- ok_head && file.size(fpath) >= 10L
    }

    if (ok_status && ok_content) return(fpath)

    # Bad attempt - clean up and (maybe) retry.
    if (file.exists(fpath)) file.remove(fpath)
    if (attempt < max_retries) Sys.sleep(2 * attempt)  # 2s, 4s backoff
  }
  return(NULL)
}

# Per-SRX BEDs are narrowPeak-style (BED8+) — observed first line:
#   chr1 10018 10209 SRX11664714.05_peak_1 416 . 20.20613 47.189
# i.e. col4 is peak name (text), col5 is score (integer), col7 is signalValue.
# We only need chr/start/end for window overlap; optional score goes into
# the returned data.frame purely for downstream convenience.
#
# `read.table` with `col.names`+`colClasses` mis-typed the 4th field as
# numeric and silently dropped every row. Line-based parse avoids that.
#' Parse a cached per-SRX BED file into a data.frame.
#' @param fpath (see function body).
#' @noRd
.read_chipatlas_srx_bed <- function(fpath) {
  if (is.null(fpath) || !file.exists(fpath) || file.size(fpath) == 0) return(NULL)
  lines <- tryCatch(readLines(fpath, warn = FALSE, encoding = "UTF-8"),
                    error = function(e) character(0))
  if (!length(lines)) return(NULL)
  # Drop UCSC track/browser header lines and blank lines.
  lines <- lines[nzchar(lines) &
                 !grepl("^(track|browser|#)", lines, ignore.case = TRUE)]
  if (!length(lines)) return(NULL)
  # Split on tab, take first 3 fields + field 5 (narrowPeak score) if present.
  split_row <- function(x) {
    v <- strsplit(x, "\t", fixed = TRUE)[[1]]
    length(v) <- 5              # NA-pad; truncate to 5
    v
  }
  mat <- do.call(rbind, lapply(lines, split_row))
  df  <- data.frame(
    chr   = mat[, 1],
    start = suppressWarnings(as.integer(mat[, 2])),
    end   = suppressWarnings(as.integer(mat[, 3])),
    score = suppressWarnings(as.numeric(mat[, 5])),
    stringsAsFactors = FALSE
  )
  # Keep only rows whose chr looks like a real chromosome name and whose
  # coordinates parsed. This also guards against any stray header that
  # slipped through the "track/browser" filter.
  keep <- !is.na(df$start) & !is.na(df$end) &
          grepl("^(chr|[0-9]|[XYM])", df$chr)
  df <- df[keep, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  df
}

# ---- window filter + coordinate conversion ----------------------------------

#' Compute genomic coordinates of the promoter window fetched by Ensembl
#'
#' Mirrors the logic in fetch_promoter_seq() so we can map ChIP-Atlas peaks
#' onto the same TSS-relative axis the rest of the pipeline uses.
#' @param gene_info (see function body).
#' @param promoter_info (see function body).
#' @param upstream (see function body).
#' @param downstream (see function body).
#' @noRd
.chipatlas_window_coords <- function(gene_info, promoter_info,
                                      upstream, downstream) {
  chr    <- gene_info$chr
  tss    <- gene_info$tss
  strand <- gene_info$strand
  if (!grepl("^chr", chr)) chr <- paste0("chr", chr)
  if (strand == 1) {
    g_start <- max(1, tss - upstream)
    g_end   <- tss + downstream
  } else {
    g_start <- max(1, tss - downstream)
    g_end   <- tss + upstream
  }
  list(chr = chr, g_start = g_start, g_end = g_end, tss = tss, strand = strand)
}

#' Filter a per-SRX BED data.frame to the promoter window; add TSS-relative cols
#' @param df (see function body).
#' @param win (see function body).
#' @noRd
.chipatlas_filter_to_window <- function(df, win) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  chr_match <- df$chr == win$chr | paste0("chr", df$chr) == win$chr |
               df$chr == sub("^chr", "", win$chr)
  df <- df[chr_match, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  mid <- (df$start + df$end) / 2
  keep <- mid >= win$g_start & mid <= win$g_end
  df <- df[keep, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  if (win$strand == 1) {
    df$pos_rel   <- (df$start + df$end) / 2 - win$tss
    df$start_rel <- df$start - win$tss
    df$end_rel   <- df$end   - win$tss
  } else {
    df$pos_rel   <- win$tss - (df$start + df$end) / 2
    df$start_rel <- win$tss - df$end
    df$end_rel   <- win$tss - df$start
  }
  df
}

# Module-level flag: has the one-shot debug trace been emitted yet?
# Flips TRUE after the first SRX fetch of the first TF with >=1 experiment.
.chipatlas_debug_done <- FALSE

# ---- public API --------------------------------------------------------------

#' Fetch ChIP-Atlas peaks for one TF, windowed + TSS-relative
#'
#' Internal worker behind run_chipatlas_scan(). For a single TF, resolves
#' SRX accessions, downloads / reads the per-SRX BEDs from the on-disk
#' cache, and returns the per-peak rows that fall inside the requested
#' window.
#'
#' @param tf HGNC symbol.
#' @param gene_info result of lookup_gene().
#' @param promoter_info result of fetch_promoter_seq().
#' @param upstream,downstream bp around TSS (must match the promoter fetch).
#' @param threshold "05" (Q<1e-5, default), "10", or "20".
#' @param max_experiments cap on SRX count per TF (default 50).
#' @param is_special_interest if TRUE, bypass `max_experiments` per-TF cap.
#' @param special_interest_cap optional numeric cap when
#'   `is_special_interest = TRUE`; NULL means scan all available SRX.
#' @param quiet logical; suppress per-SRX download messages.
#' @return data.frame with columns (srx, cell_type_class, cell_type, chr,
#'   start, end, score, pos_rel, start_rel, end_rel) or NULL if no peaks.
#' @noRd
fetch_chipatlas_peaks <- function(tf, gene_info, promoter_info,
                                  upstream = 2500, downstream = 500,
                                  threshold = "05",
                                  max_experiments = 50,
                                  is_special_interest = FALSE,
                                  special_interest_cap = NULL,
                                  quiet = TRUE) {
  srx_ids <- chipatlas_srx_for_tf(tf)
  if (!length(srx_ids)) return(NULL)
  # Cap selection logic:
  #   - Standard TF: take the newest `max_experiments` SRX.
  #   - Special-interest TF, special_interest_cap = NULL: scan ALL SRX
  #     (could be 2000+ for CTCF — slow on first run, free thereafter).
  #   - Special-interest TF, special_interest_cap = N: take newest N SRX.
  #     Useful middle-ground (e.g. 250 covers most cell-type diversity for
  #     well-studied TFs without a full 2000-BED download).
  if (is_special_interest) {
    if (!is.null(special_interest_cap) &&
        length(srx_ids) > special_interest_cap)
      srx_ids <- srx_ids[seq_len(special_interest_cap)]
  } else {
    if (length(srx_ids) > max_experiments)
      srx_ids <- srx_ids[seq_len(max_experiments)]
  }
  n_srx_scanned <- length(srx_ids)
  win <- .chipatlas_window_coords(gene_info, promoter_info, upstream, downstream)

  el  <- .chipatlas_experiment_list()
  ct_lookup <- setNames(el$cell_type, el$srx)
  cc_lookup <- setNames(el$cell_type_class, el$srx)

  # One-shot debug trace: on the first SRX fetch of the first TF that HAS
  # experiments, print the URL, the local download result, and the parsed
  # row count. Makes silent 0-peak failures visible exactly once per session.
  if (!isTRUE(.chipatlas_debug_done)) {
    first_srx <- srx_ids[[1]]
    url <- .chipatlas_srx_bed_url(first_srx, "hg38", threshold)
    message("  [chipatlas debug] First fetch: ", tf, " / ", first_srx)
    message("    URL: ", url)
    fpath_dbg <- download_chipatlas_srx_bed(first_srx, threshold = threshold,
                                            quiet = FALSE)
    if (is.null(fpath_dbg)) {
      message("    -> download FAILED (404, HTML stub, or empty)")
    } else {
      sz  <- file.size(fpath_dbg)
      fl  <- tryCatch(readLines(fpath_dbg, n = 1, warn = FALSE),
                      error = function(e) "")
      df0 <- .read_chipatlas_srx_bed(fpath_dbg)
      n   <- if (is.null(df0)) 0 else nrow(df0)
      message("    -> downloaded ", format(sz, big.mark = ","), " bytes, ",
              "first line: ", substr(fl[[1]], 1, 60))
      message("    -> parsed ", format(n, big.mark = ","), " BED rows",
              if (!is.null(df0) && nrow(df0) > 0)
                paste0(" (first chr: ", df0$chr[1], " ", df0$start[1],
                       "-", df0$end[1], ")")
              else "")
      # Also show the filter window so mismatched chromosome names are
      # obvious (e.g. "chr13" vs "13").
      message("    -> window: ", win$chr, ":",
              format(win$g_start, big.mark = ","), "-",
              format(win$g_end,   big.mark = ","))
    }
    .chipatlas_debug_done <<- TRUE
  }

  pieces <- lapply(srx_ids, function(srx) {
    fpath <- download_chipatlas_srx_bed(srx, threshold = threshold, quiet = quiet)
    df <- .read_chipatlas_srx_bed(fpath)
    df <- .chipatlas_filter_to_window(df, win)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$srx <- srx
    df$cell_type       <- ct_lookup[[srx]] %||% NA_character_
    df$cell_type_class <- cc_lookup[[srx]] %||% NA_character_
    df
  })
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (!length(pieces)) {
    # Even when no peaks survive the window filter we still want to surface
    # how many SRX were scanned, so the deck can render
    # "0 SRX / N experiments" semantics if the caller chooses.
    out <- data.frame(srx = character(0), cell_type_class = character(0),
                      cell_type = character(0), chr = character(0),
                      start = integer(0), end = integer(0),
                      score = numeric(0), pos_rel = numeric(0),
                      start_rel = numeric(0), end_rel = numeric(0),
                      stringsAsFactors = FALSE)
    attr(out, "n_srx_scanned")       <- n_srx_scanned
    attr(out, "is_special_interest") <- is_special_interest
    return(out)
  }
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out <- out[, c("srx", "cell_type_class", "cell_type", "chr", "start", "end",
                 "score", "pos_rel", "start_rel", "end_rel")]
  # Attach scan-size metadata for downstream "X SRX / Y experiments" display.
  attr(out, "n_srx_scanned")       <- n_srx_scanned
  attr(out, "is_special_interest") <- is_special_interest
  out
}

#' Batch fetch for a set of TFs
#'
#' @return named list tf -> data.frame (or NULL if no peaks in window).
#' @param tfs (see function body).
#' @param gene_info (see function body).
#' @param promoter_info (see function body).
#' @param upstream (see function body).
#' @param downstream (see function body).
#' @param threshold (see function body).
#' @param max_experiments (see function body).
#' @param special_interest_gene (see function body).
#' @param special_interest_cap (see function body).
#' @param quiet (see function body).
#' @export
run_chipatlas_scan <- function(tfs, gene_info, promoter_info,
                               upstream = 2500, downstream = 500,
                               threshold = "05",
                               max_experiments = 50,
                               special_interest_gene = NULL,
                               special_interest_cap  = NULL,
                               quiet = TRUE) {
  if (!length(tfs)) return(list())
  # Normalize special-interest list for case-insensitive matching.
  sig_set <- if (is.null(special_interest_gene)) character(0)
             else toupper(unique(trimws(special_interest_gene[
               nzchar(special_interest_gene)])))
  sig_cap_label <- if (is.null(special_interest_cap)) "no cap"
                   else paste0("cap=", special_interest_cap, " SRX/TF")
  message("  ChIP-Atlas scan: ", length(tfs), " TFs",
          " | threshold=", threshold,
          " | cap=", max_experiments, " SRX/TF",
          if (length(sig_set))
            paste0(" | special-interest (", sig_cap_label, "): ",
                   paste(sig_set, collapse = ", "))
          else "")
  # Prime the experimentList cache once so progress lines don't interleave.
  # Note: .chipatlas_experiment_list(verbose = TRUE) prints parse diagnostics
  # exactly once per R session (it caches inside a local()).
  el <- .chipatlas_experiment_list()
  # Tell the user how many SRX rows we're searching against and sanity-check
  # the first requested TF; "0 peaks" for every TF is almost always a parse
  # or filter failure, not genuinely missing data.
  message("    hg38 TF experiments loaded: ",
          format(nrow(el), big.mark = ","),
          " | unique antigens: ",
          format(length(unique(el$antigen)), big.mark = ","))
  if (length(tfs) > 0 && nrow(el) > 0) {
    first <- tfs[[1]]
    n_first <- sum(toupper(el$antigen) == toupper(first))
    message("    Probe: ", first, " -> ", n_first, " SRX in experimentList")
  }
  out <- vector("list", length(tfs)); names(out) <- tfs
  for (tf in tfs) {
    is_sig <- toupper(tf) %in% sig_set
    pk <- tryCatch(
      fetch_chipatlas_peaks(tf, gene_info, promoter_info,
                            upstream = upstream, downstream = downstream,
                            threshold = threshold,
                            max_experiments = max_experiments,
                            is_special_interest = is_sig,
                            special_interest_cap = special_interest_cap,
                            quiet = quiet),
      error = function(e) {
        warning("ChIP-Atlas fetch failed for ", tf, ": ", conditionMessage(e))
        NULL
      })
    n <- if (is.null(pk)) 0L else nrow(pk)
    n_srx <- if (is.null(pk) || nrow(pk) == 0) 0L
             else length(unique(pk$srx))
    n_scanned <- if (is.null(pk)) 0L
                 else (attr(pk, "n_srx_scanned") %||% NA_integer_)
    sig_tag <- if (is_sig) "  [SPECIAL-INTEREST: full scan]" else ""
    message(sprintf("    %s: %d peaks across %d / %s experiment(s)%s",
                    tf, n, n_srx,
                    if (is.na(n_scanned)) "?" else as.character(n_scanned),
                    sig_tag))
    out[[tf]] <- pk
  }
  out
}

# =============================================================================
# Histone-marks: per-locus fetch + cell-type filtering
# =============================================================================

#' Loose cell-type matching for ChIP-Atlas's `cell_type` field.
#'
#' ChIP-Atlas's metadata uses inconsistent spellings for the same cell line —
#' "HEK293T", "HEK 293T", "HEK-293T", "HEK293", "Hek 293/T", "293T" all refer
#' to the same line, depending on which submitter typed the metadata. This
#' helper normalises both sides (uppercase, strip non-alphanumerics) and
#' tests substring equality in BOTH directions: a row matches if the target
#' is a substring of the field OR the field is a substring of the target.
#' That way "HEK293T" matches both "HEK293T" and "HEK293" (in case some
#' submissions dropped the trailing "T").
#' @param cell_type_field (see function body).
#' @param target (see function body).
#' @noRd
.chipatlas_celltype_match <- function(cell_type_field, target) {
  if (is.null(target) || is.na(target) || !nzchar(target))
    return(rep(TRUE, length(cell_type_field)))
  norm <- function(x) gsub("[^A-Z0-9]", "", toupper(as.character(x)))
  ct_norm <- norm(cell_type_field)
  tg_norm <- norm(target)
  # Bidirectional substring: target ⊆ field OR field ⊆ target. Skips empty
  # field strings (some ChIP-Atlas rows have no cell_type recorded).
  ok_field <- nzchar(ct_norm)
  matched <- rep(FALSE, length(cell_type_field))
  matched[ok_field] <- grepl(tg_norm, ct_norm[ok_field], fixed = TRUE) |
                       vapply(ct_norm[ok_field],
                              function(s) grepl(s, tg_norm, fixed = TRUE),
                              logical(1))
  matched
}

#' Fetch peaks for a single antigen, restricted to a custom SRX subset.
#'
#' Same machinery as fetch_chipatlas_peaks() but takes the SRX list as input
#' rather than deriving it via chipatlas_srx_for_tf(). Lets the histone-mark
#' fetch pre-filter SRXs by cell type (matched bucket) or take a top-N
#' newest slice (all-cell-types bucket) without re-querying experimentList.
#' @param srx_ids (see function body).
#' @param gene_info (see function body).
#' @param promoter_info (see function body).
#' @param upstream (see function body).
#' @param downstream (see function body).
#' @param threshold (see function body).
#' @param quiet (see function body).
#' @noRd
.fetch_chipatlas_peaks_for_srxs <- function(srx_ids, gene_info, promoter_info,
                                             upstream = 2500, downstream = 500,
                                             threshold = "05",
                                             quiet = TRUE) {
  if (!length(srx_ids)) return(NULL)
  win <- .chipatlas_window_coords(gene_info, promoter_info, upstream, downstream)
  el  <- .chipatlas_experiment_list()
  ct_lookup <- setNames(el$cell_type,       el$srx)
  cc_lookup <- setNames(el$cell_type_class, el$srx)
  pieces <- lapply(srx_ids, function(srx) {
    fpath <- download_chipatlas_srx_bed(srx, threshold = threshold, quiet = quiet)
    df <- .read_chipatlas_srx_bed(fpath)
    df <- .chipatlas_filter_to_window(df, win)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$srx <- srx
    df$cell_type       <- ct_lookup[[srx]] %||% NA_character_
    df$cell_type_class <- cc_lookup[[srx]] %||% NA_character_
    df
  })
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (!length(pieces)) return(NULL)
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out[, c("srx", "cell_type_class", "cell_type", "chr", "start", "end",
          "score", "pos_rel", "start_rel", "end_rel")]
}

#' Fetch histone-mark peaks at a locus, in two parallel buckets.
#'
#' For each mark, returns peaks twice:
#'   * `matched`: SRXs whose cell_type matches `cell_type` (loose match —
#'     case-insensitive, alphanumeric-stripped substring). Caps at
#'     `max_experiments_matched`.
#'   * `all`: top-N newest SRXs across all cell types. Caps at
#'     `max_experiments_all`.
#'
#' For the histone-marks page in the epigenetic deck, one row of bars is
#' rendered per mark per bucket — six rows of cell-type-matched on top, six
#' rows of all-cell-types-aggregated on bottom. See plot_histone_marks_locus()
#' in caspex_epigenetic.R.
#'
#' @return named list `list(matched = list(mark = peaks_df, ...),
#'                          all     = list(mark = peaks_df, ...),
#'                          n_srx_matched = list(mark = n, ...),
#'                          n_srx_all     = list(mark = n, ...))`.
#' @param marks (see function body).
#' @param gene_info (see function body).
#' @param promoter_info (see function body).
#' @param upstream (see function body).
#' @param downstream (see function body).
#' @param cell_type (see function body).
#' @param max_experiments_matched (see function body).
#' @param max_experiments_all (see function body).
#' @param threshold (see function body).
#' @param quiet (see function body).
#' @noRd
fetch_histone_peaks_for_locus <- function(
    marks,
    gene_info, promoter_info,
    upstream = 2500, downstream = 500,
    cell_type = "HEK293T",
    max_experiments_matched = 50,
    max_experiments_all     = 50,
    threshold = "05",
    quiet = TRUE) {
  el <- .chipatlas_experiment_list()
  matched_peaks <- vector("list", length(marks)); names(matched_peaks) <- marks
  all_peaks     <- vector("list", length(marks)); names(all_peaks)     <- marks
  # Two counts per mark per bucket:
  #   `n_srx_*`         = effective (post-cap) — what the union bar shows
  #   `n_srx_*_total`   = available (pre-cap)  — total in ChIP-Atlas
  # The plot label uses the effective count so the (n=…) annotation
  # matches what the bar represents; the console message shows both for
  # diagnostic purposes.
  n_srx_matched       <- setNames(integer(length(marks)), marks)
  n_srx_all           <- setNames(integer(length(marks)), marks)
  n_srx_matched_total <- setNames(integer(length(marks)), marks)
  n_srx_all_total     <- setNames(integer(length(marks)), marks)

  diag_done <- FALSE
  for (mark in marks) {
    # All SRXs for this antigen, newest-first (chipatlas_srx_for_tf sorts
    # by accession-number suffix, descending).
    all_srx <- chipatlas_srx_for_tf(mark)
    if (!length(all_srx)) {
      message("  [histone] ", mark, ": 0 SRX in experimentList")
      next
    }

    # ── Cell-type-matched bucket ────────────────────────────────────────
    # Match on BOTH cell_type AND cell_type_class because ChIP-Atlas
    # submitters fill these inconsistently — some put "HEK293" in
    # cell_type and "Embryonic Kidney" in cell_type_class, others put the
    # specific line label in cell_type_class. A row matches the user's
    # target if either field passes the bidirectional substring test.
    el_idx <- match(all_srx, el$srx)
    ct_vec <- el$cell_type[el_idx]
    cc_vec <- el$cell_type_class[el_idx]
    # One-shot diagnostic: dump the top cell_type / cell_type_class values
    # seen for the first mark so the user can verify spellings against
    # what they passed in. Helps explain n=0 surprises.
    if (!diag_done) {
      ct_top <- sort(table(ct_vec[nzchar(ct_vec)]), decreasing = TRUE)
      cc_top <- sort(table(cc_vec[nzchar(cc_vec)]), decreasing = TRUE)
      n_show <- min(8, length(ct_top))
      n_show_cc <- min(5, length(cc_top))
      message("  [histone diag] First mark (", mark, ") n=",
              length(all_srx), " SRX. Top cell_type spellings: ",
              paste(sprintf("%s=%d",
                            names(ct_top)[seq_len(n_show)],
                            as.integer(ct_top[seq_len(n_show)])),
                    collapse = ", "))
      message("  [histone diag] Top cell_type_class spellings: ",
              paste(sprintf("%s=%d",
                            names(cc_top)[seq_len(n_show_cc)],
                            as.integer(cc_top[seq_len(n_show_cc)])),
                    collapse = ", "))
      diag_done <- TRUE
    }
    matched_mask <- .chipatlas_celltype_match(ct_vec, cell_type) |
                    .chipatlas_celltype_match(cc_vec, cell_type)
    matched_srx  <- all_srx[matched_mask]
    # One-shot sanity check: on the first mark, dump up to 8 unique
    # (cell_type, cell_type_class) pairs that survived the match so the
    # user can verify the bidirectional matcher isn't grabbing false
    # positives. Empty / NA rows skipped.
    if (mark == marks[1] && length(matched_srx) > 0) {
      ct_match <- ct_vec[matched_mask]
      cc_match <- cc_vec[matched_mask]
      pairs    <- paste0(ifelse(nzchar(ct_match), ct_match, "(empty)"),
                          " | ",
                          ifelse(nzchar(cc_match), cc_match, "(empty)"))
      uniq_pairs <- unique(pairs)
      n_show <- min(8, length(uniq_pairs))
      message("  [histone diag] Matched cell_type | cell_type_class pairs (",
              length(uniq_pairs), " unique, showing ", n_show, "): ",
              paste(uniq_pairs[seq_len(n_show)], collapse = "; "))
    }
    n_srx_matched_total[mark] <- length(matched_srx)
    if (length(matched_srx) > max_experiments_matched)
      matched_srx <- matched_srx[seq_len(max_experiments_matched)]
    n_srx_matched[mark] <- length(matched_srx)
    if (length(matched_srx)) {
      matched_peaks[[mark]] <- .fetch_chipatlas_peaks_for_srxs(
        matched_srx, gene_info, promoter_info,
        upstream = upstream, downstream = downstream,
        threshold = threshold, quiet = quiet)
    }

    # ── All-cell-types bucket ───────────────────────────────────────────
    n_srx_all_total[mark] <- length(all_srx)
    all_srx_capped  <- if (length(all_srx) > max_experiments_all)
                         all_srx[seq_len(max_experiments_all)]
                       else all_srx
    n_srx_all[mark] <- length(all_srx_capped)
    if (length(all_srx_capped)) {
      all_peaks[[mark]] <- .fetch_chipatlas_peaks_for_srxs(
        all_srx_capped, gene_info, promoter_info,
        upstream = upstream, downstream = downstream,
        threshold = threshold, quiet = quiet)
    }

    message(sprintf(
      "  [histone] %s: matched %s n=%d (capped to %d), all n=%d (capped to %d)",
      mark, cell_type %||% "(no filter)",
      n_srx_matched_total[mark], n_srx_matched[mark],
      n_srx_all_total[mark],     n_srx_all[mark]))
  }

  list(matched              = matched_peaks,
       all                  = all_peaks,
       # Effective (post-cap) counts — what the bar union represents
       # and what plot_histone_marks_locus puts in the row label.
       n_srx_matched        = n_srx_matched,
       n_srx_all            = n_srx_all,
       # Available (pre-cap) counts — total SRX in ChIP-Atlas for the
       # mark in this cell-type filter; useful for downstream summary
       # and for the console diagnostic.
       n_srx_matched_total  = n_srx_matched_total,
       n_srx_all_total      = n_srx_all_total,
       cell_type            = cell_type,
       marks                = marks)
}

#' Utility: clear the on-disk ChIP-Atlas cache
#'
#' @param keep_experiment_list if TRUE (default), preserve the ~300 MB
#'   experimentList.tab and only purge the per-SRX BED cache. Useful between
#'   debug runs so you don't re-download the big metadata file every time.
#' @export
clear_chipatlas_cache <- function(keep_experiment_list = TRUE) {
  dir <- .chipatlas_cache_dir()
  if (keep_experiment_list) {
    # Purge every subdirectory (e.g. hg38/bed05/*.bed) but leave
    # experimentList.tab at the top level untouched.
    kids <- list.files(dir, full.names = TRUE, include.dirs = TRUE)
    kids <- kids[basename(kids) != "experimentList.tab"]
    for (k in kids) unlink(k, recursive = TRUE)
    message("ChIP-Atlas per-experiment BED cache cleared: ", dir)
    message("  (experimentList.tab preserved; pass keep_experiment_list=FALSE ",
            "to also purge that)")
  } else {
    unlink(dir, recursive = TRUE)
    message("ChIP-Atlas cache cleared (including experimentList.tab): ", dir)
  }
  invisible(NULL)
}

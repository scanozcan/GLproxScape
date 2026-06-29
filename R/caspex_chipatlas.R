# =============================================================================
# caspex_chipatlas.R \u2014 ChIP-Atlas peak backend for CasPEX
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

# ---- species -> ChIP-Atlas genome assembly mapping --------------------------
#
# ChIP-Atlas publishes its per-SRX BEDs and the experimentList.tab under
# per-genome subtrees keyed by UCSC-style assembly codes. We expose a small
# Ensembl-species -> ChIP-Atlas genome lookup so callers can pass the same
# `species` arg they already give run_caspex (mus_musculus, homo_sapiens, ...)
# without thinking about UCSC assembly codes. Returns NULL for unknown species
# so the caller can surface a clean error rather than silently fetching the
# wrong genome's BEDs.
#
# Default mapping is the LATEST ChIP-Atlas-supported assembly per species; if
# you need a specific older assembly (hg19, mm9), pass `chipatlas_genome`
# explicitly to run_caspex() to override the auto-derivation.
.SPECIES_TO_CHIPATLAS_GENOME <- list(
  homo_sapiens               = "hg38",
  mus_musculus               = "mm10",
  rattus_norvegicus          = "rn7",
  drosophila_melanogaster    = "dm6",
  caenorhabditis_elegans     = "ce11",
  danio_rerio                = "danRer11",
  gallus_gallus              = "galGal6",
  saccharomyces_cerevisiae   = "sacCer3",
  schizosaccharomyces_pombe  = "spo2",
  arabidopsis_thaliana       = "tair10",
  oryza_sativa               = "msu7"
)

#' Map an Ensembl species string to a ChIP-Atlas genome assembly code.
#'
#' @param species Ensembl species token (e.g. "homo_sapiens", "mus_musculus").
#'   Case-insensitive; underscores or spaces tolerated. Pass NULL to get back
#'   NULL.
#' @return character genome code (e.g. "hg38") or NULL for unknown species.
#' @noRd
.species_to_chipatlas_genome <- function(species) {
  if (is.null(species) || !nzchar(species)) return(NULL)
  key <- tolower(gsub("[[:space:]]+", "_", trimws(species)))
  .SPECIES_TO_CHIPATLAS_GENOME[[key]]
}

# ---- UCSC <-> Ensembl assembly name map ------------------------------------
#
# Used to detect assembly mismatches between Ensembl REST (returns Ensembl
# assembly names like "GRCm39") and ChIP-Atlas (uses UCSC codes like "mm10").
# When a run has Ensembl-frame gene_info but is asking for ChIP-Atlas data
# in a different frame, peaks systematically miss the window. The lift-over
# helpers below resolve this.
.UCSC_TO_ENSEMBL_ASSEMBLY <- list(
  hg38     = "GRCh38",
  hg19     = "GRCh37",
  mm10     = "GRCm38",
  mm9      = "NCBIM37",
  mm39     = "GRCm39",
  rn7      = "mRatBN7.2",
  rn6      = "Rnor_6.0",
  dm6      = "BDGP6.46",
  ce11     = "WBcel235",
  danRer11 = "GRCz11",
  galGal6  = "GRCg6a",
  sacCer3  = "R64-1-1"
)

# Per-session cache of Ensembl REST /info/assembly/<species> responses.
.ENSEMBL_ASSEMBLY_CACHE <- new.env(parent = emptyenv())

#' Query Ensembl REST for the current default assembly for a species.
#'
#' @param species Ensembl species token (e.g. "mus_musculus").
#' @return character Ensembl assembly name (e.g. "GRCm39"), or NA if the
#'   query fails / species unknown.
#' @noRd
.ensembl_assembly_for_species <- function(species) {
  if (is.null(species) || !nzchar(species)) return(NA_character_)
  if (exists(species, envir = .ENSEMBL_ASSEMBLY_CACHE, inherits = FALSE))
    return(get(species, envir = .ENSEMBL_ASSEMBLY_CACHE, inherits = FALSE))
  js <- tryCatch(ensembl_get(paste0("/info/assembly/", species)),
                 error = function(e) NULL)
  asm <- if (!is.null(js) && !is.null(js$assembly_name))
           as.character(js$assembly_name) else NA_character_
  assign(species, asm, envir = .ENSEMBL_ASSEMBLY_CACHE)
  asm
}

# ---- UCSC chain-file fetch + cache (rtracklayer fallback path) -------------
#
# UCSC publishes pre-computed chain files for assembly pairs. We cache them
# in R_user_dir("caspex","cache")/liftover/ so repeated lift-overs reuse
# one disk copy. Only used when an Ensembl archive endpoint isn't known
# for the assembly pair (the archive route is the primary lift-over path).

#' @noRd
.chain_cache_dir <- function() {
  base <- tools::R_user_dir("caspex", which = "cache")
  dir  <- file.path(base, "liftover")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  dir
}

#' @noRd
.download_chain_file <- function(from_ucsc, to_ucsc, quiet = FALSE) {
  to_cap <- paste0(toupper(substr(to_ucsc, 1, 1)), substring(to_ucsc, 2))
  fname  <- sprintf("%sTo%s.over.chain.gz", from_ucsc, to_cap)
  fpath  <- file.path(.chain_cache_dir(), fname)
  if (file.exists(fpath) && file.size(fpath) > 1024) return(fpath)
  url <- sprintf("https://hgdownload.soe.ucsc.edu/goldenpath/%s/liftOver/%s",
                 from_ucsc, fname)
  if (!quiet) message("  Downloading UCSC chain file: ", url)
  res <- tryCatch(
    utils::download.file(url, fpath, mode = "wb", quiet = quiet),
    error = function(e) -1L)
  if (res != 0L || !file.exists(fpath) || file.size(fpath) < 1024) {
    if (file.exists(fpath)) file.remove(fpath)
    return(NULL)
  }
  fpath
}

#' Lift-over a gene_info list from one assembly to another via rtracklayer.
#'
#' Fallback path — used only when the Ensembl archive route is unavailable.
#' Returns gene_info unchanged on failure (missing deps, parse error, etc.).
#' @noRd
.liftover_gene_info <- function(gene_info, from_ucsc, to_ucsc) {
  if (from_ucsc == to_ucsc) return(gene_info)
  if (!requireNamespace("rtracklayer", quietly = TRUE) ||
      !requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("IRanges", quietly = TRUE)) {
    warning("Lift-over needs rtracklayer + GenomicRanges; ",
            "BiocManager::install(c('rtracklayer','GenomicRanges')). ",
            "Skipping lift-over ", from_ucsc, " -> ", to_ucsc,
            "; ChIP-Atlas peaks may not align with gene_info coords.")
    return(gene_info)
  }
  chain_gz <- .download_chain_file(from_ucsc, to_ucsc, quiet = TRUE)
  if (is.null(chain_gz)) {
    warning("Could not fetch UCSC chain file ", from_ucsc, " -> ", to_ucsc,
            "; skipping lift-over.")
    return(gene_info)
  }
  chain_path <- sub("\\.gz$", "", chain_gz)
  if (!file.exists(chain_path) || file.size(chain_path) < 1024) {
    ok <- tryCatch({
      con_in  <- gzfile(chain_gz, "rb")
      con_out <- file(chain_path, "wb")
      on.exit({ close(con_in); close(con_out) }, add = TRUE)
      while (length(buf <- readBin(con_in, "raw", 65536)) > 0)
        writeBin(buf, con_out)
      TRUE
    }, error = function(e) FALSE)
    if (!isTRUE(ok)) {
      warning("Failed to decompress chain file ", chain_gz,
              "; skipping lift-over.")
      return(gene_info)
    }
  }
  chain <- tryCatch(rtracklayer::import.chain(chain_path),
                    error = function(e) {
                      warning("rtracklayer::import.chain failed: ",
                              conditionMessage(e),
                              "; skipping lift-over.")
                      NULL
                    })
  if (is.null(chain)) return(gene_info)

  chr <- gene_info$chr
  if (!grepl("^chr", chr)) chr <- paste0("chr", chr)
  str <- if (gene_info$strand == 1) "+" else "-"
  to_lift <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges   = IRanges::IRanges(
      start = c(gene_info$tss, gene_info$start, gene_info$end),
      end   = c(gene_info$tss, gene_info$start, gene_info$end)),
    strand   = str)
  names(to_lift) <- c("tss", "start", "end")
  lifted <- tryCatch(rtracklayer::liftOver(to_lift, chain),
                     error = function(e) {
                       warning("rtracklayer::liftOver failed: ",
                               conditionMessage(e),
                               "; skipping lift-over.")
                       NULL
                     })
  if (is.null(lifted)) return(gene_info)

  pick <- function(grl, key) {
    g <- grl[[key]]
    if (length(g) == 0) return(NA_integer_)
    as.integer(GenomicRanges::start(g)[1])
  }
  new_tss   <- pick(lifted, "tss")
  new_start <- pick(lifted, "start")
  new_end   <- pick(lifted, "end")
  if (is.na(new_tss)) {
    warning("Lift-over TSS mapping empty for ", gene_info$name,
            " (", from_ucsc, " -> ", to_ucsc,
            "); keeping original coords. ChIP-Atlas overlay may misalign.")
    return(gene_info)
  }

  out <- gene_info
  out$tss   <- new_tss
  if (!is.na(new_start)) out$start <- new_start
  if (!is.na(new_end))   out$end   <- new_end
  message(sprintf(
    "  Lifted gene_info %s -> %s: TSS %d -> %d (delta %+d bp)",
    from_ucsc, to_ucsc, gene_info$tss, out$tss, out$tss - gene_info$tss))
  out
}

# ---- Ensembl archive endpoints for older assemblies (PRIMARY lift-over) ----
#
# Ensembl mirrors prior releases at versioned REST subdomains. Querying the
# archive returns coordinates in the OLD assembly's frame — for assemblies
# where an archive REST endpoint exists, this is faster and more reliable
# than UCSC chain-file lift-over. NOTE: the *.archive.ensembl.org sites are
# the WEB interface only; REST lives at e<release>.rest.ensembl.org or
# grch37.rest.ensembl.org.
.UCSC_TO_ENSEMBL_ARCHIVE <- list(
  hg19 = "https://grch37.rest.ensembl.org",  # standalone GRCh37 REST
  mm10 = "https://e102.rest.ensembl.org"     # e102 (Nov 2020) = last GRCm38
)

#' Re-fetch a gene from an Ensembl archive in a specific assembly's frame.
#'
#' Mirrors lookup_gene() but hits the archive base URL. Tries to preserve
#' the same transcript_id when possible so the TSS anchor stays consistent.
#'
#' @return gene_info list with archive-frame coords, or NULL on failure.
#' @noRd
.lookup_gene_in_archive <- function(gene_info, archive_base_url,
                                     verbose = FALSE) {
  if (!requireNamespace("httr",     quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE))
    return(NULL)
  species <- gene_info$species
  gene    <- gene_info$name
  if (is.null(species) || is.null(gene)) return(NULL)
  url <- paste0(archive_base_url, "/lookup/symbol/", species, "/", gene,
                "?expand=1")
  if (verbose) message("  [archive lookup] GET ", url)
  res <- tryCatch(
    httr::GET(url, httr::add_headers(Accept = "application/json"),
              httr::timeout(20)),
    error = function(e) {
      if (verbose) message("  [archive lookup] httr error: ",
                            conditionMessage(e))
      NULL
    })
  if (is.null(res)) return(NULL)
  status <- httr::status_code(res)
  if (status != 200) {
    if (verbose) message("  [archive lookup] HTTP ", status, " from ", url)
    return(NULL)
  }
  js <- tryCatch(httr::content(res, "parsed", simplifyVector = FALSE),
                 error = function(e) NULL)
  if (is.null(js) || is.null(js$start) || is.null(js$end)) {
    if (verbose) message("  [archive lookup] response missing start/end fields")
    return(NULL)
  }

  # Pick a transcript: same transcript_id, then is_canonical, then longest.
  tx_target <- gene_info$transcript_id
  picked_tx <- NULL
  if (!is.null(js$Transcript) && length(js$Transcript)) {
    for (tx in js$Transcript) {
      if (!is.null(tx_target) && identical(tx$id, tx_target)) {
        picked_tx <- tx; break
      }
    }
    if (is.null(picked_tx)) {
      canon_ix <- which(vapply(js$Transcript,
                                function(t) isTRUE(t$is_canonical == 1L),
                                logical(1)))
      if (length(canon_ix)) {
        picked_tx <- js$Transcript[[canon_ix[1]]]
      } else {
        spans <- vapply(js$Transcript,
                        function(t) as.integer(t$end) - as.integer(t$start),
                        integer(1))
        picked_tx <- js$Transcript[[which.max(spans)]]
        if (verbose) message(
          "  [archive lookup] no is_canonical / no matching transcript_id; ",
          "picked longest transcript: ", picked_tx$id,
          " (span ", max(spans), " bp)")
      }
    }
  }
  strand <- if (!is.null(js$strand)) as.integer(js$strand) else gene_info$strand
  out <- gene_info
  out$start <- as.integer(js$start)
  out$end   <- as.integer(js$end)
  if (!is.null(picked_tx)) {
    out$transcript_id <- picked_tx$id
    out$tss <- if (strand == 1) as.integer(picked_tx$start)
               else as.integer(picked_tx$end)
  } else {
    out$tss <- if (strand == 1) as.integer(js$start) else as.integer(js$end)
  }
  out
}

# Per-session cache for lifted gene_info. Keyed by (transcript_id|name,
# chipatlas_genome). Lift-over involves an Ensembl archive REST call;
# without this cache, every fetch_chipatlas_peaks() call would hit the
# network for the same gene — ~100 redundant calls per typical deck.
.LIFTOVER_CACHE <- new.env(parent = emptyenv())

#' Decide whether to lift over gene_info for a given chipatlas_genome,
#' and emit a clear warning when the Ensembl assembly and the
#' chipatlas_genome's expected Ensembl name disagree.
#'
#' Strategy:
#'   1. Ensembl archive re-lookup (fast, no extra deps; works for assemblies
#'      where an archive base URL is known).
#'   2. UCSC chain-file lift-over via rtracklayer (fallback).
#'   3. Original gene_info with a warning when both fail.
#'
#' Results are cached per session so the archive REST call and the
#' mismatch warning each fire exactly once per (gene, target_genome).
#'
#' @noRd
.maybe_liftover_for_chipatlas <- function(gene_info, chipatlas_genome) {
  if (is.null(chipatlas_genome) || is.null(gene_info$species))
    return(gene_info)
  ens_asm     <- .ensembl_assembly_for_species(gene_info$species)
  expected_ens <- .UCSC_TO_ENSEMBL_ASSEMBLY[[chipatlas_genome]]
  if (is.na(ens_asm) || is.null(expected_ens)) return(gene_info)
  # Strip Ensembl patch-version suffix (.p13, .p14, ...) before comparing.
  # Ensembl returns "GRCh38.p14" for current human releases but
  # .UCSC_TO_ENSEMBL_ASSEMBLY stores the un-patched assembly name "GRCh38".
  # Without this strip the comparison falsely reports a mismatch for every
  # human run, triggers an unnecessary lift-over branch, and crashes with
  # "subscript out of bounds" inside the archive / rtracklayer fallback --
  # silently swallowed upstream as a NULL return from fetch_chipatlas_peaks,
  # producing the "0 peaks across 0 / 0 experiment(s)" pattern for every TF.
  # Mouse paths are unaffected because GRCm39 / GRCm38 don't carry patch
  # suffixes.
  ens_asm_base <- sub("\\.p[0-9]+$", "", ens_asm)
  if (identical(ens_asm_base, expected_ens)) return(gene_info)

  cache_key <- paste0(gene_info$transcript_id %||% gene_info$name, "__",
                      chipatlas_genome)
  if (exists(cache_key, envir = .LIFTOVER_CACHE, inherits = FALSE))
    return(get(cache_key, envir = .LIFTOVER_CACHE, inherits = FALSE))

  warning(sprintf(
    "ChIP-Atlas coordinate frame mismatch: Ensembl returned %s coordinates ",
    ens_asm),
    sprintf("for %s, but chipatlas_genome='%s' expects %s. Lifting over via ",
            gene_info$name, chipatlas_genome, expected_ens),
    "Ensembl archive (preferred) or rtracklayer chain file (fallback).",
    call. = FALSE)

  # --- Strategy 1: Ensembl archive lookup ---------------------------------
  archive_url <- .UCSC_TO_ENSEMBL_ARCHIVE[[chipatlas_genome]]
  if (!is.null(archive_url)) {
    lifted <- .lookup_gene_in_archive(gene_info, archive_url, verbose = TRUE)
    if (!is.null(lifted) && !is.null(lifted$tss) &&
        is.finite(lifted$tss) && lifted$tss != gene_info$tss) {
      message(sprintf(
        "  Re-fetched gene_info from Ensembl archive (%s -> %s): TSS %d -> %d (delta %+d bp)",
        ens_asm, expected_ens, gene_info$tss, lifted$tss,
        lifted$tss - gene_info$tss))
      assign(cache_key, lifted, envir = .LIFTOVER_CACHE)
      return(lifted)
    }
    if (!is.null(lifted) && !is.null(lifted$tss) && lifted$tss == gene_info$tss) {
      message("  Ensembl archive returned same TSS as current release; ",
              "no lift-over needed for ", gene_info$name, ".")
      assign(cache_key, lifted, envir = .LIFTOVER_CACHE)
      return(lifted)
    }
    warning("Ensembl archive lookup at ", archive_url, " failed; ",
            "falling back to rtracklayer chain file.")
  }

  # --- Strategy 2: UCSC chain via rtracklayer (fallback) ------------------
  ens_to_ucsc <- setNames(names(.UCSC_TO_ENSEMBL_ASSEMBLY),
                           unlist(.UCSC_TO_ENSEMBL_ASSEMBLY))
  from_ucsc <- ens_to_ucsc[[ens_asm]]
  if (is.null(from_ucsc)) {
    warning("No UCSC code known for Ensembl assembly '", ens_asm,
            "'; cannot lift-over. Original coords retained.")
    assign(cache_key, gene_info, envir = .LIFTOVER_CACHE)
    return(gene_info)
  }
  result <- .liftover_gene_info(gene_info, from_ucsc, chipatlas_genome)
  assign(cache_key, result, envir = .LIFTOVER_CACHE)
  result
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
#' @examples
#' \dontrun{
#' exp <- download_chipatlas_experiment_list()
#' }
#' @export
download_chipatlas_experiment_list <- function(force = FALSE, quiet = FALSE) {
  fpath <- file.path(.chipatlas_cache_dir(), "experimentList.tab")
  # Sanity check existing cache. ChIP-Atlas experimentList.tab is ~300+ MB and
  # each line starts with an experiment accession \u2014 SRX (SRA), DRX (DDBJ), or
  # ERX (ENA) \u2014 followed by a tab. A previous failed download could leave an
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
  # Two-level cache:
  #   `cached_unfiltered` — the full parsed experimentList. Built once per
  #     R session (file is 300+ MB, parse is ~3-5 s).
  #   `cached_filtered`   — named list keyed by genome assembly, each entry
  #     the genome-filtered subset. Built lazily on first call per genome.
  # The filter message ("mm10 TF + Histone rows after filter: ...") fires
  # exactly once per (session, genome) combination, not once per call from
  # the inner loop — that previously flooded the console.
  cached_unfiltered <- NULL
  cached_filtered   <- list()
  function(force_reload = FALSE, verbose = TRUE,
           genome = "hg38") {
    if (!is.null(cached_filtered[[genome]]) && !force_reload)
      return(cached_filtered[[genome]])
    if (is.null(cached_unfiltered) || force_reload) {
      # First call this session — go fetch + parse.
      ## continues below
    } else {
      # Cache hit on the parse, miss on the genome subset. Filter without
      # printing the parse diagnostics.
      cached_filtered[[genome]] <<-
        .filter_chipatlas_for_run(cached_unfiltered, genome, verbose = verbose)
      return(cached_filtered[[genome]])
    }
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
    # Cache the full parsed table; per-genome filtered subsets are built
    # lazily by .filter_chipatlas_for_run() and cached in `cached_filtered`.
    cached_unfiltered <<- df
    cached_filtered[[genome]] <<-
      .filter_chipatlas_for_run(df, genome, verbose = verbose)
    cached_filtered[[genome]]
  }
})

#' Filter the parsed experimentList by genome assembly + antigen class.
#'
#' Keep antigen classes for TFs and Histone marks; drop DNase / ATAC /
#' Bisulfite / etc. Try several spellings because ChIP-Atlas has been known
#' to shift labels between releases.
#' @param df parsed experimentList data.frame.
#' @param genome ChIP-Atlas genome assembly code (e.g. "hg38", "mm10").
#' @param verbose whether to print the filter result.
#' @noRd
.filter_chipatlas_for_run <- function(df, genome = "hg38", verbose = FALSE) {
  keep_class_patterns <- c("TFs and others", "TFs_and_others",
                            "TF", "TFs", "Transcription factor",
                            "Histone", "Histones", "Histone modification")
  keep <- df$genome == genome & df$antigen_class %in% keep_class_patterns
  if (verbose) message("    ", genome, " TF + Histone rows after filter: ",
                       format(sum(keep), big.mark = ","))
  if (sum(keep) == 0) {
    warning("ChIP-Atlas: 0 rows matched ", genome, " + (TF | Histone) filter. ",
            "Check the antigen_class breakdown above and update ",
            "keep_class_patterns in caspex_chipatlas.R, or verify that ",
            "ChIP-Atlas has data for the '", genome, "' assembly.")
  }
  df[keep, , drop = FALSE]
}

#' Experiment IDs for a given TF (antigen symbol) on a specified genome.
#'
#' @param tf HGNC / MGI symbol; matched case-insensitively against `antigen`.
#' @param genome ChIP-Atlas genome assembly code (default "hg38").
#' @return character vector of experiment IDs (SRX/DRX/ERX), ordered
#'   newest-first by numeric suffix within each prefix.
#' @noRd
chipatlas_srx_for_tf <- function(tf, genome = "hg38") {
  el <- .chipatlas_experiment_list(genome = genome)
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
#' @param srx ChIP-Atlas experiment accession (SRX identifier).
#' @param genome ChIP-Atlas genome assembly, e.g. "hg38" or "mm10".
#' @param threshold ChIP-Atlas peak significance threshold as a string (e.g. "05" for the q < 1e-5 track).
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
#' @param srx ChIP-Atlas experiment accession (SRX identifier).
#' @param genome ChIP-Atlas genome assembly, e.g. "hg38" or "mm10".
#' @param threshold ChIP-Atlas peak significance threshold as a string (e.g. "05" for the q < 1e-5 track).
#' @param force Logical; force re-download and bypass the cache when TRUE.
#' @param quiet Logical; suppress progress messages when TRUE.
#' @param max_retries Maximum number of HTTP retry attempts on transient failures.
#' @param timeout_sec HTTP request timeout in seconds.
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

# Per-SRX BEDs are narrowPeak-style (BED8+) \u2014 observed first line:
#   chr1 10018 10209 SRX11664714.05_peak_1 416 . 20.20613 47.189
# i.e. col4 is peak name (text), col5 is score (integer), col7 is signalValue.
# We only need chr/start/end for window overlap; optional score goes into
# the returned data.frame purely for downstream convenience.
#
# `read.table` with `col.names`+`colClasses` mis-typed the 4th field as
# numeric and silently dropped every row. Line-based parse avoids that.
#' Parse a cached per-SRX BED file into a data.frame.
#' @param fpath File-system path to read or write.
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
#' @param gene_info Gene/transcript coordinate record returned by lookup_gene() (chromosome, strand, TSS, assembly).
#' @param promoter_info Promoter-window sequence and coordinate record returned by fetch_promoter_seq().
#' @param upstream Basepairs upstream of the TSS included in the analysis window.
#' @param downstream Basepairs downstream of the TSS included in the analysis window.
#' @noRd
.chipatlas_window_coords <- function(gene_info, promoter_info,
                                      upstream, downstream,
                                      chipatlas_genome = NULL) {
  # Resolve assembly-frame mismatch before computing the window. If
  # gene_info is in (say) GRCm39 but chipatlas_genome="mm10" (GRCm38),
  # the window must be lifted or peaks will miss. Returns gene_info
  # unchanged when no mismatch / no chain or archive route available.
  if (!is.null(chipatlas_genome))
    gene_info <- .maybe_liftover_for_chipatlas(gene_info, chipatlas_genome)
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
#' @param df Data.frame of ChIP-Atlas peaks to filter to the window.
#' @param win Analysis-window coordinates (TSS-relative start and end in bp).
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
                                  quiet = TRUE,
                                  genome = "hg38") {
  srx_ids <- chipatlas_srx_for_tf(tf, genome = genome)
  if (!length(srx_ids)) return(NULL)
  # Cap selection logic:
  #   - Standard TF: take the newest `max_experiments` SRX.
  #   - Special-interest TF, special_interest_cap = NULL: scan ALL SRX
  #     (could be 2000+ for CTCF \u2014 slow on first run, free thereafter).
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
  win <- .chipatlas_window_coords(gene_info, promoter_info, upstream, downstream,
                                   chipatlas_genome = genome)

  el  <- .chipatlas_experiment_list(genome = genome)
  ct_lookup <- setNames(el$cell_type, el$srx)
  cc_lookup <- setNames(el$cell_type_class, el$srx)

  pieces <- lapply(srx_ids, function(srx) {
    fpath <- download_chipatlas_srx_bed(srx, genome = genome,
                                         threshold = threshold, quiet = quiet)
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
#' @param tfs Character vector of transcription-factor symbols to include.
#' @param gene_info Gene/transcript coordinate record returned by lookup_gene() (chromosome, strand, TSS, assembly).
#' @param promoter_info Promoter-window sequence and coordinate record returned by fetch_promoter_seq().
#' @param upstream Basepairs upstream of the TSS included in the analysis window.
#' @param downstream Basepairs downstream of the TSS included in the analysis window.
#' @param threshold ChIP-Atlas peak significance threshold as a string (e.g. "05" for the q < 1e-5 track).
#' @param max_experiments Maximum number of ChIP-Atlas experiments (SRX) fetched per TF.
#' @param special_interest_gene Factor(s) exempted from the chipatlas_max_experiments cap.
#' @param special_interest_cap Per-factor SRX cap for special-interest genes; NULL means uncapped.
#' @param quiet Logical; suppress progress messages when TRUE.
#' @param genome ChIP-Atlas genome assembly code (e.g. "hg38", "mm10").
#'   Defaults to "hg38" for backward compatibility. Use the
#'   `chipatlas_genome` argument on \code{run_caspex()} to auto-derive
#'   this from the run's `species` argument (recommended).
#' @examples
#' \dontrun{
#' gi  <- lookup_gene("TERT")
#' pi  <- fetch_promoter_seq(gi)
#' res <- run_chipatlas_scan(c("CTCF", "MAZ"), gi, pi)
#' }
#' @export
run_chipatlas_scan <- function(tfs, gene_info, promoter_info,
                               upstream = 2500, downstream = 500,
                               threshold = "05",
                               max_experiments = 50,
                               special_interest_gene = NULL,
                               special_interest_cap  = NULL,
                               quiet = TRUE,
                               genome = "hg38") {
  if (!length(tfs)) return(list())
  # Normalize special-interest list for case-insensitive matching.
  sig_set <- if (is.null(special_interest_gene)) character(0)
             else toupper(unique(trimws(special_interest_gene[
               nzchar(special_interest_gene)])))
  sig_cap_label <- if (is.null(special_interest_cap)) "no cap"
                   else paste0("cap=", special_interest_cap, " SRX/TF")
  message("  ChIP-Atlas scan: ", length(tfs), " TFs",
          " | genome=", genome,
          " | threshold=", threshold,
          " | cap=", max_experiments, " SRX/TF",
          if (length(sig_set))
            paste0(" | special-interest (", sig_cap_label, "): ",
                   paste(sig_set, collapse = ", "))
          else "")
  # Prime the experimentList cache once so progress lines don't interleave.
  # Note: .chipatlas_experiment_list(verbose = TRUE) prints parse diagnostics
  # exactly once per R session (it caches inside a local()).
  el <- .chipatlas_experiment_list(genome = genome)
  # Tell the user how many SRX rows we're searching against and sanity-check
  # the first requested TF; "0 peaks" for every TF is almost always a parse
  # or filter failure, not genuinely missing data.
  message("    ", genome, " TF experiments loaded: ",
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
                            quiet = quiet,
                            genome = genome),
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
#' @param cell_type_field Name of the column holding the cell-type label in the experiment table.
#' @param target Target cell-type label to match against ChIP-Atlas metadata.
#' @noRd
.chipatlas_celltype_match <- function(cell_type_field, target) {
  if (is.null(target) || is.na(target) || !nzchar(target))
    return(rep(TRUE, length(cell_type_field)))
  norm <- function(x) gsub("[^A-Z0-9]", "", toupper(as.character(x)))
  ct_norm <- norm(cell_type_field)
  tg_norm <- norm(target)
  # Bidirectional substring: target \u2286 field OR field \u2286 target. Skips empty
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
#' @param srx_ids Character vector of ChIP-Atlas SRX accessions.
#' @param gene_info Gene/transcript coordinate record returned by lookup_gene() (chromosome, strand, TSS, assembly).
#' @param promoter_info Promoter-window sequence and coordinate record returned by fetch_promoter_seq().
#' @param upstream Basepairs upstream of the TSS included in the analysis window.
#' @param downstream Basepairs downstream of the TSS included in the analysis window.
#' @param threshold ChIP-Atlas peak significance threshold as a string (e.g. "05" for the q < 1e-5 track).
#' @param quiet Logical; suppress progress messages when TRUE.
#' @noRd
.fetch_chipatlas_peaks_for_srxs <- function(srx_ids, gene_info, promoter_info,
                                             upstream = 2500, downstream = 500,
                                             threshold = "05",
                                             quiet = TRUE,
                                             genome = "hg38") {
  if (!length(srx_ids)) return(NULL)
  win <- .chipatlas_window_coords(gene_info, promoter_info, upstream, downstream,
                                   chipatlas_genome = genome)
  el  <- .chipatlas_experiment_list(genome = genome)
  ct_lookup <- setNames(el$cell_type,       el$srx)
  cc_lookup <- setNames(el$cell_type_class, el$srx)
  pieces <- lapply(srx_ids, function(srx) {
    fpath <- download_chipatlas_srx_bed(srx, genome = genome,
                                         threshold = threshold, quiet = quiet)
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
#' @param marks Character vector of histone modifications to query.
#' @param gene_info Gene/transcript coordinate record returned by lookup_gene() (chromosome, strand, TSS, assembly).
#' @param promoter_info Promoter-window sequence and coordinate record returned by fetch_promoter_seq().
#' @param upstream Basepairs upstream of the TSS included in the analysis window.
#' @param downstream Basepairs downstream of the TSS included in the analysis window.
#' @param cell_type ChIP-Atlas cell-type label to match or filter on.
#' @param max_experiments_matched Cap on matched-cell-type ChIP-Atlas experiments fetched per TF.
#' @param max_experiments_all Cap on all-cell-type ChIP-Atlas experiments fetched per TF.
#' @param threshold ChIP-Atlas peak significance threshold as a string (e.g. "05" for the q < 1e-5 track).
#' @param quiet Logical; suppress progress messages when TRUE.
#' @noRd
fetch_histone_peaks_for_locus <- function(
    marks,
    gene_info, promoter_info,
    upstream = 2500, downstream = 500,
    cell_type = "HEK293T",
    max_experiments_matched = 50,
    max_experiments_all     = 50,
    threshold = "05",
    quiet = TRUE,
    genome = "hg38") {
  el <- .chipatlas_experiment_list(genome = genome)
  matched_peaks <- vector("list", length(marks)); names(matched_peaks) <- marks
  all_peaks     <- vector("list", length(marks)); names(all_peaks)     <- marks
  # Two counts per mark per bucket:
  #   `n_srx_*`         = effective (post-cap) \u2014 what the union bar shows
  #   `n_srx_*_total`   = available (pre-cap)  \u2014 total in ChIP-Atlas
  # The plot label uses the effective count so the (n=\u2026) annotation
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
    all_srx <- chipatlas_srx_for_tf(mark, genome = genome)
    if (!length(all_srx)) {
      message("  [histone] ", mark, ": 0 SRX in experimentList")
      next
    }

    # \u2500\u2500 Cell-type-matched bucket \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
    # Match on BOTH cell_type AND cell_type_class because ChIP-Atlas
    # submitters fill these inconsistently \u2014 some put "HEK293" in
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
        threshold = threshold, quiet = quiet, genome = genome)
    }

    # \u2500\u2500 All-cell-types bucket \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
    n_srx_all_total[mark] <- length(all_srx)
    all_srx_capped  <- if (length(all_srx) > max_experiments_all)
                         all_srx[seq_len(max_experiments_all)]
                       else all_srx
    n_srx_all[mark] <- length(all_srx_capped)
    if (length(all_srx_capped)) {
      all_peaks[[mark]] <- .fetch_chipatlas_peaks_for_srxs(
        all_srx_capped, gene_info, promoter_info,
        upstream = upstream, downstream = downstream,
        threshold = threshold, quiet = quiet, genome = genome)
    }

    message(sprintf(
      "  [histone] %s: matched %s n=%d (capped to %d), all n=%d (capped to %d)",
      mark, cell_type %||% "(no filter)",
      n_srx_matched_total[mark], n_srx_matched[mark],
      n_srx_all_total[mark],     n_srx_all[mark]))
  }

  list(matched              = matched_peaks,
       all                  = all_peaks,
       # Effective (post-cap) counts \u2014 what the bar union represents
       # and what plot_histone_marks_locus puts in the row label.
       n_srx_matched        = n_srx_matched,
       n_srx_all            = n_srx_all,
       # Available (pre-cap) counts \u2014 total SRX in ChIP-Atlas for the
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
#' @return Invisibly \code{NULL}; called for its side effect of clearing the on-disk ChIP-Atlas cache.
#' @examples
#' \dontrun{
#' clear_chipatlas_cache()
#' }
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

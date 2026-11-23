# =============================================================================
# motif_cisbp.R -- CIS-BP motif engine (third motif source)
# =============================================================================
# Used when run_glproxscape(..., motif_search_engine = "cisbp"). CIS-BP
# (Weirauch et al. 2014; Lambert et al. 2019; Build 3.10 adds the CodeBook
# motifs of Jolma et al. 2026) has no REST API: the species archive
# (TF_Information*.txt + pwms_all_motifs/) is requested from the bulk-download
# form once and cached on disk, like the HOCOMOCO bundle.
#
# One motif per TF is chosen, by source tier (cisbp_prefer = "in_vitro"):
#   1. high-throughput in vitro motifs (HT-SELEX, GHT-SELEX, PBM, SMiLE-seq)
#   2. other direct motifs not imported from JASPAR / HOCOMOCO
#   3. motifs CIS-BP imported from JASPAR / HOCOMOCO
# and, within the best available tier, the medoid motif (highest mean
# similarity to the TF's other motifs). Probability matrices are converted to
# log-odds exactly as HOCOMOCO matrices, so motif_thresh means the same thing.
# =============================================================================

.cisbp_host <- "https://cisbp.ccbr.utoronto.ca"
.cisbp_in_vitro_types <- c("SELEX", "GHT-SELEX", "PBM", "SMiLE-seq")
.cisbp_imported_sources <- c("JASPAR", "HocoMoco")

#' Map a run_glproxscape() species (Ensembl style) to a CIS-BP species name.
#' @param species e.g. "homo_sapiens", "mus_musculus", "danio_rerio".
#' @return e.g. "Homo_sapiens".
#' @noRd
.cisbp_species <- function(species) {
  s <- tolower(species %||% "homo_sapiens")
  if (s %in% c("human")) s <- "homo_sapiens"
  if (s %in% c("mouse")) s <- "mus_musculus"
  paste0(toupper(substr(s, 1, 1)), substring(s, 2))
}

#' Local cache directory for one CIS-BP species archive.
#' @noRd
.cisbp_cache_dir <- function(species = "Homo_sapiens") {
  d <- file.path(tools::R_user_dir("caspex", which = "cache"),
                 sprintf("cisbp_%s", species))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

.cisbp_ready <- function(d) {
  file.exists(file.path(d, "TF_Information_all_motifs.txt")) &&
    dir.exists(file.path(d, "pwms_all_motifs"))
}

#' Download (or import) the CIS-BP species archive
#'
#' Requests the CIS-BP species archive (TF information + PWMs) from the
#' CIS-BP bulk-download form and unzips it into the local cache
#' (\code{tools::R_user_dir("caspex", "cache")}). Called automatically on
#' first use of \code{motif_search_engine = "cisbp"}.
#'
#' @param species CIS-BP species name (e.g. \code{"Homo_sapiens"},
#'   \code{"Mus_musculus"}) or an Ensembl-style name
#'   (\code{"homo_sapiens"}).
#' @param force Re-download even if the archive is already cached.
#' @param from_file Path to a manually downloaded archive (.zip from
#'   \url{https://cisbp.ccbr.utoronto.ca/bulk.php} with "TF info" and "PWMs"
#'   ticked) or an unzipped folder, used when automatic download is blocked.
#' @return The cache directory (invisibly).
#' @examples
#' \dontrun{
#' download_cisbp_bundle("Homo_sapiens")
#' }
#' @export
download_cisbp_bundle <- function(species = "Homo_sapiens", force = FALSE,
                                  from_file = NULL) {
  if (grepl("_", species) && species == tolower(species))
    species <- .cisbp_species(species)
  dest <- .cisbp_cache_dir(species)
  if (.cisbp_ready(dest) && !force) return(invisible(dest))

  if (!is.null(from_file)) {
    if (!file.exists(from_file)) stop("from_file does not exist: ", from_file)
    if (dir.exists(from_file)) {
      file.copy(list.files(from_file, full.names = TRUE), dest,
                recursive = TRUE, overwrite = TRUE)
    } else {
      utils::unzip(from_file, exdir = dest)
    }
    if (!.cisbp_ready(dest))
      stop("from_file does not contain TF_Information_all_motifs.txt and ",
           "pwms_all_motifs/: ", from_file)
    message("Cached CIS-BP ", species, " archive from ", from_file)
    return(invisible(dest))
  }

  message("Requesting CIS-BP species archive (", species, ")...")
  res <- tryCatch(
    httr::POST(paste0(.cisbp_host, "/bulk_archive.php"),
               body = list(selSpec = species,
                           `Spec[]` = "TF_Information",
                           `Spec[]` = "PWMs",
                           submit = "Download Species Archive!"),
               encode = "form", httr::timeout(120)),
    error = function(e) e)
  link <- NULL
  if (!inherits(res, "error") && httr::status_code(res) == 200) {
    html <- httr::content(res, as = "text", encoding = "UTF-8")
    m <- regmatches(html, regexpr('href="tmp/[^"]+\\.zip"', html))
    if (length(m) == 1) link <- sub('^href="', "", sub('"$', "", m))
  }
  if (!is.null(link)) {
    url <- paste0(.cisbp_host, "/", utils::URLencode(link))
    tmp <- tempfile(fileext = ".zip")
    old_timeout <- getOption("timeout")
    on.exit(options(timeout = old_timeout), add = TRUE)
    options(timeout = max(600, old_timeout))
    ok <- tryCatch({
      utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (ok && file.exists(tmp) && file.size(tmp) > 1e5) {
      utils::unzip(tmp, exdir = dest)
      unlink(tmp)
      if (.cisbp_ready(dest)) {
        message("  cached at ", dest)
        return(invisible(dest))
      }
    }
  }
  stop(sprintf(paste0(
    "Could not download the CIS-BP %s archive automatically.\n",
    "Open https://cisbp.ccbr.utoronto.ca/bulk.php, choose '%s' under 'By ",
    "Species', tick 'TF info' and 'PWMs', download the archive, then run\n",
    "  download_cisbp_bundle(species = \"%s\", from_file = \"/path/to/archive.zip\")"),
    species, species, species), call. = FALSE)
}

# In-memory cache of parsed CIS-BP tables, PPMs and motif choices.
.cisbp_env <- new.env(parent = emptyenv())

#' Read the CIS-BP motif table of one species (direct and inferred motifs).
#' @noRd
.load_cisbp_table <- function(species = "Homo_sapiens") {
  key <- paste0("tab_", species)
  if (!is.null(.cisbp_env[[key]])) return(.cisbp_env[[key]])
  d <- .cisbp_cache_dir(species)
  if (!.cisbp_ready(d)) download_cisbp_bundle(species)
  tab <- utils::read.delim(file.path(d, "TF_Information_all_motifs.txt"),
                           stringsAsFactors = FALSE, quote = "",
                           check.names = FALSE)
  tab <- tab[tab$Motif_ID != "." & tab$TF_Status %in% c("D", "I"),
             c("TF_Name", "Motif_ID", "TF_Status", "MSource_Identifier",
               "MSource_Type")]
  tab$TF_Name <- toupper(tab$TF_Name)
  tab <- unique(tab)
  message("CIS-BP ", species, ": ", length(unique(tab$TF_Name)), " TFs, ",
          length(unique(tab$TF_Name[tab$TF_Status == "D"])),
          " with a direct (experimental) motif")
  .cisbp_env[[key]] <- tab
  tab
}

#' Read one CIS-BP PWM file (Pos A C G T probabilities) into a 4 x L matrix.
#' @noRd
.cisbp_ppm <- function(motif_id, species) {
  key <- paste0("ppm_", species, "_", motif_id)
  if (exists(key, envir = .cisbp_env, inherits = FALSE))
    return(get(key, envir = .cisbp_env, inherits = FALSE))
  path <- file.path(.cisbp_cache_dir(species), "pwms_all_motifs",
                    paste0(motif_id, ".txt"))
  ppm <- NULL
  if (file.exists(path) && file.size(path) > 0) {
    df <- tryCatch(utils::read.delim(path, check.names = FALSE),
                   error = function(e) NULL)
    if (!is.null(df) && nrow(df) > 0 && all(c("A", "C", "G", "T") %in% names(df))) {
      ppm <- t(as.matrix(df[, c("A", "C", "G", "T")]))
      storage.mode(ppm) <- "double"
      rownames(ppm) <- c("A", "C", "G", "T")
      if (anyNA(ppm)) ppm <- NULL
    }
  }
  assign(key, ppm, envir = .cisbp_env)
  ppm
}

#' Similarity of two PPMs: best mean per-column Pearson correlation over
#' ungapped offsets (overlap >= min_overlap columns) and both strands.
#' @noRd
.ppm_similarity <- function(a, b, min_overlap = 5) {
  rc <- function(m) m[4:1, ncol(m):1, drop = FALSE]
  colcor <- function(A, B) {
    A <- sweep(A, 2, colMeans(A)); B <- sweep(B, 2, colMeans(B))
    r <- colSums(A * B) / sqrt(colSums(A^2) * colSums(B^2))
    r[!is.finite(r)] <- 0
    r
  }
  best <- -Inf
  la <- ncol(a)
  for (bb in list(b, rc(b))) {
    lb <- ncol(bb); ov <- min(min_overlap, la, lb)
    for (off in seq(-(lb - ov), la - ov)) {
      ia <- seq_len(la); ib <- ia - off
      k <- ib >= 1 & ib <= lb
      if (sum(k) < ov) next
      best <- max(best, mean(colcor(a[, ia[k], drop = FALSE],
                                    bb[, ib[k], drop = FALSE])))
    }
  }
  best
}

.ppm_ic <- function(ppm) {
  p <- pmax(sweep(ppm, 2, colSums(ppm), "/"), 1e-9)
  sum(2 + colSums(p * log2(p)))
}

#' Choose the CIS-BP motif used for one TF (tiered, then medoid / IC).
#' @noRd
.select_cisbp_motif <- function(tf_name, species = "Homo_sapiens",
                                status = "direct", exclude_sources = NULL,
                                pick = "medoid", prefer = "in_vitro") {
  memo <- paste("sel", species, toupper(tf_name), status,
                paste(sort(exclude_sources), collapse = "+"), pick, prefer,
                sep = "|")
  if (exists(memo, envir = .cisbp_env, inherits = FALSE))
    return(get(memo, envir = .cisbp_env, inherits = FALSE))
  tab  <- .load_cisbp_table(species)
  rows <- tab[tab$TF_Name == toupper(tf_name), , drop = FALSE]
  if (length(exclude_sources))
    rows <- rows[!rows$MSource_Identifier %in% exclude_sources, , drop = FALSE]
  cand <- rows[rows$TF_Status == "D", , drop = FALSE]
  if (nrow(cand) == 0 && status == "inferred")
    cand <- rows[rows$TF_Status == "I", , drop = FALSE]
  cand <- cand[!duplicated(cand$Motif_ID), , drop = FALSE]
  ppms <- lapply(cand$Motif_ID, .cisbp_ppm, species = species)
  ok   <- !vapply(ppms, is.null, logical(1))
  cand <- cand[ok, , drop = FALSE]; ppms <- ppms[ok]
  out <- NULL
  if (nrow(cand) > 0) {
    cand$tier <- ifelse(cand$MSource_Type %in% .cisbp_in_vitro_types, 1L,
                 ifelse(cand$MSource_Identifier %in% .cisbp_imported_sources, 3L, 2L))
    if (prefer == "in_vitro") {
      keep <- cand$tier == min(cand$tier)
      cand <- cand[keep, , drop = FALSE]; ppms <- ppms[keep]
    }
    if (nrow(cand) == 1) {
      i <- 1L
    } else if (pick == "ic") {
      i <- order(-vapply(ppms, .ppm_ic, numeric(1)), cand$Motif_ID)[1]
    } else {
      n <- length(ppms); S <- matrix(0, n, n)
      for (a in seq_len(n - 1)) for (b in (a + 1):n)
        S[a, b] <- S[b, a] <- .ppm_similarity(ppms[[a]], ppms[[b]])
      i <- order(-rowSums(S) / (n - 1), cand$Motif_ID)[1]
    }
    out <- list(motif_id = cand$Motif_ID[i], ppm = ppms[[i]],
                status = cand$TF_Status[i], source = cand$MSource_Identifier[i],
                source_type = cand$MSource_Type[i], tier = cand$tier[i],
                n_candidates = nrow(cand))
  }
  assign(memo, out, envir = .cisbp_env)
  out
}

#' Fetch a CIS-BP PWM for a TF symbol
#'
#' Returns one CIS-BP motif per TF as a log-odds matrix, in the same shape as
#' \code{fetch_jaspar_pwm()}. When CIS-BP lists several motifs for a TF, the
#' motif is chosen by source tier (with \code{prefer = "in_vitro"}):
#' high-throughput in vitro motifs (HT-SELEX, GHT-SELEX, PBM, SMiLE-seq), then
#' other experimentally determined motifs, then motifs imported from JASPAR or
#' HOCOMOCO; within the tier, the motif most similar on average to the TF's
#' other motifs (\code{pick = "medoid"}) or the one with the highest
#' information content (\code{pick = "ic"}). Probabilities are floored at
#' 1e-6 and converted to log2-odds against a uniform background, exactly as
#' for HOCOMOCO.
#'
#' @param tf_name TF gene symbol (matched case-insensitively).
#' @param species CIS-BP species name (e.g. \code{"Homo_sapiens"}) or
#'   Ensembl-style name (\code{"homo_sapiens"}).
#' @param status \code{"direct"} (experimental motifs only, default) or
#'   \code{"inferred"} (fall back to motifs inferred from DNA-binding-domain
#'   similarity when no direct motif exists; useful for species with few
#'   measured motifs).
#' @param exclude_sources CIS-BP \code{MSource_Identifier} values to ignore,
#'   e.g. \code{c("JASPAR", "HocoMoco")}.
#' @param pick \code{"medoid"} (default) or \code{"ic"}.
#' @param prefer \code{"in_vitro"} (default, tiered choice) or \code{"none"}.
#' @return \code{list(id, name, pwm, len, source, status, tier, n_candidates)}
#'   or \code{NULL} if the TF has no usable motif.
#' @examples
#' \dontrun{
#' fetch_cisbp_pwm("MAX")
#' }
#' @export
fetch_cisbp_pwm <- function(tf_name, species = "Homo_sapiens",
                            status = c("direct", "inferred"),
                            exclude_sources = NULL,
                            pick = c("medoid", "ic"),
                            prefer = c("in_vitro", "none")) {
  status <- match.arg(status); pick <- match.arg(pick); prefer <- match.arg(prefer)
  if (species == tolower(species)) species <- .cisbp_species(species)
  sel <- .select_cisbp_motif(tf_name, species, status, exclude_sources, pick, prefer)
  if (is.null(sel)) return(NULL)
  m   <- sweep(sel$ppm, 2, colSums(sel$ppm), "/")
  pwm <- log2(pmax(m, 1e-6) / 0.25)
  list(id = sel$motif_id, name = tf_name, pwm = pwm, len = ncol(pwm),
       source = sel$source, status = sel$status, tier = sel$tier,
       n_candidates = sel$n_candidates)
}

#' List TFs with a CIS-BP motif
#'
#' @inheritParams fetch_cisbp_pwm
#' @return Sorted character vector of TF symbols (upper case).
#' @export
list_cisbp_tfs <- function(species = "Homo_sapiens",
                           status = c("direct", "inferred")) {
  status <- match.arg(status)
  if (species == tolower(species)) species <- .cisbp_species(species)
  tab <- .load_cisbp_table(species)
  keep <- if (status == "direct") tab$TF_Status == "D" else rep(TRUE, nrow(tab))
  sort(unique(tab$TF_Name[keep]))
}

#' Clear the in-memory CIS-BP cache (does not delete the on-disk archive).
#' @return \code{TRUE} (invisibly).
#' @export
clear_cisbp_memory_cache <- function() {
  rm(list = ls(.cisbp_env), envir = .cisbp_env)
  invisible(TRUE)
}

#' Human-readable label of a CIS-BP configuration.
#' @noRd
.cisbp_label <- function(species, status, prefer) {
  sprintf("CIS-BP 3.10 (%s, %s motifs%s)", species, status,
          if (identical(prefer, "in_vitro")) ", in vitro first" else "")
}

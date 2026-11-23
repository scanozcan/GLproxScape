# =============================================================================
# expression.R -- cell-line expression annotation and cell-type ChIP support
# =============================================================================
# Optional layer, run by run_glproxscape() when `cell_type` is given:
#   * human runs: TFs are annotated with RNA expression of the cell line from
#     the Human Protein Atlas (HPA) cell-line data (nTPM). TFs with motif-
#     anchored events are grouped by DNA-binding-domain family (Lambert et al.
#     2018); within a family, TFs whose events overlap (Jaccard index) share
#     sites, and the most highly expressed of these is flagged. Nothing is
#     filtered or re-weighted.
#   * any species with ChIP-Atlas: support of each binding event by ChIP-Atlas
#     experiments from the same cell type vs other cell types.
# Own expression tables (e.g. mouse) can be used through annotate_tf_expression().
# =============================================================================

.hpa_celline_url <- "https://www.proteinatlas.org/download/tsv/rna_celline.tsv.zip"

.expr_cache_dir <- function() {
  d <- file.path(tools::R_user_dir("caspex", which = "cache"), "expression")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

.download_hpa_celline <- function(force = FALSE) {
  dest <- file.path(.expr_cache_dir(), "rna_celline.tsv.zip")
  if (file.exists(dest) && file.size(dest) > 1e8 && !force) return(dest)
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(3600, old_timeout))
  message("Downloading Human Protein Atlas cell-line RNA (~200 MB, one-time)...")
  tmp <- paste0(dest, ".part")
  utils::download.file(.hpa_celline_url, tmp, mode = "wb", quiet = FALSE)
  if (!file.exists(tmp) || file.size(tmp) < 1e8) {
    unlink(tmp)
    stop("Human Protein Atlas download incomplete. Download ", .hpa_celline_url,
         " manually and copy it to:\n  ", dest, call. = FALSE)
  }
  file.rename(tmp, dest)
  dest
}

#' Cell-line RNA expression (Human Protein Atlas or own table)
#'
#' Returns gene-level RNA expression of one cell line. By default the values
#' come from the Human Protein Atlas cell-line RNA data (nTPM, about 1,200
#' human cell lines); the ~200 MB archive is downloaded once and the values
#' of each requested cell line are cached as a small file. Alternatively a
#' gene-level table can be supplied with \code{file} (e.g. for mouse or for
#' cell lines that are not in the Human Protein Atlas).
#'
#' @param cell_type Cell-line name as used by the Human Protein Atlas
#'   (case-insensitive, e.g. \code{"HEK293"}).
#' @param file Optional own table (tab- or comma-separated) with a gene
#'   column and an expression column.
#' @param gene_col,value_col Column names in \code{file}; guessed if
#'   \code{NULL} (gene: gene / symbol / name; value: ntpm / tpm / fpkm /
#'   expression).
#' @param force Re-extract the cell line from the HPA archive.
#' @return \code{data.frame(gene, expr)} with attributes \code{cell_type} and
#'   \code{source}.
#' @examples
#' \dontrun{
#' expr <- load_cell_line_expression("HEK293")
#' }
#' @export
load_cell_line_expression <- function(cell_type = NULL, file = NULL,
                                      gene_col = NULL, value_col = NULL,
                                      force = FALSE) {
  if (!is.null(file)) {
    df <- if (grepl("\\.csv$", file, ignore.case = TRUE))
      utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
    else utils::read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)
    pick <- function(pattern) {
      hit <- grep(pattern, names(df), ignore.case = TRUE, value = TRUE)
      if (!length(hit)) stop("Could not guess the column matching ", pattern,
                             "; pass gene_col / value_col.", call. = FALSE)
      hit[1]
    }
    gene_col  <- gene_col  %||% pick("^(gene.?name|gene|symbol|name)$")
    value_col <- value_col %||% pick("^(ntpm|tpm|fpkm|rpkm|expression|expr)$")
    out <- data.frame(gene = toupper(trimws(df[[gene_col]])),
                      expr = suppressWarnings(as.numeric(df[[value_col]])),
                      stringsAsFactors = FALSE)
    out <- stats::aggregate(expr ~ gene, out, max)
    attr(out, "cell_type") <- cell_type %||% basename(file)
    attr(out, "source")    <- paste0("user table (", value_col, ")")
    return(out)
  }
  if (is.null(cell_type)) stop("Give cell_type or file.", call. = FALSE)
  rds <- file.path(.expr_cache_dir(),
                   paste0("hpa_", gsub("[^A-Za-z0-9]+", "_", cell_type), ".rds"))
  if (file.exists(rds) && !force) {
    out <- readRDS(rds)
    # Caches written by earlier scripts stored the name as "cell_line".
    if (is.null(attr(out, "cell_type")))
      attr(out, "cell_type") <- attr(out, "cell_line") %||% cell_type
    return(out)
  }

  zip <- .download_hpa_celline()
  con <- unz(zip, "rna_celline.tsv")
  open(con, "r")
  on.exit(close(con), add = TRUE)
  hdr <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]
  i_gene <- match("Gene name", hdr); i_cl <- match("Cell line", hdr)
  i_val  <- match("nTPM", hdr)
  if (anyNA(c(i_gene, i_cl, i_val)))
    stop("Unexpected Human Protein Atlas header: ", paste(hdr, collapse = ", "),
         call. = FALSE)
  # The file is sorted by gene, so the first gene's rows list every cell line.
  first <- readLines(con, n = 5000)
  lines_cl <- unique(vapply(strsplit(first, "\t", fixed = TRUE), `[`, "", i_cl))
  norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  name <- lines_cl[norm(lines_cl) == norm(cell_type)]
  if (!length(name)) {
    near <- lines_cl[grepl(substr(norm(cell_type), 1, 3), norm(lines_cl), fixed = TRUE)]
    stop("Cell line '", cell_type, "' is not in the Human Protein Atlas. ",
         "Similar names: ", paste(utils::head(near, 15), collapse = ", "),
         call. = FALSE)
  }
  message("Extracting ", name[1], " from the Human Protein Atlas table (one-time)...")
  key  <- paste0("\t", name[1], "\t")
  keep <- first[grepl(key, first, fixed = TRUE)]
  repeat {
    ln <- readLines(con, n = 1e6)
    if (!length(ln)) break
    keep <- c(keep, ln[grepl(key, ln, fixed = TRUE)])
  }
  f <- strsplit(keep, "\t", fixed = TRUE)
  f <- f[vapply(f, `[`, "", i_cl) == name[1]]
  out <- data.frame(gene = toupper(vapply(f, `[`, "", i_gene)),
                    expr = as.numeric(vapply(f, `[`, "", i_val)),
                    stringsAsFactors = FALSE)
  out <- stats::aggregate(expr ~ gene, out, max)
  attr(out, "cell_type") <- name[1]
  attr(out, "source")    <- "Human Protein Atlas cell-line RNA (nTPM)"
  saveRDS(out, rds)
  out
}

#' DNA-binding-domain family of each human TF (Lambert et al. 2018)
#'
#' @return \code{data.frame(tf, family)}; shipped with the package.
#' @export
load_tf_families <- function() {
  f <- system.file("extdata", "tf_families_lambert2018.tsv.gz",
                   package = "GLproxScape")
  d <- utils::read.delim(f, stringsAsFactors = FALSE, quote = "")
  data.frame(tf = toupper(d$tf), family = d$family, stringsAsFactors = FALSE)
}

#' Jaccard index of two sets of event positions; positions within tol_bp are
#' matched one-to-one (greedy, nearest first).
#' @noRd
.jaccard_sites <- function(a, b, tol_bp = 25) {
  if (!length(a) || !length(b)) return(0)
  used <- logical(length(b)); m <- 0L
  for (x in a) {
    d <- abs(b - x); d[used] <- Inf
    k <- which.min(d)
    if (length(k) && d[k] <= tol_bp) { used[k] <- TRUE; m <- m + 1L }
  }
  m / (length(a) + length(b) - m)
}

#' ChIP-Atlas support of a run's events, unless the TF has no ChIP data.
#' @noRd
.event_chip_support <- function(ev, peaks, tol_bp = 25) {
  has_chip <- function(tf) {
    pk <- if (!is.null(peaks)) peaks[[tf]] else NULL
    !is.null(pk) && isTRUE((attr(pk, "n_srx_scanned") %||% nrow(pk)) > 0)
  }
  hit <- function(pk, pos) any(pk$start_rel <= pos + tol_bp & pk$end_rel >= pos - tol_bp)
  out <- data.frame(tf = as.character(ev$tf), position = ev$position,
                    motif_based = ev$motif_based %in% TRUE,
                    chip_data = vapply(as.character(ev$tf), has_chip, logical(1)),
                    stringsAsFactors = FALSE)
  out$supported <- out$supported_matched <- out$supported_other <- NA
  for (i in which(out$chip_data)) {
    pk <- peaks[[out$tf[i]]]
    m  <- if ("cell_type_matched" %in% names(pk)) pk$cell_type_matched %in% TRUE
          else rep(FALSE, nrow(pk))
    out$supported[i]         <- hit(pk, out$position[i])
    out$supported_matched[i] <- hit(pk[m, , drop = FALSE], out$position[i])
    out$supported_other[i]   <- hit(pk[!m, , drop = FALSE], out$position[i])
  }
  rownames(out) <- NULL
  out
}

#' ChIP-Atlas support of binding events by cell type
#'
#' For every binding event of a \code{\link{run_glproxscape}} result, reports
#' whether a ChIP-Atlas peak of the same TF lies within \code{tol_bp}, from
#' experiments of the experimental cell type (\code{cell_type} given to
#' \code{run_glproxscape}) and from other cell types separately.
#'
#' @param result A \code{\link{run_glproxscape}} result with ChIP-Atlas peaks.
#' @param tol_bp Window around each event (default 25 bp).
#' @return A list: \code{events} (one row per event) and \code{by_tf}
#'   (per TF: events, events supported by any / matched / other-cell-type
#'   experiments, and the number of matched and all experiments scanned).
#' @export
chipatlas_support_by_cell_type <- function(result, tol_bp = 25) {
  peaks <- result$chipatlas_peaks
  if (is.null(peaks)) stop("result has no ChIP-Atlas peaks (run with chipatlas = TRUE).")
  ev <- .event_chip_support(result$binding_events, peaks, tol_bp)
  tfs <- unique(ev$tf)
  by_tf <- do.call(rbind, lapply(tfs, function(tf) {
    e <- ev[ev$tf == tf, ]; pk <- peaks[[tf]]
    data.frame(tf = tf, n_events = nrow(e),
               n_supported = sum(e$supported, na.rm = TRUE),
               n_supported_matched = sum(e$supported_matched, na.rm = TRUE),
               n_supported_other = sum(e$supported_other, na.rm = TRUE),
               n_experiments_scanned = as.integer(attr(pk, "n_srx_scanned") %||% NA),
               n_experiments_matched = as.integer(attr(pk, "n_srx_matched") %||% NA),
               chip_data = e$chip_data[1], stringsAsFactors = FALSE)
  }))
  attr(ev, "cell_type") <- attr(by_tf, "cell_type") <- result$cell_type
  list(events = ev, by_tf = by_tf)
}

#' Annotate a run with cell-line expression of the detected TFs
#'
#' Adds RNA expression of the experimental cell line to the TFs of a
#' \code{\link{run_glproxscape}} result and helps to choose among TFs that
#' recognize similar motifs. TFs with motif-anchored events are grouped by
#' DNA-binding-domain family (Lambert et al. 2018; C2H2 zinc fingers and TFs
#' without an assigned family are not grouped, because their members do not
#' share motifs). Within a family, the overlap of two TFs' events is the
#' Jaccard index (events within \code{tol_bp} count as the same site); TFs
#' with a Jaccard index of at least \code{jaccard_min} share sites, and the
#' most highly expressed of these is flagged as \code{preferred}. Predictions
#' are not filtered or re-weighted.
#'
#' @param result A \code{\link{run_glproxscape}} result.
#' @param expr Output of \code{\link{load_cell_line_expression}}.
#' @param families Output of \code{\link{load_tf_families}} (default).
#' @param min_expr Expression cut-off for \code{expressed} (default 1).
#' @param tol_bp Same-site window for the Jaccard index and the ChIP-Atlas
#'   match (default 25 bp).
#' @param jaccard_min Jaccard index at which family members share sites
#'   (default 0.3).
#' @param skip_families DBD families that are not grouped.
#' @return A list: \code{tf} (one row per TF with events), \code{events},
#'   \code{pairs} (Jaccard index of family member pairs) and \code{meta}.
#' @examples
#' \dontrun{
#' ann <- annotate_tf_expression(res, load_cell_line_expression("HEK293"))
#' plot_expression_families(ann)
#' }
#' @export
annotate_tf_expression <- function(result, expr, families = load_tf_families(),
                                   min_expr = 1, tol_bp = 25, jaccard_min = 0.3,
                                   skip_families = c("C2H2 ZF", "Unknown")) {
  ev <- result$binding_events
  if (is.null(ev) || !nrow(ev)) stop("result has no binding_events.", call. = FALSE)
  ev <- data.frame(tf = as.character(ev$tf), position = ev$position,
                   motif_based = ev$motif_based %in% TRUE, stringsAsFactors = FALSE)
  tfs <- unique(ev$tf)
  tab <- data.frame(tf = tfs, stringsAsFactors = FALSE)
  e <- expr$expr[match(toupper(tfs), expr$gene)]
  tab$expr <- e
  tab$in_expression_table <- !is.na(e)
  tab$expressed <- !is.na(e) & e >= min_expr
  tab$family <- families$family[match(toupper(tfs), families$tf)]
  tab$family[is.na(tab$family)] <- "not in TF list"
  tab$motif_anchored <- tab$tf %in% ev$tf[ev$motif_based]
  sd <- result$spatial_df
  tab$composite <- if (!is.null(sd)) sd$composite[match(tab$tf, sd$protein)] else NA_real_

  ok <- tab$motif_anchored & !(tab$family %in% c(skip_families, "not in TF list")) &
        !grepl(";", tab$family)
  fam_n <- table(tab$family[ok])
  tab$family_group <- ifelse(ok & tab$family %in% names(fam_n)[fam_n >= 2],
                             tab$family, NA_character_)

  pos <- split(ev$position[ev$motif_based], ev$tf[ev$motif_based])
  pairs <- do.call(rbind, lapply(unique(stats::na.omit(tab$family_group)), function(f) {
    m <- tab$tf[tab$family_group %in% f]
    cmb <- utils::combn(m, 2)
    data.frame(family = f, tf1 = cmb[1, ], tf2 = cmb[2, ],
               jaccard = apply(cmb, 2, function(x)
                 .jaccard_sites(pos[[x[1]]], pos[[x[2]]], tol_bp)),
               stringsAsFactors = FALSE)
  }))
  if (is.null(pairs)) pairs <- data.frame(family = character(0), tf1 = character(0),
                                          tf2 = character(0), jaccard = numeric(0))
  pairs$linked <- pairs$jaccard >= jaccard_min

  # Sets of TFs sharing sites (connected by linked pairs); preferred = most expressed
  tab$sharing_set <- NA_integer_; tab$preferred <- FALSE
  lk <- pairs[pairs$linked, , drop = FALSE]; k <- 0L
  for (f in unique(lk$family)) {
    l <- lk[lk$family == f, ]
    memb <- unique(c(l$tf1, l$tf2)); comp <- stats::setNames(seq_along(memb), memb)
    repeat {
      changed <- FALSE
      for (i in seq_len(nrow(l))) {
        a <- comp[l$tf1[i]]; b <- comp[l$tf2[i]]
        if (a != b) { comp[comp == max(a, b)] <- min(a, b); changed <- TRUE }
      }
      if (!changed) break
    }
    for (cc in unique(comp)) {
      k <- k + 1L
      i <- match(names(comp)[comp == cc], tab$tf)
      tab$sharing_set[i] <- k
      x <- ifelse(is.na(tab$expr[i]), 0, tab$expr[i])
      if (max(x) > 0) tab$preferred[i[which.max(x)]] <- TRUE
    }
  }

  if (!is.null(result$chipatlas_peaks)) {
    cs <- .event_chip_support(result$binding_events, result$chipatlas_peaks, tol_bp)
    ev$supported <- cs$supported
    tab$chip_data <- vapply(tab$tf, function(tf) any(cs$chip_data[cs$tf == tf]), logical(1))
    tab$n_events  <- as.integer(table(factor(ev$tf, levels = tab$tf)))
    tab$frac_supported <- ifelse(tab$chip_data,
      vapply(tab$tf, function(tf) mean(ev$supported[ev$tf == tf] %in% TRUE), numeric(1)),
      NA_real_)
  } else {
    tab$chip_data <- FALSE
    tab$n_events  <- as.integer(table(factor(ev$tf, levels = tab$tf)))
    tab$frac_supported <- NA_real_
  }
  ev <- cbind(ev, tab[match(ev$tf, tab$tf), c("expr", "expressed", "family",
                                              "family_group", "sharing_set", "preferred")])
  tab <- tab[order(is.na(tab$family_group), tab$family_group,
                   -ifelse(is.na(tab$expr), -1, tab$expr)), ]
  rownames(tab) <- NULL
  list(tf = tab, events = ev, pairs = pairs,
       meta = list(cell_type = attr(expr, "cell_type") %||% attr(expr, "cell_line") %||% "cell line", source = attr(expr, "source"),
                   min_expr = min_expr, tol_bp = tol_bp, jaccard_min = jaccard_min,
                   skip_families = skip_families))
}

#' Plot detected TFs of the same DNA-binding-domain family by expression
#'
#' One row per DNA-binding-domain family with at least two motif-anchored
#' TFs. TFs are placed by their expression in the cell line; point size is
#' spatial enrichment (composite), fill the fraction of the TF's events with
#' a ChIP-Atlas peak (gray: no ChIP-Atlas data), and a thick outline marks
#' the most highly expressed TF among family members sharing sites.
#'
#' @param ann Output of \code{\link{annotate_tf_expression}}.
#' @return A ggplot, or \code{NULL} when no family has two motif-anchored TFs.
#' @export
plot_expression_families <- function(ann) {
  d <- ann$tf[!is.na(ann$tf$family_group), , drop = FALSE]
  if (!nrow(d)) { message("No DBD families with >= 2 motif-anchored TFs."); return(NULL) }
  d$x <- log2(ifelse(is.na(d$expr), 0, d$expr) + 1)
  ord <- names(sort(tapply(d$x, d$family_group, max)))
  d$y <- match(d$family_group, ord)
  rng <- do.call(rbind, lapply(split(d, d$family_group), function(z)
    data.frame(y = z$y[1], xmin = min(z$x), xmax = max(z$x))))
  d$stroke <- ifelse(d$preferred, 1.3, 0.3)
  m <- ann$meta
  p <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = log2(m$min_expr + 1), linetype = "dashed",
                        colour = "grey60") +
    ggplot2::geom_segment(data = rng, ggplot2::aes(x = .data$xmin, xend = .data$xmax,
                                                   y = .data$y, yend = .data$y),
                          colour = "grey75", linewidth = 0.5) +
    ggplot2::geom_point(data = d, ggplot2::aes(.data$x, .data$y, size = .data$composite,
                                               fill = .data$frac_supported,
                                               stroke = .data$stroke),
                        shape = 21, colour = "black") +
    ggplot2::scale_fill_gradient(low = "#A8DADC", high = "#1D3557", limits = c(0, 1),
                                 na.value = "grey85", name = "Events with\nChIP-Atlas peak") +
    ggplot2::scale_size_continuous(range = c(1.5, 6), name = "Spatial\nenrichment") +
    ggplot2::scale_y_continuous(breaks = seq_along(ord), labels = ord,
                                expand = ggplot2::expansion(add = c(0.6, 0.4))) +
    ggplot2::labs(x = sprintf("Expression in %s, log2(nTPM + 1)", m$cell_type), y = NULL,
                  title = "Detected TFs of the same DNA-binding-domain family",
                  subtitle = sprintf(paste0("families: Lambert et al. 2018 (C2H2 ZF and Unknown not grouped)\n",
                                            "thick outline = most-expressed TF among members sharing sites ",
                                            "(Jaccard >= %.1f, events within %d bp)\n",
                                            "dashed line = %g nTPM"),
                                     m$jaccard_min, m$tol_bp, m$min_expr)) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  if (requireNamespace("ggrepel", quietly = TRUE))
    p + ggrepel::geom_text_repel(data = d, ggplot2::aes(.data$x, .data$y, label = .data$tf),
                                 size = 2.6, nudge_y = -0.2, direction = "x",
                                 min.segment.length = grid::unit(1.5, "lines"),
                                 segment.colour = "grey50", seed = 1)
  else p + ggplot2::geom_text(data = d, ggplot2::aes(.data$x, .data$y, label = .data$tf),
                              size = 2.6, vjust = 2)
}

#' Run the cell-type steps at the end of run_glproxscape().
#' @noRd
.run_cell_type_layer <- function(result, out_dir, save_plots = TRUE) {
  ct <- result$cell_type
  if (is.null(ct) || !nzchar(ct)) return(result)
  if (!is.null(result$chipatlas_peaks)) {
    cs <- tryCatch(chipatlas_support_by_cell_type(result), error = function(e) {
      message("  Cell-type ChIP-Atlas support skipped: ", conditionMessage(e)); NULL })
    if (!is.null(cs)) {
      result$chipatlas_support <- cs
      utils::write.csv(cs$by_tf, file.path(out_dir, "chipatlas_support_by_cell_type.csv"),
                       row.names = FALSE)
      message(sprintf("  ChIP-Atlas support (%s): %d events supported by %s experiments, %d by other cell types",
                      ct, sum(cs$events$supported_matched, na.rm = TRUE), ct,
                      sum(cs$events$supported_other, na.rm = TRUE)))
    }
  }
  if (!identical(tolower(result$species %||% "homo_sapiens"), "homo_sapiens")) {
    message("  Expression annotation skipped: Human Protein Atlas data are human only ",
            "(use annotate_tf_expression() with load_cell_line_expression(file = ...)).")
    return(result)
  }
  ann <- tryCatch({
    expr <- load_cell_line_expression(ct)
    annotate_tf_expression(result, expr)
  }, error = function(e) {
    message("  Expression annotation skipped: ", conditionMessage(e)); NULL })
  if (is.null(ann)) return(result)
  result$expression <- ann
  utils::write.csv(ann$tf, file.path(out_dir, "tf_expression_annotation.csv"), row.names = FALSE)
  utils::write.csv(ann$pairs, file.path(out_dir, "tf_family_shared_sites.csv"), row.names = FALSE)
  if (isTRUE(save_plots)) {
    p <- plot_expression_families(ann)
    if (!is.null(p)) {
      n <- length(unique(stats::na.omit(ann$tf$family_group)))
      grDevices::pdf(file.path(out_dir, "tf_expression_families.pdf"), width = 9,
                     height = max(4, 1.1 * n + 2))
      print(p)
      grDevices::dev.off()
    }
  }
  message(sprintf("  Expression (%s, Human Protein Atlas): %d of %d TFs with events expressed",
                  ann$meta$cell_type, sum(ann$tf$expressed), nrow(ann$tf)))
  result
}

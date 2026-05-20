#' GLproxScape: Spatial deconvolution of dCas9-APEX2 proximity proteomics
#'
#' GLproxScape recovers transcription-factor and chromatin-regulator binding
#' predictions from dCas9-APEX2 promoter-tiling proximity-proteomics data.
#' Each guide is modelled as a Gaussian biotinylation cone (default
#' \eqn{\sigma = 300} bp); per-region enrichment values are forward-smeared into
#' a continuous spatial track \eqn{s(x)}; coverage normalisation
#' \eqn{\beta(x) = s(x) / C(x)} removes the bias toward dense-tiling
#' interiors; and motif-anchored binding events are recovered via
#' non-negative least squares against JASPAR position-weight matrices, with a
#' separate zone-based path for sequence-specific-motif-less chromatin
#' factors.
#'
#' @section Main entry points:
#' \describe{
#'   \item{\code{\link{run_caspex}}}{Top-level pipeline: gene lookup, gRNA
#'     matching, signal modelling, motif scan, deconvolution, plotting, and
#'     ChIP-Atlas overlay.}
#'   \item{\code{\link{run_caspex_extras}}}{Diagnostic plot pack
#'     (A1-A6, B1-B3, C1-C3, D1-D3) for a completed run.}
#'   \item{\code{\link{run_caspex_epigenetic}}}{Zone-based deconvolution
#'     for chromatin readers/writers/erasers/remodellers without a
#'     sequence-specific motif.}
#'   \item{\code{\link{load_caspex_inputs}}}{Read the
#'     \code{grnas.tsv} + per-region \code{Region*.txt} input layout.}
#' }
#'
#' @section Bundled example data:
#' A small example sgRNA manifest + per-region differential file is
#' shipped under \code{inst/extdata/examples/} for end-to-end smoke
#' testing and demonstration of the input file layout. Resolve a path
#' with \code{system.file("extdata/examples", package = "GLproxScape")}.
#'
#' @keywords internal
#' @import ggplot2
#' @import patchwork
#' @importFrom httr GET modify_url add_headers content status_code timeout
#' @importFrom nnls nnls
#' @importFrom scales comma percent
#' @importFrom grDevices dev.list dev.off pdf
#' @importFrom stats aggregate approx complete.cases cor embed median p.adjust qnorm quantile reshape sd setNames var weighted.mean
#' @importFrom tools R_user_dir
#' @importFrom utils head read.csv read.delim tail write.csv
"_PACKAGE"

# \u2500\u2500 R CMD check NSE-binding declarations \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
# ggplot2 evaluates aes() arguments as data-frame column names; R CMD check
# can't tell that and reports them as "no visible binding for global variable".
# Declaring them here is the standard fix; same for column names produced by
# our reshape() / aggregate() / merge() steps. Anything that's not actually
# a column name (e.g. mistyped helper) belongs in code, not on this list \u2014
# keep this list trimmed to the real NSE bindings.
utils::globalVariables(c(
  # primary event / signal columns
  "tf", "position", "weight", "motif_based", "n_motifs_merged",
  "lfc", "pval", "protein", "region", "pos", "centroid", "composite",
  "specificity", "sig", "stars", "frac", "label", "cell_id", "panel",
  # spatial-track / motif-tick geometry
  "x", "y", "y_lo", "y_hi", "y_mid", "y_top", "y_bottom", "y_base",
  "y_signal", "y0", "y1", "yc", "xs", "xe", "xc", "ymin", "ymax",
  "xmin", "xmax", "angle", "hjust",
  # epigenetic deck columns
  "zone_start", "zone_end", "core_start", "core_end", "centroid_pos",
  "peak_beta", "is_det",
  # ChIP-Atlas / histone columns
  "start_rel", "end_rel", "n_srx", "fill",
  # diagnostic-plot columns
  "n_surv", "n_drops", "logp", "kind", "lab",
  "null_mean", "null_sd", "perm_p", "perm_fdr",
  "rescued", "local_coverage", "distance_to_nearest_grna",
  "survived", "beta_norm", "surv", "confidence",
  "cov_floor", "sigma", "family", "n", "status", "w", "r", "R1", "R2",
  # internal field names referenced as column accessors after reshape()
  ".vj", ".hj"
))

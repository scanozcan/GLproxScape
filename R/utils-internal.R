# =============================================================================
# Package-internal utilities used by multiple R/ files.
# Nothing in here is exported; everything is documented @noRd so roxygen2
# leaves the user-facing manual clean.
# =============================================================================

#' Null-coalescing operator
#'
#' Returns `b` if `a` is NULL; otherwise returns `a`. Used pervasively
#' for default-fallback patterns such as `result$weight_mode %||% "z"`
#' so callers can pass a sparse result list and individual fields
#' default cleanly.
#'
#' IMPORTANT: this used to also fall back when `a` was length-0 or
#' NA in its first element, but `is.na(a[[1]])` errors in R 4.3+ when
#' `a` is a list whose first element is itself a multi-field list (e.g.
#' `result$motif_results_extra`), because `is.na()` then returns a
#' multi-element logical vector and R 4.3+ disallows non-scalars on
#' either side of `||`. That bug silently knocked out the B2 / D2
#' sensitivity sweeps via the `motif_results_extra %||% list()` call in
#' `run_sigma_sensitivity` / `run_covfloor_sensitivity`. Keeping the
#' operator scalar-safe by checking NULL only.
#'
#' Defined module-level (not inside a function) so every R/ file in the
#' package has it visible at parse time. Replaces the three duplicate
#' definitions previously in caspex_analysis.R, caspex_chipatlas.R, and
#' caspex_epigenetic.R from when those files needed to be sourceable in
#' isolation.
#'
#' @param a left-hand operand.
#' @param b fallback value.
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' GLproxScape colour palette
#'
#' Internal palette used by every plotting helper so figures stay
#' visually coherent across decks. Slots:
#' \describe{
#'   \item{tss / high}{Punchy red `#E63946` — TSS line and motif-anchored
#'     event bubbles.}
#'   \item{guide}{Mid blue `#457B9D` — gRNA cut sites and no-motif events.}
#'   \item{guide_ns}{Pale teal `#A8DADC` — gRNAs without a matched sequence.}
#'   \item{mid}{Warm orange `#F4A261` — JASPAR motif tick marks.}
#'   \item{low}{Forest teal `#2A9D8F` — coverage curves, repressive marks.}
#'   \item{neutral}{Slate `#264653` — text annotations, baseline strokes.}
#' }
#' @noRd
COLS <- list(
  tss      = "#E63946",
  guide    = "#457B9D",
  guide_ns = "#A8DADC",
  mid      = "#F4A261",
  high     = "#E63946",
  low      = "#2A9D8F",
  neutral  = "#264653"
)

#' Package-internal JASPAR PWM cache
#'
#' Keyed by `toupper(tf_name)` + host. Lives in the package namespace
#' (parent = emptyenv() so cached PWM lists don't accidentally inherit
#' bindings from the package or global env). One cache per loaded
#' package session. Replaces the historical \code{globalenv()}-anchored
#' cache, which only existed to survive re-sourcing of the script-style
#' codebase — re-sourcing isn't a thing once we're a package.
#'
#' Clear interactively with
#' \code{rm(list = ls(GLproxScape:::.caspex_pwm_cache),
#'        envir = GLproxScape:::.caspex_pwm_cache)}.
#' @noRd
.caspex_pwm_cache <- new.env(parent = emptyenv())

#' Safe wrapper around the `pdf() / for-print / dev.off()` idiom for
#' multi-page PDFs.
#'
#' Common failure mode is `Error in dev.off() : write failed` which fires
#' when the target PDF is held open by another process (Preview / Acrobat
#' on macOS still showing the previous version of the file). Without this
#' wrapper, a single such failure aborts the entire run partway through
#' plot generation; with it, the offending PDF is skipped (with a clear
#' warning telling the user to close any open viewer) and the rest of the
#' deck still renders. `on.exit` ensures the device is closed even when
#' `print()` throws, so subsequent `pdf()` calls don't end up writing to
#' a stale device.
#'
#' Used by run_caspex (analysis-side) and run_caspex_epigenetic
#' (epigenetic-side); kept here so both paths share the same retry policy.
#'
#' @param path output PDF path.
#' @param width,height device dimensions in inches.
#' @param plots list of grobs / ggplot objects to print, one per page.
#' @param label short label for the warning message; defaults to the
#'   PDF basename.
#' @return TRUE on success, FALSE if the device couldn't be opened or any
#'   page failed to print. Always returns invisibly.
#' @noRd
.safe_pdf <- function(path, width, height, plots,
                      label = basename(path)) {
  ok <- tryCatch({
    grDevices::pdf(path, width = width, height = height)
    on.exit({
      while (length(grDevices::dev.list()) > 0)
        try(invisible(grDevices::dev.off()), silent = TRUE)
    }, add = TRUE)
    for (pl in plots) print(pl)
    TRUE
  }, error = function(e) {
    message("  WARNING: ", label, " failed to write: ",
            conditionMessage(e))
    message("    (typical cause: PDF is already open in another ",
            "app \u2014 close it and re-run, or delete the stale file)")
    while (length(grDevices::dev.list()) > 0)
      try(invisible(grDevices::dev.off()), silent = TRUE)
    FALSE
  })
  invisible(ok)
}

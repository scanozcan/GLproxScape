# =============================================================================
# GLproxScape-deprecated.R -- former function names (version < 0.2.0)
# =============================================================================
# The package API was renamed from the CasPEX-era names to GLproxScape names
# in version 0.2.0. The old names keep working and point to the new ones.

#' Deprecated function names
#'
#' These functions were renamed in GLproxScape 0.2.0. The old names still
#' work and call the new functions, with a deprecation message.
#'
#' \tabular{ll}{
#'   Old name \tab New name \cr
#'   \code{run_caspex()} \tab \code{\link{run_glproxscape}()} \cr
#'   \code{run_caspex_epigenetic()} \tab \code{\link{run_glproxscape_epigenetic}()} \cr
#'   \code{run_caspex_extras()} \tab \code{\link{run_glproxscape_extras}()} \cr
#'   \code{load_caspex_inputs()} \tab \code{\link{load_glproxscape_inputs}()} \cr
#'   \code{plot_binding_deconvolution()} \tab \code{\link{plot_binding_events}()} \cr
#' }
#'
#' @param ... Arguments passed to the new function.
#' @return The value of the new function.
#' @name GLproxScape-deprecated
#' @keywords internal
NULL

#' @rdname GLproxScape-deprecated
#' @export
run_caspex <- function(...) {
  .Deprecated("run_glproxscape", package = "GLproxScape")
  run_glproxscape(...)
}

#' @rdname GLproxScape-deprecated
#' @export
run_caspex_epigenetic <- function(...) {
  .Deprecated("run_glproxscape_epigenetic", package = "GLproxScape")
  run_glproxscape_epigenetic(...)
}

#' @rdname GLproxScape-deprecated
#' @export
run_caspex_extras <- function(...) {
  .Deprecated("run_glproxscape_extras", package = "GLproxScape")
  run_glproxscape_extras(...)
}

#' @rdname GLproxScape-deprecated
#' @export
load_caspex_inputs <- function(...) {
  .Deprecated("load_glproxscape_inputs", package = "GLproxScape")
  load_glproxscape_inputs(...)
}

#' @rdname GLproxScape-deprecated
#' @export
plot_binding_deconvolution <- function(...) {
  .Deprecated("plot_binding_events", package = "GLproxScape")
  plot_binding_events(...)
}

# Internal theme kept under its old name for any user code that called it
# via GLproxScape:::theme_caspex().
theme_caspex <- function(...) theme_glproxscape(...)

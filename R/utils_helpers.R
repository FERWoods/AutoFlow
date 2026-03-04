# Null coalescing operator
# Returns `a` if it is not NULL, otherwise returns `b`.
# Useful for setting defaults.
#' @noRd
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

# Canonicalise strings for matching
# Lowercases and removes non-alphanumeric characters.
# Helps when matching channel names or markers.
#' @param x Character vector
#' @return Canonicalised character vector
#' @noRdß
canon <- function(x) {
  gsub("[^a-z0-9]+", "", tolower(x))
}

#' Add a density curve or fallback histogram to a plotly object
#'
#' Helper for QC plotting: given a numeric vector, this adds a kernel density
#' estimate as a line trace when there is enough variability; otherwise it
#' falls back to a 1-bin histogram for a single observation, or does nothing
#' for empty/invalid inputs.
#'
#' Non-finite values are dropped after coercing to numeric.
#'
#' @param p A `plotly` object to add a trace to.
#' @param x A numeric vector (or coercible to numeric) of observations.
#' @param nm Character scalar used as the trace name.
#'
#' @return A `plotly` object. If `x` has >=2 finite values and non-zero variance,
#'   a density line trace is added. If `x` has exactly 1 finite value, a 1-bin
#'   histogram trace is added. Otherwise `p` is returned unchanged.
#'
#' @noRd
dens_or_hist <- function(p, x, nm) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) == 0) return(p)

  if (length(x) == 1) {
    return(plotly::add_trace(
      p,
      x = x,
      type = "histogram",
      nbinsx = 1,
      name = paste0(nm, " (1 obs)")
    ))
  }

  v <- suppressWarnings(stats::var(x))
  if (!is.finite(v) || v <= 0) return(p)

  d <- stats::density(x)
  plotly::add_trace(
    p,
    x = d$x,
    y = d$y,
    type = "scatter",
    mode = "lines",
    name = nm
  )
}

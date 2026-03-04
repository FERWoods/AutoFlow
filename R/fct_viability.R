# R/fct_viability.R
# Helpers for extracting marker vectors and estimating viability thresholds

#' Extract a marker vector across a list of flowFrames
#'
#' Pulls values for a single channel (by NAME) from each flowFrame and concatenates
#' them into one numeric vector. Frames that do not contain the channel contribute
#' zero values (are skipped).
#'
#' @param frames List of `flowCore::flowFrame`.
#' @param name_col Character scalar; channel NAME to extract (must match `flowCore::colnames(ff)`).
#'
#' @return Numeric vector containing concatenated marker values across frames.
#'   Returns `numeric(0)` if `frames` is NULL/empty, `name_col` is NULL, or the channel
#'   is not found in any frame.
#' @noRd

get_marker_vector <- function(frames, name_col) {
  if (is.null(frames) || !length(frames) || is.null(name_col)) return(numeric(0))
  unlist(lapply(frames, function(ff) {
    cn <- flowCore::colnames(ff)
    if (!(name_col %in% cn)) return(numeric(0))
    as.numeric(flowCore::exprs(ff)[, name_col])
  }), use.names = FALSE)
}
#' Min-max normalise a numeric vector to [0, 1]
#'
#' Computes `(x - min(x)) / (max(x) - min(x))` using `na.rm = TRUE`.
#' If the range is non-finite or zero (all values equal), returns a vector of zeros.
#'
#' @param x Numeric vector.
#'
#' @return Numeric vector of the same length as `x`, scaled to [0, 1] where possible.
#'   If scaling cannot be performed, returns `rep(0, length(x))`.
#' @noRd
norm_minmax <- function(x) {
  if (!any(is.finite(x))) return(rep(0, length(x)))
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r)) || diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / (r[2] - r[1])
}

#' Estimate a bimodal threshold using a 2-component Gaussian mixture model
#'
#' Fits a 2-component Gaussian mixture (via `mclust::Mclust`) to a numeric vector,
#' and estimates the decision threshold as the intersection point of the two
#' weighted Gaussian component densities.
#'
#' This is used to propose an automatic viability threshold for a marker.
#'
#' Failure modes return `NA_real_`:
#' - fewer than 200 finite values
#' - mixture model fit fails
#' - intersection root finding fails
#'
#' @param x Numeric vector (e.g., marker intensities). Non-finite values are removed.
#'
#' @return Numeric scalar threshold (intersection of weighted component densities),
#'   or `NA_real_` if a robust threshold cannot be found.
#' @noRd
auto_threshold_gmm <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 200) return(NA_real_)

  m <- try(mclust::Mclust(x, G = 2, verbose = FALSE), silent = TRUE)
  if (inherits(m, "try-error")) return(NA_real_)

  # Use the fitted classification to get robust component stats (avoids mclust variance-shape quirks)
  cl <- m$classification
  if (is.null(cl) || length(unique(cl)) < 2) return(NA_real_)

  mu <- as.numeric(tapply(x, cl, mean))
  sd <- as.numeric(tapply(x, cl, stats::sd))
  pr <- as.numeric(table(cl)) / length(cl)

  # order low -> high mean
  o  <- order(mu)
  mu <- mu[o]; sd <- sd[o]; pr <- pr[o]

  if (any(!is.finite(mu)) || any(!is.finite(sd)) || any(sd <= 0) || any(!is.finite(pr)) || any(pr <= 0)) {
    return(NA_real_)
  }

  m1 <- mu[1]; m2 <- mu[2]
  s1 <- sd[1]; s2 <- sd[2]
  w1 <- pr[1]; w2 <- pr[2]

  # Solve w1*N(m1,s1) = w2*N(m2,s2)
  # => quadratic in t
  if (abs(s1 - s2) < 1e-12) {
    # equal-variance case -> linear solution
    s <- s1
    thr <- (m1 + m2) / 2 + (s^2 / (m2 - m1)) * log(w2 / w1)
  } else {
    a <- 1/(2*s2^2) - 1/(2*s1^2)
    b <- m1/(s1^2) - m2/(s2^2)
    c <- (m2^2)/(2*s2^2) - (m1^2)/(2*s1^2) + log((w2*s1)/(w1*s2))

    disc <- b^2 - 4*a*c
    if (!is.finite(disc) || disc < 0) return(NA_real_)

    r1 <- (-b + sqrt(disc)) / (2*a)
    r2 <- (-b - sqrt(disc)) / (2*a)

    # choose the root between the component means if possible
    lo <- min(m1, m2); hi <- max(m1, m2)
    cand <- c(r1, r2)
    cand <- cand[is.finite(cand)]
    between <- cand[cand > lo & cand < hi]

    thr <- if (length(between)) between[1] else {
      # fallback: closest root to midpoint
      mid <- (m1 + m2) / 2
      cand[which.min(abs(cand - mid))]
    }
  }

  if (!is.finite(thr)) NA_real_ else thr
}

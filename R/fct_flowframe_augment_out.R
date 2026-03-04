
#' Appends extra columns to exprs and extends the parameters block while
#' preserving the internal flowCore "parameters" class invariants.
#'
#' @param ff flowCore::flowFrame
#' @param add numeric/integer matrix or data.frame with nrow == nrow(exprs(ff))
#' @param extra_desc optional character vector length ncol(add) for parameter desc
#' @return flowCore::flowFrame with extra columns appended
#' @noRd
ff_add_cols_safe <- function(ff, add, extra_desc = NULL) {
  stopifnot(inherits(ff, "flowFrame"))

  x <- flowCore::exprs(ff)
  if (!is.matrix(x)) x <- as.matrix(x)

  addm <- add
  if (is.data.frame(addm)) addm <- as.matrix(addm)
  if (!is.matrix(addm)) stop("add must be a matrix or data.frame")

  if (nrow(addm) != nrow(x)) {
    stop(sprintf("Row mismatch: exprs has %d rows; add has %d rows", nrow(x), nrow(addm)))
  }
  if (is.null(colnames(addm)) || any(!nzchar(colnames(addm)))) {
    stop("add must have non-empty colnames")
  }

  # prevent name collisions
  new_names <- colnames(addm)
  collide <- intersect(new_names, colnames(x))
  if (length(collide)) {
    new_names <- make.unique(c(colnames(x), new_names))
    new_names <- tail(new_names, ncol(addm))
    colnames(addm) <- new_names
  }

  # append expression matrix
  x2 <- cbind(x, addm)

  # --- parameter stats for new cols ---
  mins <- apply(addm, 2, function(v) suppressWarnings(min(as.numeric(v), na.rm = TRUE)))
  maxs <- apply(addm, 2, function(v) suppressWarnings(max(as.numeric(v), na.rm = TRUE)))
  mins[!is.finite(mins)] <- 0
  maxs[!is.finite(maxs)] <- 1
  span <- pmax(1, ceiling(maxs - mins))

  if (is.null(extra_desc)) extra_desc <- colnames(addm)
  if (length(extra_desc) != ncol(addm)) extra_desc <- rep(extra_desc[1], ncol(addm))

  new_pd <- data.frame(
    name     = colnames(addm),
    desc     = as.character(extra_desc),
    range    = as.numeric(span),
    minRange = as.numeric(mins),
    maxRange = as.numeric(maxs),
    stringsAsFactors = FALSE
  )

  # --- extend parameters slot ---
  params <- flowCore::parameters(ff)  # S4 "parameters"
  pd     <- flowCore::pData(params)

  need_cols <- c("name", "desc", "range", "minRange", "maxRange")
  for (cc in need_cols) if (!cc %in% names(pd)) pd[[cc]] <- NA

  pd2 <- rbind(pd[, need_cols, drop = FALSE], new_pd[, need_cols, drop = FALSE])

  # invariants
  pd2$name <- colnames(x2)
  rownames(pd2) <- paste0("$P", seq_len(nrow(pd2)))

  Biobase::pData(params) <- pd2

  ff2 <- ff
  methods::slot(ff2, "parameters") <- params
  flowCore::exprs(ff2) <- x2

  # --- FIXED: keep/extend $Pn* keywords instead of deleting ---
  kw <- tryCatch(flowCore::keyword(ff2), error = function(e) NULL)
  if (is.null(kw)) kw <- list()

  # Ensure $PAR/$TOT correct
  kw[["$PAR"]] <- as.character(ncol(x2))
  kw[["$TOT"]] <- as.character(nrow(x2))

  # Helper to set per-parameter keywords safely
  set_pn <- function(k, N = NULL, S = NULL, R = NULL, B = NULL, E = NULL) {
    kk <- as.integer(k)
    if (!is.null(N)) kw[[sprintf("$P%dN", kk)]] <- as.character(N)
    if (!is.null(S)) kw[[sprintf("$P%dS", kk)]] <- as.character(S)
    if (!is.null(R)) kw[[sprintf("$P%dR", kk)]] <- as.character(R)
    if (!is.null(B)) kw[[sprintf("$P%dB", kk)]] <- as.character(B)
    if (!is.null(E)) kw[[sprintf("$P%dE", kk)]] <- as.character(E)
  }

  # For ALL params, ensure $PnN at least matches name
  for (k in seq_len(ncol(x2))) {
    set_pn(k, N = colnames(x2)[k])
  }

  # For the NEW appended cols, set sensible defaults
  old_p <- ncol(x)
  for (j in seq_len(ncol(addm))) {
    k <- old_p + j
    set_pn(
      k,
      S = extra_desc[j],
      R = max(1L, as.integer(span[j])),
      B = 32L,      # 32-bit storage is safe default
      E = "0,0"     # linear scale
    )
  }

  # Write back keywords
  flowCore::keyword(ff2) <- kw

  ff2
}

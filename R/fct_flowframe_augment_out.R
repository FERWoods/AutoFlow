
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

  # prevent name collisions (against existing exprs colnames)
  new_names <- colnames(addm)
  collide <- intersect(new_names, colnames(x))
  if (length(collide)) {
    new_names <- make.unique(c(colnames(x), new_names))
    new_names <- tail(new_names, ncol(addm))
    colnames(addm) <- new_names
  }

  # append expression matrix
  x2 <- cbind(x, addm)

  # --- build new parameter rows ---
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

  # --- extend existing parameters object IN PLACE (preserve class) ---
  params <- flowCore::parameters(ff)             # S4 "parameters"
  pd     <- flowCore::pData(params)

  need_cols <- c("name", "desc", "range", "minRange", "maxRange")
  for (cc in need_cols) if (!cc %in% names(pd)) pd[[cc]] <- NA

  pd2 <- rbind(pd[, need_cols, drop = FALSE], new_pd[, need_cols, drop = FALSE])

  # Critical invariants:
  # 1) name must match exprs colnames
  pd2$name <- colnames(x2)

  # 2) rownames must be valid and sequential ($P1..$Pn) for many read.FCS-derived frames
  rownames(pd2) <- paste0("$P", seq_len(nrow(pd2)))

  # assign back into the SAME parameters object (keeps class)
  Biobase::pData(params) <- pd2

  # --- update the flowFrame in place ---
  ff2 <- ff

  # set expanded parameters first (so exprs replacement validates)
  methods::slot(ff2, "parameters") <- params

  # now replace exprs (colnames must match parameters$name)
  flowCore::exprs(ff2) <- x2

  # --- sanitize keywords for export ---
  kw <- flowCore::keyword(ff2)
  if (is.null(kw)) kw <- list()

  # remove all parameter-specific keywords ($P1N, $P1S, $P1R, $P1B, $P1E, etc.)
  drop_idx <- grepl("^\\$P\\d+", names(kw))
  if (any(drop_idx)) kw <- kw[!drop_idx]

  # set key header fields consistently
  kw[["$PAR"]] <- as.character(ncol(x2))
  kw[["$TOT"]] <- as.character(nrow(x2))

  flowCore::keyword(ff2) <- kw

  ff2
}

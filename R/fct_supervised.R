#' Scale a dataset using mean/SD parameters stored in a model bundle
#'
#' Applies z-scoring using per-feature means and standard deviations (typically
#' stored alongside a trained supervised model). Only the columns named in
#' `means` are used (and must exist in `df_mat`).
#'
#' Any missing/zero standard deviations are replaced with 1 to avoid division by
#' zero.
#'
#' @param df_mat A data.frame or matrix of observations (rows) by features (columns).
#' @param means Named numeric vector of feature means. Names must be columns in `df_mat`.
#' @param sds Named numeric vector of feature standard deviations. Names must be columns in `df_mat`.
#'
#' @return A numeric matrix with the same number of rows as `df_mat` and columns
#'   ordered to match `names(means)`. Values are scaled as:
#'   \deqn{(x - mean) / sd}.
#'
#' @noRd
scale_with_bundle <- function(df_mat, means, sds) {
  stopifnot(all(names(means) %in% colnames(df_mat)), all(names(sds) %in% colnames(df_mat)))
  X <- as.matrix(df_mat[, names(means), drop = FALSE])
  sds[is.na(sds) | sds == 0] <- 1
  X <- sweep(X, 2, means[names(means)], "-")
  X <- sweep(X, 2, sds[names(sds)],   "/")
  X
}

#' Auto-map model feature names to dataset column names
#'
#' Attempts to match each model feature (as stored in a trained model bundle) to
#' a column in an uploaded dataset. Matching is performed in stages:
#' \enumerate{
#'   \item Case-insensitive exact match on raw names.
#'   \item Exact match after applying `make.names()` to both.
#'   \item Exact match after canonicalising (lowercase, remove non-alphanumeric).
#'   \item Soft match where canonical strings start with each other.
#' }
#'
#' The result is a named character vector with one entry per model feature,
#' containing the matched dataset column name, or `NA` if no match was found.
#'
#' @param model_feats Character vector of model feature names (e.g. markers used in training).
#' @param dataset_cols Character vector of available dataset column names.
#'
#' @return Named character vector of length `length(model_feats)`.
#'   Names are `model_feats`; values are matched dataset columns or `NA_character_`.
#'
#' @examples
#' model_feats <- c("CD3", "CD45RA", "Live/Dead")
#' dataset_cols <- c("cd3", "CD45RA", "Live.Dead")
#' auto_map_features(model_feats, dataset_cols)
#' @noRd
auto_map_features <- function(model_feats, dataset_cols) {
  stopifnot(is.character(model_feats), is.character(dataset_cols))
  if (!length(model_feats) || !length(dataset_cols)) {
    return(setNames(rep(NA_character_, length(model_feats)), model_feats))
  }

  out <- setNames(rep(NA_character_, length(model_feats)), model_feats)

  ds_raw   <- dataset_cols
  ds_low   <- tolower(dataset_cols)
  ds_mn    <- make.names(dataset_cols)
  ds_canon <- canon(dataset_cols)

  mf_raw   <- model_feats
  mf_low   <- tolower(model_feats)
  mf_mn    <- make.names(model_feats)
  mf_canon <- canon(model_feats)

  # 1) case-insensitive exact match on raw names
  for (i in seq_along(mf_raw)) {
    hit <- which(ds_low == mf_low[i])
    if (length(hit)) { out[i] <- ds_raw[hit[1]]; next }
  }

  # 2) exact match after make.names on both sides
  for (i in seq_along(mf_raw)) if (is.na(out[i])) {
    hit <- which(ds_mn == mf_mn[i])
    if (length(hit)) { out[i] <- ds_raw[hit[1]]; next }
  }

  # 3) exact match on canonicalised strings (strip punctuation/case)
  for (i in seq_along(mf_raw)) if (is.na(out[i])) {
    hit <- which(ds_canon == mf_canon[i])
    if (length(hit)) { out[i] <- ds_raw[hit[1]]; next }
  }

  # 4) Controlled fallback: match "base + token" where token is A/H/W (Area/Height/Width)
  #    Example: "SSC-A" matches "SSC.A" and also "SSC-A" etc.
  #    But it will NOT match "SSC.A.Width" unless you asked for SSC-W.
  parse_base_token <- function(x) {
    x0 <- toupper(gsub("\\s+", "", x))
    parts <- unlist(strsplit(x0, "[^A-Z0-9]+"))
    parts <- parts[nzchar(parts)]

    # Normalize long tokens
    norm_tok <- function(tok) {
      tok <- toupper(tok)
      if (tok %in% c("AREA"))   return("A")
      if (tok %in% c("HEIGHT")) return("H")
      if (tok %in% c("WIDTH"))  return("W")
      tok
    }

    token <- NA_character_

    if (length(parts) >= 2) {
      last  <- norm_tok(parts[length(parts)])
      prev  <- norm_tok(parts[length(parts) - 1])

      # Case: "SSC.A.Width" -> parts = SSC, A, WIDTH
      # Interpret as base=SSC token=A (width is a suffix/stat; ignore it for mapping)
      if (last %in% c("W","H","A") && prev %in% c("A","H","W")) {
        token <- prev
        parts <- parts[seq_len(length(parts) - 2)]
      } else if (last %in% c("A","H","W")) {
        # Case: "SSC-A" or "SSC.A" -> base=SSC token=A
        token <- last
        parts <- parts[seq_len(length(parts) - 1)]
      }
    } else if (length(parts) == 1) {
      # no token
    }

    base <- paste(parts, collapse = "")
    list(base = base, token = token)
  }
  ds_bt <- lapply(ds_raw, parse_base_token)
  ds_base  <- vapply(ds_bt, `[[`, character(1), "base")
  ds_token <- vapply(ds_bt, function(z) z$token %||% NA_character_, character(1))

  for (i in seq_along(mf_raw)) if (is.na(out[i])) {
    q <- parse_base_token(mf_raw[i])
    q_base  <- q$base
    q_token <- q$token

    if (!nzchar(q_base)) next

    # candidates share base
    cand <- which(ds_base == q_base)
    if (!length(cand)) next

    if (!is.na(q_token)) {
      # prefer exact token match (A/H/W), but allow missing token as fallback
      cand2 <- cand[which(ds_token[cand] == q_token)]
      if (length(cand2)) { out[i] <- ds_raw[cand2[1]]; next }

      cand3 <- cand[which(is.na(ds_token[cand]))]
      if (length(cand3)) { out[i] <- ds_raw[cand3[1]]; next }
    } else {
      # no token requested: prefer missing token first, else take first base match
      cand0 <- cand[which(is.na(ds_token[cand]))]
      out[i] <- ds_raw[ if (length(cand0)) cand0[1] else cand[1] ]
      next
    }
  }

  out
}

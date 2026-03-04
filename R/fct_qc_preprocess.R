#' Flag debris events using low-FSC clustering
#'
#' Heuristic debris detector that assumes debris has low forward scatter.
#' If both FSC-A and FSC-H channels are available, it computes a 2D radial
#' magnitude score; otherwise it falls back to 1D on the available FSC channel.
#'
#' It then fits a 2-component Gaussian mixture model (mclust) on the score and
#' labels the lowest-mean component as debris. If mclust fails, falls back to
#' k-means (k=2) and labels the lowest-mean cluster as debris.
#'
#' Returned vector is aligned to the original `flowFrame` rows.
#'
#' @param ff A `flowCore::flowFrame`.
#' @return Integer vector of length `nrow(ff)`, where 1 indicates debris and 0 indicates keep.
#'   If channels are missing or there are too few finite events (<50), returns a vector of 0s
#'   (or integer(0) if the frame has 0 rows).
#' @noRd
detect_debris_flags <- function(ff) {
  n <- nrow(ff@exprs); if (!n) return(integer(0))
  cn <- flowCore::colnames(ff)
  fsc_a <- pick_fsc(cn, "A")
  fsc_h <- pick_fsc(cn, "H")

  # build scalar to cluster on
  score <- rep(NA_real_, n)
  if (!is.na(fsc_a) && !is.na(fsc_h)) {
    A <- suppressWarnings(as.numeric(ff@exprs[, fsc_a]))
    H <- suppressWarnings(as.numeric(ff@exprs[, fsc_h]))
    ok <- is.finite(A) & is.finite(H)
    if (sum(ok) < 50) return(integer(n))
    score[ok] <- log1p(sqrt(pmax(A[ok], 0)^2 + pmax(H[ok], 0)^2))
  } else {
    # if only one channel found use that
    one <- if (!is.na(fsc_a)) fsc_a else if (!is.na(fsc_h)) fsc_h else NA
    if (is.na(one)) return(integer(n))
    x <- suppressWarnings(as.numeric(ff@exprs[, one]))
    ok <- is.finite(x)
    if (sum(ok) < 50) return(integer(n))
    score[ok] <- log1p(pmax(x[ok], 0))
  }

  s_ok <- is.finite(score)
  s    <- score[s_ok]
  if (length(s) < 50) return(integer(n))

  # Try Mclust 2 Gaussians
  flg_ok <- try({
    fit <- mclust::Mclust(s, G = 2, verbose = FALSE)
    debris_comp <- which.min(fit$parameters$mean) # low FSC cluster
    as.integer(fit$classification == debris_comp)
  }, silent = TRUE)

  if (inherits(flg_ok, "try-error")) {
    # Fallback: kmeans (k=2) on scalar
    km <- stats::kmeans(s, centers = 2, iter.max = 50)
    cmeans <- tapply(s, km$cluster, mean)
    debris_label <- as.integer(names(which.min(cmeans)))
    flg_ok <- as.integer(km$cluster == debris_label)
  }

  flg <- integer(n)
  flg[s_ok] <- flg_ok
  # 1 = debris, 0 = keep
  flg
}

#' Remove debris events (if detected)
#'
#' Applies `detect_debris_flags()` and subsets the flowFrame to keep only non-debris events.
#' If no debris are detected (or detection is not possible), returns the input unchanged.
#'
#' @param ff A `flowCore::flowFrame`.
#' @return A `flowCore::flowFrame` with debris removed where possible.
#' @noRd
remove_debris_if_any <- function(ff) {
  flg <- detect_debris_flags(ff)
  if (!length(flg) || !any(flg == 1L)) return(ff)
  ff[which(flg == 0L), ]
}

#' Robust preprocessing: compensate (if needed) and logicle-transform fluorescence
#'
#' - Compensation: if a spillover matrix exists in the FCS keywords (SPILL/$SPILLOVER/spillover)
#'   and the frame does not appear already compensated, apply compensation on the common channels.
#' - Transformation: applies `estimateLogicle()` to fluorescence-like channels (non-FSC/SSC and non-Time),
#'   plus any channel whose NAME/DESC suggests viability (matches "viab|live|livedead").
#'
#' All operations are wrapped so failures return the original frame unchanged.
#'
#' @param ff A `flowCore::flowFrame`.
#' @return A `flowCore::flowFrame` (possibly compensated and/or transformed).
#' @noRd

preprocess_robust <- function(ff) {
  pd <- flowCore::pData(flowCore::parameters(ff))
  nm <- as.character(pd$name)
  already_comp <- any(grepl("^Comp", nm, ignore.case = TRUE))

  kw <- flowCore::keyword(ff)
  spill <- kw[["SPILL"]] %||% kw[["$SPILLOVER"]] %||% kw[["spillover"]] %||% NULL
  spill_ok <- is.matrix(spill) && length(colnames(spill)) > 1

  if (!already_comp && spill_ok) {
    common <- intersect(colnames(spill), flowCore::colnames(ff))
    if (length(common) > 1) {
      spill_use <- spill[common, common, drop = FALSE]
      ff <- tryCatch(flowCore::compensate(ff, spill_use), error = function(e) ff)
    }
  }

  # choose channels for logicle (fluorescence + viability-like even if not in spill)
  is_fsc_ssc <- grepl("^(FSC|SSC)", nm, ignore.case = TRUE)
  is_time    <- grepl("^Time$", nm,  ignore.case = TRUE)
  trans_channels <- if (spill_ok) intersect(colnames(spill), flowCore::colnames(ff))
  else flowCore::colnames(ff)[!(is_fsc_ssc | is_time)]

  pd_all <- flowCore::pData(flowCore::parameters(ff))
  desc_fb <- ifelse(!is.na(pd_all$desc) & nzchar(pd_all$desc), pd_all$desc, pd_all$name)
  pool    <- unique(c(desc_fb, pd_all$name))
  vi_idx  <- grepl("viab|live|livedead", canon(pool)) # find viability channel
  vi_names<- intersect(pool[vi_idx], flowCore::colnames(ff))
  trans_channels <- unique(c(trans_channels, vi_names))
  trans_channels <- trans_channels[trans_channels %in% flowCore::colnames(ff)]

  if (length(trans_channels)) {
    ff <- tryCatch(
      flowCore::transform(ff, flowCore::estimateLogicle(ff, channels = trans_channels)),
      error = function(e) ff
    )
  }
  ff
}

#' Remove margin events using PeacoQC
#'
#' Wrapper around `PeacoQC::RemoveMargins()` that returns the original flowFrame if it errors.
#'
#' @param ff A `flowCore::flowFrame`.
#' @return A `flowCore::flowFrame` with margins removed where possible.
#' @noRd
run_qc_remove_margins <- function(ff) {
  tryCatch(PeacoQC::RemoveMargins(ff, channels = 1:ncol(ff), save_fcs = FALSE, report = FALSE),
           error = function(e) ff)
}


#' Flag doublets using PeacoQC FSC-A vs FSC-H
#'
#' Returns an integer vector aligned to the input flowFrame rows:
#' 1 = doublet, 0 = keep.
#'
#' @param ff flowCore::flowFrame
#' @param channelA Optional FSC-A name
#' @param channelH Optional FSC-H name
#' @return Integer vector length nrow(ff) (or integer(0) if empty)
#' @noRd
detect_doublet_flags <- function(ff, channelA = NULL, channelH = NULL) {
  n <- nrow(ff@exprs); if (!n) return(integer(0))

  cn <- flowCore::colnames(ff)
  chA <- channelA %||% (pick_fsc(cn, "A") %||% NA_character_)
  chH <- channelH %||% (pick_fsc(cn, "H") %||% NA_character_)

  if (is.na(chA) || is.na(chH)) return(integer(n))
  if (!(chA %in% cn) || !(chH %in% cn)) return(integer(n))

  out <- tryCatch(
    PeacoQC::RemoveDoublets(ff, channel1 = chA, channel2 = chH, output = "full"),
    error = function(e) NULL
  )
  if (is.null(out)) return(integer(n))

  dbl <- out$indices_doublets
  if (is.null(dbl) || !length(dbl)) return(integer(n))

  dbl <- dbl[dbl >= 1 & dbl <= n]
  flg <- integer(n)
  flg[dbl] <- 1L
  flg
}

#' Flag bad QC events using PeacoQC::PeacoQC
#'
#' Returns integer vector aligned to input rows:
#' 1 = bad QC, 0 = keep.
#'
#' @param ff flowCore::flowFrame
#' @param min_events Minimum events required to run PeacoQC
#' @return Integer vector length nrow(ff) (or integer(0) if empty)
#' @noRd
detect_badqc_flags <- function(ff, min_events = 5000) {
  n <- nrow(ff@exprs); if (!n) return(integer(0))
  if (n < min_events) return(integer(n))  # skip => mark none bad

  cn <- flowCore::colnames(ff)
  sig_channels <- cn[!grepl("^Time$", cn, ignore.case = TRUE)]

  out <- tryCatch(
    PeacoQC::PeacoQC(ff, channels = sig_channels, report = FALSE, save_fcs = FALSE),
    error = function(e) NULL
  )
  if (is.null(out) || is.null(out$GoodCells)) return(integer(n))

  keep <- as.logical(out$GoodCells)
  keep <- keep[seq_len(n)]
  # Safety: if too aggressive, skip QC entirely
  if (mean(keep, na.rm = TRUE) < 0.2) return(integer(n))

  as.integer(!keep)  # 1 = bad qc
}

preprocess_flowframes <- function(frames, min_events_peacoqc = 5000) {
  if (is.null(frames) || !length(frames)) {
    return(list(frames = list(), flags = list(), counts = data.frame()))
  }

  n_files <- length(frames)

  out_frames <- vector("list", n_files)
  out_flags  <- vector("list", n_files)

  counts <- data.frame(
    file_i = seq_len(n_files),
    file_ok = FALSE,
    n_before = NA_integer_,
    n_after_transform = NA_integer_,
    n_debris = NA_integer_,
    n_doublet = NA_integer_,
    n_badqc = NA_integer_,
    err = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_files)) {
    ff0 <- frames[[i]]

    if (!inherits(ff0, "flowFrame")) {
      out_frames[i] <- list(NULL)  # <-- IMPORTANT: preserve length
      out_flags[i]  <- list(NULL)  # <-- IMPORTANT: preserve length
      counts$err[i] <- "not a flowFrame"
      next
    }

    res_i <- tryCatch(
      withCallingHandlers({
        n0 <- nrow(ff0@exprs)

        ff <- preprocess_robust(ff0)
        n1 <- nrow(ff@exprs)

        debris <- detect_debris_flags(ff)
        dbl    <- detect_doublet_flags(ff)
        badqc  <- detect_badqc_flags(ff, min_events = min_events_peacoqc)

        if (length(debris) != n1) debris <- integer(n1)
        if (length(dbl)    != n1) dbl    <- integer(n1)
        if (length(badqc)  != n1) badqc  <- integer(n1)

        flags_df <- data.frame(
          AutoFlow_debris  = as.integer(debris),
          AutoFlow_doublet = as.integer(dbl),
          AutoFlow_badqc   = as.integer(badqc),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )

        list(
          ff = ff,
          flags = flags_df,
          counts = list(
            n_before = n0,
            n_after_transform = n1,
            n_debris  = sum(debris == 1L, na.rm = TRUE),
            n_doublet = sum(dbl == 1L, na.rm = TRUE),
            n_badqc   = sum(badqc == 1L, na.rm = TRUE)
          )
        )
      }, warning = function(w) invokeRestart("muffleWarning")),
      error = function(e) {
        counts$err[i] <<- conditionMessage(e)
        NULL
      }
    )

    if (is.null(res_i)) {
      out_frames[i] <- list(NULL)  # <-- preserve length
      out_flags[i]  <- list(NULL)  # <-- preserve length
      next
    }

    out_frames[i] <- list(res_i$ff)     # <-- preserve length
    out_flags[i]  <- list(res_i$flags)  # <-- preserve length

    counts$file_ok[i]           <- TRUE
    counts$n_before[i]          <- as.integer(res_i$counts$n_before)
    counts$n_after_transform[i] <- as.integer(res_i$counts$n_after_transform)
    counts$n_debris[i]          <- as.integer(res_i$counts$n_debris)
    counts$n_doublet[i]         <- as.integer(res_i$counts$n_doublet)
    counts$n_badqc[i]           <- as.integer(res_i$counts$n_badqc)
    if (is.na(counts$err[i])) counts$err[i] <- ""
  }

  list(frames = out_frames, flags = out_flags, counts = counts)
}

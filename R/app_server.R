# app_server.R
#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#' @import shiny
#' @import Seurat
#' @import flowCore
#' @importFrom magrittr %>%
#' @noRd
app_server <- function(input, output, session) {

  ## Utilities/helpers
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # Read parameter table from fcs file
  param_table <- function(ff) {
    pd <- flowCore::pData(flowCore::parameters(ff))
    data.frame(
      name = as.character(pd$name),
      desc = {x <- as.character(pd$desc); x[is.na(x)] <- ""; x},
      stringsAsFactors = FALSE
    )
  }

  # Build NAME -> DESC maps across a list of frames (merge by NAME, label by majority DESC)
  build_name_label_maps <- function(frames) {
    tabs <- lapply(frames, param_table)
    union_names <- Reduce(union, lapply(tabs, `[[`, "name"))

    # collect candidate descs per name across files
    desc_map <- lapply(union_names, function(nm) {
      v <- unlist(lapply(tabs, function(tb) tb$desc[tb$name == nm]), use.names = FALSE)
      v[!is.na(v) & nzchar(v)]
    })
    names(desc_map) <- union_names

    # majority label (case-insensitive), fallback to name
    pick_desc <- function(v) {
      if (!length(v)) return(NA_character_)
      key <- tolower(v)
      v[ match(names(sort(table(key), decreasing = TRUE))[1], key) ]
    }
    chosen  <- vapply(desc_map, pick_desc, character(1))
    display <- ifelse(is.na(chosen) | !nzchar(chosen), union_names, chosen)

    # disambiguate duplicate labels by appending [NAME] (i.e. if raw file contains labels for fluorescent channels)
    dups <- duplicated(display) | duplicated(display, fromLast = TRUE)
    display[dups] <- paste0(display[dups], " [", union_names[dups], "]")

    list(
      union_names   = union_names,
      name_to_label = setNames(display, union_names),
      label_to_name = setNames(union_names, display)
    )
  }

  # for quick searches on strings later on
  canon <- function(x) gsub("[^a-z0-9]+", "", tolower(x))

  # find forward scatter area and height channels
  .pick_fsc <- function(cn, which = c("A","H")) {
    which <- match.arg(which)
    pat <- if (which == "A") "^FSC[._-]?(A|Area)$" else "^FSC[._-]?(H|Height)$"
    hit <- cn[grepl(pat, cn, ignore.case = TRUE)]
    if (length(hit)) hit[1] else NA_character_
  }

  # Debris detection: robust GMM (2D radial if A+H, else 1D; kmeans fallback)
  # assumption here that debris = low FSC-H/A
  detect_debris_flags <- function(ff) {
    n <- nrow(ff@exprs); if (!n) return(integer(0))
    cn <- flowCore::colnames(ff)
    fsc_a <- .pick_fsc(cn, "A")
    fsc_h <- .pick_fsc(cn, "H")

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

  remove_debris_if_any <- function(ff) {
    flg <- detect_debris_flags(ff)
    if (!length(flg) || !any(flg == 1L)) return(ff)
    ff[which(flg == 0L), ]
  }

  # Robust preprocessing (comp if spill exists & not already comped; logicle on fluoresc.)
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

  # QC helpers
  run_qc_remove_margins <- function(ff) {
    tryCatch(PeacoQC::RemoveMargins(ff, channels = 1:ncol(ff), save_fcs = FALSE, report = FALSE),
             error = function(e) ff)
  }

  # Doublet removal using FSC-A vs FSC-H .
  # Returns a subset flowFrame (keeps singlets). If channels are missing or it fails, returns input unchanged.
  run_qc_remove_doublets <- function(ff, channelA = NULL, channelH = NULL) {
    cn <- flowCore::colnames(ff)

    # Prefer typical names if present; otherwise try to infer via .pick_fsc()
    chA <- channelA %||% (.pick_fsc(cn, "A") %||% NA_character_)
    chH <- channelH %||% (.pick_fsc(cn, "H") %||% NA_character_)

    if (is.na(chA) || is.na(chH)) return(ff)
    if (!(chA %in% cn) || !(chH %in% cn)) return(ff)

    out <- tryCatch(
      PeacoQC::RemoveDoublets(ff, channel1 = chA, channel2 = chH, output = "full"),
      error = function(e) NULL
    )
    if (is.null(out)) return(ff)

    # 'indices_doublets' are the event indices flagged as doublets
    dbl <- out$indices_doublets
    if (is.null(dbl) || !length(dbl)) return(ff)

    keep <- setdiff(seq_len(nrow(ff@exprs)), dbl)
    if (!length(keep)) return(ff)

    ff[keep, ]
  }

  run_qc_peacoqc <- function(ff, min_events = 5000) {

    n <- nrow(ff@exprs)

    # Skip if too small
    if (n < min_events) {
      message("Skipping PeacoQC: only ", n, " events (< ", min_events, ")")
      return(ff)
    }

    # Exclude Time from QC channels
    cn <- flowCore::colnames(ff)
    sig_channels <- cn[!grepl("^Time$", cn, ignore.case = TRUE)]

    out <- tryCatch(
      PeacoQC::PeacoQC(
        ff,
        channels = sig_channels,
        report = FALSE,
        save_fcs = FALSE
      ),
      error = function(e) {
        warning("PeacoQC failed: ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(out)) {
      return(ff)
    }

    # Safety: don't allow pathological removals
    keep <- out$GoodCells
    if (!is.null(keep) && mean(keep, na.rm = TRUE) < 0.2) {
      warning("PeacoQC would remove >80% of events - skipping for this file")
      return(ff)
    }

    out$FinalFF %||% ff
  }


  # Marker vector across frames
  get_marker_vector <- function(frames, name_col) {
    if (is.null(frames) || !length(frames) || is.null(name_col)) return(numeric(0))
    unlist(lapply(frames, function(ff) {
      cn <- flowCore::colnames(ff)
      if (!(name_col %in% cn)) return(numeric(0))
      as.numeric(ff@exprs[, name_col])
    }), use.names = FALSE)
  }
  # min-max normalise func
  norm_minmax <- function(x){
    r <- range(x, na.rm = TRUE)
    if (!all(is.finite(r)) || diff(r) == 0) return(rep(0, length(x)))
    (x - r[1]) / (r[2] - r[1])
  }

  # auto threshold for viability marker -- auto find pos neg populations and fit as two gaussians
  # fall back to nothing if fails (<200 events, GMM fit failure, or component intersection can't be found)
  auto_threshold_gmm <- function(x) {
    x <- x[is.finite(x)]
    # Not enough data to infer a bimodal split
    if (length(x) < 200)
      return(NA_real_)
    m <- try(mclust::Mclust(x, G = 2), silent = TRUE)
    if (inherits(m, "try-error"))
      return(NA_real_)
    # Order components by increasing mean (low -> high)
    o  <- order(m$parameters$mean)
    mu <- m$parameters$mean[o]
    sd <- sqrt(m$parameters$variance$sigmasq)[o]
    pr <- m$parameters$pro[o]
    # Weighted component densities
    f1 <- function(t) pr[1] * dnorm(t, mu[1], sd[1])
    f2 <- function(t) pr[2] * dnorm(t, mu[2], sd[2])
    # Search within robust range (avoid outliers)
    rng <- range(quantile(x, c(0.001, 0.999), na.rm = TRUE))
    ur <- try(
      uniroot(function(t) f1(t) - f2(t), interval = rng)$root,
      silent = TRUE
    )
    if (inherits(ur, "try-error")) NA_real_ else ur
  }

  # for downstream -- add debris/viability/single-cell labels as columns to the original
  # flowframe. Allows writing pre-processed to .fsc files.
  ff_with_extra_cols <- function(ff, extra_mat, extra_desc = NULL) {
    stopifnot(is.matrix(extra_mat) || is.data.frame(extra_mat))
    extra_mat <- as.matrix(extra_mat)
    old_expr <- flowCore::exprs(ff)
    new_expr <- cbind(old_expr, extra_mat)

    old_par <- flowCore::pData(flowCore::parameters(ff))
    tpl <- old_par[1, , drop = FALSE]
    mk_row <- function(nm, ds, vals) {
      row <- tpl
      cols <- colnames(row)
      if ("name"      %in% cols) row[,"name"]      <- nm
      if ("desc"      %in% cols) row[,"desc"]      <- ds %||% nm
      rng <- range(vals, na.rm = TRUE); if (!all(is.finite(rng))) rng <- c(0, 1)
      span <- max(1, ceiling(rng[2] - rng[1]))
      if ("range"     %in% cols) row[,"range"]     <- span
      if ("minRange"  %in% cols) row[,"minRange"]  <- rng[1]
      if ("maxRange"  %in% cols) row[,"maxRange"]  <- rng[2]
      row
    }
    extra_names <- colnames(extra_mat)
    extra_descs <- extra_desc %||% extra_names
    add_par <- do.call(rbind, Map(function(nm, ds) mk_row(nm, ds, new_expr[, nm]), extra_names, extra_descs))

    par2 <- Biobase::AnnotatedDataFrame(
      data.frame(
        rbind(as.data.frame(old_par, check.names = FALSE, stringsAsFactors = FALSE),
              as.data.frame(add_par, check.names = FALSE, stringsAsFactors = FALSE)),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    )
    kw <- flowCore::keyword(ff)
    kw[["AUTOFLOW_EXTRA_CHANNELS"]] <- paste(extra_names, collapse = ",")
    flowCore::flowFrame(exprs = new_expr, parameters = par2, description = kw)
  }

  # helper for supervised - allow scaling test data with the bundled model info
  scale_with_bundle <- function(df_mat, means, sds) {
    stopifnot(all(names(means) %in% colnames(df_mat)), all(names(sds) %in% colnames(df_mat)))
    X <- as.matrix(df_mat[, names(means), drop = FALSE])
    sds[is.na(sds) | sds == 0] <- 1
    X <- sweep(X, 2, means[names(means)], "-")
    X <- sweep(X, 2, sds[names(sds)],   "/")
    X
  }

  # Auto-map model features to dataset columns
  # Returns a named character vector: names = model features, values = matched
  # dataset column names (or NA if no match).
  auto_map_features <- function(model_feats, dataset_cols) {
    stopifnot(is.character(model_feats), is.character(dataset_cols))
    if (!length(model_feats) || !length(dataset_cols)) {
      return(setNames(rep(NA_character_, length(model_feats)), model_feats))
    }

    canon <- function(x) tolower(gsub("[^a-z0-9]+", "", x))

    out <- setNames(rep(NA_character_, length(model_feats)), model_feats)

    ds_raw   <- dataset_cols
    ds_low   <- tolower(dataset_cols)
    ds_mn    <- make.names(dataset_cols)
    ds_canon <- canon(dataset_cols)

    mf_raw   <- model_feats
    mf_low   <- tolower(model_feats)
    mf_mn    <- make.names(model_feats)
    mf_canon <- canon(model_feats)

    # case-insensitive exact match on raw names
    for (i in seq_along(mf_raw)) {
      hit <- which(ds_low == mf_low[i])
      if (length(hit)) { out[i] <- ds_raw[hit[1]]; next }
    }

    # exact match after make.names on both sides
    for (i in seq_along(mf_raw)) if (is.na(out[i])) {
      hit <- which(ds_mn == mf_mn[i])
      if (length(hit)) { out[i] <- ds_raw[hit[1]]; next }
    }

    # exact match on canonicalised strings (strip punctuation/case)
    for (i in seq_along(mf_raw)) if (is.na(out[i])) {
      hit <- which(ds_canon == mf_canon[i])
      if (length(hit)) { out[i] <- ds_raw[hit[1]]; next }
    }

    # soft match: canonical startsWith in either direction (take first hit)
    for (i in seq_along(mf_raw)) if (is.na(out[i])) {
      hit <- which(startsWith(ds_canon, mf_canon[i]) | startsWith(mf_canon[i], ds_canon))
      if (length(hit)) out[i] <- ds_raw[hit[1]]
    }

    out
  }

  ## Reactive state
  volumes <- c(Home = fs::path_home(), shinyFiles::getVolumes()())
  shinyFiles::shinyFileChoose(input, "files", roots = volumes, session = session)

  raw_ff   <- reactiveVal(NULL)
  proc_ff  <- reactiveVal(NULL)
  use_ff   <- reactive({ proc_ff() %||% raw_ff() })

  # reactive val placeholders
  name_to_label <- reactiveVal(NULL)
  label_to_name <- reactiveVal(NULL)

  default_viab_name <- reactiveVal(NULL)
  auto_viab_thr_val <- reactiveVal(NULL)

  seurat_dat_comb <- reactiveVal(NULL)
  seurat_meta_comb<- reactiveVal(NULL)

  unsup_obj      <- reactiveVal(NULL)
  supervised_obj <- reactiveVal(NULL)

  model_bundle  <- reactiveVal(NULL)
  model_features<- reactiveVal(NULL)
  scaling_means <- reactiveVal(NULL)
  scaling_sds   <- reactiveVal(NULL)
  class_levels  <- reactiveVal(NULL)
  feature_map   <- reactiveVal(NULL)
  viab_thr_rv  <- reactiveVal(NULL)     # the actual threshold used everywhere
  viab_thr_src <- reactiveVal("auto")   # "auto" or "manual"
  # Debounced version of the slider input (wait 1s after last change)
  viab_thr_deb <- debounce(
    reactive(as.numeric(input$viability_threshold)),
    millis = 1000
  )
  ## File read in (RAW only)
  output$files <- renderPrint({
    if (is.integer(input$files)) {
      cat("No directory has been selected")
    } else {
      tmp <- shinyFiles::parseFilePaths(volumes, input$files)
      cat(paste(nrow(tmp), "files selected"))
    }
  })

  # prep for reading in .fcs and initialise
  files_all <- reactive({
    if (is.integer(input$files)) return(NULL)
    paths <- shinyFiles::parseFilePaths(volumes, input$files)$datapath
    paths <- as.character(paths)
    if (!length(paths)) return(NULL)

    R <- lapply(paths, function(f) {
      tryCatch(flowCore::read.FCS(f, alter.names = TRUE, transformation = NULL),
               error = function(e) { message("Error reading: ", f); NULL })
    })
    R <- Filter(Negate(is.null), R)
    if (!length(R)) return(NULL)

    raw_ff(R)
    proc_ff(NULL)
    name_to_label(NULL); label_to_name(NULL)
    seurat_dat_comb(NULL); seurat_meta_comb(NULL)
    unsup_obj(NULL); supervised_obj(NULL)
    default_viab_name(NULL); auto_viab_thr_val(NULL)
    TRUE
  })

  ## Preprocess toggle
  observeEvent(list(raw_ff(), input$preprocess), {
    req(raw_ff())

    if (!identical(input$preprocess, "Yes")) {
      # No preprocessing: just (re)build label maps from raw
      proc_ff(NULL)
      maps <- build_name_label_maps(raw_ff())
      name_to_label(maps$name_to_label)
      label_to_name(maps$label_to_name)

      ntl <- maps$name_to_label
      vi_idx_all  <- grep("viab|live|livedead", canon(unname(ntl)))# search viable column
      vi_idx_A <- vi_idx_all[ grepl("\\.A\\]", unname(ntl)[vi_idx_all]) | grepl("\\.A$", unname(ntl)[vi_idx_all]) ] #choose .A not .H

      default_viab_name(if (length(vi_idx_A)) names(ntl)[vi_idx_A[1]] else names(ntl)[1])

      auto_viab_thr_val(NULL)
      return(invisible())
    }

    R <- raw_ff()
    if (!length(R)) return(invisible())

    # Safe wrappers that never throw
    safe_RemoveMargins <- function(ff) {
      tryCatch(run_qc_remove_margins(ff), error = function(e) ff)
    }
    safe_remove_debris <- function(ff) {
      tryCatch(remove_debris_if_any(ff), error = function(e) ff)
    }
    safe_preprocess    <- function(ff) {
      tryCatch(preprocess_robust(ff), error = function(e) ff)
    }
    safe_PeacoQC       <- function(ff) {
      tryCatch({
        out <- run_qc_peacoqc(ff)
        if (is.null(out)) ff else out
      }, error = function(e) ff)
    }
    safe_remove_doublets <- function(ff) {
      tryCatch(run_qc_remove_doublets(ff), error = function(e) ff)
    }

    # Apply stepwise per file; track debris counts but don't fail
    n_before   <- vapply(R, function(ff) nrow(ff@exprs), numeric(1))
    deb_counts <- numeric(length(R))
    dbl_counts <- numeric(length(R))

    P <- lapply(seq_along(R), function(i){
      ff0 <- R[[i]]
      ff1 <- safe_RemoveMargins(ff0)
      # debris detection
      dc  <- tryCatch(sum(detect_debris_flags(ff1) == 1L, na.rm = TRUE), error = function(e) 0)
      deb_counts[i] <<- dc
      ff2 <- safe_remove_debris(ff1)
      n2   <- nrow(ff2@exprs)
      ff2b <- safe_remove_doublets(ff2)
      dbl_counts[i] <<- max(0, n2 - nrow(ff2b@exprs))
      ff3 <- safe_preprocess(ff2b)
      ff4 <- safe_PeacoQC(ff3)
      ff4
    })

    # If something catastrophic produced zero-length or invalid frames, fallback to raw
    ok_lengths <- vapply(P, function(ff) inherits(ff, "flowFrame") && nrow(ff@exprs) > 0, logical(1))
    if (!length(P) || !all(ok_lengths)) {
      proc_ff(NULL)  # ensure downstream uses raw
      showNotification("Pre-processing failed or produced empty output - proceeding with RAW data.", type = "warning", duration = 6)
    } else {
      proc_ff(P)

      if (sum(deb_counts) > 0) {
        pct <- 100 * sum(deb_counts) / sum(n_before)
        for (i in seq_along(R)) {
          if (deb_counts[i] > 0) {
            showNotification(sprintf("File %02d: removed %s debris events (%.2f%%).",
                                     i, format(deb_counts[i]), 100*deb_counts[i]/n_before[i]),
                             type = "message", duration = 4)
          }
        }
        showNotification(sprintf("Total debris removed: %s (%.2f%%).",
                                 format(sum(deb_counts)), pct),
                         type = "message", duration = 5)
      }
      if (sum(dbl_counts) > 0) {
        pct_dbl <- 100 * sum(dbl_counts) / sum(n_before)
        for (i in seq_along(R)) {
          if (dbl_counts[i] > 0) {
            showNotification(sprintf("File %02d: removed %s doublets (%.2f%%).",
                                     i, format(dbl_counts[i]), 100*dbl_counts[i]/n_before[i]),
                             type = "message", duration = 4)
          }
        }
        showNotification(sprintf("Total doublets removed: %s (%.2f%%).",
                                 format(sum(dbl_counts)), pct_dbl),
                         type = "message", duration = 5)
      }
    }

    # Build maps & default viability from whatever we ended up with (processed OR raw)
    FF_use <- proc_ff() %||% raw_ff()
    maps <- build_name_label_maps(FF_use)
    name_to_label(maps$name_to_label)
    label_to_name(maps$label_to_name)

    ntl <- maps$name_to_label
    vi_idx_all  <- grep("viab|live|livedead", canon(unname(ntl)))# search viable column
    vi_idx_A <- vi_idx_all[ grepl("\\.A\\]", unname(ntl)[vi_idx_all]) | grepl("\\.A$", unname(ntl)[vi_idx_all]) ] #choose .A not .H

    default_viab_name(if (length(vi_idx_A)) names(ntl)[vi_idx_A[1]] else names(ntl)[1])

    # Auto-threshold viability from the chosen set
    vv <- get_marker_vector(FF_use, default_viab_name())
    auto_viab_thr_val(if (length(vv)) auto_threshold_gmm(vv) else 0.5)
  }, ignoreInit = TRUE)

  ## Rebuild combined matrix - allowed combine event-by-channel matrix from loaded flowframes
  # harmonises channel names, remove duplicates within a file, row binds all event for one matrix + per event metadata
  rebuild_combined <- function() {
    CF <- use_ff(); req(CF)

    per_file_cols <- lapply(CF, function(ff) {
      m <- flowCore::exprs(ff)
      cn <- colnames(m) %||% character(0)
      if (!length(cn)) cn else make.unique(cn, sep = "__dupcol")
    })
    # Union schema across all files (allows combining files with different panels/derived cols)
    name_union <- Reduce(union, per_file_cols)
    if (!length(name_union)) stop("No columns found across files (name_union empty).")

    mats <- lapply(seq_along(CF), function(i) {
      ff <- CF[[i]]
      m  <- flowCore::exprs(ff)
      if (!is.matrix(m)) m <- as.matrix(m)
      cn <- make.unique(colnames(m) %||% character(0), sep = "__dupcol")
      colnames(m) <- cn
      # Add missing columns as NA so every file conforms to the same schema before row-binding
      miss <- setdiff(name_union, cn)
      if (length(miss)) {
        add <- matrix(NA_real_, nrow = nrow(m), ncol = length(miss))
        colnames(add) <- miss
        m <- cbind(m, add)
      }
      m <- m[, name_union, drop = FALSE]
      storage.mode(m) <- "double"
      as.data.frame(m, check.names = FALSE, stringsAsFactors = FALSE)
    })

    dat_comb <- as.data.frame(
      data.table::rbindlist(mats, use.names = TRUE, fill = TRUE),
      check.names = FALSE, stringsAsFactors = FALSE
    )

    paths <- shinyFiles::parseFilePaths(volumes, input$files)$datapath
    paths <- as.character(paths)
    meta_list <- lapply(seq_along(CF), function(i) {
      nc  <- nrow(CF[[i]]@exprs)
      data.frame(
        filename1 = rep(basename(paths[i] %||% sprintf("file_%d", i)), nc),
        filename2 = rep(NA_character_, nc),
        stringsAsFactors = FALSE, check.names = FALSE
      )
    })
    meta_comb <- as.data.frame(
      data.table::rbindlist(meta_list, use.names = TRUE, fill = TRUE),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    # Each event gets a metadata row with its source file; must stay row-aligned with dat_comb
    if (nrow(dat_comb) != nrow(meta_comb)) {
      stop(sprintf("Row mismatch after bind: data=%d, meta=%d", nrow(dat_comb), nrow(meta_comb)))
    }

    rownames(dat_comb)  <- paste0("cell_", seq_len(nrow(dat_comb)))
    rownames(meta_comb) <- rownames(dat_comb)

    maps <- build_name_label_maps(CF)
    seurat_dat_comb(dat_comb)
    seurat_meta_comb(meta_comb)
    name_to_label(maps$name_to_label)
    label_to_name(maps$label_to_name)
    # if input changed -> invalidate downstream objects so they are recomputed from the new combined matrix
    unsup_obj(NULL); supervised_obj(NULL)
  }
  # Rebuild combined table whenever the active flowFrames change
  observeEvent(use_ff(), { rebuild_combined() })

  ## Sidebar UI
  desc_choice_map <- reactive({
    ntl <- name_to_label(); req(ntl, seurat_dat_comb())
    present <- intersect(names(ntl), colnames(seurat_dat_comb()))
    setNames(present, ntl[present])  # labels=pretty, values=NAME
  })

  final_to_base <- reactive({
    name_to_label() %||% setNames(colnames(seurat_dat_comb() %||% data.frame()),
                                  colnames(seurat_dat_comb() %||% data.frame()))
  })

  output$columnSelector <- renderUI({
    req(input$model_type == "Unsupervised", files_all(), seurat_dat_comb(), desc_choice_map())
    choices <- desc_choice_map()
    checkboxGroupInput("columns", "Select columns (markers):",
                       choices = choices, selected = names(choices))
  })

  ## Viability
  output$viability_controls <- renderUI({
    req(files_all())
    ntl <- name_to_label(); req(ntl)

    # choose a sensible starting threshold from whatever data we have, runs on render
    FF_use <- use_ff(); req(FF_use)
    start_thr <- {
      vv0 <- get_marker_vector(FF_use, default_viab_name() %||% names(ntl)[1])
      thr0 <- if (length(vv0)) auto_threshold_gmm(vv0) else 0.5
      if (!is.finite(thr0)) 0.5 else thr0
    }
    if (is.null(viab_thr_rv())) {
      viab_thr_rv(start_thr)
      viab_thr_src("auto")
    }
    # for display in ui
    tagList(
      h4("Marker / Viability QC"),
      selectInput(
        "viability_marker_name",
        "Marker (by channel NAME; labels shown):",
        choices  = setNames(names(ntl), ntl),
        selected = default_viab_name() %||% names(ntl)[1]
      ),
      sliderInput(
        "viability_threshold",
        "Threshold (marker < threshold -> viable):",
        min = -1, max = 6,
        value = round(if (!is.finite(viab_thr_rv())) start_thr else viab_thr_rv(), 3),
        step = 0.001
      ),
      # Force re-computation of an automatic threshold for the currently selected marker
      actionButton("recalc_auto_thr", "Auto-threshold"),
      tags$div(style="margin-top:6px;",
               tags$small(textOutput("viability_summary", inline = TRUE)))
    )
  })
  # if the user moves the threshold use this, or recalc if they click auto-threshold
  observeEvent(viab_thr_deb(), {
    req(!is.null(viab_thr_deb()))
    viab_thr_rv(as.numeric(viab_thr_deb()))
    viab_thr_src("manual")
  }, ignoreInit = TRUE)

  # Keep the marker selector synced to the app-wide default viability marker
  observeEvent(default_viab_name(), {
    req(default_viab_name(), name_to_label())
    updateSelectInput(session, "viability_marker_name",
                      selected = default_viab_name(),
                      choices  = setNames(names(name_to_label()), name_to_label()))
  }, ignoreInit = TRUE)

  # When the selected marker changes (or th flowFrames change),
  # update the slider bounds to a robust range, and (optionally) recompute
  # the auto-threshold if we are still in "auto" mode.
  observeEvent(list(input$viability_marker_name, use_ff()), {
    req(use_ff(), input$viability_marker_name)
    vv <- get_marker_vector(use_ff(), input$viability_marker_name)
    vv <- vv[is.finite(vv)]
    if (!length(vv)) return()
    # Robust slider bounds (ignore extreme outliers)
    rng <- quantile(vv, c(0.001, 0.999), na.rm = TRUE)
    updateSliderInput(session, "viability_threshold",
                      min = round(rng[1], 3), max = round(rng[2], 3))
    # Only auto-update threshold if:
    #  it hasn't been set yet, or
    # the current threshold came from automation (not a user override)
    if (is.null(viab_thr_rv()) || identical(viab_thr_src(), "auto")) {
      # Prefer cached auto threshold if available; otherwise compute from current marker values.
      thr <- auto_viab_thr_val() %||% auto_threshold_gmm(vv)
      viab_thr_rv(thr)
      viab_thr_src("auto")
      updateSliderInput(session, "viability_threshold", value = round(thr, 3))
    }
  }, ignoreInit = TRUE)

  # recalculate the autothresholding if button clicked
  observeEvent(input$recalc_auto_thr, {
    req(use_ff(), input$viability_marker_name)
    vv <- get_marker_vector(use_ff(), input$viability_marker_name)
    vv <- vv[is.finite(vv)]
    if (!length(vv)) return()

    thr <- auto_threshold_gmm(vv)
    auto_viab_thr_val(thr)

    viab_thr_rv(thr)
    viab_thr_src("auto")

    updateSliderInput(session, "viability_threshold", value = round(thr, 3))
    showNotification(sprintf("Auto-threshold set to %.3f", thr), type = "message", duration = 3)
  })

  # summary test for thresholding
  output$viability_summary <- renderText({
    req(use_ff(), input$viability_marker_name, !is.null(input$viability_threshold))
    vv <- get_marker_vector(use_ff(), input$viability_marker_name)
    vv <- vv[is.finite(vv)]
    if (!length(vv)) return("No values available for this marker.")
    # Use stored threshold (single source of truth). Fall back to slider if needed.
    thr <- viab_thr_rv() %||% as.numeric(input$viability_threshold)
    pct_non_viable <- mean(vv >= thr, na.rm = TRUE) * 100
    total <- length(vv)
    sprintf("Non-viable at current threshold: %.1f%% (n=%s of %s)",
            pct_non_viable,
            format(round(total * pct_non_viable/100)),
            format(total))
  })

  # Current viability settings to be used by downstream analysis (not just UI)
  viability_settings <- reactive({
    marker <- input$viability_marker_name %||% default_viab_name()
    thr    <- viab_thr_rv()
    list(marker = marker, thr = thr)
  })

  ## Unsupervised
  # allow marker density plots of chosed cluster, or pairs of clusters for visual comparison
  output$unsup_qc_controls <- renderUI({
    req(unsup_obj(), seurat_dat_comb(), label_to_name())
    asn <- sort(unique(as.character(unsup_obj()@meta.data$assignment)))
    labels <- names(label_to_name())
    if (!length(labels)) return(NULL)
    tagList(
      h4("Unsupervised marker densities by assignment"),
      selectInput("unsup_marker_desc", "Marker (DESC):",
                  choices = labels, selected = labels[1]),
      radioButtons("unsup_compare", "Compare",
                   c("Selected vs all others" = "one_vs_all", "Pick two assignments" = "pair"),
                   selected = "one_vs_all"),
      conditionalPanel("input.unsup_compare == 'one_vs_all'",
                       selectInput("unsup_assign_a", "Selected assignment:", choices = asn)),
      conditionalPanel("input.unsup_compare == 'pair'",
                       fluidRow(
                         column(6, selectInput("unsup_assign_a", "Assignment A:", choices = asn)),
                         column(6, selectInput("unsup_assign_b", "Assignment B:", choices = asn))
                       )
      )
    )
  })

  observeEvent(input$runClustering, {
    req(input$model_type == "Unsupervised")
    req(seurat_dat_comb(), seurat_meta_comb())

    use_cols <- input$columns
    if (is.null(use_cols) || !length(use_cols)) use_cols <- colnames(seurat_dat_comb())
    use_cols <- intersect(use_cols, colnames(seurat_dat_comb()))
    validate(need(length(use_cols) > 0, "No valid columns selected for clustering."))

    # Build the matrix that becomes Seurat counts (cells x features)
    X <- seurat_dat_comb()[, use_cols, drop = FALSE]
    X <- apply(X, 2, norm_minmax) # apply minmax per col to full dataset

    ntl <- name_to_label()
    disp_names <- ntl[use_cols]
    na_idx <- which(is.na(disp_names)); if (length(na_idx)) disp_names[na_idx] <- use_cols[na_idx]
    disp_names <- make.unique(as.character(disp_names), sep = "__dup")
    colnames(X) <- disp_names

    meta <- seurat_meta_comb()

    # If subsampling is needed, KEEP row names in lockstep with X
    #if (nrow(X) > 1e5) {
    #  set.seed(1L)
    #  idx  <- sample(nrow(X), 1e5)
    #  X    <- X[idx, , drop = FALSE]
    #  meta <- meta[idx, , drop = FALSE]
    #}

    # Make 100% sure meta rownames match cell ids (rownames of X)
    stopifnot(identical(rownames(X), rownames(meta)))
    row_ids <- rownames(X)

    # Build Seurat; Seurat uses columns as cells
    seur <- SeuratObject::CreateSeuratObject(counts = t(X), meta.data = meta)

    # Store the original row ids we used to build this Seurat object
    seur@meta.data$autoflow_row_id <- row_ids
    # run_unsup_fun is a helper taking the users inputs for res and logfold
    seur <- run_unsupervised_func(
      seur,
      res     = as.numeric(input$res_umap),
      logfold = as.numeric(input$lf_umap)
    )

    unsup_obj(seur)
    showNotification("Unsupervised clustering complete.", type = "message")
  })

  ## Marker QC plots
  # plot for viability qc (or any marker, user can choose)
  output$marker_qc <- plotly::renderPlotly({
    req(files_all())
    validate(need(!is.null(use_ff()), "Load data to view viability densities."))
    req(input$viability_marker_name, input$viability_threshold, use_ff())

    vv <- get_marker_vector(use_ff(), input$viability_marker_name)
    vv <- vv[is.finite(vv)]
    validate(need(length(vv) >= 2, "Not enough values to draw a density."))

    d <- density(vv, n = 1024)
    ntl <- name_to_label() %||% c()
    desc_label <- ntl[[input$viability_marker_name]] %||% input$viability_marker_name
    # Base density plot for the marker
    p <- plotly::plot_ly()
    p <- plotly::add_lines(p, x = d$x, y = d$y, name = paste0("Density: ", desc_label))
    # Add vertical line for threshold
    p <- plotly::add_lines(p,
                           x = c(input$viability_threshold, input$viability_threshold),
                           y = c(0, max(d$y)),
                           name = paste0("Threshold = ", signif(input$viability_threshold, 4)))
    plotly::layout(p, xaxis = list(title = "Marker intensity"), yaxis = list(title = "Density"))
  })

  dens_or_hist <- function(p, x, nm) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) >= 2 && stats::var(x) > 0) {
      d <- stats::density(x)
      plotly::add_lines(p, x = d$x, y = d$y, name = nm)
    } else if (length(x) == 1) {
      plotly::add_histogram(p, x = x, name = paste0(nm, " (1 obs)"), nbinsx = 1)
    } else {
      p
    }
  }

  # render the unsup cluster expression plots
  output$unsup_marker_qc <- plotly::renderPlotly({
    req(unsup_obj(), seurat_dat_comb(), label_to_name())

    labels <- names(label_to_name()); req(length(labels) > 0)
    ui_label <- input$unsup_marker_desc %||% labels[1]
    validate(need(ui_label %in% labels, paste0("Marker '", ui_label, "' not present.")))

    mk_name <- label_to_name()[[ui_label]]
    validate(need(mk_name %in% colnames(seurat_dat_comb()), "Selected marker name not present in data."))

    # Expression values from the matrix (cells x NAME)
    vals_whole <- suppressWarnings(as.numeric(seurat_dat_comb()[, mk_name]))

    # Use the exact row ids that built the Seurat object (stored at clustering time)
    md  <- unsup_obj()@meta.data
    ids <- md$autoflow_row_id
    if (is.null(ids) || all(is.na(ids))) ids <- rownames(md)  # fallback

    # Intersect against the current combined matrix and align all vectors
    common <- intersect(ids, rownames(seurat_dat_comb()))
    validate(need(length(common) > 0, "No overlapping cells between clustering and matrix."))

    vals_common <- vals_whole[match(common, rownames(seurat_dat_comb()))]
    asg_common  <- as.character(md[match(common, ids), "assignment", drop = TRUE])

    df <- data.frame(
      val        = vals_common,
      assignment = asg_common,
      stringsAsFactors = FALSE, check.names = FALSE
    )
    df <- df[is.finite(df$val), , drop = FALSE]
    validate(need(nrow(df) > 0, "No finite expression values available for the selected marker."))

    # Build the plot
    mode <- input$unsup_compare %||% "one_vs_all"
    p <- plotly::plot_ly()
    n_layers <- 0L
    # choice for pairs of single cluster plot
    if (mode == "pair") {
      req(input$unsup_assign_a, input$unsup_assign_b)
      keep <- df$assignment %in% c(input$unsup_assign_a, input$unsup_assign_b)
      d2 <- df[keep, , drop = FALSE]
      for (g in unique(d2$assignment)) {
        xg <- d2$val[d2$assignment == g]
        p  <- dens_or_hist(p, xg, g)
        if (length(xg)) n_layers <- n_layers + 1L
      }
    } else {
      req(input$unsup_assign_a)
      d2 <- transform(df, group = ifelse(assignment == input$unsup_assign_a,
                                         input$unsup_assign_a, "All other"))
      for (g in unique(d2$group)) {
        xg <- d2$val[d2$group == g]
        p  <- dens_or_hist(p, xg, g)
        if (length(xg)) n_layers <- n_layers + 1L
      }
    }

    validate(need(n_layers > 0, "Nothing to plot: all selected groups had no finite values."))
    plotly::layout(
      p,
      xaxis = list(title = ui_label),
      yaxis = list(title = "Density / count"),
      legend = list(orientation = "h")
    )
  })


  ## Supervised
  output$supervised_controls <- renderUI({
    req(input$model_type == "Supervised")
    if (is.null(input$model_file)) return(helpText("Awaiting model bundle upload..."))
    req(model_bundle())
    b <- model_bundle()
    tags$div(
      tags$p(HTML(sprintf("<b>Bundle loaded.</b> %s",
                          if (!is.null(b$meta$dataset)) paste0("Dataset: ", b$meta$dataset) else ""))),
      tags$p(sprintf("Features in bundle: %d", length(model_features() %||% character(0))))
    )
  })

  # user uploads the model bundle RDS: read it in and extract the relevant info for later use (features, scaling, levels)
  observeEvent(input$model_file, {
    req(input$model_type == "Supervised", input$model_file$datapath)
    b <- tryCatch(readRDS(input$model_file$datapath), error = function(e) NULL)
    validate(need(!is.null(b), "Could not read the model bundle RDS."))

    feats <- b$features %||% NULL
    means <- b$scaling$means %||% NULL
    sds   <- b$scaling$sds   %||% NULL
    lvls  <- b$levels        %||% NULL
    validate(need(is.character(feats) && length(feats) > 0, "Bundle missing `features`."))
    validate(need(is.numeric(means) && is.numeric(sds), "Bundle missing `scaling$means`/`$sds`."))
    model_bundle(b); model_features(feats); scaling_means(means); scaling_sds(sds); class_levels(lvls)

    if (!is.null(seurat_dat_comb())) {
      feature_map(auto_map_features(feats, colnames(seurat_dat_comb())))
    } else {
      feature_map(setNames(rep(NA_character_, length(feats)), feats))
    }
  })

  # map features from model bundle to columns in the uploaded dataset; allow user override if auto-mapping fails or is incomplete
  output$feature_mapper_ui <- renderUI({
    req(input$model_type == "Supervised", model_features(), !is.null(seurat_dat_comb()))
    feats <- model_features()
    current_map <- feature_map()
    choices <- colnames(seurat_dat_comb())
    # If the current map is complete and matches the model features, show a success message instead of the mapping UI
    if (!is.null(current_map) && all(!is.na(current_map)) && identical(unname(current_map), feats)) {
      return(tags$p(HTML("<b>All model features were matched automatically.</b>")))
    }
    # setup options for user to manually map each model feature to a dataset column; pre-populate with any existing mapping if available
    rows <- lapply(seq_along(feats), function(i) {
      mf <- feats[i]
      sel <- current_map[[mf]] %||% NA_character_
      fluidRow(
        column(6, tags$code(mf)),
        column(6, selectizeInput(
          inputId = paste0("map_", i),
          label = NULL,
          choices = c("", choices),
          selected = if (!is.na(sel) && sel %in% choices) sel else "",
          options = list(placeholder = "Choose dataset column")
        ))
      )
    })
    # Whenever any of the mapping inputs change, update the reactive feature_map with the current selections
    observe({
      req(length(feats) > 0)
      new_map <- setNames(rep(NA_character_, length(feats)), feats)
      for (i in seq_along(feats)) {
        val <- input[[paste0("map_", i)]]
        new_map[i] <- if (!is.null(val) && nzchar(val)) val else NA_character_
      }
      feature_map(new_map)
    })

    tagList(
      tags$hr(),
      tags$h5("Feature mapping (model -> uploaded dataset)"),
      actionButton("autoMap", "Auto-map by name"),
      tags$small("  (case-insensitive exact; then make.names heuristic)"),
      tags$br(), tags$br(),
      div(style = "max-height: 300px; overflow-y: auto;", rows),
      tags$br(),
      uiOutput("mapping_summary")
    )
  })

  observeEvent(input$autoMap, {
    req(model_features(), !is.null(seurat_dat_comb()))
    feature_map(auto_map_features(model_features(), colnames(seurat_dat_comb())))
  })
  # summary of mapping status: how many features mapped, and which (if any) are missing
  output$mapping_summary <- renderUI({
    req(feature_map())
    fm <- feature_map(); n_ok <- sum(!is.na(fm) & nzchar(fm)); n_tot <- length(fm)
    missing <- names(fm)[is.na(fm) | !nzchar(fm)]
    tagList(
      tags$p(sprintf("Mapped %d / %d features.", n_ok, n_tot)),
      if (length(missing))
        tags$p(style = "color:#a94442;", paste("Missing mappings:", paste(missing, collapse = ", ")))
    )
  })
  # run the supervised prediction using the uploaded model bundle and the mapped features from the combined dataset
  # add results to the combined dataset and make a new Seurat object for downstream output
  observeEvent(input$runSupervised, {
    req(input$model_type == "Supervised")
    req(model_bundle(), model_features(), scaling_means(), scaling_sds())
    req(!is.null(seurat_dat_comb()), !is.null(seurat_meta_comb()))
    # Check that all model features have a mapped dataset column, and that those columns exist in the combined dataset
    feats <- model_features()
    fm <- feature_map()
    if (is.null(fm) || any(is.na(fm) | !nzchar(fm))) {
      feature_map(auto_map_features(feats, colnames(seurat_dat_comb())))
      fm <- feature_map()
    }
    missing <- names(fm)[is.na(fm) | !nzchar(fm)]
    validate(need(!length(missing), paste("Please map all features before running. Missing:", paste(missing, collapse = ", "))))
    # Now that we have a complete mapping, check that the mapped columns actually exist in the combined dataset
    data_cols <- unname(fm[feats])
    validate(need(all(data_cols %in% colnames(seurat_dat_comb())), "Mapped columns not found in uploaded data."))
    # Build the feature matrix for prediction (cells x features), applying the user-provided mapping to get from model features to dataset columns
    pred_df <- seurat_dat_comb()[, data_cols, drop = FALSE]
    colnames(pred_df) <- feats
    # Apply the scaling from the model bundle to the feature matrix (using the mapped columns from the combined dataset)
    Xs <- scale_with_bundle(pred_df, means = scaling_means(), sds = scaling_sds())

    b <- model_bundle()
    yhat <- NULL; probs <- NULL
    # Try to predict using $model if available; otherwise fall back to $caret_fit.
    # Handle different prediction output formats for ranger vs randomForest, and factor vs numeric predictions.
    if (!is.null(b$model)) {
      if (inherits(b$model, "ranger")) {
        pr <- predict(b$model, data = as.data.frame(Xs))
        if (!is.null(pr$predictions) && is.matrix(pr$predictions)) {
          probs <- pr$predictions
          cls <- colnames(probs)[max.col(probs, ties.method = "first")]
          yhat <- factor(cls, levels = colnames(probs))
        } else {
          yhat <- pr$predictions
        }
      } else if (inherits(b$model, "randomForest")) {
        pr <- predict(b$model, newdata = as.data.frame(Xs), type = "prob")
        if (is.matrix(pr) || is.data.frame(pr)) {
          probs <- as.matrix(pr)
          cls <- colnames(probs)[max.col(probs, ties.method = "first")]
          yhat <- factor(cls, levels = colnames(probs))
        } else {
          yhat <- predict(b$model, newdata = as.data.frame(Xs), type = "response")
        }
      } else {
        yhat <- predict(b$model, newdata = as.data.frame(Xs))
      }
    } else if (!is.null(b$caret_fit)) {
      yhat <- predict(b$caret_fit, newdata = as.data.frame(Xs))
      p_try <- try(predict(b$caret_fit, newdata = as.data.frame(Xs), type = "prob"), silent = TRUE)
      if (!inherits(p_try, "try-error")) probs <- as.matrix(p_try)
    } else {
      showNotification("Bundle missing $model or $caret_fit.", type = "error")
      return(invisible(NULL))
    }

    if (!is.null(class_levels())) {
      yhat <- factor(as.character(yhat), levels = class_levels())
    }

    dat <- seurat_dat_comb()
    dat$predicted_class <- as.character(yhat)
    if (!is.null(probs)) {
      for (cc in colnames(probs)) {
        dat[[paste0("pred_prob_", make.names(cc))]] <- as.numeric(probs[, cc])
      }
    }
    seurat_dat_comb(dat)
    # todo update naming convention so not seurat for supervised
    cm <- t(as.matrix(pred_df))
    stopifnot(ncol(cm) == nrow(seurat_meta_comb()))
    colnames(cm) <- rownames(seurat_meta_comb())
    seur_sup <- SeuratObject::CreateSeuratObject(counts = cm, meta.data = seurat_meta_comb())
    seur_sup@meta.data$assignment <- seurat_dat_comb()$predicted_class
    supervised_obj(seur_sup)

    showNotification("Supervised predictions completed.", type = "message")
  })

  ## Outputs
  out_dat_reactive <- reactive({
    if (identical(input$model_type, "Supervised")) {
      req(supervised_obj()); supervised_obj()
    } else {
      req(unsup_obj()); unsup_obj()
    }
  })
  # UMAP plot if available, otherwise bar plot of cluster counts; color by cluster/assignment
  output$plotdata <- plotly::renderPlotly({
    req(out_dat_reactive())
    out_dat <- out_dat_reactive()
    if ("umap" %in% names(out_dat@reductions)) {
      umap_data <- out_dat[["umap"]]@cell.embeddings
      rownames(out_dat@meta.data) <- rownames(umap_data)
      df <- cbind(
        data.frame(
          umap1 = umap_data[, 1],
          umap2 = umap_data[, 2],
          umap3 = if (ncol(umap_data) >= 3) umap_data[, 3] else 0
        ),
        out_dat@meta.data
      )
      plotly::plot_ly(
        data = df,
        x = ~umap1, y = ~umap2, z = ~umap3,
        color = ~assignment, type = "scatter3d", mode = "markers",
        marker = list(size = 2, width = 2),
        hoverinfo = "text", text = ~assignment, showlegend = FALSE
      )
    } else {
      md <- out_dat@meta.data
      agg <- md %>% dplyr::count(assignment, name = "count")
      plotly::plot_ly(data = agg, x = ~assignment, y = ~count, type = "bar")
    }
  })

  # setup summary tables for output and downloads; group by cluster/assignment and proliferation if available, count number of cells in each group
  summary_tab <- reactive({
    req(out_dat_reactive())
    out_dat <- out_dat_reactive()

    seurat_metadata <- out_dat@meta.data |>
      dplyr::select(-dplyr::matches("^n(Count|Feature)_"))

    filename_cols <- grep("^filename", colnames(seurat_metadata), value = TRUE)

    seurat_metadata <- seurat_metadata |>
      dplyr::mutate(
        Sample = if (length(filename_cols))
          apply(dplyr::select(seurat_metadata, dplyr::all_of(filename_cols)),
                1, paste, collapse = "/")
        else
          rownames(seurat_metadata)
      )

    if ("proliferation" %in% colnames(seurat_metadata)) {
      seurat_metadata |>
        dplyr::group_by(Sample, assignment, proliferation) |>
        dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
        dplyr::arrange(Sample, assignment)
    } else {
      seurat_metadata |>
        dplyr::group_by(Sample, assignment) |>
        dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
        dplyr::arrange(Sample, assignment)
    }
  })

  output$tablecounts <- DT::renderDataTable({
    req(summary_tab())
    summary_tab()
  })

  ## Downloads
  output$downloadcounts <- downloadHandler(
    filename = function() paste(Sys.Date(), "AutoFlow_longformat_counts.csv", sep = "_"),
    content = function(fname) {
      req(summary_tab())
      tab <- summary_tab() |>
        dplyr::mutate(count = as.integer(count)) |>
        dplyr::group_by(dplyr::across(-count)) |>
        dplyr::summarise(count = sum(count), .groups = "drop")
      utils::write.csv(data.table::data.table(tab), fname, row.names = FALSE)
    }
  )

  output$downloadcountsdelta <- downloadHandler(
    filename = function() paste(Sys.Date(), "AutoFlow_wideformat_counts.csv", sep = "_"),
    content = function(fname) {
      req(out_dat_reactive(), summary_tab())

      wide <- summary_tab() |>
        dplyr::mutate(
          assignment = paste0(as.character(assignment), " Count"),
          count = as.numeric(count)
        ) |>
        tidyr::pivot_wider(
          names_from  = assignment,
          values_from = count,
          values_fill = list(count = 0)
        )

      utils::write.csv(wide, fname, row.names = FALSE)
    }
  )
  observe({
    shinyjs::toggleState(
      "downloadprocessed",
      condition = !is.null(proc_ff()) && length(proc_ff()) > 0
    )
  })
  # Export processed FCS
  output$downloadprocessed <- downloadHandler(
    filename = function() paste0(Sys.Date(), "_processed_fcs.zip"),
    content = function(zipfile) {
      req(raw_ff())
      validate(need(identical(input$preprocess, "Yes"),
                    "Turn preprocessing on to export processed FCS."))
      P <- proc_ff(); req(P)

      vl_name <- input$viability_marker_name %||% default_viab_name()
      thr     <- input$viability_threshold %||% auto_viab_thr_val() %||% 0

      td <- tempfile("autoflow_fcs_"); dir.create(td, showWarnings = FALSE)

      for (i in seq_along(P)) {
        ff <- P[[i]]

        debris <- detect_debris_flags(ff)          # annotate
        if (!length(debris)) debris <- rep(0L, nrow(ff@exprs))

        viable <- rep(1L, nrow(ff@exprs))
        if (!is.null(vl_name) && vl_name %in% flowCore::colnames(ff)) {
          viable <- as.integer(ff@exprs[, vl_name] < thr)
        }

        single <- rep(1L, nrow(ff@exprs))
        if (all(c("FSC.A","FSC.H") %in% flowCore::colnames(ff))) {
          dres <- try(PeacoQC::RemoveDoublets(ff, channel1 = "FSC.A", channel2 = "FSC.H", output = "full"),
                      silent = TRUE)
          if (!inherits(dres, "try-error") && length(dres$indices_doublets)) {
            single[dres$indices_doublets] <- 0L
          }
        }

        add <- cbind(
          AutoFlow_debris     = debris,
          AutoFlow_viable     = viable,
          AutoFlow_singlecell = single
        )

        new_names <- colnames(add)
        collide <- intersect(new_names, flowCore::colnames(ff))
        if (length(collide)) new_names <- paste0(new_names, "_af")
        colnames(add) <- new_names

        ff2 <- ff_with_extra_cols(ff, add, extra_desc = paste0(new_names, " (AutoFlow)"))
        out_path <- file.path(td, sprintf("AutoFlow_processed_%02d.fcs", i))
        flowCore::write.FCS(ff2, filename = out_path)
      }

      old_wd <- setwd(td); on.exit(setwd(old_wd), add = TRUE)
      files <- list.files(td, pattern = "\\.fcs$", full.names = FALSE)
      if (length(files) == 0) stop("No processed FCS files generated.")
      if (requireNamespace("zip", quietly = TRUE)) {
        zip::zip(zipfile, files = files)
      } else {
        utils::zip(zipfile, files = files)
      }
    }
  )
  output$downloadseurat <- downloadHandler(
    filename = function() paste0(Sys.Date(), "_AutoFlow_seurat.rds"),
    content = function(file) {
      req(out_dat_reactive())
      saveRDS(out_dat_reactive(), file)
    }
  )
}

## Unsupervised helper
run_unsupervised_func <- function(flow_data, res = 0.5, logfold = 0.25, percentage_cells = 0.25, batch_correct = FALSE) {

  # Ensure Timepoint column exists
  if (!("Timepoint" %in% colnames(flow_data@meta.data))) {
    flow_data@meta.data$Timepoint <- NA
  }

  ref <- flow_data
  ref <- Seurat::NormalizeData(ref)
  ref <- Seurat::FindVariableFeatures(ref)
  ref <- Seurat::ScaleData(ref)
  ref <- Seurat::RunPCA(ref, features = VariableFeatures(ref))

  # Run UMAP in 3D
  ref <- Seurat::RunUMAP(ref, dims = 1:ncol(ref$pca), n.components = 3L, return.model = TRUE)
  ref <- Seurat::FindNeighbors(ref, dims = 1:ncol(ref$pca), reduction = "pca")

  # Prefer Leiden (algorithm = 4) if supported, else Louvain (algorithm = 1)
  algo_to_use <- if ("algorithm" %in% names(formals(Seurat::FindClusters))) 4 else 1
  message(sprintf("Running clustering with algorithm = %s (%s)",
                  algo_to_use, ifelse(algo_to_use == 4, "Leiden", "Louvain")))
  ref <- Seurat::FindClusters(ref, resolution = res, algorithm = algo_to_use)

  # Marker detection
  suppressWarnings({
    markers.pos     <- Seurat::FindAllMarkers(ref, only.pos = TRUE,  min.pct = percentage_cells, logfc.threshold = logfold)
    markers.pos_neg <- Seurat::FindAllMarkers(ref, only.pos = FALSE, min.pct = percentage_cells, logfc.threshold = logfold)
  })
  all_markers <- merge(markers.pos_neg, markers.pos, by = c("cluster", "gene"), all = TRUE)

  # Clean up marker labels
  if (nrow(all_markers)) {
    all_markers$gene <- sub("__dup.*$", "", all_markers$gene)
  }

  # Generate cluster assignments with +/- marker labels
  if (nrow(all_markers)) {
    all_table_neg <- table(all_markers[all_markers$avg_log2FC.x < 0, ]$cluster,
                           all_markers[all_markers$avg_log2FC.x < 0, ]$gene)
    all_table_pos <- table(all_markers[all_markers$avg_log2FC.x > 0, ]$cluster,
                           all_markers[all_markers$avg_log2FC.x > 0, ]$gene)

    all_labels <- lapply(1:nrow(all_table_neg), function(i) {
      labs_neg <- paste0(colnames(all_table_neg)[all_table_neg[i, ] == 1], "-")
      labs_pos <- paste0(colnames(all_table_pos)[all_table_pos[i, ] == 1], "+")
      c("cluster" = rownames(all_table_pos)[i],
        "assignment" = paste(paste(labs_pos, collapse = ""), paste(labs_neg, collapse = ""), collapse = ""))
    })

    if (length(all_labels)) {
      manual <- data.frame(do.call("rbind", all_labels))
      manual <- manual[order(as.numeric(manual$cluster)), ]
      ref@meta.data <- ref@meta.data %>%
        dplyr::left_join(manual, by = c("seurat_clusters" = "cluster"))
    }
  }

  # Ensure assignment column exists
  rownames(ref@meta.data) <- colnames(ref)
  if (!"assignment" %in% colnames(ref@meta.data)) {
    ref@meta.data$assignment <- as.character(Seurat::Idents(ref))
  }

  ref
}

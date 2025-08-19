# app_server.R
#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#' @import shiny
#' @import Seurat
#' @import dplyr
#' @import flowCore
#' @noRd
app_server <- function(input, output, session) {

  ## ───────────────────────────── Utilities / helpers ──────────────────────────

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # Rename matrix columns to DESC where available (fallback to NAME); make unique with NAME
  processFCS_desc <- function(ff) {
    pd <- flowCore::pData(flowCore::parameters(ff))
    nm <- as.character(pd$name)
    ds <- as.character(pd$desc); ds[is.na(ds)] <- ""
    final <- ifelse(nzchar(ds), ds, nm)
    dup <- duplicated(final) | duplicated(final, fromLast = TRUE)
    final[dup] <- paste0(final[dup], "_", nm[dup])
    colnames(ff) <- final
    ff
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

    is_fsc_ssc <- grepl("^(FSC|SSC)", nm, ignore.case = TRUE)
    is_time    <- grepl("^Time$", nm,  ignore.case = TRUE)
    trans_channels <- if (spill_ok) intersect(colnames(spill), flowCore::colnames(ff))
    else flowCore::colnames(ff)[!(is_fsc_ssc | is_time)]

    # include viability-ish channels even if not in spill
    pd_all <- flowCore::pData(flowCore::parameters(ff))
    desc_fb <- ifelse(!is.na(pd_all$desc) & nzchar(pd_all$desc), pd_all$desc, pd_all$name)
    canon   <- function(x) tolower(gsub("[^a-z0-9]+", "", x))
    pool    <- unique(c(desc_fb, pd_all$name))
    vi_idx  <- grepl("viab|live|livedead|7aad|amine|dapi|v525", canon(pool))
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
    tryCatch(PeacoQC::RemoveMargins(ff, channels = 1:ncol(ff)), error = function(e) ff)
  }
  run_qc_peacoqc <- function(ff) {
    out <- tryCatch(
      PeacoQC::PeacoQC(ff, channels = c("Time"), report = FALSE, save_fcs = FALSE),
      error = function(e) NULL
    )
    if (is.null(out)) ff else (out$FinalFF %||% ff)
  }

  get_marker_vector <- function(frames, name_col) {
    if (is.null(frames) || !length(frames) || is.null(name_col)) return(numeric(0))
    unlist(lapply(frames, function(ff) {
      cn <- flowCore::colnames(ff)
      if (!(name_col %in% cn)) return(numeric(0))
      as.numeric(ff@exprs[, name_col])
    }), use.names = FALSE)
  }

  norm_minmax <- function(x){
    r <- range(x, na.rm = TRUE)
    if (!all(is.finite(r)) || diff(r) == 0) return(rep(0, length(x)))
    (x - r[1]) / (r[2] - r[1])
  }

  auto_threshold_gmm <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 200) return(median(x, na.rm = TRUE))
    m <- try(mclust::Mclust(x, G = 2), silent = TRUE)
    if (inherits(m, "try-error")) return(median(x, na.rm = TRUE))
    o  <- order(m$parameters$mean)
    mu <- m$parameters$mean[o]
    sd <- sqrt(m$parameters$variance$sigmasq)[o]
    pr <- m$parameters$pro[o]
    f1 <- function(t) pr[1] * dnorm(t, mu[1], sd[1])
    f2 <- function(t) pr[2] * dnorm(t, mu[2], sd[2])
    rng <- range(quantile(x, c(0.001, 0.999), na.rm = TRUE))
    ur  <- try(uniroot(function(t) f1(t) - f2(t), interval = rng)$root, silent = TRUE)
    if (inherits(ur, "try-error")) median(x, na.rm = TRUE) else ur
  }

  scale_with_bundle <- function(df_mat, means, sds) {
    stopifnot(all(names(means) %in% colnames(df_mat)), all(names(sds) %in% colnames(df_mat)))
    X <- as.matrix(df_mat[, names(means), drop = FALSE])
    sds[is.na(sds) | sds == 0] <- 1
    X <- sweep(X, 2, means[names(means)], "-")
    X <- sweep(X, 2, sds[names(sds)],   "/")
    X
  }

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

    par2 <- Biobase::AnnotatedDataFrame(rbind(old_par, add_par))
    kw <- flowCore::keyword(ff)
    kw[["AUTOFLOW_EXTRA_CHANNELS"]] <- paste(extra_names, collapse = ",")
    flowCore::flowFrame(exprs = new_expr, parameters = par2, description = kw)
  }

  ## ───────────────────────────── Reactive state ──────────────────────────────

  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())
  shinyFiles::shinyFileChoose(input, "files", roots = volumes, session = session)

  raw_ff   <- reactiveVal(NULL)
  proc_ff  <- reactiveVal(NULL)
  use_ff   <- reactive({ proc_ff() %||% raw_ff() })

  raw_map  <- reactiveVal(NULL)  # names=DESC, values=name  (from RAW)
  proc_map <- reactiveVal(NULL)  # names=DESC, values=name  (from PROC)
  cur_map  <- reactive({ proc_map() %||% raw_map() })

  default_viab_name <- reactiveVal(NULL)
  auto_viab_thr     <- reactiveVal(NULL)

  seurat_dat_comb <- reactiveVal(NULL)  # DESC columns
  seurat_meta_comb<- reactiveVal(NULL)

  unsup_obj      <- reactiveVal(NULL)
  supervised_obj <- reactiveVal(NULL)

  model_bundle  <- reactiveVal(NULL)
  model_features<- reactiveVal(NULL)
  scaling_means <- reactiveVal(NULL)
  scaling_sds   <- reactiveVal(NULL)
  class_levels  <- reactiveVal(NULL)
  feature_map   <- reactiveVal(NULL)

  preproc_labels <- reactiveVal(NULL)

  ## ───────────────────────────── File ingest (RAW only) ──────────────────────

  output$files <- renderPrint({
    if (is.integer(input$files)) {
      cat("No directory has been selected")
    } else {
      tmp <- shinyFiles::parseFilePaths(volumes, input$files)
      cat(paste(nrow(tmp), "files selected"))
    }
  })

  files_all <- reactive({
    if (is.integer(input$files)) return(NULL)
    paths <- shinyFiles::parseFilePaths(volumes, input$files)$datapath
    paths <- as.character(paths)
    if (!length(paths)) return(NULL)

    # Read RAW (only)
    R <- lapply(paths, function(f) {
      tryCatch(flowCore::read.FCS(f, alter.names = TRUE, transformation = NULL),
               error = function(e) { message("Error reading: ", f); NULL })
    })
    R <- Filter(Negate(is.null), R)
    if (!length(R)) return(NULL)
    raw_ff(R)

    # Build RAW DESC map for UI
    pd <- flowCore::pData(flowCore::parameters(R[[1]]))
    nm <- as.character(pd$name)
    ds <- as.character(pd$desc); ds[is.na(ds)] <- ""
    desc_final <- ifelse(nzchar(ds), ds, nm)
    raw_map(setNames(nm, desc_final))

    # Default viability from RAW
    vi_idx <- grep("viab|live|livedead|7aad|amine|dapi|v525", tolower(names(raw_map())))
    default_viab_name(if (length(vi_idx)) unname(raw_map()[vi_idx[1]]) else unname(raw_map()[1]))

    # clear any previous processing; dedicated observer will run if needed
    proc_ff(NULL); proc_map(NULL); auto_viab_thr(NULL)

    TRUE
  })

  ## ───────────────────────────── Preprocess toggle (single place) ────────────

  observeEvent(list(raw_ff(), input$preprocess), {
    req(raw_ff())
    if (identical(input$preprocess, "Yes")) {
      # do preprocessing once here (not in files_all)
      Rm <- lapply(raw_ff(), run_qc_remove_margins)
      P0 <- lapply(Rm, preprocess_robust)
      P  <- lapply(P0, run_qc_peacoqc)
      proc_ff(P)

      # build PROC map
      pdp <- flowCore::pData(flowCore::parameters(P[[1]]))
      nmp <- as.character(pdp$name)
      dsp <- as.character(pdp$desc); dsp[is.na(dsp)] <- ""
      descp <- ifelse(nzchar(dsp), dsp, nmp)
      proc_map(setNames(nmp, descp))

      # default viability from PROC + auto thr
      vi_idx <- grep("viab|live|livedead|7aad|amine|dapi|v525", tolower(names(proc_map())))
      default_viab_name(if (length(vi_idx)) unname(proc_map()[vi_idx[1]]) else unname(proc_map()[1]))
      vv <- get_marker_vector(P, default_viab_name()); auto_viab_thr(auto_threshold_gmm(vv))

    } else {
      proc_ff(NULL); proc_map(NULL); auto_viab_thr(NULL)
      vi_idx <- grep("viab|live|livedead|7aad|amine|dapi|v525", tolower(names(raw_map())))
      default_viab_name(if (length(vi_idx)) unname(raw_map()[vi_idx[1]]) else unname(raw_map()[1]))
    }
  }, ignoreInit = TRUE)

  ## ───────────────── Rebuild DESC matrices whenever chosen frames change ────

  rebuild_combined <- function() {
    CF <- use_ff()
    req(CF)
    # ensure DESC cols
    CF_desc <- lapply(CF, processFCS_desc)
    mats <- lapply(CF_desc, flowCore::exprs)

    # filenames metadata
    paths <- shinyFiles::parseFilePaths(volumes, input$files)$datapath
    paths <- as.character(paths)
    fn_meta <- lapply(seq_along(mats), function(i) paste(basename(paths[i]), NA, sep = "/"))
    meta_list <- lapply(seq_along(mats), function(i) {
      nc  <- nrow(mats[[i]]); vec <- fn_meta[[i]]
      df  <- as.data.frame(t(matrix(vec, nrow = length(vec), ncol = nc)))
      colnames(df) <- paste0("filename", seq_len(ncol(df)))
      df
    })

    dat_comb  <- as.data.frame(do.call("rbind", mats))
    meta_comb <- do.call("rbind", meta_list)
    rownames(dat_comb) <- rownames(meta_comb)
    seurat_dat_comb(dat_comb)
    seurat_meta_comb(meta_comb)
    unsup_obj(NULL); supervised_obj(NULL)
  }

  observeEvent(use_ff(), { rebuild_combined() })

  ## ───────────────────────────── Sidebar UI (dynamic) ────────────────────────

  # Pretty label→final column map (for checkbox selector)
  desc_choice_map <- reactive({
    lst <- use_ff(); req(lst, length(lst) > 0)
    pd <- flowCore::pData(flowCore::parameters(lst[[1]]))
    nm <- as.character(pd$name)
    ds <- as.character(pd$desc); ds[is.na(ds)] <- ""
    base  <- ifelse(nzchar(ds), ds, nm)
    final <- base
    dup   <- duplicated(base) | duplicated(base, fromLast = TRUE)
    final[dup] <- paste0(base[dup], "_", nm[dup])
    lab <- base
    lab[dup] <- paste0(base[dup], " [", nm[dup], "]")
    setNames(final, lab)
  })

  # final → clean DESC
  final_to_base <- reactive({
    lst <- use_ff(); req(lst, length(lst) > 0)
    pd <- flowCore::pData(flowCore::parameters(lst[[1]]))
    nm <- as.character(pd$name)
    ds <- as.character(pd$desc); ds[is.na(ds)] <- ""
    base  <- ifelse(nzchar(ds), ds, nm)
    final <- base
    dup   <- duplicated(base) | duplicated(base, fromLast = TRUE)
    final[dup] <- paste0(base[dup], "_", nm[dup])
    setNames(base, final)
  })

  output$columnSelector <- renderUI({
    req(input$model_type == "Unsupervised", files_all(), seurat_dat_comb(), desc_choice_map())
    choices <- desc_choice_map()
    present  <- intersect(unname(choices), colnames(seurat_dat_comb()))
    selected <- if (length(input$columns %||% character())) intersect(input$columns, present) else present
    checkboxGroupInput("columns", "Select columns (markers):", choices = choices, selected = selected)
  })

  ## ───────────────────── Viability controls + live non-viable % ─────────────

  output$viability_controls <- renderUI({
    req(files_all())
    if (!identical(input$preprocess, "Yes")) return(NULL)
    tagList(
      h4("Marker / Viability QC"),
      selectInput("viability_marker_name",
                  "Marker (matrix 'name') for density:",
                  choices  = cur_map() %||% c(),
                  selected = default_viab_name()),
      sliderInput("viability_threshold",
                  "Threshold (marker < threshold → viable):",
                  min = -1, max = 5, value = auto_viab_thr() %||% 0.5, step = 0.001),
      actionButton("recalc_auto_thr", "Auto-threshold"),
      tags$div(style="margin-top:6px;", tags$small(textOutput("viability_summary", inline = TRUE)))
    )
  })

  # NEW: keep the viability picker in sync with the current map
  observeEvent(cur_map(), {
    req(cur_map())
    sel <- default_viab_name() %||% unname(cur_map()[1])
    updateSelectInput(session, "viability_marker_name",
                      choices = cur_map(), selected = sel)
  }, ignoreInit = FALSE)

  # Update slider bounds & default when marker or proc data change
  observeEvent(list(input$viability_marker_name, proc_ff()), {
    req(proc_ff(), input$viability_marker_name)
    vv <- get_marker_vector(proc_ff(), input$viability_marker_name)
    if (!length(vv)) return()
    rng <- stats::quantile(vv, c(0.001, 0.999), na.rm = TRUE)
    updateSliderInput(session, "viability_threshold",
                      min = round(rng[1], 3), max = round(rng[2], 3),
                      value = round(auto_viab_thr() %||% auto_threshold_gmm(vv), 3))
  })

  observeEvent(input$recalc_auto_thr, {
    req(proc_ff(), input$viability_marker_name)
    vv  <- get_marker_vector(proc_ff(), input$viability_marker_name)
    thr <- auto_threshold_gmm(vv)
    auto_viab_thr(thr)
    updateSliderInput(session, "viability_threshold", value = round(thr, 3))
  })

  output$viability_summary <- renderText({
    req(proc_ff(), input$viability_marker_name, !is.null(input$viability_threshold))
    vv <- get_marker_vector(proc_ff(), input$viability_marker_name)
    vv <- vv[is.finite(vv)]
    if (!length(vv)) return("No values available for this marker.")
    thr <- as.numeric(input$viability_threshold)
    pct_non_viable <- mean(vv >= thr, na.rm = TRUE) * 100
    total <- length(vv)
    sprintf("Non-viable at current threshold: %.1f%% (n=%s of %s)",
            pct_non_viable,
            format(round(total * pct_non_viable/100)),
            format(total))
  })

  ## ───────────────────────────── Unsupervised ────────────────────────────────

  # Controls UI (restored) + defaults
  output$unsup_qc_controls <- renderUI({
    req(unsup_obj(), seurat_dat_comb())
    asn <- sort(unique(as.character(unsup_obj()@meta.data$assignment)))
    tagList(
      h4("Unsupervised marker densities by assignment"),
      selectInput("unsup_marker_desc", "Marker (DESC):",
                  choices = colnames(seurat_dat_comb()), selected = colnames(seurat_dat_comb())[1]),
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

  # After clustering: seed sensible defaults for assignment inputs
  observeEvent(unsup_obj(), {
    req(unsup_obj())
    md <- unsup_obj()@meta.data
    tab <- sort(table(as.character(md$assignment)), decreasing = TRUE)
    if (length(tab)) {
      top1 <- names(tab)[1]
      top2 <- names(tab)[min(2, length(tab))]
      updateRadioButtons(session, "unsup_compare", selected = "one_vs_all")
      updateSelectInput(session, "unsup_assign_a", choices = names(tab), selected = top1)
      updateSelectInput(session, "unsup_assign_b", choices = names(tab), selected = top2)
    }
    # default marker: first column
    if (!is.null(seurat_dat_comb())) {
      updateSelectInput(session, "unsup_marker_desc",
                        choices = colnames(seurat_dat_comb()),
                        selected = colnames(seurat_dat_comb())[1])
    }
  })

  observeEvent(input$runClustering, {
    req(input$model_type == "Unsupervised")
    req(seurat_dat_comb(), seurat_meta_comb())

    use_cols <- input$columns
    if (is.null(use_cols) || !length(use_cols)) use_cols <- colnames(seurat_dat_comb())
    use_cols <- intersect(use_cols, colnames(seurat_dat_comb()))
    validate(need(length(use_cols) > 0, "No valid columns selected for clustering."))

    X <- seurat_dat_comb()[, use_cols, drop = FALSE]
    X <- apply(X, 2, norm_minmax)

    # rename to clean DESC (unique internally)
    ftb <- final_to_base()
    disp_names <- unname(ftb[use_cols])
    na_idx <- which(is.na(disp_names)); if (length(na_idx)) disp_names[na_idx] <- use_cols[na_idx]
    disp_names <- make.unique(disp_names, sep = "__dup")
    colnames(X) <- disp_names

    meta <- seurat_meta_comb()
    if (nrow(X) > 1e5) {
      idx <- sample(nrow(X), 1e5)
      X   <- X[idx, , drop = FALSE]
      meta<- meta[idx, , drop = FALSE]
    }

    seur <- SeuratObject::CreateSeuratObject(counts = t(X), meta.data = meta)
    seur <- run_unsupervised_func(seur, res = as.numeric(input$res_umap), logfold = as.numeric(input$lf_umap))
    unsup_obj(seur)
    showNotification("Unsupervised clustering complete.", type = "message")
  })

  ## ───────────────────────────── Marker QC plots ────────────────────────────

  output$marker_qc <- plotly::renderPlotly({
    req(files_all())
    validate(need(identical(input$preprocess, "Yes"),
                  "Enable pre-processing to view viability densities."))
    req(input$viability_marker_name, input$viability_threshold, proc_ff())
    vv <- get_marker_vector(proc_ff(), input$viability_marker_name)
    vv <- vv[is.finite(vv)]
    validate(need(length(vv) >= 2, "Not enough values to draw a density."))

    d <- density(vv, n = 1024)
    desc_label <- names(cur_map())[cur_map() == input$viability_marker_name][1] %||% input$viability_marker_name

    p <- plotly::plot_ly()
    p <- plotly::add_lines(p, x = d$x, y = d$y, name = paste0("Density: ", desc_label))
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

  output$unsup_marker_qc <- plotly::renderPlotly({
    req(unsup_obj(), seurat_dat_comb())

    mk_desc <- input$unsup_marker_desc %||% colnames(seurat_dat_comb())[1]
    validate(need(mk_desc %in% colnames(seurat_dat_comb()),
                  paste0("Marker '", mk_desc, "' not present.")))

    vals <- suppressWarnings(as.numeric(seurat_dat_comb()[, mk_desc]))
    md   <- unsup_obj()@meta.data

    common <- intersect(rownames(md), rownames(seurat_dat_comb()))
    validate(need(length(common) > 0, "No overlapping cells between matrix and clustering."))

    idx <- match(common, rownames(seurat_dat_comb()))
    df <- data.frame(
      val        = vals[idx],
      assignment = as.character(md[common, "assignment", drop = TRUE]),
      stringsAsFactors = FALSE, check.names = FALSE
    )

    df <- df[is.finite(df$val), , drop = FALSE]
    validate(need(nrow(df) > 0, "No finite expression values available for the selected marker."))

    mode <- input$unsup_compare %||% "one_vs_all"
    p <- plotly::plot_ly()
    n_layers <- 0L

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
      xaxis = list(title = (final_to_base()[[mk_desc]] %||% mk_desc)),
      yaxis = list(title = "Density / count")
    )
  })

  ## ───────────────────────────── Supervised (unchanged core) ─────────────────

  output$supervised_controls <- renderUI({
    req(input$model_type == "Supervised")
    if (is.null(input$model_file)) return(helpText("Awaiting model bundle upload…"))
    req(model_bundle())
    b <- model_bundle()
    tags$div(
      tags$p(HTML(sprintf("<b>Bundle loaded.</b> %s",
                          if (!is.null(b$meta$dataset)) paste0("Dataset: ", b$meta$dataset) else ""))),
      tags$p(sprintf("Features in bundle: %d", length(model_features() %||% character(0))))
    )
  })

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

  output$feature_mapper_ui <- renderUI({
    req(input$model_type == "Supervised", model_features(), !is.null(seurat_dat_comb()))
    feats <- model_features()
    current_map <- feature_map()
    choices <- colnames(seurat_dat_comb())

    if (!is.null(current_map) && all(!is.na(current_map)) && identical(unname(current_map), feats)) {
      return(tags$p(HTML("<b>All model features were matched automatically.</b>")))
    }

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
      tags$h5("Feature mapping (model → uploaded dataset)"),
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

  observeEvent(input$runSupervised, {
    req(input$model_type == "Supervised")
    req(model_bundle(), model_features(), scaling_means(), scaling_sds())
    req(!is.null(seurat_dat_comb()), !is.null(seurat_meta_comb()))

    feats <- model_features()
    fm <- feature_map()
    if (is.null(fm) || any(is.na(fm) | !nzchar(fm))) {
      feature_map(auto_map_features(feats, colnames(seurat_dat_comb())))
      fm <- feature_map()
    }
    missing <- names(fm)[is.na(fm) | !nzchar(fm)]
    validate(need(!length(missing), paste("Please map all features before running. Missing:", paste(missing, collapse = ", "))))

    data_cols <- unname(fm[feats])
    validate(need(all(data_cols %in% colnames(seurat_dat_comb())), "Mapped columns not found in uploaded data."))

    pred_df <- seurat_dat_comb()[, data_cols, drop = FALSE]
    colnames(pred_df) <- feats
    Xs <- scale_with_bundle(pred_df, means = scaling_means(), sds = scaling_sds())

    b <- model_bundle()
    yhat <- NULL; probs <- NULL

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

    cm <- t(as.matrix(pred_df))
    stopifnot(ncol(cm) == nrow(seurat_meta_comb()))
    colnames(cm) <- rownames(seurat_meta_comb())
    seur_sup <- SeuratObject::CreateSeuratObject(counts = cm, meta.data = seurat_meta_comb())
    seur_sup@meta.data$assignment <- seurat_dat_comb()$predicted_class
    supervised_obj(seur_sup)

    showNotification("Supervised predictions completed.", type = "message")
  })

  ## ───────────────────────────── Outputs ─────────────────────────────────────

  out_dat_reactive <- reactive({
    if (identical(input$model_type, "Supervised")) {
      req(supervised_obj()); supervised_obj()
    } else {
      req(unsup_obj()); unsup_obj()
    }
  })

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

  output$tablecounts <- DT::renderDataTable({
    req(out_dat_reactive())
    out_dat <- out_dat_reactive()
    seurat_metadata <- out_dat@meta.data

    filename_cols <- grep("^filename", colnames(seurat_metadata), value = TRUE)
    seurat_metadata <- seurat_metadata %>%
      dplyr::mutate(Sample = if (length(filename_cols))
        apply(dplyr::select(., dplyr::all_of(filename_cols)), 1, paste, collapse = "/")
        else rownames(seurat_metadata))

    if ("proliferation" %in% colnames(seurat_metadata)) {
      summary_counts <- seurat_metadata %>%
        dplyr::group_by(Sample, assignment, proliferation) %>%
        dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
        setNames(c("Sample", "assignment", "proliferation", "count"))

      summary_tab <<- seurat_metadata %>%
        dplyr::left_join(summary_counts, by = c("Sample", "assignment", "proliferation"), keep = FALSE) %>%
        dplyr::select(-dplyr::any_of("nCount_RNA"))
    } else {
      summary_counts <- seurat_metadata %>%
        dplyr::group_by(Sample, assignment) %>%
        dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
        setNames(c("Sample", "assignment", "count"))

      summary_tab <<- seurat_metadata %>%
        dplyr::left_join(summary_counts, by = c("Sample", "assignment"), keep = FALSE) %>%
        dplyr::select(-dplyr::any_of("nCount_RNA")) %>%
        dplyr::distinct()
    }
    summary_tab
  })

  output$plottreatment <- plotly::renderPlotly({
    req(out_dat_reactive())
    plotly::plot_ly(mtcars, x = ~wt, y = ~mpg, type = "scatter", mode = "markers")
  })

  ## ───────────────────────────── Downloads ──────────────────────────────────

  output$downloadcounts <- downloadHandler(
    filename = function() paste(Sys.Date(), "AutoFlow_counts.csv", sep = "_"),
    content = function(fname) {
      req(summary_tab)
      utils::write.csv(data.table::data.table(summary_tab), fname, row.names = FALSE)
    }
  )

  output$downloadcountsdelta <- downloadHandler(
    filename = function() paste(Sys.Date(), "DeltaFlow_counts.csv", sep = "_"),
    content = function(fname) {
      req(out_dat_reactive(), summary_tab)
      library(tidyr)
      wide <- summary_tab %>%
        mutate(assignment = paste0(assignment, " Count"),
               count = as.numeric(count),
               assignment = as.character(assignment)) %>%
        tidyr::pivot_wider(names_from = assignment, values_from = count, values_fill = list(count = 0))
      utils::write.csv(wide, fname, row.names = FALSE)
    }
  )

  # Export processed FCS (one per file) as a ZIP
  output$downloadprocessed <- downloadHandler(
    filename = function() paste0(Sys.Date(), "_processed_fcs.zip"),
    content = function(zipfile) {
      req(raw_ff())
      validate(need(identical(input$preprocess, "Yes"),
                    "Turn preprocessing on to export processed FCS."))
      P <- proc_ff(); req(P)

      vl_name <- input$viability_marker_name %||% default_viab_name()
      thr     <- input$viability_threshold %||% auto_viab_thr() %||% 0

      td <- tempfile("autoflow_fcs_"); dir.create(td, showWarnings = FALSE)

      for (i in seq_along(P)) {
        ff <- P[[i]]
        debris <- rep(0L, nrow(ff@exprs))
        if ("FSC.A" %in% flowCore::colnames(ff)) {
          m <- try(mclust::Mclust(ff@exprs[, "FSC.A"], G = 3), silent = TRUE)
          if (!inherits(m, "try-error")) {
            debris_clus <- which.min(m$parameters$mean)
            debris <- as.integer(m$classification == debris_clus)
          }
        }
        viable <- rep(1L, nrow(ff@exprs))
        if (!is.null(vl_name) && vl_name %in% flowCore::colnames(ff)) {
          viable <- as.integer(ff@exprs[, vl_name] < thr)
        }
        single <- rep(1L, nrow(ff@exprs))
        if (all(c("FSC.A","FSC.H") %in% flowCore::colnames(ff))) {
          dres <- try(PeacoQC::RemoveDoublets(ff, channel1 = "FSC.A", channel2 = "FSC.H", output = "full"), silent = TRUE)
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
}


## ───────────────────────── Unsupervised helper ─────────────────────────────

run_unsupervised_func <- function(flow_data, res = 0.5, logfold = 0.25, percentage_cells = 0.25, batch_correct = FALSE) {
  library(Seurat); library(dplyr)

  if (!("Timepoint" %in% colnames(flow_data@meta.data))) {
    flow_data@meta.data$Timepoint <- NA
  }

  ref <- flow_data
  ref <- NormalizeData(ref)
  ref <- FindVariableFeatures(ref)
  ref <- ScaleData(ref)
  ref <- RunPCA(ref, "RNA", features = VariableFeatures(ref))
  ref <- RunUMAP(ref, dims = 1:ncol(ref$pca), n.components = 3L, return.model = TRUE)
  ref <- FindNeighbors(ref, dims = 1:ncol(ref$pca), reduction = "pca")
  ref <- FindClusters(ref, resolution = res)

  suppressWarnings({
    markers.pos     <- FindAllMarkers(ref, only.pos = TRUE,  min.pct = percentage_cells, logfc.threshold = logfold)
    markers.pos_neg <- FindAllMarkers(ref, only.pos = FALSE, min.pct = percentage_cells, logfc.threshold = logfold)
  })
  all_markers <- merge(markers.pos_neg, markers.pos, by = c("cluster", "gene"), all = TRUE)

  # keep plain DESC in labels (remove internal __dup suffix)
  if (nrow(all_markers)) {
    all_markers$gene <- sub("__dup.*$", "", all_markers$gene)
  }

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
  rownames(ref@meta.data) <- colnames(ref)
  if (!"assignment" %in% colnames(ref@meta.data)) ref@meta.data$assignment <- as.character(Seurat::Idents(ref))
  ref
}

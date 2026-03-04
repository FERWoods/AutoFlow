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
  qc_flags  <- reactiveVal(NULL)  # list of per-file data.frames (nrow = events)
  qc_counts <- reactiveVal(NULL)  # per-file summary (optional)
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
    qc_flags(NULL); qc_counts(NULL)
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
      qc_flags(NULL); qc_counts(NULL)
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

    res <- tryCatch(
      preprocess_flowframes(R, min_events_peacoqc = 5000),
      error = function(e) NULL
    )

    # catastrophic fallback
    if (is.null(res) || is.null(res$frames) || !length(res$frames)) {
      proc_ff(NULL)
      qc_flags(NULL); qc_counts(NULL)
      showNotification("Pre-processing failed - proceeding with RAW data.", type = "warning", duration = 6)
    } else {
      ok_lengths <- vapply(res$frames, function(ff) inherits(ff, "flowFrame") && nrow(ff@exprs) > 0, logical(1))
      if (!all(ok_lengths)) {
        proc_ff(NULL)
        qc_flags(NULL); qc_counts(NULL)
        showNotification("Pre-processing produced empty output - proceeding with RAW data.", type = "warning", duration = 6)
      } else {
        proc_ff(res$frames)
        qc_flags(res$flags)
        qc_counts(res$counts)

        # notifications from counts (optional)
        if (!is.null(res$counts) && nrow(res$counts)) {
          # debris
          if (sum(res$counts$n_debris, na.rm = TRUE) > 0) {
            for (i in seq_len(nrow(res$counts))) {
              if (res$counts$n_debris[i] > 0) {
                showNotification(
                  sprintf("File %02d: flagged %s debris events.",
                          i, format(res$counts$n_debris[i])),
                  type = "message", duration = 4
                )
              }
            }
          }
          # doublets
          if (sum(res$counts$n_doublet, na.rm = TRUE) > 0) {
            for (i in seq_len(nrow(res$counts))) {
              if (res$counts$n_doublet[i] > 0) {
                showNotification(
                  sprintf("File %02d: flagged %s doublets.",
                          i, format(res$counts$n_doublet[i])),
                  type = "message", duration = 4
                )
              }
            }
          }
          # badqc
          if (sum(res$counts$n_badqc, na.rm = TRUE) > 0) {
            for (i in seq_len(nrow(res$counts))) {
              if (res$counts$n_badqc[i] > 0) {
                showNotification(
                  sprintf("File %02d: flagged %s bad-QC events.",
                          i, format(res$counts$n_badqc[i])),
                  type = "message", duration = 4
                )
              }
            }
          }
        }
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

    flags_list <- qc_flags()  # may be NULL

    meta_list <- lapply(seq_along(CF), function(i) {
      nc <- nrow(CF[[i]]@exprs)

      base <- data.frame(
        filename1 = rep(basename(paths[i] %||% sprintf("file_%d", i)), nc),
        #filename2 = rep(NA_character_, nc),
        stringsAsFactors = FALSE, check.names = FALSE
      )

      need_cols <- c("AutoFlow_debris","AutoFlow_doublet","AutoFlow_badqc")

      if (!is.null(flags_list) && length(flags_list) >= i && !is.null(flags_list[[i]])) {
        fl <- flags_list[[i]]
        if (is.data.frame(fl) && nrow(fl) == nc && all(need_cols %in% names(fl))) {
          fl <- fl[, need_cols, drop = FALSE]
          cbind(base, fl)
        } else {
          cbind(
            base,
            data.frame(
              AutoFlow_debris  = integer(nc),
              AutoFlow_doublet = integer(nc),
              AutoFlow_badqc   = integer(nc),
              stringsAsFactors = FALSE, check.names = FALSE
            )
          )
        }
      } else {
        cbind(
          base,
          data.frame(
            AutoFlow_debris  = integer(nc),
            AutoFlow_doublet = integer(nc),
            AutoFlow_badqc   = integer(nc),
            stringsAsFactors = FALSE, check.names = FALSE
          )
        )
      }
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
  output$run_unsup_button <- renderUI({
    lab <- if (isTRUE(input$use_projection)) {
      "Run Unsupervised (fast: train+project+predict)"
    } else {
      "Run Unsupervised (full)"
    }
    actionButton("runClustering", lab)
  })
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

  observe({
    req(seurat_dat_comb())
    n <- nrow(seurat_dat_comb())
    shinyjs::toggleState("use_projection", condition = n > 120000)
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

    # Optional filtering when preprocessing was chosen
    if (identical(input$preprocess, "Yes")) {
      keep <- rep(TRUE, nrow(meta))

      # drop debris/doublets/badqc
      if ("AutoFlow_debris" %in% names(meta))  keep <- keep & (meta$AutoFlow_debris  == 0L)
      if ("AutoFlow_doublet" %in% names(meta)) keep <- keep & (meta$AutoFlow_doublet == 0L)
      if ("AutoFlow_badqc" %in% names(meta))   keep <- keep & (meta$AutoFlow_badqc   == 0L)

      # viability filter (computed on the fly; still traceable)
      vl_name <- input$viability_marker_name %||% default_viab_name()
      thr     <- viab_thr_rv() %||% as.numeric(input$viability_threshold)

      if (!is.null(vl_name) && vl_name %in% colnames(seurat_dat_comb()) && is.finite(thr)) {
        v <- suppressWarnings(as.numeric(seurat_dat_comb()[, vl_name]))
        viable_flag <- as.integer(v < thr) # 1 = viable, 0 = non-viable
        meta$AutoFlow_viable <- viable_flag
        keep <- keep & (meta$AutoFlow_viable == 1L)
      }

      # apply filter to X + meta together
      X    <- X[keep, , drop = FALSE]
      meta <- meta[keep, , drop = FALSE]

      validate(need(nrow(X) > 10, "Too few events left after preprocessing filters."))
    }
    stopifnot(identical(rownames(X), rownames(meta)))
    row_ids <- rownames(X)

    # Build Seurat; Seurat uses columns as cells
    seur <- SeuratObject::CreateSeuratObject(counts = t(X), meta.data = meta)

    # Store the original row ids we used to build this Seurat object
    seur@meta.data$autoflow_row_id <- row_ids
    # run_unsup_fun is a helper taking the users inputs for res and logfold
    seur <- seur <- run_unsupervised_func(
      seur,
      res     = as.numeric(input$res_umap),
      logfold = as.numeric(input$lf_umap),
      percentage_cells = 0.25,  # or expose in UI later if you want
      max_umap_train = if (isTRUE(input$use_projection)) as.integer(input$max_umap_train) else NULL,
      seed = 1L,  # or expose
      k_transfer = if (isTRUE(input$use_projection)) as.integer(input$k_transfer) else 15L,
      min_transfer_conf = if (isTRUE(input$use_projection)) as.numeric(input$min_transfer_conf) else NULL,
      low_conf_label = if (isTRUE(input$use_projection)) input$low_conf_label else "Uncertain"
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
  qc_file_summary <- reactive({
    req(raw_ff())
    R <- raw_ff()
    P <- proc_ff() %||% R

    paths <- shinyFiles::parseFilePaths(volumes, input$files)$datapath
    paths <- as.character(paths)

    # counts from preprocess_flowframes()
    cc <- qc_counts()

    # viability settings
    vl_name <- input$viability_marker_name %||% default_viab_name()
    thr     <- input$viability_threshold %||% auto_viab_thr_val() %||% 0

    # compute viable + totals per file using proc frames (row space matches flags)
    flags_list <- qc_flags()

    out <- lapply(seq_along(P), function(i) {
      ff <- P[[i]]
      if (!inherits(ff, "flowFrame")) return(NULL)
      n <- nrow(ff@exprs)

      fl <- NULL
      if (!is.null(flags_list) && length(flags_list) >= i) fl <- flags_list[[i]]
      if (is.null(fl) || !is.data.frame(fl) || nrow(fl) != n) {
        fl <- data.frame(
          AutoFlow_debris  = integer(n),
          AutoFlow_doublet = integer(n),
          AutoFlow_badqc   = integer(n),
          stringsAsFactors = FALSE
        )
      }

      keep <- (fl$AutoFlow_debris == 0L) &
        (fl$AutoFlow_doublet == 0L) &
        (fl$AutoFlow_badqc == 0L)

      # viable computed on *all* events (or you can do viable among keep only)
      viable <- rep(NA_integer_, n)
      if (!is.null(vl_name) && nzchar(vl_name) && vl_name %in% flowCore::colnames(ff)) {
        v <- suppressWarnings(as.numeric(ff@exprs[, vl_name]))
        viable <- as.integer(is.finite(v) & (v < thr))
      }

      n_viable_all  <- if (all(is.na(viable))) NA_integer_ else sum(viable == 1L, na.rm = TRUE)
      n_viable_keep <- if (all(is.na(viable))) NA_integer_ else sum(viable[keep] == 1L, na.rm = TRUE)

      # per-file counts: prefer qc_counts if present, else compute from flags
      n_deb <- if (!is.null(cc) && nrow(cc) >= i) cc$n_debris[i] else sum(fl$AutoFlow_debris == 1L)
      n_dbl <- if (!is.null(cc) && nrow(cc) >= i) cc$n_doublet[i] else sum(fl$AutoFlow_doublet == 1L)
      n_bad <- if (!is.null(cc) && nrow(cc) >= i) cc$n_badqc[i] else sum(fl$AutoFlow_badqc == 1L)

      data.frame(
        file_i = i,
        filename = basename(paths[i] %||% sprintf("file_%02d.fcs", i)),
        preprocess_on = identical(input$preprocess, "Yes"),
        n_events = n,
        n_debris = as.integer(n_deb),
        pct_debris = 100 * as.numeric(n_deb) / n,
        n_doublet = as.integer(n_dbl),
        pct_doublet = 100 * as.numeric(n_dbl) / n,
        n_badqc = as.integer(n_bad),
        pct_badqc = 100 * as.numeric(n_bad) / n,
        n_keep = sum(keep),
        pct_keep = 100 * sum(keep) / n,
        viability_marker = if (!is.null(vl_name)) vl_name else NA_character_,
        viability_threshold = as.numeric(thr),
        n_viable_all = as.integer(n_viable_all),
        pct_viable_all = if (is.na(n_viable_all)) NA_real_ else 100 * n_viable_all / n,
        n_viable_keep = as.integer(n_viable_keep),
        pct_viable_keep = if (is.na(n_viable_keep)) NA_real_ else 100 * n_viable_keep / sum(keep),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    })

    out <- Filter(Negate(is.null), out)
    if (!length(out)) return(data.frame())
    dplyr::bind_rows(out)
  })
  output$table_qcfiles <- DT::renderDataTable({
    req(qc_file_summary())
    qc_file_summary()
  })
  ## Downloads
  output$download_qcfiles <- downloadHandler(
    filename = function() paste(Sys.Date(), "AutoFlow_file_QC_summary.csv", sep = "_"),
    content = function(fname) {
      req(qc_file_summary())
      utils::write.csv(qc_file_summary(), fname, row.names = FALSE)
    }
  )
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
  output$downloadqcout <- downloadHandler(
    filename = function() paste(Sys.Date(), "AutoFlow_QC_flags.csv", sep = "_"),
    content = function(fname) {
      req(qc_file_summary())
      utils::write.csv(qc_file_summary(), fname, row.names = FALSE)
    }
  )
  observe({
    shinyjs::toggleState(
      "downloadprocessed",
      condition = !is.null(proc_ff()) && length(proc_ff()) > 0
    )
  })
  # Export PROCESSED FCS annotated with already-computed flags (light debug)
  output$downloadprocessed <- downloadHandler(
    filename = function() paste0(Sys.Date(), "_AutoFlow_processed_annotated_fcs.zip"),
    content = function(zipfile) {

      logd <- function(...) {
        message(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = ""))
      }

      # Compact per-file snapshot only when something goes wrong
      snap_ff <- function(i, ff, add = NULL, fl = NULL, tag = "") {
        logd("---- SNAPSHOT i=", i, " ", tag, " ----")

        expr <- tryCatch(flowCore::exprs(ff), error = function(e) NULL)
        if (!is.null(expr)) {
          logd("exprs: nrow=", nrow(expr), " ncol=", ncol(expr))
          cn <- colnames(expr) %||% character()
          logd("colnames head: ", paste(utils::head(cn, 8), collapse = ", "))
        } else {
          logd("exprs: <FAILED>")
        }

        pd <- tryCatch(flowCore::pData(flowCore::parameters(ff)), error = function(e) NULL)
        if (!is.null(pd)) {
          logd("pData: nrow=", nrow(pd), " ncol=", ncol(pd))
          if ("name" %in% names(pd)) {
            logd("pd$name head: ", paste(utils::head(pd$name, 8), collapse = ", "))
          }
          logd("rownames(pd) head: ", paste(utils::head(rownames(pd), 8), collapse = ", "))
        } else {
          logd("pData: <FAILED>")
        }

        kw <- tryCatch(flowCore::keyword(ff), error = function(e) NULL)
        if (!is.null(kw)) {
          has_pn <- any(grepl("^\\$P\\d+", names(kw) %||% character()))
          logd("keyword: $PAR=", kw[["$PAR"]], " $TOT=", kw[["$TOT"]], " has $Pn*=", has_pn)
        } else {
          logd("keyword: <FAILED>")
        }

        if (!is.null(fl)) {
          logd("flags: class=", paste(class(fl), collapse = "/"),
               if (is.data.frame(fl)) paste0(" nrow=", nrow(fl), " ncol=", ncol(fl)) else "")
          if (is.data.frame(fl)) logd("flags cols: ", paste(colnames(fl), collapse = ", "))
        }

        if (!is.null(add)) {
          logd("add: nrow=", nrow(add), " ncol=", ncol(add))
          logd("add cols: ", paste(colnames(add), collapse = ", "))
        }

        logd("---- END SNAPSHOT ----")
      }

      # If rlang exists, last_trace() is nicer than base traceback
      log_trace <- function() {
        if (requireNamespace("rlang", quietly = TRUE)) {
          tr <- tryCatch(rlang::last_trace(), error = function(e) NULL)
          if (!is.null(tr)) {
            logd("rlang::last_trace():")
            logd(paste(utils::capture.output(print(tr)), collapse = "\n"))
            return(invisible())
          }
        }
        logd("base traceback():")
        logd(paste(utils::capture.output(traceback()), collapse = "\n"))
      }

      # Preconditions
      validate(need(identical(input$preprocess, "Yes"),
                    "Turn preprocessing on to export processed annotated FCS."))

      P <- proc_ff()
      validate(need(!is.null(P) && length(P) > 0, "No processed flowFrames available to export."))

      flags_list <- qc_flags()
      validate(need(!is.null(flags_list) && length(flags_list) == length(P),
                    "QC flags not available (or wrong length). Re-run preprocessing."))

      paths <- shinyFiles::parseFilePaths(volumes, input$files)$datapath
      paths <- as.character(paths)
      if (length(paths) < length(P)) paths <- c(paths, rep(NA_character_, length(P) - length(paths)))

      vl_name <- input$viability_marker_name %||% default_viab_name()
      thr     <- input$viability_threshold %||% auto_viab_thr_val() %||% 0

      td <- tempfile("autoflow_fcs_")
      dir.create(td, showWarnings = FALSE, recursive = TRUE)

      logd("downloadprocessed: n_files=", length(P),
           "  vl_name=", vl_name, " thr=", thr,
           "  td=", td)

      out_files <- character(0)
      need_cols <- c("AutoFlow_debris", "AutoFlow_doublet", "AutoFlow_badqc")

      # throttle progress logging
      every_k <- 5L

      for (i in seq_along(P)) {
        ff <- P[[i]]
        fl <- flags_list[[i]]

        if (!inherits(ff, "flowFrame")) {
          logd("SKIP i=", i, ": not flowFrame")
          next
        }

        expr <- tryCatch(flowCore::exprs(ff), error = function(e) {
          logd("ERROR i=", i, ": exprs(ff) failed: ", conditionMessage(e))
          snap_ff(i, ff, fl = fl, tag = "[exprs failed]")
          log_trace()
          NULL
        })
        if (is.null(expr)) next

        n <- nrow(expr)
        if (!is.finite(n) || n <= 0) {
          logd("SKIP i=", i, ": n<=0")
          next
        }

        # build required vectors from precomputed flags
        if (!(is.data.frame(fl) && nrow(fl) == n && all(need_cols %in% names(fl)))) {
          logd("FATAL i=", i, ": flags misaligned (flags_rows=",
               if (is.data.frame(fl)) nrow(fl) else NA, " expected=", n, ")")
          snap_ff(i, ff, fl = fl, tag = "[flags misaligned]")
          stop(sprintf("Flags missing/misaligned for file %d.", i))
        }

        debris  <- as.integer(fl$AutoFlow_debris)
        doublet <- as.integer(fl$AutoFlow_doublet)
        badqc   <- as.integer(fl$AutoFlow_badqc)

        debris[!is.finite(debris)   | !(debris %in% c(0L,1L))]   <- 0L
        doublet[!is.finite(doublet) | !(doublet %in% c(0L,1L))] <- 0L
        badqc[!is.finite(badqc)     | !(badqc %in% c(0L,1L))]   <- 0L

        # viability on processed data
        viable <- rep.int(1L, n)
        if (!is.null(vl_name) && (vl_name %in% flowCore::colnames(ff))) {
          v <- suppressWarnings(as.numeric(expr[, vl_name]))
          ok <- is.finite(v)
          viable[ok] <- as.integer(v[ok] < thr)
        }

        add <- data.frame(
          AutoFlow_debris     = debris,
          AutoFlow_doublet    = doublet,
          AutoFlow_badqc      = badqc,
          AutoFlow_viable     = as.integer(viable),
          AutoFlow_singlecell = as.integer(doublet == 0L),
          check.names = FALSE
        )

        if (i %% every_k == 1L || i == length(P)) {
          logd("i=", i, "/", length(P), "  n=", n, "  p=", ncol(expr))
        }

        ff2 <- tryCatch(
          ff_add_cols_safe(ff, add, extra_desc = paste0(colnames(add), " (AutoFlow)")),
          error = function(e) {
            logd("ERROR i=", i, ": ff_add_cols_safe failed: ", conditionMessage(e))
            snap_ff(i, ff, add = add, fl = fl, tag = "[ff_add_cols_safe failed]")
            log_trace()
            NULL
          }
        )
        if (is.null(ff2)) next

        in_base <- basename(paths[i] %||% sprintf("file_%02d.fcs", i))
        in_base <- sub("\\.fcs$", "", in_base, ignore.case = TRUE)
        out_name <- sprintf("AutoFlow_processed__%02d__%s.fcs", i, in_base)
        out_path <- file.path(td, out_name)

        ok_write <- tryCatch({
          flowCore::write.FCS(ff2, filename = out_path)
          TRUE
        }, error = function(e) {
          logd("ERROR i=", i, ": write.FCS failed: ", conditionMessage(e))
          # snapshot after building ff2 is useful too:
          snap_ff(i, ff2, add = add, fl = fl, tag = "[write.FCS failed; snapshot ff2]")
          log_trace()
          FALSE
        })

        if (ok_write && file.exists(out_path)) {
          out_files <- c(out_files, out_path)
        }
      }

      logd("written=", length(out_files), "/", length(P))
      validate(need(length(out_files) > 0,
                    "No processed annotated FCS files were generated (see console logs)."))

      if (requireNamespace("zip", quietly = TRUE)) {
        zip::zipr(zipfile, files = out_files, root = td)
      } else {
        old <- getwd(); on.exit(setwd(old), add = TRUE)
        setwd(td)
        utils::zip(zipfile, files = basename(out_files))
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


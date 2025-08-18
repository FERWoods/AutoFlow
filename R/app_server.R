#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#' @import shiny
#' @import Seurat
#' @import dplyr
#' @import flowCore
#' @export
#' @noRd
app_server <- function(input, output, session) {
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())

  observe({ print(input$files) })
  shinyFiles::shinyFileChoose(input, "files", roots = volumes, session = session, restrictions = system.file(package = "base"))
  observe({ print(input$files) })

  output$files <- renderPrint({
    if (is.integer(input$files)) {
      cat("No directory has been selected")
    } else {
      tmp <- shinyFiles::parseFilePaths(volumes, input$files)
      cat(paste(nrow(tmp), "files selected"))
    }
  })

  shinyFiles::shinyFileChoose(input, "directory_meta", roots = volumes, session = session, restrictions = system.file(package = "base"))
  observe({ print(input$directory_meta) })

  output$meta_directory <- renderPrint({
    if (is.integer(input$directory_meta) & is.integer(input$directory)) {
      cat("Select metadata file and directory of .fcs files")
    } else {
      fileselected <- shinyFiles::parseFilePaths(volumes, input$directory_meta)
      renderText(as.character(fileselected$datapath))
    }
  })

  ## =========================== Robust preprocessing =========================
  .align_spill_to_ff <- function(ff, spill) {
    if (is.null(spill)) return(NULL)
    if (is.list(spill) && !is.null(spill$SPILL)) spill <- spill$SPILL
    if (!is.matrix(spill)) return(NULL)

    ff_cols <- colnames(ff)
    sp_cols <- colnames(spill)

    .match_once <- function(sp_name) {
      cands <- c(
        sp_name,
        gsub("-", ".", sp_name, fixed = TRUE),
        gsub("\\.", "-", sp_name),
        paste0("Comp-", gsub("\\.", "-", sp_name)),
        paste0("Comp-", gsub("-", ".", sp_name))
      )
      hit <- cands[cands %in% ff_cols]
      if (length(hit)) hit[1] else NA_character_
    }

    target <- vapply(sp_cols, .match_once, character(1))
    ok <- !is.na(target)

    # If poor mapping AND already looks compensated, skip
    looks_comp <- any(grepl("^Comp-", ff_cols))
    mapped_ratio <- mean(ok)
    if (mapped_ratio < 0.6 && looks_comp) {
      message(sprintf(
        "Compensation skipped: only %d/%d spill channels matched (%.0f%%) and data appears already compensated.",
        sum(ok), length(ok), 100 * mapped_ratio
      ))
      return(NULL)
    }

    if (!all(ok)) {
      message(sprintf(
        "Compensation with partial mapping: %d/%d spill channels matched (%.0f%%).",
        sum(ok), length(ok), 100 * mapped_ratio
      ))
    }

    sp_use <- spill[ok, ok, drop = FALSE]
    colnames(sp_use) <- rownames(sp_use) <- target[ok]
    sp_use
  }

  safe_preprocess <- function(ff) {
    cn <- colnames(ff)

    # If looks already compensated (Comp- channels), skip compensation
    looks_comp <- any(grepl("^Comp-", cn))
    if (!looks_comp) {
      spill <- flowCore::keyword(ff)$SPILL
      if (is.null(spill)) spill <- flowCore::keyword(ff)$`$SPILLOVER`
      spill_aligned <- .align_spill_to_ff(ff, spill)

      if (is.matrix(spill_aligned) && nrow(spill_aligned) > 0) {
        ff <- tryCatch(
          flowCore::compensate(ff, spill_aligned),
          error = function(e) {
            message("Compensate failed (using aligned spill): ", conditionMessage(e), " — proceeding without compensation.")
            ff
          }
        )
      } else {
        message("No usable SPILL mapping; skipping compensation.")
      }
    } else {
      message("Data appears compensated (Comp-* channels found); skipping compensation.")
    }

    # Choose fluorescence channels for logicle transform
    cn <- colnames(ff)
    trans_cols <- cn[ grepl("-A$", cn) & !grepl("FSC|SSC|Time", cn, ignore.case = TRUE) ]
    if (length(trans_cols)) {
      ff <- tryCatch(
        flowCore::transform(ff, flowCore::estimateLogicle(ff, trans_cols)),
        error = function(e) {
          message("Logicle transform failed: ", conditionMessage(e), " — leaving data untransformed.")
          ff
        }
      )
    } else {
      message("No fluorescence channels (-A) found for transformation; skipping logicle.")
    }
    ff
  }

  ## ============================== File ingest ===============================
  files_all <- reactive({
    if (!is.integer(input$files) && !is.integer(input$metadata_file)) {
      files <- shinyFiles::parseFilePaths(volumes, input$files)$datapath
      files <- as.character(files)

      all_read <- lapply(files, function(file) {
        tryCatch(
          flowCore::read.FCS(file, alter.names = TRUE, transformation = NULL),
          error = function(e) { message("Error reading file: ", file); NULL }
        )
      })
      all_read <- Filter(Negate(is.null), all_read)

      # Simple filename metadata
      fn_metadata <<- lapply(seq_along(all_read), function(i) {
        fcs_obj <- all_read[[i]]
        meta <- as.data.frame(fcs_obj@description)[1, , drop = FALSE]
        well <- if ("WELL.ID" %in% colnames(meta) && "FILENAME" %in% colnames(meta)) {
          ifelse(!is.na(meta$WELL.ID), meta$WELL.ID, basename(files[i]))
        } else { NA }
        paste(basename(files[i]), well, sep = "/")
      })

      if (input$preprocess == "Yes") {
        cat("\nRunning compensation (if needed) and transformation\n")

        # Robust preprocessing (comp + logicle) BEFORE any renaming
        all_ff <- lapply(all_read, safe_preprocess)

        # channel subset (unchanged)
        all_ff_a <- lapply(all_ff, channel_select)

        # QC stages
        all_fc_margin <- lapply(all_ff_a, function(fcs_obj) {
          tryCatch(PeacoQC::RemoveMargins(fcs_obj, channels = 1:ncol(fcs_obj)), error = function(e) fcs_obj)
        })
        all_fc_qc <- lapply(all_fc_margin, function(fcs_obj) {
          tryCatch(PeacoQC::PeacoQC(fcs_obj, report = FALSE, channels = c("Time"), save_fcs = FALSE)$FinalFF, error = function(e) fcs_obj)
        })

        # Debris removal
        all_fc_int <- lapply(all_fc_qc, function(fcs_obj) {
          tryCatch({
            mdl <- mclust::Mclust(fcs_obj@exprs[, "FSC.A"], G = 2)
            debris_clus <- which.min(mdl$parameters$mean)
            fcs_obj[mdl$classification != debris_clus, ]
          }, error = function(e) fcs_obj)
        })

        # Single-cell refinement + viability threshold (if present)
        sc_only <- lapply(seq_along(all_fc_int), function(i) {
          fit <- tryCatch(lm(all_fc_int[[i]]@exprs[, "FSC.H"] ~ all_fc_int[[i]]@exprs[, "FSC.A"]), error = function(e) NULL)
          if (!is.null(fit)) {
            slope <- coef(fit)[2]; intercept <- coef(fit)[1]
            remove_these <- remove_outliers(all_fc_int[[i]]@exprs[, "FSC.H"], all_fc_int[[i]]@exprs[, "FSC.A"], slope = slope, intercept = intercept)
            fcs_obj <- all_fc_int[[i]][-remove_these, ]
            if ("Viability" %in% colnames(fcs_obj@exprs)) {
              fcs_obj <- fcs_obj[fcs_obj@exprs[, "Viability"] < 2, ]
            }
            fcs_obj
          } else {
            all_fc_int[[i]]
          }
        })

        # Rename to desc after preprocessing
        sc_only_rename <<- lapply(sc_only, processFCS)

        # Bind expression + filename metadata
        seurat_comb <<- lapply(seq_along(sc_only_rename), function(i) {
          expr_data <- sc_only_rename[[i]]@exprs
          if (nrow(expr_data) > 0) {
            metadata_vector <- fn_metadata[[i]]
            metadata_df <- as.data.frame(t(matrix(metadata_vector, nrow = length(metadata_vector), ncol = nrow(expr_data))))
            colnames(metadata_df) <- paste0("filename", seq_len(ncol(metadata_df)))
            cbind(metadata_df, expr_data)
          } else {
            data.frame()
          }
        })
        seurat_comb_dat <- do.call(rbind, seurat_comb)

        seurat_meta_comb <<- seurat_comb_dat[, grepl("filename", colnames(seurat_comb_dat)), drop = FALSE]
        seurat_dat_comb  <<- seurat_comb_dat[, !grepl("filename", colnames(seurat_comb_dat)), drop = FALSE]

        # Merge user metadata if provided
        if (!is.null(input$metadata_file)) {
          print("User uploaded metadata file - merging!")
          meta_file <- read_file(input$metadata_file$datapath)
          common_cols <- find_common_columns(seurat_meta_comb, meta_file)
          merged_data <- merge(seurat_meta_comb, meta_file, by = common_cols, all.x = TRUE)
          seurat_meta_comb <<- merged_data
        }

        rownames(seurat_dat_comb) <<- rownames(seurat_meta_comb)
        return(seurat_dat_comb)

      } else {  # No preprocessing
        no_pp_rename <<- lapply(all_read, processFCS)
        all_exprs <- lapply(no_pp_rename, flowCore::exprs)
        meta_file <- lapply(seq_along(all_exprs), function(i) {
          metadata_vector <- fn_metadata[[i]]
          expr_data <- all_exprs[[i]]
          metadata_df <- as.data.frame(t(matrix(metadata_vector, nrow = length(metadata_vector), ncol = nrow(expr_data))))
          colnames(metadata_df) <- paste0("filename", seq_len(ncol(metadata_df)))
          metadata_df
        })
        seurat_dat_comb <<- as.data.frame(do.call("rbind", all_exprs))
        seurat_meta_comb <<- do.call("rbind", meta_file)
        rownames(seurat_dat_comb) <<- rownames(seurat_meta_comb)
        return(seurat_dat_comb)
      }
    }
  })

  ## ======================= Column selector (hide for Supervised) ============
  output$columnSelector <- renderUI({
    req(files_all())
    if (identical(input$model_type, "Supervised")) return(NULL)
    colnames <- colnames(seurat_dat_comb)
    checkboxGroupInput("columns", "Select columns for analysis:", choices = colnames, selected = colnames)
  })

  ## ======================== UNSUPERVISED (unchanged) ========================
  cluster_dat <- eventReactive(input$runClustering, {
    if (!is.integer(input$files) & input$model_type == "Unsupervised") {
      req(input$columns)
      seurat_dat_comb <- seurat_dat_comb[, (colnames(seurat_dat_comb) %in% input$columns)]
      seurat_dat_comb <- apply(seurat_dat_comb, 2, norm_minmax)

      # optional downsample
      if (nrow(seurat_dat_comb) > 1e5) {
        sample_rows <- sample(nrow(seurat_dat_comb), 1e5)
        seurat_dat_comb <- seurat_dat_comb[sample_rows, ]
        seurat_meta_comb <- seurat_meta_comb[sample_rows, ]
      }

      seurat_obj <- SeuratObject::CreateSeuratObject(counts = t(seurat_dat_comb), meta.data = seurat_meta_comb)
      cluster_dat <- run_unsupervised_func(seurat_obj, res = as.numeric(input$res_umap), logfold = as.numeric(input$lf_umap))
      return(cluster_dat)
    }
  })

  ## ======================== SUPERVISED ======================================
  model_bundle   <- reactiveVal(NULL)
  model_features <- reactiveVal(NULL)
  scaling_means  <- reactiveVal(NULL)
  scaling_sds    <- reactiveVal(NULL)
  class_levels   <- reactiveVal(NULL)
  feature_map    <- reactiveVal(NULL)
  supervised_obj <- reactiveVal(NULL)

  scale_with_bundle <- function(df_mat, means, sds) {
    stopifnot(all(names(means) %in% colnames(df_mat)), all(names(sds) %in% colnames(df_mat)))
    X <- as.matrix(df_mat[, names(means), drop = FALSE])
    sds[is.na(sds) | sds == 0] <- 1
    X <- sweep(X, 2, means[names(means)], "-")
    X <- sweep(X, 2, sds[names(sds)],   "/")
    X
  }
  auto_map_features <- function(model_feats, data_cols) {
    dc_lower <- tolower(data_cols)
    mf_lower <- tolower(model_feats)
    exact <- match(mf_lower, dc_lower)
    out <- rep(NA_character_, length(model_feats))
    out[!is.na(exact)] <- data_cols[exact[!is.na(exact)]]
    still_na <- which(is.na(out))
    if (length(still_na)) {
      dc_made <- make.names(data_cols, unique = FALSE)
      mf_made <- make.names(model_feats, unique = FALSE)
      mm <- match(mf_made[still_na], dc_made)
      idx <- which(!is.na(mm))
      if (length(idx)) out[still_na[idx]] <- data_cols[mm[idx]]
    }
    names(out) <- model_feats
    out
  }

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

    model_bundle(b)
    model_features(feats)
    scaling_means(means)
    scaling_sds(sds)
    class_levels(lvls)

    if (exists("seurat_dat_comb", inherits = TRUE) && is.data.frame(seurat_dat_comb)) {
      feature_map(auto_map_features(feats, colnames(seurat_dat_comb)))
    } else {
      feature_map(setNames(rep(NA_character_, length(feats)), feats))
    }
  })

  output$feature_mapper_ui <- renderUI({
    req(input$model_type == "Supervised", model_features())
    req(exists("seurat_dat_comb", inherits = TRUE) && is.data.frame(seurat_dat_comb))

    feats <- model_features()
    current_map <- feature_map()
    choices <- colnames(seurat_dat_comb)

    # Hide UI if perfectly matched (1:1 and identical order)
    if (!is.null(current_map) &&
        all(!is.na(current_map)) &&
        all(nzchar(current_map)) &&
        all(names(current_map) == feats) &&
        all(current_map == feats)) {
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
    req(model_features())
    req(exists("seurat_dat_comb", inherits = TRUE) && is.data.frame(seurat_dat_comb))
    feature_map(auto_map_features(model_features(), colnames(seurat_dat_comb)))
  })

  output$mapping_summary <- renderUI({
    req(feature_map())
    fm <- feature_map()
    n_ok <- sum(!is.na(fm) & nzchar(fm))
    n_tot <- length(fm)
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
    req(exists("seurat_dat_comb", inherits = TRUE) && is.data.frame(seurat_dat_comb))

    feats <- model_features()
    fm <- feature_map()

    if (is.null(fm) || any(is.na(fm) | !nzchar(fm))) {
      feature_map(auto_map_features(feats, colnames(seurat_dat_comb)))
      fm <- feature_map()
    }

    missing <- names(fm)[is.na(fm) | !nzchar(fm)]
    validate(need(!length(missing), paste("Please map all features before running. Missing:", paste(missing, collapse = ", "))))

    data_cols <- unname(fm[feats])
    validate(need(all(data_cols %in% colnames(seurat_dat_comb)), "Mapped columns not found in uploaded data."))

    pred_df <- seurat_dat_comb[, data_cols, drop = FALSE]
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

    # Attach predictions back to the wide table used for Seurat meta
    seurat_dat_comb$predicted_class <- as.character(yhat)
    if (!is.null(probs)) {
      for (cc in colnames(probs)) {
        seurat_dat_comb[[paste0("pred_prob_", make.names(cc))]] <- as.numeric(probs[, cc])
      }
    }

    # Build a minimal Seurat object (no UMAP) just to carry metadata and support downstream code
    counts_mat <- t(as.matrix(seurat_dat_comb[, feats, drop = FALSE]))
    stopifnot(ncol(counts_mat) == nrow(seurat_meta_comb))
    colnames(counts_mat) <- rownames(seurat_meta_comb)

    seur_sup <- SeuratObject::CreateSeuratObject(
      counts = counts_mat,
      meta.data = seurat_meta_comb
    )
    seur_sup@meta.data$assignment <- seurat_dat_comb$predicted_class
    supervised_obj(seur_sup)

    showNotification("Supervised predictions added to dataset.", type = "message")
  })

  ## ======================== Downstream consumers ============================
  out_dat_reactive <- reactive({
    if (identical(input$model_type, "Supervised")) {
      req(supervised_obj()); supervised_obj()
    } else {
      cluster_dat()
    }
  })

  observe({ out_dat <- out_dat_reactive() })

  output$plotdata <- plotly::renderPlotly({
    req(out_dat_reactive)
    out_dat <- out_dat_reactive()
    if ("umap" %in% names(out_dat@reductions)) {
      umap_data <- out_dat[["umap"]]@cell.embeddings
      rownames(out_dat@meta.data) <- rownames(umap_data)
      seurat_to_df <<- Seurat::FetchData(out_dat, vars = c("umap_1", "umap_2", "umap_3", colnames(out_dat@meta.data)))
      plotly::plot_ly(
        data = seurat_to_df,
        x = ~umap_1, y = ~umap_2, z = ~umap_3,
        color = ~assignment,
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 2, width = 2),
        hoverinfo = "text",
        size = 10, alpha = I(1), text = ~assignment, showlegend = FALSE
      )
    } else {
      md <- out_dat@meta.data
      agg <- md %>% dplyr::count(assignment, name = "count")
      plotly::plot_ly(data = agg, x = ~assignment, y = ~count, type = "bar")
    }
  })

  ## ======================== Counts table (fixed) ============================
  output$tablecounts <- DT::renderDataTable({
    req(out_dat_reactive)
    out_dat <- out_dat_reactive()
    if (!is.null(out_dat)) {
      seurat_metadata <- out_dat@meta.data

      # Build per-sample key from all filename* columns
      filename_cols <- grep("^filename", colnames(seurat_metadata), value = TRUE)
      if (!length(filename_cols)) {
        # Fallback: create a single file tag if none present
        seurat_metadata$filename1 <- "sample"
        filename_cols <- "filename1"
      }

      # A single "Sample" label per event
      seurat_metadata <- seurat_metadata %>%
        mutate(Sample = apply(dplyr::select(., dplyr::all_of(filename_cols)), 1, paste, collapse = "/"))

      # Pure aggregated table (no left_join back to per-event rows!)
      summary_tab <- seurat_metadata %>%
        group_by(Sample, assignment) %>%
        summarise(count = dplyr::n(), .groups = "drop") %>%
        arrange(Sample, desc(count))

      # expose globally for download handlers
      assign("summary_tab", summary_tab, envir = .GlobalEnv)

      summary_tab
    }
  })

  ## ============================== Downloads =================================
  output$downloadcounts <- downloadHandler(
    filename = function() paste(Sys.Date(), "AutoFlow_counts.csv", sep = "_"),
    content = function(fname) {
      req(exists("summary_tab", envir = .GlobalEnv))
      write.csv(get("summary_tab", envir = .GlobalEnv), fname, row.names = FALSE)
    }
  )

  output$downloadcountsdelta <- downloadHandler(
    filename = function() paste(Sys.Date(), "DeltaFlow_counts.csv", sep = "_"),
    content = function(fname) {
      req(exists("summary_tab", envir = .GlobalEnv))
      library(tidyr)
      summary_tab <- get("summary_tab", envir = .GlobalEnv)
      seurat_metadata_wide <- summary_tab %>%
        mutate(assignment = paste0(assignment, " Count")) %>%
        tidyr::pivot_wider(names_from = assignment, values_from = count, values_fill = list(count = 0))
      write.csv(seurat_metadata_wide, fname, row.names = FALSE)
    }
  )

  output$downloadseurat <- downloadHandler(
    filename = function() paste(Sys.Date(), "AutoFlow_seurat_obj.rds", sep = "_"),
    content = function(fname) {
      req(out_dat_reactive)
      out_dat <- out_dat_reactive()
      saveRDS(out_dat, fname)
    }
  )

  ## ===================== Plot selectors (unchanged) =========================
  output$select_color_column <- renderUI({
    req(out_dat_reactive)
    out_dat <- out_dat_reactive()
    selectInput("color_column", "Select Treatment for Coloring:",
                choices = colnames(out_dat@meta.data),
                selected = colnames(out_dat@meta.data)[1])
  })

  output$select_x_column <- renderUI({
    req(out_dat_reactive)
    out_dat <- out_dat_reactive()
    selectInput("x_column", "Select Timepoint Column (X-axis):",
                choices = colnames(out_dat@meta.data),
                selected = colnames(out_dat@meta.data)[2])
  })

  output$select_cells <- renderUI({
    req(out_dat_reactive)
    out_dat <- out_dat_reactive()
    selectInput("cell_assignment", "Select cell assignment for plot",
                choices = unique(out_dat@meta.data$assignment),
                selected = unique(out_dat@meta.data$assignment)[min(2, length(unique(out_dat@meta.data$assignment)))])
  })

  observe({
    req(out_dat_reactive())
    updateSelectInput(session, "color_column", choices = colnames(out_dat_reactive()@meta.data))
    updateSelectInput(session, "x_column", choices = colnames(out_dat_reactive()@meta.data))
    updateSelectInput(session, "cell_assignment", choices = unique(out_dat_reactive()@meta.data$assignment))
  })

  observeEvent(input$show_treatment_plot, {
    if (isTruthy(input$metadata_file) && isTruthy(input$show_treatment_plot)) {
      showTab(inputId = "tabs", target = "Treatment Plot")
    } else {
      hideTab(inputId = "tabs", target = "Treatment Plot")
    }
  })

  output$plottreatment <- plotly::renderPlotly({
    req(out_dat_reactive())
    plotly::plot_ly(mtcars, x = ~wt, y = ~mpg, type = "scatter", mode = "markers")
  })
}


# helper functions

# pre-process functions
preprocess_2 <- function(ff){
  ff <- flowCore::compensate(ff, flowCore::keyword(ff)$SPILL)
  ff <- flowCore::transform(ff,
                            flowCore::estimateLogicle(ff,
                                                      colnames(flowCore::keyword(ff)$SPILL)))
  #scale(ff)
}
preprocess_3 <- function(ff){
  ff <- flowCore::compensate(ff, flowCore::keyword(ff)$SPILL)
  ff <- flowCore::transform(ff,
                            flowCore::estimateLogicle(ff,
                                                      colnames(flowCore::keyword(ff)$SPILL)))
  ff <- flowStats::gaussNorm(flowSet(ff), channel.names = colnames(flowCore::keyword(ff)$SPILL))
  #scale(ff)
}
channel_select <- function(ff){
  library(data.table)
  ff = ff[,c((flowCore::colnames(ff) %like% ".A") | flowCore::colnames(ff) %like% "Time") | flowCore::colnames(ff) %like% "FSC.H"]

  ff
}
# Functions for .fsc data importing and mapping to metadata

# Below sections are to allow for cbind.fill -- package origin had depreciated
cbind_fill<-function(...,fill=NULL)
{
  inputs<-list(...)
  inputs<-lapply(inputs,vert)
  maxlength<-max(unlist(lapply(inputs,len)))
  bufferedInputs<-lapply(inputs,buffer,length.out=maxlength,fill,preserveClass=FALSE)
  return(Reduce(cbind.data.frame,bufferedInputs))
}

vert<-function(object)
{
  #result<-as.data.frame(cbind(as.matrix(object)))
  if(is.list(object))
    object<-cbind(object)
  object<-data.frame(object)

  return(object)
}

len <- function(data)
{
  result<-ifelse(is.null(nrow(data)),length(data),nrow(data))
  return(result)
}

buffer<-function(x,length.out=len(x),fill=NULL,preserveClass=TRUE)
{
  xclass<-class(x)
  input<-lapply(vert(x),unlist)
  results<-as.data.frame(lapply(input,rep,length.out=length.out))
  if(length.out>len(x) && !is.null(fill))
  {
    results<-t(results)
    results[(length(unlist(x))+1):length(unlist(results))]<-fill
    results<-t(results)
  }
  if(preserveClass)
    results<-as2(results,xclass)
  return(results)
}

# Now process each with the preprocess function which applies the set spillover file and normalises with biexponential logicle function
preprocess <- function(ff) {
  ff <- compensate(ff, ff@description$SPILL)
  ff <- transform(ff, flowCore::transformList(colnames(ff@description$SPILL), flowCore::logicleTransform()))
  ff@exprs <- scale(ff@exprs)
  return(ff)
}

# function to convert wellIDs to match those attached to flowjo objects
checkWellNames = function(wellNames) {
  #' Check and convert well names to the appropriate format
  #' E.g. A1 -> A01
  #'
  #' @param wellNames String vector (e.g. c("A1", "A2", "B3", "B4"))
  #'
  #' @return String vector
  #'
  o = wellNames
  rows = substr(wellNames, 1, 1)
  stopifnot(all(rows %in% toupper(letters)[1:8]))
  columns = as.integer(substr(wellNames, 2, 10))
  stopifnot(all(columns >= 1 & columns <= 12))
  columns = as.character(columns)
  columns = sapply(columns, function(x) {
    if (nchar(x) == 1) {
      return(paste0("0", x))
    } else {
      return(x)
    }
  })
  return(paste0(rows, columns))
}

# Function to extract labels from fcs file paired with workspace --
workspace_cell_labels <- function(flowset, workspace, cell_types=cell_types_fixed, ws_group){
  groups = fj_ws_get_sample_groups(open_flowjo_xml(workspace))#, execute = FALSE)
  group_pos = groups[match(ws_group, groups$groupName), 2] + 1 # returns the group ID for the matched group +1 as the numbering begins at 0
  gating_manual = GetFlowJoLabels(flowset,
                                  workspace,
                                  group = group_pos,
                                  cellTypes = cell_types,
                                  get_data = TRUE)
  manual_labels = do.call("c",lapply(gating_manual,function(x){x$manual}))
  rm(gating_manual)
  return(manual_labels)
}


# wiht this the condition is a on points that deviate from FSC.H ~ FSC.A linearly, plus removing debrit (v. small vals in both)

# Function to remove data 5% from the line
remove_outliers <- function(x, y, p = 0.05, slope, intercept) {
  # Get the predicted values
  predicted <- slope * x + intercept

  # Calculate the residuals
  residuals <- y - predicted

  # Calculate the threshold for outliers
  threshold <- quantile(abs(residuals), 1 - p)

  # Identify the outliers
  outliers <- which(abs(residuals) > threshold)

  # Return the new data without outliers
  return(outliers)
}

log_transform_with_shift <- function(data) {
  min_value <- min(data)

  # Calculate the shift based on the minimum value
  shift <- ifelse(min_value < 0, abs(min_value) + 1, 0)

  # Log-transform the adjusted data
  log_transformed_data <- log(data + shift)

  return(log_transformed_data)
}


cell_labelling_bm <- function(seur) {
  library(stringr)
  seur$cell_type <- ifelse(str_detect(seur$assignment, "CD34\\+"), "HSC",
                           ifelse(str_detect(seur$assignment, "CD38\\+"), "LDPs",
                                  ifelse(str_detect(seur$assignment, "CD13\\+") & str_detect(seur$assignment, "CD16\\+"), "Late Granulocytes",
                                         ifelse(str_detect(seur$assignment, "CD13\\+") & str_detect(seur$assignment, "CD36\\+"), "Late Monocytes",
                                                ifelse(str_detect(seur$assignment, "CD13\\+"), "Early Granulocytes",
                                                       ifelse(str_detect(seur$assignment, "CD36\\+") & str_detect(seur$assignment, "CD71\\+"), "Early Erythroid",
                                                              ifelse(str_detect(seur$assignment, "CD41\\+") | str_detect(seur$assignment, "CD42b\\+"), "Megakaryocytes",
                                                                     ifelse(str_detect(seur$assignment, "CD235a\\+"), "Late Erythroid", seur$assignment))))))))

  return(seur)
}

norm_minmax <- function(x){
  tmp = (x - min(x)) / (max(x) - min(x))
  return(tmp)
}

#' @export

run_unsupervised_func <- function(flow_data, res = 0.5, logfold = 0.25, percentage_cells = 0.25, batch_correct = FALSE) {
  library(Seurat)
  library(dplyr)
  if (!("Timepoint" %in% colnames(flow_data@meta.data))) {
    cat("Creating 'Timepoint' column with NULL values\n")
    flow_data@meta.data$Timepoint <- NA
  }

  ref_seurat <- flow_data

  if (batch_correct && length(unique(ref_seurat$Timepoint)) > 1) {
    cat("Running batch correction based on timepoint\n")
    day.list <- SplitObject(ref_seurat, split.by = "Timepoint")
    features <- SelectIntegrationFeatures(object.list = day.list, nfeatures = 3000)
    anchors <- FindIntegrationAnchors(object.list = day.list, anchor.features = features,
                                      normalization.method = "LogNormalize", reduction = "cca",
                                      dims = 1:(nrow(ref_seurat) - 1), verbose = TRUE, k.filter = 7)
    day.integrated <- IntegrateData(anchorset = anchors,
                                    normalization.method = "LogNormalize",
                                    dims = 1:(nrow(ref_seurat) - 1), features = features,
                                    features.to.integrate = features,
                                    k.weight = 5)
    day.integrated <- ScaleData(day.integrated)
    day.integrated <- RunPCA(day.integrated, verbose = TRUE, features = features)
    day.integrated <- RunUMAP(day.integrated, dims = 1:ncol(ref_seurat$pca), n.components = 3L)
    day.integrated <- FindNeighbors(day.integrated, dims = 1:ncol(ref_seurat$pca), reduction = "pca", verbose = TRUE)
    day.integrated <- FindClusters(day.integrated, resolution = res, verbose = TRUE)
    ref_seurat <- day.integrated
  } else if (!batch_correct || (batch_correct && length(unique(ref_seurat$Timepoint)) == 1)) {
    if (batch_correct && length(unique(ref_seurat$Timepoint)) == 1) {
      cat("Cannot run batch correction as only 1 study day found")
    }
    ref_seurat <- NormalizeData(ref_seurat)
    ref_seurat <- FindVariableFeatures(ref_seurat)#, loess.span = 10)
    ref_seurat <- ScaleData(ref_seurat)
    ref_seurat <- RunPCA(ref_seurat, "RNA", features = VariableFeatures(ref_seurat))
    ref_seurat <- RunUMAP(ref_seurat,dims=1:ncol(ref_seurat$pca), n.components = 3L, return.model = TRUE)
    ref_seurat <- FindNeighbors(ref_seurat, dims = 1:ncol(ref_seurat$pca), reduction = "pca", verbose = TRUE)
    ref_seurat <- FindClusters(ref_seurat, resolution = res, verbose = TRUE)
  }
  markers.pos <- FindAllMarkers(ref_seurat, only.pos = TRUE, min.pct = percentage_cells, logfc.threshold = logfold)
  markers.pos_neg <- FindAllMarkers(ref_seurat, only.pos = FALSE, min.pct = percentage_cells, logfc.threshold = logfold)

  all_markers <- merge(markers.pos_neg, markers.pos, by = c("cluster", "gene"), all = TRUE)
  all_table_neg <- table(all_markers[all_markers$avg_log2FC.x < 0, ]$cluster, all_markers[all_markers$avg_log2FC.x < 0, ]$gene)
  all_table_pos <- table(all_markers[all_markers$avg_log2FC.x > 0, ]$cluster, all_markers[all_markers$avg_log2FC.x > 0, ]$gene)

  all_labels <- lapply(1:nrow(all_table_neg), function(i) {
    labels_tmp_neg <- paste0(colnames(all_table_neg)[all_table_neg[i, ] == 1], "-")
    labels_tmp_pos <- paste0(colnames(all_table_pos)[all_table_pos[i, ] == 1], "+")

    c("cluster" = rownames(all_table_pos)[i],
      "assignment" = paste(paste(labels_tmp_pos, collapse = ""), paste(labels_tmp_neg, collapse = ""), collapse = "")
    )
  })

  manual_labelling <- data.frame(do.call("rbind", all_labels))
  manual_labelling_sorted <- manual_labelling[order(as.numeric(manual_labelling$cluster)), ]

  new.cluster.ids <- manual_labelling_sorted$assignment
  names(new.cluster.ids) <- levels(ref_seurat)

  ref_seurat@meta.data <- ref_seurat@meta.data %>%
    left_join(manual_labelling_sorted, by = c("seurat_clusters" = "cluster"))

  rownames(ref_seurat@meta.data) <- colnames(ref_seurat)

  return(ref_seurat)
}

# Function to read either .xlsx or .csv based on the file extension
read_file <- function(file_path) {
  library(readxl)  # For reading Excel files
  library(readr)   # For reading CSV files
  # Detect the file extension
  file_ext <- tools::file_ext(file_path)

  # Read file based on the extension
  if (file_ext == "xlsx") {
    data <- read_excel(file_path)  # Read Excel file
  } else if (file_ext == "csv") {
    data <- read_csv(file_path)    # Read CSV file
  } else {
    stop("Unsupported file format. Please provide a .xlsx or .csv file.")
  }

  return(data)
}

find_common_columns <- function(df1, df2) {
  common_cols <- list()  # Initialize an empty list to store matched columns
  matched_cols <- c()    # Keep track of already matched columns in df2

  for (col1 in colnames(df1)) {
    for (col2 in setdiff(colnames(df2), matched_cols)) {  # Exclude already matched columns
      # Check for common values between columns
      if (length(intersect(unique(df1[[col1]]), unique(df2[[col2]]))) > 0) {
        common_cols[[col1]] <- col2  # Map the match
        matched_cols <- c(matched_cols, col2)  # Mark the column in df2 as matched
        break  # Stop looking for matches for this col1 once a match is found
      }
    }
  }

  return(common_cols)
}

# Check each fcsObject's column names after applying processFCS
processFCS <- function(fcsObject) {
  param_names <- flowCore::pData(flowCore::parameters(fcsObject))[,"name"]  # e.g., "FITC.A"
  param_desc <- flowCore::pData(flowCore::parameters(fcsObject))[,"desc"]   # e.g., "EdU"

  # Use desc where available, fallback to name
  final_colnames <- ifelse(!is.na(param_desc) & param_desc != "", param_desc, param_names)

  # If there are duplicates, append the original param name to make it unique
  dupes <- duplicated(final_colnames) | duplicated(final_colnames, fromLast = TRUE)
  final_colnames[dupes] <- paste0(final_colnames[dupes], "_", param_names[dupes])

  colnames(fcsObject) <- final_colnames
  return(fcsObject)
}


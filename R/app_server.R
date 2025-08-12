#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @import Seurat
#' @import dplyr
#' @import flowCore
#' @export
#' @noRd
app_server <- function(input, output, session) {
  #source("R/helper_functions.R")
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())
  #shinyFiles::shinyDirChoose(input, "directory", roots = volumes, session = session, restrictions = system.file(package = "base"))

  ## button changes after clicking and selection
  observe({
    cat("\ninput$directory value:\n\n")
    print(input$files)
  })
  ## button changes after clicking and selection
  shinyFiles::shinyFileChoose(input, "files", roots = volumes, session = session, restrictions = system.file(package = "base"))
  observe({
    cat("\ninput$files value:\n\n")
    print(input$files)
  })

  output$files <- renderPrint({
    if (is.integer(input$files)) {
      cat("No directory has been selected")
    } else {
      tmp = shinyFiles::parseFilePaths(volumes, input$files)
      cat(paste(nrow(tmp), "files selected"))
    }
  })
  ## button changes after clicking and selection
  shinyFiles::shinyFileChoose(input, "directory_meta", roots = volumes, session = session, restrictions = system.file(package = "base"))
  observe({
    cat("\ninput$directory_meta value:\n\n")
    print(input$directory_meta)
  })

  #model_selected <- reactive({
  #  if(input$model_type == "Unsupervised"){
  #    return("Unsupervised")
  #  }
  #  if(input$model_type == "Supervised"){
  #    return("Supervised")
  #  }
  #})

  output$meta_directory <- renderPrint({
    if (is.integer(input$directory_meta) &
        is.integer(input$directory)) {
      cat("Select metadata file and directory of .fcs files")
    } else {
      fileselected <- shinyFiles::parseFilePaths(volumes, input$directory_meta)
      renderText(as.character(fileselected$datapath))
    }
  })

  files_all <- reactive({
    if (!is.integer(input$files) && !is.integer(input$metadata_file)) {
      # Parse file paths
      files <- shinyFiles::parseFilePaths(volumes, input$files)$datapath
      files <- as.character(files)

      # Read FCS files using flowCore
      all_read <- lapply(files, function(file) {
        tryCatch(
          flowCore::read.FCS(file, alter.names = TRUE, transformation = NULL),
          error = function(e) {
            message("Error reading file: ", file)
            NULL
          }
        )
      })
      all_read <- Filter(Negate(is.null), all_read)  # Remove NULL entries

      # Extract metadata from FCS files
      fn_metadata <<- lapply(seq_along(all_read), function(i) {
        fcs_obj <- all_read[[i]]
        meta <- as.data.frame(fcs_obj@description)[1, , drop = FALSE]

        well <- if ("WELL.ID" %in% colnames(meta) && "FILENAME" %in% colnames(meta)) {
          ifelse(!is.na(meta$WELL.ID), meta$WELL.ID, basename(files[i]))
        } else {
          NA
        }

        # Pair the specific filename with metadata
        paste(basename(files[i]), well, sep = "/")
      })
#      print(paste0("basename_files:", basename(files)))
      if (input$preprocess == "Yes") {
        cat("\nRunning compensation and transformation\n")

        # Preprocess data
        all_ff <- lapply(all_read, preprocess_2)
        all_ff_a <- lapply(all_ff, channel_select)

        # Remove margin events and perform QC
        all_fc_margin <- lapply(all_ff_a, function(fcs_obj) {
          tryCatch(PeacoQC::RemoveMargins(fcs_obj, channels = 1:ncol(fcs_obj)), error = function(e) fcs_obj)
        })
        all_fc_qc <- lapply(all_fc_margin, function(fcs_obj) {
          tryCatch(PeacoQC::PeacoQC(fcs_obj, report = FALSE, channels = c("Time"), save_fcs = FALSE)$FinalFF, error = function(e) fcs_obj)
        })

        # Remove debris using Gaussian Mixture Model clustering
        all_fc_int <- lapply(all_fc_qc, function(fcs_obj) {
          tryCatch({
            mdl <- mclust::Mclust(fcs_obj@exprs[, "FSC.A"], G = 2)
            debris_clus <- which.min(mdl$parameters$mean)
            fcs_obj[mdl$classification != debris_clus, ]
          }, error = function(e) fcs_obj)
        })

        # Remove outliers and check cell viability
        sc_only <- lapply(seq_along(all_fc_int), function(i) {
          fit <- tryCatch(lm(all_fc_int[[i]]@exprs[, "FSC.H"] ~ all_fc_int[[i]]@exprs[, "FSC.A"]), error = function(e) NULL)
          if (!is.null(fit)) {
            slope <- coef(fit)[2]
            intercept <- coef(fit)[1]
            remove_these <- remove_outliers(all_fc_int[[i]]@exprs[, "FSC.H"], all_fc_int[[i]]@exprs[, "FSC.A"], slope=slope, intercept=intercept)
            fcs_obj <- all_fc_int[[i]][-remove_these, ]
            if ("Viability" %in% colnames(fcs_obj@exprs)) {
              fcs_obj <- fcs_obj[fcs_obj@exprs[, "Viability"] < 2, ]
            }
            return(fcs_obj)
          } else {
            return(all_fc_int[[i]])
          }
        })

        # Standardize column names and combine metadata with expression data
        sc_only_rename <<- lapply(sc_only, processFCS)
        seurat_comb <<- lapply(seq_along(sc_only_rename), function(i) {
          expr_data <- sc_only_rename[[i]]@exprs
          if (nrow(expr_data) > 0) {
            metadata_vector <- fn_metadata[[i]]

            # Transpose the metadata to match nrow with expr_data
            metadata_df <- as.data.frame(t(matrix(metadata_vector, nrow = length(metadata_vector), ncol = nrow(expr_data))))
            colnames(metadata_df) <- paste0("filename", seq_len(ncol(metadata_df)))
            cbind(metadata_df, expr_data)
          } else {
            data.frame()
          }
        })
        seurat_comb_dat <- do.call(rbind, seurat_comb)
        #print(colnames(seurat_comb_dat))

        seurat_meta_comb <<- seurat_comb_dat[, grepl("filename", colnames(seurat_comb_dat)),drop=FALSE]
        seurat_dat_comb <<- seurat_comb_dat[, !grepl("filename", colnames(seurat_comb_dat)),drop=FALSE]

        # Merge user-provided metadata
        if (!is.null(input$metadata_file)) {
          print("User uploaded metadata file - merging!")
          meta_file <- read_file(input$metadata_file$datapath)
          common_cols <- find_common_columns(seurat_meta_comb, meta_file)
          merged_data <- merge(seurat_meta_comb, meta_file, by = common_cols, all.x = TRUE)
          seurat_meta_comb <<- merged_data
        }
        #print(paste0("meta_check:", seurat_meta_comb))

        rownames(seurat_dat_comb) <<- rownames(seurat_meta_comb)
        return(seurat_dat_comb)
      } else if (input$preprocess == "No") {
        no_pp_rename <<- lapply(all_read, processFCS) #rename using labels from flow cytometry
        all_exprs <- lapply(no_pp_rename, flowCore::exprs)
        meta_file <- lapply(seq_along(all_exprs), function(i) {
          metadata_vector <- fn_metadata[[i]]
          expr_data <- all_exprs[[i]]

          # Transpose the metadata to match nrow with expr_data
          metadata_df <- as.data.frame(t(matrix(metadata_vector, nrow = length(metadata_vector), ncol = nrow(expr_data))))
          colnames(metadata_df) <- paste0("filename", seq_len(ncol(metadata_df)))
          #cbind(metadata_df, expr_data)
          metadata_df
        })
        seurat_dat_comb <<- as.data.frame(do.call("rbind",
                                                  all_exprs))
        seurat_meta_comb <<- do.call("rbind", meta_file)
        rownames(seurat_dat_comb) <<- rownames(seurat_meta_comb)
        return(seurat_dat_comb)
      }
    }
  })



  # Render column selector dynamically based on uploaded data
  output$columnSelector <- renderUI({
    req(files_all())
    colnames <- colnames(seurat_dat_comb)
    checkboxGroupInput("columns", "Select columns for analysis:", choices = colnames, selected = colnames)
  })

  cluster_dat <- eventReactive(input$runClustering, {
    if(!is.integer(input$files) & input$model_type == "Unsupervised") {
      req(input$columns)
      # select columns
      seurat_dat_comb <- seurat_dat_comb[, (colnames(seurat_dat_comb) %in% input$columns)]
      seurat_dat_comb <- apply(seurat_dat_comb, 2, norm_minmax)
      # dimensionality reduction and clustering
      # convert to seurat first
      # condition to sample to reduce computation
      sample_rows <- NULL  # Initialize sample_rows
      if (nrow(seurat_dat_comb) > 1e5) {
        sample_rows <- sample(nrow(seurat_dat_comb), 1e5)
        seurat_dat_comb <- seurat_dat_comb[sample_rows,]
        seurat_meta_comb <- seurat_meta_comb[sample_rows,]
      }
      seurat_meta_comb
      print(colnames(seurat_meta_comb))
      seurat_obj <- SeuratObject::CreateSeuratObject(counts=t(seurat_dat_comb), meta.data = seurat_meta_comb)
      cluster_dat <- run_unsupervised_func(seurat_obj, res=as.numeric(input$res_umap), logfold=as.numeric(input$lf_umap))
      return(cluster_dat)
    } else if (!is.integer(input$files) & input$model_type == "Supervised" & !is.integer(input$model_file)) {
      # supervised model code here
      model_path <- as.character(shinyFiles::parseFilePaths(volumes, input$input$model_file)$datapath)
      model_single_cell <- readRDS(model_path)

      # check model built with same column headers
      # Extract the column names (features) used in the model
      if("randomForest" %in% class(model_single_cell)){
        features <- model_single_cell$forest$xlevels
        feature_names <- names(features)
      } else if ("svm" %in% model_single_cell){
        feature_names <- colnames(model_single_cell$SV)
      }

      # Check if all features are present in seurat_dat_comb
      missing_features <- setdiff(feature_names, colnames(seurat_dat_comb))

      if (length(missing_features) > 0) {
        stop(paste("Missing features in seurat_dat_comb:", paste(missing_features, collapse = ", ")))
      } else {
        print("All required features are present.")
      }
      # Scale the data
      scale_data <- function(data, means, sds) {
        scaled_data <- sweep(data, 2, means, FUN = "-")
        scaled_data <- sweep(scaled_data, 2, sds, FUN = "/")
        return(scaled_data)
      }
      # Check if scaling_parameters are available
      if (exists("scaling_parameters")) {
        # Select only the relevant columns from seurat_dat_comb
        relevant_data <- seurat_dat_comb[, feature_names, drop = FALSE]

        # Scale the data
        scaled_seurat_data <- scale_data(relevant_data, scaling_parameters$means, scaling_parameters$sds)
        scaled_seurat_data <- apply(scaled_seurat_data, 2, norm_minmax)

        # Make predictions using the model
        if ("randomForest" %in% class(model_single_cell)) {
          predictions <- predict(model_single_cell, newdata = scaled_seurat_data)
        } else if ("svm" %in% class(model_single_cell)) {
          predictions <- predict(model_single_cell, newdata = scaled_seurat_data)
        } else {
          stop("Unsupported model type.")
        }

        # Add the predictions as a new column to the original dataset
        seurat_dat_comb$predicted_class <- predictions

        # Print the updated dataset with predictions
      } else {
        stop("Scaling parameters are not available.")
      }

    }


  })

  # Define a reactive expression to handle out_dat
  out_dat_reactive <- reactive({
    cluster_dat()
  })

  # Observe changes in out_dat_reactive
  observe({
    out_dat <- out_dat_reactive()
    #print(out_dat)
  })
  output$plotdata <- plotly::renderPlotly({
    req(out_dat_reactive)  # Ensure out_dat is not NULL before proceeding
    out_dat <- out_dat_reactive()  # Retrieve the value
    #req(input$cell_assignment)
    #label <- input$cell_assignment
    # Check if the UMAP data is present in the Seurat object
    if ("umap" %in% names(out_dat@reductions)) {
      umap_data <- out_dat[["umap"]]@cell.embeddings
      print(head(umap_data))

    } else {
      cat("UMAP data is not available in the Seurat object.")
    }
    rownames(out_dat@meta.data) <- rownames(out_dat[["umap"]]@cell.embeddings)
    if (!is.null(out_dat)) {

      seurat_to_df <<- Seurat::FetchData(object = out_dat, vars = c("umap_1", "umap_2", "umap_3", colnames(out_dat@meta.data)))

      plotly::plot_ly(data = seurat_to_df,#[plot.data$ID != "",],
                      x = ~umap_1, y = ~umap_2, z = ~umap_3,
                      color = ~assignment,
                      type = "scatter3d",
                      mode = "markers",
                      marker = list(size = 2, width=2), #This is that extra column we made earlier for which we will use for cell ID
                      hoverinfo="text",
                      size = 10, alpha = I(1), text =~assignment, showlegend = FALSE) #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names

    }
  })
  output$tablecounts <- DT::renderDataTable({
    req(out_dat_reactive)  # Ensure out_dat is not NULL before proceeding
    out_dat <- out_dat_reactive()  # Retrieve the value
    if(!is.null(out_dat)){
      print(out_dat@meta.data)
      library(dplyr)
      seurat_metadata <<- out_dat@meta.data

      # Dynamically get all column names with the prefix "filename"
      filename_cols <- grep("^filename", colnames(seurat_metadata), value = TRUE)

      # Concatenate all filename columns to create a unique full file path column
      seurat_metadata <- seurat_metadata %>%
        mutate(Sample = apply(select(., all_of(filename_cols)), 1, paste, collapse = "/"))

      # Summarize counts and merge back with the original metadata
      if ("proliferation" %in% colnames(seurat_metadata)) {
        # Summarize counts
        summary_counts <- seurat_metadata %>%
          group_by(Sample, assignment, proliferation) %>%
          summarise(count = n(), .groups = "drop") %>%
          setNames(c("Sample", "assignment", "proliferation", "count"))

        # Merge with original metadata
        summary_tab <<- seurat_metadata %>%
          left_join(summary_counts, by = c("Sample", "assignment", "proliferation"), keep = FALSE) %>%
          select(-nCount_RNA)
      } else {
        # Summarize counts
        summary_counts <- seurat_metadata %>%
          group_by(Sample, assignment) %>%
          summarise(count = n(), .groups = "drop") %>%
          setNames(c("Sample", "assignment", "count"))

        # Merge with original metadata
        summary_tab <<- seurat_metadata %>%
          left_join(summary_counts, by = c("Sample", "assignment"), keep = FALSE) %>%
          select(-nCount_RNA)
        summary_tab <<- summary_tab %>%
          distinct()
      }

      summary_tab
    }
  })


  # Downloadable zip of processed FCS files
  output$downloadcounts <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "processed_fcs.zip", sep = "_")
    },
    content = function(fname) {
      req(summary_tab)
      write.csv(data.table(summary_tab), fname, row.names = FALSE)
    }
  )
  # Downloadable csv of count data per file
  output$downloadcounts <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "AutoFlow_counts.csv", sep = "_")
    },
    content = function(fname) {
      req(summary_tab)
      write.csv(data.table(summary_tab), fname, row.names = FALSE)
    }
  )
  output$downloadcountsdelta <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "DeltaFlow_counts.csv", sep = "_")
    },
    content = function(fname) {
      req(out_dat_reactive)  # Ensure out_dat is not NULL before proceeding
      req(summary_tab)
      library(tidyr)
      out_dat <- out_dat_reactive()  # Retrieve the value

      if (!is.null(out_dat)) {
        #print(out_dat@meta.data)
        seurat_metadata <- out_dat@meta.data
        summary_tab_prep <- summary_tab %>%
          mutate(assignment = paste0(assignment, " Count")) %>%
          mutate(count = as.numeric(count),
                 assignment = as.character(assignment))
        #print(summary_tab_prep)
        seurat_metadata_wide <- summary_tab_prep %>%
          pivot_wider(names_from = assignment,
                      values_from = count,
                      values_fill = list(count = 0))
        # Write the wide-format data to a CSV file
        write.csv(seurat_metadata_wide, fname, row.names = FALSE)
      }
    }
  )
  # Downloadable seurat object as RDS object
  output$downloadseurat <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "AutoFlow_seurat_obj.rds", sep = "_")
    },
    content = function(fname) {
      req(out_dat_reactive)
      out_dat <- out_dat_reactive()  # Retrieve the value
      saveRDS(out_dat, fname)
    }
  )
  # Dynamically create dropdowns for user-selected columns
  output$select_color_column <- renderUI({
    req(out_dat_reactive)  # Ensure out_dat is not NULL before proceeding
    out_dat <- out_dat_reactive()  # Retrieve the value
    selectInput(
      "color_column",
      "Select Treatment for Coloring:",
      choices = colnames(out_dat@meta.data),
      selected = colnames(out_dat@meta.data)[1]  # Default to the first column
    )
  })

  output$select_x_column <- renderUI({
    req(out_dat_reactive)  # Ensure out_dat is not NULL before proceeding
    out_dat <- out_dat_reactive()  # Retrieve the value
    selectInput(
      "x_column",
      "Select Timepoint Column (X-axis):",
      choices = colnames(out_dat@meta.data),
      selected = colnames(out_dat@meta.data)[2]  # Default to the second column
    )
  })

  output$select_cells <- renderUI({
    req(out_dat_reactive)  # Ensure out_dat is not NULL before proceeding
    out_dat <- out_dat_reactive()  # Retrieve the value
    selectInput(
      "cell_assignment",
      "Select cell assignment for plot",
      choices = unique(out_dat@meta.data$assignment),
      selected = unique(out_dat@meta.data$assignment)[2]  # Default to the second
    )
  })

  observe({
    req(out_dat_reactive())
    updateSelectInput(
      session, "color_column",
      choices = colnames(out_dat_reactive()@meta.data)
    )
    updateSelectInput(
      session, "x_column",
      choices = colnames(out_dat_reactive()@meta.data)
    )
    updateSelectInput(
      session, "cell_assignment",
      choices = unique(out_dat_reactive()@meta.data$assignment)
    )
  })
  observe({
    if (!is.null(out_dat_reactive())) {
      #print("Debug: out_dat_reactive is valid")
      # Only print inputs if out_dat_reactive is ready
      #print(paste("Debug: input$color_column =", input$color_column))
      #print(paste("Debug: input$x_column =", input$x_column))
      #print(paste("Debug: input$cell_assignment =", input$cell_assignment))
    } else {
      #print("Debug: Skipping observe logic because out_dat_reactive is NULL")
    }
  })

  #
  #   # Dynamically render the treatment plot tab based on user input
  #   output$treatment_plot_tab <- renderUI({
  #     req(out_dat_reactive)  # Ensure out_dat is not NULL before proceeding
  #     out_dat <- out_dat_reactive()  # Retrieve the value
  #     if (isTruthy(input$show_treatment_plot) && isTruthy(input$metadata_file)) {
  #       tabPanel("Treatment Plot", plotly::plotlyOutput(outputId = "plottreatment"))
  #     } else {
  #       NULL
  #     }
  #   })
  #
  #   observeEvent(input$show_treatment_plot, {
  #     print(paste0("checking treatment plot logic, isTruthy(input$show_treatment_plot):", isTruthy(input$show_treatment_plot),
  #                  "isTruthy(input$metadata_file)", isTruthy(input$metadata_file)))
  #     if (isTruthy(input$metadata_file) && input$show_treatment_plot) {
  #       showTab("tabs", "Treatment Plot")
  #     } else {
  #       hideTab("tabs", "Treatment Plot")
  #     }
  #   })

  # Observe changes in the show_treatment_plot input and metadata_file
  observeEvent(input$show_treatment_plot, {
    # Check and print the logic for debugging
    print(paste0("checking treatment plot logic, isTruthy(input$show_treatment_plot):", isTruthy(input$show_treatment_plot),
                 " isTruthy(input$metadata_file):", isTruthy(input$metadata_file)))

    # Show or hide the "Treatment Plot" tab based on conditions
    if (isTruthy(input$metadata_file) && isTruthy(input$show_treatment_plot)) {
      showTab(inputId = "tabs", target = "Treatment Plot")
    } else {
      hideTab(inputId = "tabs", target = "Treatment Plot")
    }
  })

  output$plottreatment <- plotly::renderPlotly({
    req(out_dat_reactive())  # Ensure out_dat is available

    print("Debug: Entered output$plottreatment")  # Check if rendering is being triggered

    # Simple Plotly plot for testing
    p <- plot_ly(mtcars, x = ~wt, y = ~mpg, type = "scatter", mode = "markers")

    return(p)
  })



  # # Render the Plotly Plot
  # output$plottreatment <- plotly::renderPlotly({
  #   print("Debug: Entered output$plottreatment")
  #   req(out_dat_reactive)  # Ensure out_dat is not NULL before proceeding
  #   out_dat <- out_dat_reactive()  # Retrieve the value
  #   req(input$color_column, input$x_column, input$cell_assignment)  # Ensure input is ready
  #
  #   # Extract user-selected columns
  #   treatment_column <- input$color_column
  #   x_column <- input$x_column
  #   cell_assignment <- input$cell_assignment
  #
  #   # Example preprocessing of out_dat@meta.data
  #   summary_data <- out_dat@meta.data %>%
  #     # Group by the desired columns
  #     group_by(x_column,  assignment, treatment_column) %>%
  #     # Summarize the counts for each group
  #     summarise(count = n(), .groups = "drop")
  #
  #   # View summary data (for debugging)
  #   print(summary_data)
  #
  #   # Test ggplot to visualize the counts
  #   p <- ggplot(summary_data %>%
  #                 filter(assignment == cell_assignment),
  #               aes(x = as.factor(x_column),
  #                   y = count, color = treatment_column,
  #                   group = treatment_column)) +
  #     geom_line(size = 1) +                # Add line for each Treatment code
  #     geom_point(size = 3) +               # Add points to highlight data points
  #     labs(
  #       x = "Timepoint",
  #       y = "Count",
  #       color = "Treatment Code",
  #       title = "Counts by Timepoint and Treatment (Line Plot)"
  #     ) +
  #     theme_minimal()
  #     plotly::ggplotly(p)
  # })


}

# functions

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

# Check each fcsObject's column names after applying `processFCS`
processFCS <- function(fcsObject) {
  param_names <- pData(parameters(fcsObject))[,"name"]  # e.g., "FITC.A"
  param_desc <- pData(parameters(fcsObject))[,"desc"]   # e.g., "EdU"

  # Use desc where available, fallback to name
  final_colnames <- ifelse(!is.na(param_desc) & param_desc != "", param_desc, param_names)

  # If there are duplicates, append the original param name to make it unique
  dupes <- duplicated(final_colnames) | duplicated(final_colnames, fromLast = TRUE)
  final_colnames[dupes] <- paste0(final_colnames[dupes], "_", param_names[dupes])

  colnames(fcsObject) <- final_colnames
  return(fcsObject)
}


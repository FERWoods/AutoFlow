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
  #volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())
  #shinyFiles::shinyDirChoose(input, "directory", roots = volumes, session = session, restrictions = system.file(package = "base"))

  # Define the directories you want to include
  volumes <- c("Home (~)" = '~', "Current directory" = getwd(), "Root (/)" = '/')

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
  ## button changes after clicking and selection
  shinyFiles::shinyFileChoose(input, "model_supervised", roots = volumes, session = session, restrictions = system.file(package = "base"))
  observe({
    cat("\ninput$model_supervised value:\n\n")
    print(input$model_supervised)
  })


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
    if(!is.integer(input$files)){
      source
      print(getwd())
      marker_panel <<- read.csv("data-raw/panel.csv")
      # list filenames select
      files <- as.character(shinyFiles::parseFilePaths(volumes, input$files)$datapath)

      #
      #files <- list.files(path="~/2024-1/D0 05022024/Sample plate/", full.names=TRUE, pattern=".fcs")

      #files <- files[files %like% "Samples "]

      #files <- list.files("~/Miltenyi Immunoprofiling plates/CD14+ Day0 Miltenyi immunoprofiling plates/Raw/Plate_001/Group_001/", pattern="fcs", full.names = TRUE)

      #read.FCS(files[[1]])
      # Attempt to read in with flowCore
      all_read <- sapply(files, function(x){tryCatch(flowCore::read.FCS(x, alter.names = TRUE, transformation = NULL))})

      # extract file name folders for metadata
      # add wells
      wells <- lapply(all_read, function(x) x@description$"WELL ID")
      files_add <- paste(files, wells, sep="/")
      fn_metadata <- stringr::str_split(files_add, "/")
      #####

      if(input$preprocess == "Yes"){
        cat("\n Running compensation and transformation \n")
        #####
        # # pre-process data
        all_ff <- lapply(all_read, preprocess_2)

        # rename with marker_panel
        all_ff_rename <- all_ff
        for(j in 1:length(all_ff_rename)){
          all_ff_rename[[j]] <- all_ff[[j]]
          for(i in 1:length(marker_panel$Marker)){
            flowCore::colnames(all_ff_rename[[j]])[which(grepl(marker_panel$Fluorochrome[i], flowCore::colnames(all_ff_rename[[j]])))] <- marker_panel$Marker[i]

          }
          # subset to get only cols we care about:
          #fs_to_mat_rename[[j]] <- fs_to_mat_rename[[j]][,colnames(fs_to_mat_rename[[j]]) %in% marker_list$Marker |
          #                                                colnames(fs_to_mat_rename[[j]]) %in% c("FSC.A", "SSC.A", "FITC.A")]
        }

        # remove debris with gauss mix clustering (MCLUST) on V525.50
        # select channels
        all_ff_a <- lapply(all_ff_rename, channel_select)

       # Remove margin events using PeacoQC's RemoveMargins()
        all_fc_margin <- lapply(all_ff_a, function(x) {
          tryCatch({
            PeacoQC::RemoveMargins(x, channels=1:ncol(x))
          }, error = function(e) {
            message("RemoveMargins failed, continuing without")
            x  # Return the original object on error
          })
        })

        # QC in time using PeacoQC selecting only time parameter
        all_fc_qc <- lapply(all_fc_margin, function(x) {
          tryCatch({
            PeacoQC::PeacoQC(x, report=FALSE, channels=c("Time"), save_fcs=FALSE)$FinalFF
          }, error = function(e) {
            message("PeacoQC failed, continuing without")
            x  # Return the original object on error
          })
        })
        library(mclust)
        # debris removed with threshold <FSC.H/A 0.3*10^5
        #all_fc_int <- lapply(all_fc_qc, function(x, threshold = 0.3*10^5){x[x@exprs[,"FSC.A"] > threshold & x@exprs[,"FSC.H"] > threshold,]})
        all_fc_int <- lapply(all_fc_qc, function(x) {
          tryCatch({
            mdl_tmp = mclust::Mclust(x@exprs[,c("FSC.A")], G=2)
            debris_clus = which.min(mdl_tmp$parameters$mean)
            x = x[mdl_tmp$classification != names(debris_clus),]
            x
          }, error = function(e) {
          message("Debris removal failed, continuing without")
          #print(x)

          x
          })
        })

        # Fit a linear regression line with error handling
        fits <- lapply(all_fc_int, function(x) {
          tryCatch({
            lm(x@exprs[, "FSC.H"] ~ x@exprs[, "FSC.A"])
          }, error = function(e) {
            return(NULL)  # Return NULL if lm fails
          })
        })

        # Remove outliers and select for cell viability if 'Viability' column exists
        sc_only <- lapply(seq_along(all_fc_int), function(i) {
          fit <- fits[[i]]
          if (is.null(fit)) {
            # If the lm model is NULL, return the original flowFrame (which could be empty)
            return(all_fc_int[[i]])
          } else {
            # Proceed with outlier removal and viability check
            slope <- coef(fit)[2]
            intercept <- coef(fit)[1]
            remove_these <- remove_outliers(all_fc_int[[i]]@exprs[, "FSC.H"], all_fc_int[[i]]@exprs[, "FSC.A"],
                                            slope = slope, intercept = intercept)
            dat_tmp <- all_fc_int[[i]][-remove_these,]

            # Check if 'Viability' column exists before applying the check
            if ("Viability" %in% colnames(dat_tmp@exprs)) {
              dat_tmp <- dat_tmp[dat_tmp@exprs[, "Viability"] < 2,]
            }

            return(dat_tmp)
          }
        })

        # Assuming fn_metadata is defined elsewhere and is a list with metadata for each flowFrame
        # Create a combined metadata and data frame, ensuring matching rows
        seurat_comb <- lapply(seq_along(sc_only), function(i) {
          # Extract expression data
          expr_data <- if (nrow(sc_only[[i]]@exprs) > 0) sc_only[[i]]@exprs else NULL

          # Create metadata only if expression data is not NULL
          if (!is.null(expr_data)) {
            # Generate metadata dataframe with the correct number of rows
            metadata_matrix <- t(replicate(nrow(expr_data), fn_metadata[[i]]))
            metadata_df <- data.frame(metadata_matrix)
            colnames(metadata_df) <- paste0("filename", seq_len(ncol(metadata_df)))

            # Combine expression data and metadata
            combined_df <- cbind(metadata_df, expr_data)
            return(combined_df)
          } else {
            # Return an empty data frame for empty flowFrames
            return(data.frame())
          }
        })

        # Bind all the combined data frames together
        seurat_comb_dat <- do.call(rbind, seurat_comb)

        # You can now extract the metadata and data into separate data frames if needed
        seurat_meta_comb <- seurat_comb_dat[, grep("filename", colnames(seurat_comb_dat))]
        seurat_dat_comb <- seurat_comb_dat[, !grepl("filename", colnames(seurat_comb_dat))]
        #seurat_meta_comb$proliferation <- edu_binary
        rownames(seurat_dat_comb) <- rownames(seurat_meta_comb,)

        return(seurat_dat_comb)

      } else if(input$preprocess == "No"){
        all_exprs <- lapply(all_read, function(x){flowCore::exprs(x)})
        meta_file <- list()
        for(i in 1:length(all_exprs)){
          meta_file[[i]] <- data.frame(t(replicate(nrow(all_exprs[[i]]),c(files[[i]], fn_metadata[[i]]))))
          colnames(meta_file[[i]]) <- paste0("filename", 1:ncol(meta_file[[i]]))

        }
        seurat_dat_comb <<- as.data.frame(do.call("rbind", all_exprs))
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
    if(!is.integer(input$files) & input$model_type == "Unsupervised"){
    req(input$columns)
    # select columns
    seurat_dat_comb <- seurat_dat_comb[, (colnames(seurat_dat_comb) %in% input$columns)]
    seurat_dat_comb <- apply(seurat_dat_comb, 2, norm_minmax)
    # dimensionality reduction and clustering
    # convert to seurat first
    # condition to sample to reduce computation
    if (nrow(seurat_dat_comb) > 1e5){
      sample_rows = sample(nrow(seurat_dat_comb), 1e5)
      seurat_dat_comb <- seurat_dat_comb[sample_rows,]
      seurat_meta_comb <- seurat_meta_comb[sample_rows,]
    }
    #sample_rows
    seurat_obj <- SeuratObject::CreateSeuratObject(counts=t(seurat_dat_comb), meta.data = seurat_meta_comb)
    #test <- flowCore::read.FCS("~/App_benchmarking_autoflow/Mosmann_rare.fcs")
    #test_dat <- data.frame(test@exprs)
    #test_dat_select <- test_dat[,!(colnames(test_dat) %in% c("FSC.A", "FSC.H", "FSC.W", "SSC.A", "SSC.H",
    #                                                         "SSC.W", "Live_Dead", "Time", "label"))]
    #test_cluster <- run_unsupervised(test_dat_select)

    cluster_dat <- run_unsupervised(seurat_obj, res=as.numeric(input$res_umap), logfold=as.numeric(input$lf_umap))

    #out_dat <<- files_all()
    return(cluster_dat)
  } else if (!is.integer(input$files) & input$model_type == "Supervised" & !is.integer(input$model_supervised)){
    #############################################
    # Apply Cell type model
    #############################################
    # read in model
    celltype_mdl <- readRDS(as.character(shinyFiles::parseFilePaths(volumes, input$model_supervised)$datapath))

    # features for model
    features_mdl <- names(celltype_mdl$forest$xlevels)
    cat(paste0("features for model: ", feature_mdl))

    # predict on data
    seurat_dat_comb <- seurat_dat_comb[, (colnames(seurat_dat_comb) %in% input$columns)]
    seurat_dat_comb <- apply(seurat_dat_comb, 2, norm_minmax)

    cat("\n Running cell type model...")
    class_pred <- tryCatch({
      predict(celltype_mdl, seurat_dat_comb[, colnames(seurat_dat_comb) %in% features_mdl])
    }, error = function(e) {
      # Return a vector of NAs with the same length as the number of rows in seurat_dat_comb#
      cat("supervised model failed - returning NA. Check features align with those in model, or try Unsupervised.")
      rep(NA, nrow(seurat_dat_comb))
    })

    seurat_meta_comb$assignment <- class_pred

    #  Define the custom S4 class - ease of access later
    setClass(
      "ClusterData",
      slots = list(
        data = "data.frame",
        meta.data = "data.frame"
      )
    )

    # Step 2: Create an instance of the class
    cluster_dat <- new("ClusterData",
                       data = seurat_dat_comb,
                       meta.data = seurat_meta_comb)

    return(cluster_dat)
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
    if (!is.null(out_dat) & input$model_type == "Unsupervised") {
      # all_files <- list.files("~/R_Code/2019-06/2019-6 Data/D46 04-11-19/", recursive=TRUE, pattern="fcs", full.names = TRUE)
      #all_read <- sapply(all_files, function(x){tryCatch(read.flowSet(x, alter.names=FALSE, transform=NULL), error=function(e){})})
      #autoplot(all_read[[1]], "AmCyan.A", "SSC.A")
      seurat_to_df <<- Seurat::FetchData(object = out_dat, vars = c("umap_1", "umap_2", "umap_3", colnames(out_dat@meta.data)))

      cat(colnames(seurat_to_df))
      #tryCatch({p <- plotly::plot_ly(data = seurat_to_df,#[plot.data$ID != "",],
      #                             x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
      #                            color = ~assignment,
      #                           type = "scatter3d",
      #                          mode = "markers",
      #                         marker = list(size = 1, width=2), #This is that extra column we made earlier for which we will use for cell ID
      #                        hoverinfo="text",
      #                       size = 10, alpha = I(1), text =~assignment) #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names
      # p},
      #error=function(e){cat("plot_error")})
      plotly::plot_ly(data = seurat_to_df,#[plot.data$ID != "",],
                      x = ~umap_1, y = ~umap_2, z = ~umap_3,
                      color = ~assignment,
                      type = "scatter3d",
                      mode = "markers",
                      marker = list(size = 1, width=2), #This is that extra column we made earlier for which we will use for cell ID
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
      summary_tab <<- if("proliferation" %in% seurat_metadata){
        seurat_metadata %>%
        group_by(eval(parse(text=paste0("filename", 1))), assignment, proliferation) %>%
        summarise("count" = n()) %>%
        setNames(c("filename1", "assignment", "proliferation", "count"))
      } else {
        seurat_metadata %>%
          group_by(eval(parse(text=paste0("filename", 1))), assignment) %>%
          summarise("count" = n()) %>%
          setNames(c("filename1", "assignment", "count"))
      }
        summary_tab
    }
  })


  observe({
    if(exists("out_dat")){
      tryCatch({
        #updateSelectInput(session, "var1", choices =  unique(csv_data()$Type), selected = mytype)
        choices <- unique(out_dat@meta.data$assignment)
        updateSelectInput(session, inputId = "cell_assignment", choices = choices)
      },
      error=function(e){cat("")}
      )
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
}

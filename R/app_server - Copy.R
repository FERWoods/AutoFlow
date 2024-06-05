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
    if(!is.integer(input$files)){
      source
      print(getwd())
      marker_panel <<- read.csv("data-raw/panel.csv")
      # list filenames select
      files <- as.character(shinyFiles::parseFilePaths(volumes, input$files)$datapath)

      #
      #files <- list.files(path="~/2023-1 BM MPS_Menin/5. Results/All Raw FCS Files/D0/Sample plate/", full.names=TRUE, pattern=".fcs")
      #files <- files[85]

      # Attempt to read in with flowCore
      all_read <- sapply(files, function(x){tryCatch(flowCore::read.FCS(x, alter.names = TRUE, transformation = NULL))})

      # extract file name folders for metadata
      fn_metadata <- stringr::str_split(files, "/")
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

        # remove margin events using PeacoQC's RemoveMargins()
        all_fc_margin <- lapply(all_ff_a, function(x){PeacoQC::RemoveMargins(x, channels=1:ncol(x))})

        # QC in time using PeacoQC selecting only time parameter
        all_fc_qc <- lapply(all_fc_margin, function(x){tryCatch(PeacoQC::PeacoQC(x, report=FALSE,channels=c("Time"), save_fcs = FALSE)$FinalFF)})
        library(mclust)
        # debris removed with threshold <FSC.H/A 0.3*10^5
        all_fc_int <- lapply(all_fc_qc, function(x, threshold = 0.3*10^5){x[x@exprs[,"FSC.A"] > threshold & x@exprs[,"FSC.H"] > threshold,]})

        # single cell identification with linear regression model
        # Fit a linear regression line
        fits <- lapply(all_fc_int, function(x){
          lm(x@exprs[, "FSC.H"] ~ x@exprs[, "FSC.A"])})

        # Get the slope and intercept of the line
        slopes <- lapply(fits, function(x){coef(x)[2]})
        intercepts <- lapply(fits, function(x){coef(x)[1]})

        #  remove_outliers(x$"FSC-H", x$"FSC-A")})

        sc_only <- list()
        for (i in 1:length(all_fc_int)){
          remove_these <- remove_outliers(all_fc_int[[i]]@exprs[,"FSC.H"], all_fc_int[[i]]@exprs[, "FSC.A"],
                                          slope = slopes[[i]],
                                          intercept = intercepts[[i]])
          dat_tmp <- all_fc_int[[i]]
          sc_only[[i]] <- dat_tmp[-remove_these,]
          sc_only[[i]] <- sc_only[[i]][sc_only[[i]]@exprs[,"Viability"] < 2,] # cell viability
          seurat_dat <- lapply(sc_only, flowCore::exprs)

          meta_file <- list()
          for(i in 1:length(seurat_dat)){
            meta_file[[i]] <- data.frame(t(replicate(nrow(seurat_dat[[i]]),fn_metadata[[i]])))
            colnames(meta_file[[i]]) <- paste0("filename", 1:ncol(meta_file[[i]]))
          }

          seurat_dat_comb <- as.data.frame(do.call("rbind", seurat_dat))
          #edu_mdl <- mclust::Mclust(seurat_dat_comb[,c("EdU")], 2)
          #edu_binary <- ifelse(edu_mdl$classification == which.min(edu_mdl$parameters$mean), 0, 1)

          seurat_meta_comb <<- do.call("rbind", meta_file)
          #seurat_meta_comb$proliferation <- edu_binary
          rownames(seurat_dat_comb) <- rownames(seurat_meta_comb)

        }
      } else if(input$preprocess == "No"){
        all_exprs <- lapply(all_read, function(x){flowCore::exprs(x)})
        meta_file <- list()
        for(i in 1:length(all_exprs)){
          meta_file[[i]] <- data.frame(t(replicate(nrow(all_exprs[[i]]),fn_metadata[[i]])))
          colnames(meta_file[[i]]) <- paste0("filename", 1:ncol(meta_file[[i]]))

        }
        seurat_dat_comb <<- as.data.frame(do.call("rbind", all_exprs))
        seurat_meta_comb <<- do.call("rbind", meta_file)
        rownames(seurat_dat_comb) <<- rownames(seurat_meta_comb)

      }


      # Render column selector dynamically based on uploaded data
      output$columnSelector <- renderUI({
        if (!is.null(seurat_dat_comb)) {
          req(seurat_dat_comb)
          colnames <- colnames(seurat_dat_comb)
          checkboxGroupInput("columns", "Select columns for analysis:", choices = colnames, selected = colnames)
          print(colnames)
        }
      })

      if(!is.integer(input$files) & input$model_type == "Unsupervised"){

      # select columns
      seurat_dat_comb <- seurat_dat_comb[, !(colnames(seurat_dat_comb) %in% input$colnames)]
      # dimensionality reduction and clustering
      # convert to seurat first

      seurat_obj <- SeuratObject::CreateSeuratObject(counts=t(seurat_dat_comb), meta.data = seurat_meta_comb)
      cluster_dat <- run_unsupervised(seurat_obj, res=as.numeric(input$res_umap), logfold=as.numeric(input$lf_umap))

      #out_dat <<- files_all()
      return(cluster_dat)
    } else if (!is.integer(input$files) & input$model_type == "Supervised"){

    }
    }
    })

  # Define a reactive expression to handle out_dat
  out_dat_reactive <- reactive({
    files_all()
  })

  # Observe changes in out_dat_reactive
  observe({
    out_dat <- out_dat_reactive()
    print(out_dat)
  })
  output$plotdata <- plotly::renderPlotly({
    req(out_dat_reactive)  # Ensure out_dat is not NULL before proceeding
    out_dat <- out_dat_reactive()  # Retrieve the value
    if (!is.null(out_dat)) {
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
                      size = 10, alpha = I(1), text =~assignment) #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names

    }
  })
  output$tablecounts <- DT::renderDataTable({
    req(out_dat_reactive)  # Ensure out_dat is not NULL before proceeding
    out_dat <- out_dat_reactive()  # Retrieve the value
    if(!is.null(out_dat)){
      print(out_dat@meta.data)
      library(dplyr)
      seurat_metadata <<- out_dat@meta.data
      summary_tab <<- seurat_metadata %>%
        group_by(eval(parse(text=paste0("filename", ncol(seurat_meta_comb)))), assignment, proliferation) %>%
        summarise("count" = n()) %>%
        setNames(c("filename", "assignment", "proliferation", "count"))
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

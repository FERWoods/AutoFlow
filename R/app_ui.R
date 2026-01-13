#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    # golem_add_external_resources(),

    fluidPage(
      theme = shinythemes::shinytheme("paper"),
      tags$style("
        body {
          -moz-transform: scale(0.7, 0.7);
          zoom: 0.7;
          zoom: 70%;
        }
      "),

      titlePanel("AutoFlow"),

      sidebarPanel(
        tags$p("Select .FCS files for processing"),
        shinyFiles::shinyFilesButton(
          "files", "Select FCS Files", "Upload",
          filter = "FCS files (*.fcs)", multiple = TRUE
        ),
        verbatimTextOutput("files", placeholder = TRUE),
        tags$hr(),

        # Metadata upload
        #fileInput("metadata_file", "Upload Metadata", accept = c(".csv", ".xlsx")),
        #tags$hr(),

        # Preprocess toggle (optional)
        radioButtons("preprocess", "Pre-process files",
                     choices = c("No", "Yes"), selected = "No"),
        # Viability + marker QC controls (server renders; only relevant when preprocess == Yes)
        conditionalPanel(
          condition = "input.preprocess == 'Yes'",
          uiOutput("viability_controls")
        ),
        tags$hr(),

        # Model type
        radioButtons("model_type", "Model Type",
                     choices = c("Unsupervised", "Supervised"),
                     selected = "Unsupervised"),

        # ---------------- Unsupervised controls (independent of preprocessing) ----------------
        conditionalPanel(
          condition = "input.model_type == 'Unsupervised'",
          tags$p("Select UMAP clustering parameters (unsupervised)"),
          selectInput("res_umap", "Resolution (higher number = more clusters)",
                      choices = c(0.001, 0.01, 0.1, 1, 10, 20, 50), selected = 1),
          selectInput("lf_umap", "Log-fold change threshold",
                      choices = c(0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5), selected = 0.2),

          # Column selection for clustering (always available in Unsupervised)
          uiOutput("columnSelector"),

          # Run button (always available in Unsupervised)
          actionButton("runClustering", "Run Clustering (unsupervised)"),
          tags$hr(),

          # Extra controls for density-by-assignment (appear after clustering is available)
          uiOutput("unsup_qc_controls")
        ),

        # Supervised controls
        conditionalPanel(
          condition = "input.model_type == 'Supervised'",
          tags$p("Upload pre-trained model bundle (.rds) with: features, scaling$means/sds, model, (optional) levels"),
          fileInput("model_file", "Upload Pre-trained Model", accept = c(".RDS", ".rds")),
          uiOutput("supervised_controls"),
          uiOutput("feature_mapper_ui"),
          actionButton("runSupervised", "Apply Supervised Model"),
          tags$hr()
        ),

        # Optional treatment plot
        #checkboxInput("show_treatment_plot", "Show Treatment Plot", value = FALSE),
        #uiOutput("select_color_column"),
        #uiOutput("select_x_column"),
        #uiOutput("select_cells"),
        #tags$hr(),

        # Downloads
        downloadButton("downloadprocessed", "Download pre-processed .fcs"),
        downloadButton("downloadcounts", "Download cell counts (Long format)"),
        downloadButton("downloadcountsdelta", "Download cell counts (Wide format)"),
        downloadButton("downloadseurat", "Download Seurat Object")
      ),

      mainPanel(
        tabsetPanel(
          id = "tabs",
          type = "tabs",
          selected = "Welcome",

          tabPanel(
            "Welcome",
            includeMarkdown(
              system.file("README.md", package = "autoflow") # adjust pkg name if needed
            )
          ),

          # Marker QC tab appears whenever preprocessing is enabled
          tabPanel(
            "Marker QC",
            conditionalPanel(
              condition = "input.preprocess == 'Yes'",
              plotly::plotlyOutput("marker_qc", height = "400px")
            )
          ),

          #  Unsupervised marker density tab (works after clustering; preprocess optional)
          tabPanel(
            "Unsupervised Marker Density",
            plotly::plotlyOutput("unsup_marker_qc", height = "400px")
          ),

          tabPanel("Cell count table", DT::dataTableOutput("tablecounts")),
          tabPanel("UMAP plot data", plotly::plotlyOutput("plotdata", height = "600px", width = "800px"))#,

          #tabPanel(
          # "Treatment Plot",
          #  conditionalPanel(
          #    condition = "input.show_treatment_plot && input.metadata_file != null",
          #    plotly::plotlyOutput("plottreatment", height = "400px")
          #)
          #)
        )
      )
    )
  )
}

#' Add external Resources to the Application
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "AutoflowV2"
    )
  )
}

#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    #golem_add_external_resources(),

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
        shinyFiles::shinyFilesButton("files", "Select FCS Files", "Upload", filter = "FCS files (*.fcs)", multiple = TRUE),
        verbatimTextOutput("files", placeholder = TRUE),
        tags$p(),

        # Metadata upload button, always visible
        fileInput("metadata_file", "Upload Metadata", accept = c(".csv", ".xlsx")),

        radioButtons("preprocess", "Pre-process files", choices = c("No", "Yes")),
        tags$p(),

        radioButtons("model_type", "Model Type", choices = c("Unsupervised", "Supervised")),
        uiOutput("columnSelector"),

        conditionalPanel(
          condition = "input.model_type == 'Unsupervised'",
          tags$p("Select UMAP clustering parameters (unsupervised)"),
          selectInput("res_umap", "Resolution (higher number = more clusters)", choices = c(0.001, 0.01, 0.1, 1, 10, 20, 50), selected = 1),
          selectInput("lf_umap", "Log-fold change threshold", choices = c(0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5), selected = 0.2),
          actionButton("runClustering", "Run Clustering (unsupervised)"),
        ),
        conditionalPanel(
          condition = "input.model_type == 'Supervised'",
          tags$p("Upload pre-trained model: ensure column names match those after pre-processing (see doc.)"),
          fileInput("model_file", "Upload Pre-trained Model", accept = c(".RDS", ".rds"))
        ),

        # Show Treatment Plot checkbox, independent of metadata
        checkboxInput("show_treatment_plot", "Show Treatment Plot", value = FALSE),
        uiOutput("select_color_column"),  # Dropdown for selecting the coloring variable
        uiOutput("select_x_column"),
        uiOutput("select_cells"), # Dropdown for selecting the timepoint variable
        tags$p(),
        downloadButton("downloadprocessed", "Download pre-processed .fcs"),
        downloadButton("downloadcounts", "Download cell counts"),
        downloadButton("downloadcountsdelta", "Download cell counts (DeltaFlow format)"),
        downloadButton("downloadseurat", "Download Seurat Object")

      ),

      mainPanel(
        tabsetPanel(
          id = "tabs",
          type = "tabs",
          tabPanel("Cell count table", DT::dataTableOutput("tablecounts")),
          tabPanel("UMAP plot data", plotly::plotlyOutput("plotdata", height = "600px", width = "800px")),

          # Treatment plot only shows if checkbox is checked and a file is uploaded
          tabPanel("Treatment Plot",
                   conditionalPanel(
                     condition = "input.show_treatment_plot && input.metadata_file != null"
                     ,
                     plotly::plotlyOutput("plottreatment")
                   )
          )
        )
      )
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
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
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}

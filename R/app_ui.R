#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic
    fluidPage(theme = shinythemes::shinytheme("paper"),
              tags$style("
              body {
    -moz-transform: scale(0.7, 0.7); /* Moz-browsers */
    zoom: 0.7; /* Other non-webkit browsers */
    zoom: 70%; /* Webkit browsers */
}
              "),

              titlePanel(
                "AutoFlow",
              ),
              sidebarPanel(
                tags$p("Select .FCS files for processing"),
                shinyFiles::shinyFilesButton("files", "Select FCS Files", "Upload", filter = "FCS files (*.fcs)", multiple = TRUE),
                verbatimTextOutput("files", placeholder = TRUE),
                tags$p(),
                radioButtons(
                  inputId = "preprocess",
                  label = "Pre-process files",
                  choices = c("No","Yes")
                ),
                tags$p(),
                radioButtons(
                  inputId = "model_type",
                  label = "Model Type",
                  choices = c("Unsupervised","Supervised")
                ),
                uiOutput("columnSelector"),  # Dynamic UI for selecting columns


                downloadButton("downloadprocessed", "Download pre-processed .fcs"),
                downloadButton("downloadcounts", "Download cell counts"),
                downloadButton("downloadseurat", "Download Seurat Object"),
                conditionalPanel(
                  condition = "input.model_type == 'Unsupervised'",
                  tags$p("Select UMAP clustering parameters (unsupervised)"),
                  selectInput("res_umap", "Resolution (higher number = more clusters)", choices = c(0.001, 0.01, 0.1, 1, 10, 20, 50), selected=1),
                  selectInput("lf_umap", "Log-fold change threshold", choices = c(0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5), selected=0.2)
                )

              ),

              mainPanel(
                tabsetPanel(type = "tabs",
                            tabPanel("Cell count table", DT::dataTableOutput(outputId = "tablecounts")),
                            tabPanel("UMAP plot data", plotly::plotlyOutput(outputId = "plotdata")),
                            #tabPanel("Treatment plot", plotly::plotlyOutput(outputId = "treatmentplot"))
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
      app_title = "AutoFlowApp"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}

#' markers UI Function
#'
#' @description A shiny Module for marker gene analysis.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_markers_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      box(
        title = "Marker Gene Analysis Settings",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        fluidRow(
          column(3,
                 selectInput(ns("marker_test"), "Statistical test:",
                             choices = c("Wilcoxon Rank Sum" = "wilcox",
                                         "Student's t-test" = "t",
                                         "DESeq2" = "DESeq2",
                                         "MAST" = "MAST"),
                             selected = "wilcox")
          ),
          column(3,
                 radioButtons(ns("marker_comparison"), "Comparison type:",
                              choices = c("One vs All" = "one_vs_all",
                                          "Pairwise" = "pairwise"),
                              selected = "one_vs_all")
          ),
          column(3,
                 numericInput(ns("marker_logfc"), "Min log fold change:",
                              value = app_config$default_marker_logfc, 
                              min = 0, max = 2, step = 0.1)
          ),
          column(3,
                 numericInput(ns("marker_minpct"), "Min percent expressed:",
                              value = app_config$default_marker_minpct, 
                              min = 0, max = 1, step = 0.05)
          )
        ),
        
        fluidRow(
          column(3,
                 checkboxInput(ns("marker_onlypos"), 
                               "Only positive markers", 
                               value = TRUE)
          ),
          column(3,
                 numericInput(ns("marker_topn"), "Top N genes per cluster:",
                              value = app_config$default_marker_topn, 
                              min = 1, max = 50)
          ),
          column(3,
                 numericInput(ns("marker_display_genes"), "Max genes to display:",
                              value = app_config$default_marker_display_genes, 
                              min = 10, max = 100, step = 5)
          ),
          column(3,
                 conditionalPanel(
                   condition = paste0("input['", ns("marker_comparison"), "'] == 'pairwise'"),
                   selectInput(ns("marker_cluster1"), "Cluster 1:", choices = NULL),
                   selectInput(ns("marker_cluster2"), "Cluster 2:", choices = NULL)
                 )
          )
        ),
        
        actionButton(ns("calculate_markers"), "Find Marker Genes", 
                     class = "btn-success"),
        
        hr(),
        
        conditionalPanel(
          condition = paste0("output['", ns("markersCalculated"), "']"),
          h4("Quick Actions:"),
          fluidRow(
            column(4,
                   downloadButton(ns("downloadMarkerTable"), 
                                  "Download All Markers (Excel)", 
                                  class = "btn-info")
            ),
            column(4,
                   downloadButton(ns("downloadTopMarkers"), 
                                  "Download Top Markers (CSV)", 
                                  class = "btn-info")
            ),
            column(4,
                   actionButton(ns("clearMarkers"), "Clear Results", 
                                class = "btn-warning")
            )
          )
        )
      )
    ),
    
    fluidRow(
      box(
        title = "Marker Gene Table",
        status = "info",
        solidHeader = TRUE,
        width = 12,
        conditionalPanel(
          condition = paste0("output['", ns("markersCalculated"), "']"),
          DTOutput(ns("markerTable"))
        ),
        conditionalPanel(
          condition = paste0("!output['", ns("markersCalculated"), "']"),
          h4("No markers calculated yet. Click 'Find Marker Genes' to start analysis.",
             style = "color: #999; text-align: center; padding: 50px;")
        )
      )
    ),
    
    box(
      title = "Top Markers Heatmap",
      status = "success",
      solidHeader = TRUE,
      width = 12,
      collapsible = TRUE,
      uiOutput(ns("dynamic_heatmap_ui"))
    ),
    
    box(
      title = "Marker Dot Plot",
      status = "success",
      solidHeader = TRUE,
      width = 12,
      collapsible = TRUE,
      uiOutput(ns("dynamic_dotplot_ui"))
    ),
    
    fluidRow(
      box(
        title = "Selected Marker Visualization",
        status = "warning",
        solidHeader = TRUE,
        width = 12,
        conditionalPanel(
          condition = paste0("output['", ns("markersCalculated"), "']"),
          fluidRow(
            column(4,
                   selectInput(ns("visualize_marker"), "Select marker gene:",
                               choices = NULL)
            ),
            column(4,
                   radioButtons(ns("marker_plot_type"), "Plot type:",
                                choices = c("Feature Plot" = "feature",
                                            "Violin Plot" = "violin",
                                            "Ridge Plot" = "ridge"),
                                selected = "feature",
                                inline = TRUE)
            ),
            column(4,
                   br(),
                   actionButton(ns("update_marker_plot"), "Update Plot", 
                                class = "btn-primary")
            )
          ),
          
          plotOutput(ns("selectedMarkerPlot"), height = "500px")
        )
      )
    )
  )
}
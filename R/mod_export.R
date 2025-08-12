#' export UI Function
#'
#' @description A shiny Module for data and plot export.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_export_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      box(
        title = "Plot Export Settings",
        status = "primary",
        solidHeader = TRUE,
        width = 6,
        h4("Export Configuration"),
        selectInput(ns("export_plot_type"), "Select plot to export:",
                    choices = c(
                      "Current Cluster Plot" = "cluster",
                      "Cell Cycle Plot" = "cellcycle",
                      "Cell Cycle Distribution" = "cellcycle_distribution",
                      "CancerSEA Feature Plot" = "cancersea_feature",
                      "CancerSEA Heatmap" = "cancersea_heatmap",
                      "Pathway Correlation Matrix" = "correlation_matrix",
                      "Pathway by Cluster" = "pathway_by_cluster",
                      "Pathway Comparison Heatmap" = "comparison_heatmap",
                      "Cluster Tree" = "tree",
                      "Custom Gene Expression" = "gene",
                      "Marker Heatmap" = "marker_heatmap",
                      "Main Trajectory" = "trajectory_main",
                      "Stratified Trajectory" = "trajectory_stratified",  # ADD THIS LINE
                      "Pseudotime" = "pseudotime"
                    )),
        
        conditionalPanel(
          condition = paste0("input['", ns("export_plot_type"), "'] == 'gene'"),
          selectInput(ns("export_gene"), "Select gene:", choices = c("none"))
        ),
        
        conditionalPanel(
          condition = paste0("input['", ns("export_plot_type"), "'] == 'marker_heatmap'"),
          p("Note: Marker genes must be calculated in the 'Marker Genes' tab before exporting these plots.",
            style = "color: #666; font-style: italic;")
        ),
        
        conditionalPanel(
          condition = paste0("(input['", ns("export_plot_type"), "'] == 'trajectory_main' || ",
                             "input['", ns("export_plot_type"), "'] == 'trajectory_stratified' || ",
                             "input['", ns("export_plot_type"), "'] == 'pseudotime')"),
          p("Note: Trajectory analysis must be completed in the 'Trajectory Analysis' tab before exporting these plots.",
            style = "color: #666; font-style: italic;")
        ),
        
        conditionalPanel(
          condition = paste0("(input['", ns("export_plot_type"), "'] == 'cancersea_heatmap' || ",
                             "input['", ns("export_plot_type"), "'] == 'correlation_matrix' || ",
                             "input['", ns("export_plot_type"), "'] == 'pathway_by_cluster' || ",
                             "input['", ns("export_plot_type"), "'] == 'comparison_heatmap')"),
          p("Note: At least 2 CancerSEA pathways must be calculated before exporting comparison plots.",
            style = "color: #666; font-style: italic;")
        ),
        
        fluidRow(
          column(6,
                 selectInput(ns("export_format"), "Format:",
                             choices = c("PNG" = "png", 
                                         "PDF" = "pdf", 
                                         "SVG" = "svg",
                                         "TIFF" = "tiff"))
          ),
          column(6,
                 selectInput(ns("export_dpi"), "Resolution (DPI):",
                             choices = c("Screen (72)" = 72,
                                         "Print (300)" = 300,
                                         "Publication (600)" = 600))
          )
        ),
        
        fluidRow(
          column(6,
                 numericInput(ns("export_width"), "Width (inches):", 
                              value = 8, 
                              min = 2, max = 20)
          ),
          column(6,
                 numericInput(ns("export_height"), "Height (inches):", 
                              value = 6, 
                              min = 2, max = 20)
          )
        ),
        
        checkboxInput(ns("export_white_bg"), "White background", value = TRUE),
        checkboxInput(ns("export_no_axes"), "Remove axes labels", value = FALSE),
        
        textInput(ns("export_title"), "Custom title (optional):", value = ""),
        
        downloadButton(ns("downloadPlot"), "Download Plot", class = "btn-success"),
        
        hr(),
        h4("Batch Export"),
        p("Export all plots at once with current settings:"),
        downloadButton(ns("downloadAllPlots"), "Download All Plots (ZIP)", class = "btn-warning")
      ),
      
      box(
        title = "Data Export",
        status = "info",
        solidHeader = TRUE,
        width = 6,
        h4("Export Data Tables"),
        
        checkboxGroupInput(ns("export_data_options"), "Select data to export:",
                           choices = c(
                             "Cell metadata" = "metadata",
                             "Cluster statistics" = "cluster_stats",
                             "CancerSEA scores" = "cancersea_scores",
                             "Average pathway scores per cluster" = "pathway_avg",
                             "Cell cycle proportions" = "cellcycle_props",
                             "Cluster marker genes" = "markers",
                             "Trajectory results" = "trajectory"
                           ),
                           selected = c("metadata", "cluster_stats")
        ),
        
        downloadButton(ns("downloadData"), "Download Excel File", class = "btn-success"),
        
        hr(),
        h4("Export Seurat Object"),
        p("Save the processed Seurat object for future use:"),
        downloadButton(ns("downloadSeurat"), "Download Seurat Object (.rds)", class = "btn-info")
      )
    ),
    
    fluidRow(
      box(
        title = "Citation & Methods",
        status = "success",
        solidHeader = TRUE,
        width = 12,
        h4("Automated Methods Text"),
        p("Copy this text for your methods section:"),
        verbatimTextOutput(ns("methodsText")),
        actionButton(ns("copyMethods"), "Copy to Clipboard", class = "btn-primary"),
        
        hr(),
        h4("Software Citations"),
        p("Please cite the following tools used in this analysis:"),
        verbatimTextOutput(ns("citationText")),
        
        hr(),
        h4("Figure Legend Template"),
        p("Template for your figure legend:"),
        verbatimTextOutput(ns("legendText")),
        
        hr(),
        h4("Analysis Summary"),
        verbatimTextOutput(ns("analysisSummary"))
      )
    )
  )
}
#' cluster_analysis UI Function
#'
#' @description A shiny Module for cluster analysis.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
#' @import plotly
#' @import DT
mod_cluster_analysis_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      box(
        title = "Cluster Visualization",
        status = "primary",
        solidHeader = TRUE,
        width = 8,
        selectInput(ns("reduction_cluster"), "Select reduction:",
                    choices = c("umap", "tsne", "pca")),
        plotlyOutput(ns("clusterPlot"), height = "600px")
      ),
      
      box(
        title = "Cluster Tree",
        status = "primary",
        solidHeader = TRUE,
        width = 4,
        plotOutput(ns("clusterTree"), height = "600px")
      )
    ),
    
    fluidRow(
      box(
        title = "Cluster Statistics",
        status = "info",
        solidHeader = TRUE,
        width = 12,
        DTOutput(ns("clusterStats"))
      )
    )
  )
}

#' cluster_analysis Server Functions
#'
#' @import Seurat
#' @import ggplot2
#' @noRd 
mod_cluster_analysis_server <- function(id, seurat_obj, processed){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
    # Update reduction choices when object is loaded
    observe({
      req(seurat_obj())
      
      available_reductions <- names(seurat_obj()@reductions)
      if (length(available_reductions) > 0) {
        updateSelectInput(session, "reduction_cluster", 
                          choices = available_reductions)
      }
    })
    
    # Output: Cluster plot
    output$clusterPlot <- renderPlotly({
      req(processed())
      req(seurat_obj())
      
      p <- DimPlot(seurat_obj(), 
                   reduction = input$reduction_cluster,
                   label = TRUE,
                   label.size = 4) + 
        theme_minimal() + 
        ggtitle("Cluster Distribution")
      
      ggplotly(p)
    })
    
    # Output: Cluster tree
    output$clusterTree <- renderPlot({
      req(processed())
      req(seurat_obj())
      
      # Ensure we have clusters
      if (length(unique(Idents(seurat_obj()))) > 1) {
        tryCatch({
          PlotClusterTree(seurat_obj())
        }, error = function(e) {
          plot.new()
          text(0.5, 0.5, "Unable to create cluster tree", cex = 1.2)
          if (grepl("plot.new", e$message)) {
            text(0.5, 0.4, "Tree may not be built yet", cex = 1)
          }
        })
      } else {
        plot.new()
        text(0.5, 0.5, "Need at least 2 clusters for tree", cex = 1.2)
      }
    })
    
    # Output: Cluster statistics
    output$clusterStats <- renderDT({
      req(processed())
      req(seurat_obj())
      
      cluster_counts <- table(Idents(seurat_obj()))
      cluster_props <- prop.table(cluster_counts) * 100
      
      stats_df <- data.frame(
        Cluster = names(cluster_counts),
        Cell_Count = as.numeric(cluster_counts),
        Percentage = round(as.numeric(cluster_props), 2)
      )
      
      datatable(stats_df, options = list(pageLength = 15))
    })
    
  })
}
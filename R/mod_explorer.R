#' explorer UI Function
#'
#' @description A shiny Module for interactive cell exploration.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_explorer_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      box(
        title = "Interactive Cell Explorer",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        fluidRow(
          column(4,
                 selectInput(ns("explorer_color"), "Color by:",
                             choices = c("seurat_clusters", "Phase"))
          ),
          column(4,
                 selectInput(ns("explorer_reduction"), "Reduction:",
                             choices = c("umap", "tsne"))
          ),
          column(4,
                 selectInput(ns("explorer_gene"), "Gene expression:",
                             choices = c("none"))
          )
        ),
        plotlyOutput(ns("explorerPlot"), height = "600px"),
        hr(),
        h4("Selected Cell Information:"),
        verbatimTextOutput(ns("cellInfo"))
      )
    )
  )
}

#' explorer Server Functions
#'
#' @import Seurat
#' @import plotly
#' @noRd 
mod_explorer_server <- function(id, seurat_obj, processed){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
    # Update choices when Seurat object changes
    observe({
      req(seurat_obj())
      
      # Update metadata choices for coloring
      meta_cols <- colnames(seurat_obj()@meta.data)
      color_choices <- c("seurat_clusters", meta_cols)
      updateSelectInput(session, "explorer_color", 
                        choices = unique(color_choices))
      
      # Update reduction choices
      available_reductions <- names(seurat_obj()@reductions)
      if (length(available_reductions) > 0) {
        updateSelectInput(session, "explorer_reduction", 
                          choices = available_reductions)
      }
      
      # Update gene choices
      gene_choices <- c("none", rownames(seurat_obj()))
      updateSelectizeInput(session, "explorer_gene", 
                           choices = gene_choices,
                           server = TRUE)
    })
    
    # Store selected cell information
    selected_cell <- reactiveVal(NULL)
    
    # Output: Explorer plot
    output$explorerPlot <- renderPlotly({
      req(processed())
      req(seurat_obj())
      
      # Get data for plotting
      plot_data <- FetchData(seurat_obj(), 
                             vars = c(input$explorer_color),
                             slot = "data")
      
      # Get reduction coordinates
      reduction_coords <- Embeddings(seurat_obj()[[input$explorer_reduction]])[, 1:2]
      colnames(reduction_coords) <- c("Dim1", "Dim2")
      
      plot_data <- cbind(plot_data, reduction_coords)
      plot_data$cell_id <- rownames(plot_data)
      
      # Create plot based on whether gene is selected
      if (input$explorer_gene != "none") {
        # Add gene expression data
        gene_expr <- FetchData(seurat_obj(), 
                               vars = input$explorer_gene,
                               slot = "data")[,1]
        plot_data$gene_expr <- gene_expr
        
        p <- plot_ly(plot_data, 
                     x = ~Dim1, 
                     y = ~Dim2,
                     color = ~gene_expr,
                     colors = viridis(100),
                     type = 'scatter',
                     mode = 'markers',
                     marker = list(size = 3),
                     text = ~paste("Cell:", cell_id,
                                   "<br>Cluster:", get(input$explorer_color),
                                   "<br>", input$explorer_gene, ":", 
                                   round(gene_expr, 2)),
                     hoverinfo = 'text',
                     customdata = ~cell_id,
                     source = "explorer_plot"
        ) %>%
          layout(
            xaxis = list(title = paste0(toupper(input$explorer_reduction), "_1")),
            yaxis = list(title = paste0(toupper(input$explorer_reduction), "_2")),
            title = paste(input$explorer_gene, "Expression"),
            colorbar = list(title = "Expression")
          )
      } else {
        # Color by metadata
        p <- plot_ly(plot_data, 
                     x = ~Dim1, 
                     y = ~Dim2,
                     color = ~get(input$explorer_color),
                     type = 'scatter',
                     mode = 'markers',
                     marker = list(size = 3),
                     text = ~paste("Cell:", cell_id,
                                   "<br>", input$explorer_color, ":", 
                                   get(input$explorer_color)),
                     hoverinfo = 'text',
                     customdata = ~cell_id,
                     source = "explorer_plot"
        ) %>%
          layout(
            xaxis = list(title = paste0(toupper(input$explorer_reduction), "_1")),
            yaxis = list(title = paste0(toupper(input$explorer_reduction), "_2")),
            title = paste("Colored by", input$explorer_color)
          )
      }
      
      # Add click event
      p %>% event_register("plotly_click")
    })
    
    # Handle click events
    observeEvent(event_data("plotly_click", source = "explorer_plot"), {
      click_data <- event_data("plotly_click", source = "explorer_plot")
      if (!is.null(click_data$customdata)) {
        selected_cell(click_data$customdata)
      }
    })
    
    # Output: Cell information
    output$cellInfo <- renderPrint({
      if (is.null(selected_cell())) {
        cat("Click on a cell to see its information")
      } else {
        req(seurat_obj())
        cell_id <- selected_cell()
        
        cat("Cell ID:", cell_id, "\n\n")
        
        # Get all metadata for this cell
        cell_meta <- seurat_obj()@meta.data[cell_id, , drop = FALSE]
        
        # Display metadata
        cat("=== Cell Metadata ===\n")
        for (col in colnames(cell_meta)) {
          cat(col, ":", cell_meta[[col]], "\n")
        }
        
        # Display top expressed genes if requested
        if (input$explorer_gene != "none") {
          cat("\n=== Gene Expression ===\n")
          cat(input$explorer_gene, ":", 
              FetchData(seurat_obj(), 
                        vars = input$explorer_gene,
                        cells = cell_id)[1,1], "\n")
        }
        
        # Get top 10 expressed genes for this cell
        cat("\n=== Top 10 Expressed Genes ===\n")
        expr_data <- GetAssayData(seurat_obj(), slot = "data")[, cell_id]
        top_genes <- sort(expr_data, decreasing = TRUE)[1:10]
        for (i in 1:10) {
          cat(names(top_genes)[i], ":", round(top_genes[i], 2), "\n")
        }
      }
    })
    
  })
}
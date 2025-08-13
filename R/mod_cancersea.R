#' cancersea UI Function
#'
#' @description A shiny Module for CancerSEA functional state analysis.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_cancersea_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      box(
        title = "Functional State Selection",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        selectInput(ns("cancersea_pathway"), "Select CancerSEA pathway:",
                    choices = NULL), # Will be updated in server
        actionButton(ns("calculate_pathway"), "Calculate/Update", class = "btn-info")
      )
    ),
    
    fluidRow(
      box(
        title = "Functional State Visualization",
        status = "success",
        solidHeader = TRUE,
        width = 6,
        plotlyOutput(ns("cancerseaFeature"), height = "500px")
      ),
      
      box(
        title = "Violin Plot by Cluster",
        status = "success",
        solidHeader = TRUE,
        width = 6,
        plotOutput(ns("cancerseaViolin"), height = "500px")
      )
    ),
    
    fluidRow(
      box(
        title = "Functional State Heatmap",
        status = "warning",
        solidHeader = TRUE,
        width = 12,
        plotOutput(ns("cancerseaHeatmap"), height = "600px")
      )
    )
  )
}

#' Helper function to prepare assay for scoring (Seurat v5 compatible)
#'
#' @param seurat_obj Seurat object
#' @return Prepared Seurat object
#' @noRd
check_and_prepare_assay <- function(seurat_obj) {
  # Ensure we're using the right assay for scoring
  if("SCT" %in% names(seurat_obj@assays)) {
    # If SCT is available, use RNA assay for module scoring
    DefaultAssay(seurat_obj) <- "RNA"
    
    # Check if we need to normalize - compatible with both v4 and v5
    if(inherits(seurat_obj@assays$RNA, "Assay5")) {
      # Seurat v5 - check layers (and handle NULL cases)
      data_layer <- seurat_obj@assays$RNA@layers$data
      counts_layer <- seurat_obj@assays$RNA@layers$counts
      
      if(is.null(data_layer) || 
         (is.matrix(data_layer) && is.matrix(counts_layer) && 
          max(data_layer) == max(counts_layer))) {
        # Data hasn't been normalized, normalize it
        seurat_obj <- NormalizeData(seurat_obj, assay = "RNA")
      }
    } else {
      # Seurat v4 - check slots
      if(max(seurat_obj@assays$RNA@data) == max(seurat_obj@assays$RNA@counts)) {
        # Data hasn't been normalized, normalize it
        seurat_obj <- NormalizeData(seurat_obj, assay = "RNA")
      }
    }
  } else {
    # Use RNA assay
    DefaultAssay(seurat_obj) <- "RNA"
  }
  
  return(seurat_obj)
}

#' cancersea Server Functions
#'
#' @import Seurat
#' @import viridis
#' @noRd 
mod_cancersea_server <- function(id, seurat_obj, processed, main_values){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
    # Load CancerSEA pathways when module initializes
    observe({
      if (!requireNamespace("cancersea", quietly = TRUE)) {
        showNotification("cancersea package required", type = "error")
        return()
      }
      
      tryCatch({
        # Try to load available_pathways data
        data('available_pathways', package = 'cancersea')
        if (exists("available_pathways")) {
          updateSelectInput(session, "cancersea_pathway", 
                            choices = available_pathways)
        } else {
          # Fallback: use common pathway names
          pathway_choices <- c(
            "Angiogenesis", "Apoptosis", "Cell.Cycle", "Differentiation",
            "DNA.Damage", "EMT", "Hypoxia", "Inflammation", "Invasion",
            "Metastasis", "Proliferation", "Quiescence", "Stemness"
          )
          updateSelectInput(session, "cancersea_pathway", 
                            choices = pathway_choices)
        }
      }, error = function(e) {
        showNotification(paste("Error loading pathways:", e$message), type = "warning")
      })
    })
    
    # Calculate CancerSEA pathway - UPDATED VERSION
    observeEvent(input$calculate_pathway, {
      req(processed())
      req(seurat_obj())
      req(input$cancersea_pathway)
      
      withProgress(message = paste('Calculating', input$cancersea_pathway, 'score...'), 
                   value = 0.5, {
                     
                     # ADDED: Prepare the assay for better scoring
                     temp_obj <- check_and_prepare_assay(seurat_obj())
                     
                     all_genes <- rownames(temp_obj)
                     #gene_list <- get(input$cancersea_pathway)$symbol
                     if (!requireNamespace("cancersea", quietly = TRUE)) {
                       showNotification("cancersea package required", type = "error")
                       return()
                     }
                     pathway_data <- get(input$cancersea_pathway, envir = asNamespace("cancersea"))
                     gene_list <- pathway_data$symbol
                     gene_list_filtered <- gene_list[gene_list %in% all_genes]
                     
                     if (length(gene_list_filtered) > 0) {
                       
                       # ULTRA SIMPLE: Use the exact same parameters as your working old app
                       tryCatch({
                         # Try with the old app's exact parameters first
                         main_values$seurat_obj <- AddModuleScore(
                           object = temp_obj,
                           features = list(gene_list_filtered),
                           name = paste0(input$cancersea_pathway, "_"),
                           ctrl = 20  # Same as old app
                           # No nbin parameter = uses default (safer)
                         )
                         
                         # Store the score column name
                         score_name <- paste0(input$cancersea_pathway, "_1")
                         main_values$cancersea_scores[[input$cancersea_pathway]] <- score_name
                         
                         showNotification(
                           paste("Successfully calculated", input$cancersea_pathway, "score"),
                           type = "message",
                           duration = 3
                         )
                         
                       }, error = function(e) {
                         # If that still fails, try with absolute minimal parameters
                         if (grepl("Insufficient data", e$message) || grepl("bin", e$message)) {
                           
                           tryCatch({
                             # FALLBACK: Absolutely minimal parameters
                             main_values$seurat_obj <- AddModuleScore(
                               object = temp_obj,
                               features = list(gene_list_filtered),
                               name = paste0(input$cancersea_pathway, "_"),
                               ctrl = 5,     # Minimal control genes
                               nbin = 5      # Minimal bins
                             )
                             
                             score_name <- paste0(input$cancersea_pathway, "_1")
                             main_values$cancersea_scores[[input$cancersea_pathway]] <- score_name
                             
                             showNotification(
                               paste("Calculated", input$cancersea_pathway, "score with minimal parameters"),
                               type = "warning",
                               duration = 5
                             )
                             
                           }, error = function(e2) {
                             # If even that fails, try without any binning at all
                             tryCatch({
                               # LAST RESORT: Disable sophisticated control matching
                               main_values$seurat_obj <- AddModuleScore(
                                 object = temp_obj,
                                 features = list(gene_list_filtered),
                                 name = paste0(input$cancersea_pathway, "_"),
                                 ctrl = 1,      # Just 1 control gene per feature
                                 nbin = 1       # No binning
                               )
                               
                               score_name <- paste0(input$cancersea_pathway, "_1")
                               main_values$cancersea_scores[[input$cancersea_pathway]] <- score_name
                               
                               showNotification(
                                 paste("Calculated", input$cancersea_pathway, "score with basic parameters (no control matching)"),
                                 type = "warning",
                                 duration = 7
                               )
                               
                             }, error = function(e3) {
                               showNotification(
                                 paste("Failed to calculate", input$cancersea_pathway, "score:", e3$message),
                                 type = "error",
                                 duration = 10
                               )
                             })
                           })
                           
                         } else {
                           showNotification(
                             paste("Error calculating", input$cancersea_pathway, "score:", e$message),
                             type = "error",
                             duration = 10
                           )
                         }
                       })
                       
                     } else {
                       showNotification(
                         paste("No genes from", input$cancersea_pathway, "pathway found in dataset"),
                         type = "warning",
                         duration = 5
                       )
                     }
                   })
    })
    
    # Get current pathway score name
    current_score <- reactive({
      req(input$cancersea_pathway)
      main_values$cancersea_scores[[input$cancersea_pathway]]
    })
    
    # Output: CancerSEA feature plot
    output$cancerseaFeature <- renderPlotly({
      req(current_score())
      req(seurat_obj())
      
      p <- FeaturePlot(seurat_obj(), 
                       features = current_score(),
                       reduction = "umap") + 
        scale_color_viridis() +
        theme_minimal() +
        ggtitle(paste(input$cancersea_pathway, "Score"))
      
      ggplotly(p)
    })
    
    # Output: CancerSEA violin plot
    output$cancerseaViolin <- renderPlot({
      req(current_score())
      req(seurat_obj())
      
      VlnPlot(seurat_obj(), 
              features = current_score(),
              pt.size = 0) + 
        theme_minimal() +
        ggtitle(paste(input$cancersea_pathway, "Score by Cluster"))
    })
    
    # Output: CancerSEA heatmap
    # Output: CancerSEA heatmap (fixed version)
    output$cancerseaHeatmap <- renderPlot({
      req(length(main_values$cancersea_scores) > 0)
      req(seurat_obj())
      
      if (length(main_values$cancersea_scores) >= 2) {
        # Create a matrix of average scores per cluster
        score_matrix <- matrix(
          nrow = length(unique(Idents(seurat_obj()))), 
          ncol = length(main_values$cancersea_scores)
        )
        
        for (i in 1:length(main_values$cancersea_scores)) {
          score_col <- main_values$cancersea_scores[[i]]
          avg_scores <- aggregate(
            seurat_obj()@meta.data[[score_col]], 
            by = list(Idents(seurat_obj())), 
            FUN = mean
          )
          score_matrix[, i] <- avg_scores$x
        }
        
        rownames(score_matrix) <- paste("Cluster", levels(Idents(seurat_obj())))
        colnames(score_matrix) <- names(main_values$cancersea_scores)
        
        # Scale and plot
        score_matrix_scaled <- t(scale(t(score_matrix)))
        
        # Handle any NA values from scaling
        score_matrix_scaled[is.na(score_matrix_scaled)] <- 0
        
        # Fixed color palette - use viridis to avoid RColorBrewer issues
        heatmap(score_matrix_scaled, 
                col = viridis::viridis(100),  # Use viridis instead
                scale = "none",
                margins = c(10, 8),
                main = "CancerSEA Pathway Scores by Cluster")
      } else {
        plot.new()
        text(0.5, 0.5, "Calculate at least 2 pathways to show heatmap", 
             cex = 1.5, col = "gray40")
      }
    })
    
  })
}
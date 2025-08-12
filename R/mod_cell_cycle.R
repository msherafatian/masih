#' cell_cycle UI Function
#'
#' @description A shiny Module for cell cycle analysis.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_cell_cycle_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      box(
        title = "Cell Cycle Phase Distribution",
        status = "primary",
        solidHeader = TRUE,
        width = 6,
        plotlyOutput(ns("cellCyclePlot"), height = "500px")
      ),
      
      box(
        title = "Cell Cycle Scores",
        status = "primary",
        solidHeader = TRUE,
        width = 6,
        radioButtons(ns("cycle_score"), "Select score:",
                     choices = c("S.Score", "G2M.Score")),
        plotlyOutput(ns("cycleScorePlot"), height = "450px")
      )
    ),
    
    fluidRow(
      box(
        title = "Phase Distribution by Cluster",
        status = "info",
        solidHeader = TRUE,
        width = 12,
        plotOutput(ns("phaseByCluster"), height = "400px")
      )
    )
  )
}

#' cell_cycle Server Functions
#'
#' @import Seurat
#' @import ggplot2
#' @import dplyr
#' @noRd 
mod_cell_cycle_server <- function(id, seurat_obj, processed){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
    # Check if cell cycle scores are available
    has_cell_cycle <- reactive({
      req(seurat_obj())
      all(c("S.Score", "G2M.Score", "Phase") %in% 
            colnames(seurat_obj()@meta.data))
    })
    
    # Output: Cell cycle plot
    output$cellCyclePlot <- renderPlotly({
      req(processed())
      req(seurat_obj())
      req(has_cell_cycle())
      
      p <- DimPlot(seurat_obj(), 
                   reduction = "umap",
                   group.by = "Phase",
                   label = TRUE,
                   label.size = 4) + 
        theme_minimal() + 
        ggtitle("Cell Cycle Phases") +
        scale_color_manual(values = c("G1" = "#E41A1C", 
                                      "G2M" = "#377EB8", 
                                      "S" = "#4DAF4A"))
      
      ggplotly(p)
    })
    
    # Output: Cycle score plot
    output$cycleScorePlot <- renderPlotly({
      req(processed())
      req(seurat_obj())
      req(has_cell_cycle())
      
      p <- FeaturePlot(seurat_obj(), 
                       features = input$cycle_score,
                       reduction = "umap") + 
        scale_color_viridis() +
        theme_minimal() +
        ggtitle(paste(input$cycle_score))
      
      ggplotly(p)
    })
    
    # Output: Phase by cluster
    output$phaseByCluster <- renderPlot({
      req(processed())
      req(seurat_obj())
      
      # Check if we have cell cycle data
      if (!has_cell_cycle()) {
        create_empty_plot("Cell cycle scores not available")
        return()
      }
      
      tryCatch({
        phase_data <- data.frame(
          Cluster = Idents(seurat_obj()),
          Phase = seurat_obj()$Phase
        )
        
        phase_props <- phase_data %>%
          group_by(Cluster, Phase) %>%
          summarise(count = n(), .groups = "drop") %>%
          group_by(Cluster) %>%
          mutate(prop = count / sum(count))
        
        ggplot(phase_props, aes(x = Cluster, y = prop, fill = Phase)) +
          geom_bar(stat = "identity") +
          theme_minimal() +
          labs(y = "Proportion", 
               title = "Cell Cycle Phase Distribution by Cluster") +
          scale_fill_manual(values = c("G1" = "#E41A1C", 
                                       "G2M" = "#377EB8", 
                                       "S" = "#4DAF4A")) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }, error = function(e) {
        create_empty_plot(paste("Error creating plot:", e$message))
      })
    })
    
    # Show notification if cell cycle scores are missing
    observe({
      if (processed() && !has_cell_cycle()) {
        showNotification(
          "Cell cycle scores not found. Please run cell cycle scoring during data processing.",
          type = "warning",
          duration = NULL,
          id = "cell_cycle_warning"
        )
      }
    })
    
  })
}
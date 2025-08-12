#' trajectory UI Function
#'
#' @description A shiny Module for trajectory analysis.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_trajectory_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      box(
        title = "Trajectory Analysis Settings",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        
        h4("Select Cell Types for Trajectory"),
        p("Choose which cell types to include in the trajectory analysis (based on your annotation column)."),
        
        fluidRow(
          column(6,
                 selectInput(ns("trajectory_annotation"), 
                             "Annotation column:",
                             choices = NULL,
                             selected = "annotate")
          ),
          column(6,
                 selectInput(ns("trajectory_celltypes"), 
                             "Select cell identities/types:", 
                             choices = NULL,
                             multiple = TRUE)
          )
        ),
        
        fluidRow(
          column(6,
                 selectInput(ns("trajectory_neotype"), 
                             "NeoType column (optional):",
                             choices = c("None" = "none"))
          ),
          column(6,
                 selectInput(ns("trajectory_neotypes"), 
                             "Select NeoTypes (optional):",
                             choices = NULL,
                             multiple = TRUE)
          )
        ),
        
        hr(),
        h4("Trajectory Parameters"),
        
        fluidRow(
          column(4,
                 numericInput(ns("trajectory_dims"), 
                              "Number of PCs to use:",
                              value = 10,
                              min = 2, 
                              max = 50)
          ),
          column(4,
                 selectInput(ns("trajectory_start"), 
                             "Starting cluster/cell type (optional):",
                             choices = c("Auto-detect" = "auto"))
          ),
          column(4,
                 selectInput(ns("trajectory_end"), 
                             "Ending cluster/cell type (optional):",
                             choices = c("Auto-detect" = "auto"))
          )
        ),
        
        fluidRow(
          column(12,
                 p("💡 Tip: Leave start/end as 'Auto-detect' to let Slingshot determine the trajectory automatically. Select specific clusters to define the trajectory direction.", 
                   style = "font-size: 12px; color: #666; font-style: italic;")
          )
        ),
        
        hr(),
        h5("Filter Preview:"),
        verbatimTextOutput(ns("trajectory_filter_preview")),
        
        hr(),
        actionButton(ns("run_trajectory"), 
                     "Run Trajectory Analysis", 
                     class = "btn-success"),
        
        conditionalPanel(
          condition = paste0("output['", ns("trajectoryCalculated"), "']"),
          hr(),
          h4("Export Options:"),
          downloadButton(ns("downloadTrajectory"), 
                         "Download Trajectory Data", 
                         class = "btn-info")
        )
      )
    ),
    
    fluidRow(
      box(
        title = "Trajectory Visualization",
        status = "success",
        solidHeader = TRUE,
        width = 12,
        conditionalPanel(
          condition = paste0("output['", ns("trajectoryCalculated"), "']"),
          tabsetPanel(
            tabPanel("Main Trajectory", 
                     plotOutput(ns("trajectoryPlot"), height = "600px")),
            tabPanel("Stratified Trajectory", 
                     plotOutput(ns("trajectorySeuratStyle"), height = "600px")),
            tabPanel("Pseudotime", 
                     plotOutput(ns("pseudotimePlot"), height = "600px"))
          )
        ),
        conditionalPanel(
          condition = paste0("!output['", ns("trajectoryCalculated"), "']"),
          h4("No trajectory calculated yet. Select cell types and click 'Run Trajectory Analysis'.",
             style = "color: #999; text-align: center; padding: 50px;")
        )
      )
    ),
    
    fluidRow(
      box(
        title = "Trajectory Statistics",
        status = "info",
        solidHeader = TRUE,
        width = 12,
        conditionalPanel(
          condition = paste0("output['", ns("trajectoryCalculated"), "']"),
          verbatimTextOutput(ns("trajectoryStats"))
        )
      )
    )
  )
}
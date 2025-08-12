#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinydashboard
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    
    # Your application UI logic
    shinydashboard::dashboardPage(
      shinydashboard::dashboardHeader(title = "MASIH"),
      
      shinydashboard::dashboardSidebar(
        shinydashboard::sidebarMenu(
          shinydashboard::menuItem("Data Upload & Processing", 
                                   tabName = "upload", 
                                   icon = icon("upload")),
          shinydashboard::menuItem("Cluster Analysis", 
                                   tabName = "clusters", 
                                   icon = icon("project-diagram")),
          shinydashboard::menuItem("CancerSEA States", 
                                   tabName = "cancersea", 
                                   icon = icon("dna")),
          shinydashboard::menuItem("Cell Cycle Analysis", 
                                   tabName = "cellcycle", 
                                   icon = icon("sync")),
          shinydashboard::menuItem("Interactive Explorer", 
                                   tabName = "explorer", 
                                   icon = icon("search")),
          shinydashboard::menuItem("Marker Genes", 
                                   tabName = "markers", 
                                   icon = icon("dna")),
          shinydashboard::menuItem("Trajectory Analysis", 
                                   tabName = "trajectory", 
                                   icon = icon("route")),
          shinydashboard::menuItem("Comparative Analysis", 
                                   tabName = "compare", 
                                   icon = icon("chart-bar")),
          shinydashboard::menuItem("Export & Citation", 
                                   tabName = "export", 
                                   icon = icon("download"))
        )
      ),
      
      shinydashboard::dashboardBody(
        tags$head(
          tags$style(HTML("
            .content-wrapper, .right-side {
              background-color: #f4f4f4;
              min-height: 100vh !important;
              height: auto !important;
              overflow-y: visible !important;
            }
            .wrapper {
              overflow-y: visible !important;
            }
            .box {
              border-radius: 8px;
              height: auto !important;
            }
            .box-body {
              height: auto !important;
              overflow-y: visible !important;
            }
          "))
        ),
        
        shinydashboard::tabItems(
          shinydashboard::tabItem(tabName = "upload", mod_upload_ui("upload")),
          shinydashboard::tabItem(tabName = "clusters", mod_cluster_analysis_ui("cluster_analysis")),
          shinydashboard::tabItem(tabName = "cancersea", mod_cancersea_ui("cancersea")),
          shinydashboard::tabItem(tabName = "cellcycle", mod_cell_cycle_ui("cell_cycle")),
          shinydashboard::tabItem(tabName = "explorer", mod_explorer_ui("explorer")),
          shinydashboard::tabItem(tabName = "markers", mod_markers_ui("markers")),
          shinydashboard::tabItem(tabName = "trajectory", mod_trajectory_ui("trajectory")),
          shinydashboard::tabItem(tabName = "compare", mod_compare_ui("compare")),
          shinydashboard::tabItem(tabName = "export", mod_export_ui("export"))
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
      app_title = "MASIH"
    ),
    # Add here other external resources
    # Include JavaScript for clipboard functionality
    tags$script(HTML("
      Shiny.addCustomMessageHandler('copyToClipboard', function(message) {
        navigator.clipboard.writeText(message.text).then(function() {
          Shiny.setInputValue(message.id, 'copied', {priority: 'event'});
        });
      });
    "))
  )
}
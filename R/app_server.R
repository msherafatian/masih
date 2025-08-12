#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  
  # Set max upload size
  options(shiny.maxRequestSize = app_config$max_upload_size)
  
  # Initialize CancerSEA pathways
  data('available_pathways', package = 'cancersea')
  
  # Main reactive values that will be shared across modules
  main_values <- reactiveValues(
    seurat_obj = NULL,
    processed = FALSE,
    cancersea_scores = list(),
    markers = NULL,
    trajectory = NULL
  )
  
  # Module: Upload and Processing
  upload_result <- mod_upload_server(
    "upload",
    main_values = main_values
  )
  
  # Module: Cluster Analysis
  mod_cluster_analysis_server(
    "cluster_analysis",
    seurat_obj = reactive(main_values$seurat_obj),
    processed = reactive(main_values$processed)
  )
  
  # Module: CancerSEA
  mod_cancersea_server(
    "cancersea",
    seurat_obj = reactive(main_values$seurat_obj),
    processed = reactive(main_values$processed),
    main_values = main_values
  )
  
  # Module: Cell Cycle
  mod_cell_cycle_server(
    "cell_cycle",
    seurat_obj = reactive(main_values$seurat_obj),
    processed = reactive(main_values$processed)
  )
  
  # Module: Explorer
  mod_explorer_server(
    "explorer",
    seurat_obj = reactive(main_values$seurat_obj),
    processed = reactive(main_values$processed)
  )
  
  # Module: Markers
  mod_markers_server(
    "markers",
    seurat_obj = reactive(main_values$seurat_obj),
    processed = reactive(main_values$processed),
    main_values = main_values
  )
  
  # Module: Trajectory
  mod_trajectory_server(
    "trajectory",
    seurat_obj = reactive(main_values$seurat_obj),
    processed = reactive(main_values$processed),
    main_values = main_values
  )
  
  # Module: Compare
  mod_compare_server(
    "compare",
    seurat_obj = reactive(main_values$seurat_obj),
    processed = reactive(main_values$processed),
    cancersea_scores = reactive(main_values$cancersea_scores),
    main_values = main_values
  )
  
  # Module: Export
  mod_export_server(
    "export",
    main_values = main_values
  )
}
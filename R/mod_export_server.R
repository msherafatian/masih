#' export Server Functions
#'
#' @import openxlsx
#' @noRd 
mod_export_server <- function(id, main_values){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
    # Update gene choices when Seurat object changes
    observe({
      req(main_values$seurat_obj)
      
      gene_choices <- rownames(main_values$seurat_obj)
      updateSelectInput(session, "export_gene", choices = gene_choices)
    })
    
    # Download single plot
    output$downloadPlot <- downloadHandler(
      filename = function() {
        paste0("masih_plot_", Sys.Date(), ".", input$export_format)
      },
      content = function(file) {
        req(main_values$seurat_obj)
        
        # Create the selected plot
        p <- create_export_plot(
          plot_type = input$export_plot_type,
          seurat_obj = main_values$seurat_obj,
          main_values = main_values,
          gene = input$export_gene,
          custom_title = input$export_title,
          remove_axes = input$export_no_axes,
          white_bg = input$export_white_bg
        )
        
        # Define which plots need base R export vs ggplot export
        base_plots <- c("tree", "cancersea_heatmap", "correlation_matrix", 
                        "comparison_heatmap", "trajectory_main", "pseudotime")
        
        # Handle different plot types
        if (input$export_plot_type %in% base_plots || is.function(p)) {
          # Base R plots or functions
          if (is.function(p)) {
            export_base_plot(
              plot_func = p,
              file = file,
              format = input$export_format,
              width = input$export_width,
              height = input$export_height,
              dpi = as.numeric(input$export_dpi)
            )
          } else {
            export_base_plot(
              plot_func = function() print(p),
              file = file,
              format = input$export_format,
              width = input$export_width,
              height = input$export_height,
              dpi = as.numeric(input$export_dpi)
            )
          }
        } else if (!is.null(p) && inherits(p, "ggplot")) {
          # ggplot objects
          ggsave(
            file, 
            plot = p, 
            width = input$export_width, 
            height = input$export_height,
            dpi = as.numeric(input$export_dpi),
            device = input$export_format
          )
        } else {
          showNotification("Plot not available for export", type = "error")
        }
      }
    )
    
    # Download all plots
    output$downloadAllPlots <- downloadHandler(
      filename = function() {
        paste0("masih_all_plots_", Sys.Date(), ".zip")
      },
      content = function(file) {
        req(main_values$seurat_obj)
        
        # Create temporary directory
        temp_dir <- tempdir()
        plots_dir <- file.path(temp_dir, "masih_plots")
        dir.create(plots_dir, showWarnings = FALSE)
        
        withProgress(message = 'Generating plots...', value = 0, {
          
          # Generate all available plots
          plot_types <- get_available_plot_types(main_values)
          
          # Define which plots need base R export vs ggplot export
          base_plots <- c("tree", "cancersea_heatmap", "correlation_matrix", 
                          "comparison_heatmap", "trajectory_main", "trajectory_stratified", "pseudotime") 
          
          for (i in seq_along(plot_types)) {
            plot_type <- plot_types[i]
            
            incProgress(1/length(plot_types), 
                        detail = paste("Creating", plot_type))
            
            tryCatch({
              file_name <- paste0(plot_type, ".", input$export_format)
              file_path <- file.path(plots_dir, file_name)
              
              p <- create_export_plot(
                plot_type = plot_type,
                seurat_obj = main_values$seurat_obj,
                main_values = main_values,
                white_bg = input$export_white_bg
              )
              
              if (!is.null(p)) {
                if (plot_type %in% base_plots || is.function(p)) {
                  # Handle base R plots (functions)
                  if (is.function(p)) {
                    export_base_plot(
                      plot_func = p,
                      file = file_path,
                      format = input$export_format,
                      width = input$export_width,
                      height = input$export_height,
                      dpi = as.numeric(input$export_dpi)
                    )
                  } else {
                    # Handle base plots that return NULL but need special handling
                    export_base_plot(
                      plot_func = function() print(p),
                      file = file_path,
                      format = input$export_format,
                      width = input$export_width,
                      height = input$export_height,
                      dpi = as.numeric(input$export_dpi)
                    )
                  }
                } else {
                  # Handle ggplot objects
                  if (inherits(p, "ggplot")) {
                    ggsave(
                      file_path, 
                      plot = p, 
                      width = input$export_width, 
                      height = input$export_height,
                      dpi = as.numeric(input$export_dpi),
                      device = input$export_format
                    )
                  }
                }
              }
            }, error = function(e) {
              message("Error creating ", plot_type, ": ", e$message)
            })
          }
        })
        
        # Create zip file
        zip(file, plots_dir)
      }
    )
    
    # Download data
    output$downloadData <- downloadHandler(
      filename = function() {
        paste0("masih_data_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        req(main_values$seurat_obj)
        
        wb <- createWorkbook()
        
        withProgress(message = 'Preparing data export...', value = 0, {
          
          # Cell metadata
          if ("metadata" %in% input$export_data_options) {
            incProgress(1/7, detail = "Adding metadata...")
            addWorksheet(wb, "Cell_Metadata")
            writeData(wb, "Cell_Metadata", main_values$seurat_obj@meta.data)
          }
          
          # Cluster statistics
          if ("cluster_stats" %in% input$export_data_options) {
            incProgress(1/7, detail = "Calculating cluster stats...")
            stats_df <- calculate_cluster_stats(main_values$seurat_obj)
            addWorksheet(wb, "Cluster_Statistics")
            writeData(wb, "Cluster_Statistics", stats_df)
          }
          
          # CancerSEA scores
          if ("cancersea_scores" %in% input$export_data_options && 
              length(main_values$cancersea_scores) > 0) {
            incProgress(1/7, detail = "Adding CancerSEA scores...")
            score_cols <- unlist(main_values$cancersea_scores)
            score_data <- main_values$seurat_obj@meta.data[, score_cols, drop = FALSE]
            colnames(score_data) <- names(main_values$cancersea_scores)
            addWorksheet(wb, "CancerSEA_Scores")
            writeData(wb, "CancerSEA_Scores", score_data)
          }
          
          # Average pathway scores
          if ("pathway_avg" %in% input$export_data_options && 
              length(main_values$cancersea_scores) > 0) {
            incProgress(1/7, detail = "Calculating pathway averages...")
            avg_scores <- calculate_pathway_averages(
              main_values$seurat_obj,
              main_values$cancersea_scores
            )
            addWorksheet(wb, "Pathway_Averages")
            writeData(wb, "Pathway_Averages", avg_scores)
          }
          
          # Cell cycle proportions
          if ("cellcycle_props" %in% input$export_data_options) {
            incProgress(1/7, detail = "Calculating cell cycle proportions...")
            if ("Phase" %in% colnames(main_values$seurat_obj@meta.data)) {
              cc_props <- calculate_cellcycle_proportions(main_values$seurat_obj)
              addWorksheet(wb, "Cell_Cycle_Proportions")
              writeData(wb, "Cell_Cycle_Proportions", cc_props)
            }
          }
          
          # Marker genes
          if ("markers" %in% input$export_data_options && 
              !is.null(main_values$markers)) {
            incProgress(1/7, detail = "Adding marker genes...")
            addWorksheet(wb, "Marker_Genes")
            writeData(wb, "Marker_Genes", main_values$markers)
          }
          
          # Trajectory results
          if ("trajectory" %in% input$export_data_options && 
              !is.null(main_values$trajectory)) {
            incProgress(1/7, detail = "Adding trajectory results...")
            traj_summary <- summarize_trajectory(main_values$trajectory)
            addWorksheet(wb, "Trajectory_Summary")
            writeData(wb, "Trajectory_Summary", traj_summary)
          }
        })
        
        saveWorkbook(wb, file, overwrite = TRUE)
      }
    )
    
    # Download Seurat object
    output$downloadSeurat <- downloadHandler(
      filename = function() {
        paste0("masih_seurat_object_", Sys.Date(), ".rds")
      },
      content = function(file) {
        req(main_values$seurat_obj)
        saveRDS(main_values$seurat_obj, file)
      }
    )
    
    # Generate methods text
    output$methodsText <- renderText({
      req(main_values$processed)
      req(main_values$seurat_obj)
      
      generate_methods_text(main_values)
    })
    
    # Copy methods to clipboard
    observeEvent(input$copyMethods, {
      methods_text <- generate_methods_text(main_values)
      
      # Use session$sendCustomMessage instead of runjs
      session$sendCustomMessage('copyToClipboard', list(
        text = methods_text,
        id = ns("clipboard_status")
      ))
    })
    
    # Show notification when copied
    observeEvent(input$clipboard_status, {
      showNotification("Methods text copied to clipboard!", 
                       type = "success", duration = 2)
    })
    
    # Generate citations
    output$citationText <- renderText({
      generate_citations()
    })
    
    # Generate figure legend
    output$legendText <- renderText({
      req(main_values$seurat_obj)
      generate_figure_legend(main_values$seurat_obj)
    })
    
    # Analysis summary
    output$analysisSummary <- renderPrint({
      req(main_values$seurat_obj)
      print_analysis_summary(main_values)
    })
    
  })
}
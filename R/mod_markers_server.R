#' markers Server Functions
#'
#' @import Seurat
#' @import dplyr
#' @import openxlsx
#' @noRd 
mod_markers_server <- function(id, seurat_obj, processed, main_values){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
    # Local reactive values
    local_values <- reactiveValues(
      markers = NULL,
      top_markers = NULL,
      heatmap_genes = NULL
    )
    
    # Update cluster choices when Seurat object changes
    observe({
      req(processed())
      req(seurat_obj())
      
      clusters <- levels(Idents(seurat_obj()))
      updateSelectInput(session, "marker_cluster1", choices = clusters)
      updateSelectInput(session, "marker_cluster2", choices = clusters)
    })
    
    # Marker calculation status
    output$markersCalculated <- reactive({
      return(!is.null(local_values$markers))
    })
    outputOptions(output, "markersCalculated", suspendWhenHidden = FALSE)
    
    # Calculate marker genes
    observeEvent(input$calculate_markers, {
      req(processed())
      req(seurat_obj())
      
      withProgress(message = 'Finding marker genes...', value = 0, {
        
        tryCatch({
          if (input$marker_comparison == "one_vs_all") {
            # Find markers for all clusters
            incProgress(0.3, detail = "Running differential expression...")
            
            markers <- FindAllMarkers(
              seurat_obj(),
              test.use = input$marker_test,
              min.pct = input$marker_minpct,
              logfc.threshold = input$marker_logfc,
              only.pos = input$marker_onlypos,
              verbose = FALSE
            )
            
            incProgress(0.7, detail = "Processing results...")
            
            # Add gene column if it doesn't exist
            if (!"gene" %in% colnames(markers)) {
              markers$gene <- rownames(markers)
            }
            
            # Get top markers per cluster
            top_markers <- markers %>%
              group_by(cluster) %>%
              arrange(desc(avg_log2FC)) %>%
              slice_head(n = input$marker_topn) %>%
              ungroup()
            
          } else {
            # Pairwise comparison
            incProgress(0.5, detail = "Running pairwise comparison...")
            
            markers <- FindMarkers(
              seurat_obj(),
              ident.1 = input$marker_cluster1,
              ident.2 = input$marker_cluster2,
              test.use = input$marker_test,
              min.pct = input$marker_minpct,
              logfc.threshold = input$marker_logfc,
              only.pos = input$marker_onlypos,
              verbose = FALSE
            )
            
            markers$gene <- rownames(markers)
            markers$cluster <- input$marker_cluster1
            top_markers <- head(markers, input$marker_topn)
          }
          
          incProgress(1, detail = "Complete!")
          
          # Store results
          local_values$markers <- markers
          local_values$top_markers <- top_markers
          
          # Update main values for export module
          main_values$markers <- markers
          
          # Update gene choices for visualization
          marker_genes <- unique(markers$gene)
          updateSelectInput(session, "visualize_marker", 
                            choices = marker_genes,
                            selected = marker_genes[1])
          
          showNotification("Marker genes identified successfully!", 
                           type = "message", duration = 3)
          
        }, error = function(e) {
          showNotification(paste("Error finding markers:", e$message), 
                           type = "error", duration = 5)
        })
      })
    })
    
    # Clear markers
    observeEvent(input$clearMarkers, {
      local_values$markers <- NULL
      local_values$top_markers <- NULL
      local_values$heatmap_genes <- NULL
      main_values$markers <- NULL
      updateSelectInput(session, "visualize_marker", choices = NULL)
      showNotification("Marker results cleared", type = "message", duration = 2)
    })
    
    # Output: Marker table
    output$markerTable <- renderDT({
      req(local_values$markers)
      
      # Format the table
      markers_display <- local_values$markers %>%
        mutate(across(where(is.numeric), ~round(., 3))) %>%
        arrange(cluster, desc(avg_log2FC))
      
      datatable(markers_display,
                options = list(
                  pageLength = 25,
                  scrollX = TRUE,
                  searchHighlight = TRUE,
                  dom = 'Bfrtip',
                  buttons = c('copy', 'csv', 'excel')
                ),
                filter = 'top',
                rownames = FALSE) %>%
        formatStyle('avg_log2FC',
                    background = styleColorBar(markers_display$avg_log2FC, 'lightblue'),
                    backgroundSize = '100% 90%',
                    backgroundRepeat = 'no-repeat',
                    backgroundPosition = 'center')
    })
    
    # Output: Marker heatmap
    output$markerHeatmap <- renderPlot({
      # Early return with empty plot if no data
      if (is.null(local_values$top_markers)) {
        create_empty_plot("No markers calculated yet. Click 'Find Marker Genes' to start.")
        return()
      }
      
      # Your actual plotting code wrapped in tryCatch
      tryCatch({
        # Get top genes - limit number based on user input
        top_genes <- unique(local_values$top_markers$gene)
        
        if (length(top_genes) > input$marker_display_genes) {
          top_genes <- top_genes[1:input$marker_display_genes]
        }
        
        # Store for export
        local_values$heatmap_genes <- top_genes
        
        # Create heatmap
        DoHeatmap(seurat_obj(), 
                  features = top_genes,
                  size = 3,
                  angle = 0,
                  hjust = 0.5,
                  label = FALSE,
                  group.bar.height = 0.02) + 
          theme(axis.text.y = element_text(size = max(6, min(10, 300/length(top_genes)))),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
      }, error = function(e) {
        create_empty_plot(paste("Error creating heatmap:", e$message))
      })
    })
    
    # Dynamic height for heatmap
    output$dynamic_heatmap_ui <- renderUI({
      req(local_values$top_markers)
      n_genes <- min(length(unique(local_values$top_markers$gene)), 
                     input$marker_display_genes)
      height_px <- max(600, n_genes * 35 + 200)
      
      plotOutput(ns("markerHeatmap"), height = paste0(height_px, "px"))
    })
    
    # Output: Marker dot plot
    output$markerDotPlot <- renderPlot({
      req(local_values$top_markers)
      req(seurat_obj())
      
      # Get top genes per cluster
      if (input$marker_comparison == "one_vs_all") {
        top_genes <- local_values$top_markers %>%
          group_by(cluster) %>%
          arrange(desc(avg_log2FC)) %>%
          slice_head(n = min(3, input$marker_topn)) %>%
          pull(gene) %>%
          unique()
      } else {
        top_genes <- head(local_values$top_markers$gene, 
                          min(10, nrow(local_values$top_markers)))
      }
      
      DotPlot(seurat_obj(), 
              features = top_genes,
              dot.scale = 8,
              col.min = -2,
              col.max = 2) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
              axis.text.y = element_text(size = 10)) +
        coord_flip()
    })
    
    # Dynamic height for dot plot
    output$dynamic_dotplot_ui <- renderUI({
      req(local_values$top_markers)
      
      n_genes <- if (input$marker_comparison == "one_vs_all") {
        length(unique(local_values$top_markers %>%
                        group_by(cluster) %>%
                        arrange(desc(avg_log2FC)) %>%
                        slice_head(n = min(3, input$marker_topn)) %>%
                        pull(gene)))
      } else {
        min(10, nrow(local_values$top_markers))
      }
      
      height_px <- max(600, n_genes * 35 + 200)
      
      plotOutput(ns("markerDotPlot"), height = paste0(height_px, "px"))
    })
    
    # Output: Selected marker plot
    output$selectedMarkerPlot <- renderPlot({
      req(input$update_marker_plot)
      req(input$visualize_marker)
      req(seurat_obj())
      
      isolate({
        if (input$marker_plot_type == "feature") {
          create_feature_plot(seurat_obj(), input$visualize_marker)
        } else if (input$marker_plot_type == "violin") {
          create_violin_plot(seurat_obj(), input$visualize_marker)
        } else if (input$marker_plot_type == "ridge") {
          RidgePlot(seurat_obj(),
                    features = input$visualize_marker) +
            theme_minimal()
        }
      })
    })
    
    # Download all markers
    output$downloadMarkerTable <- downloadHandler(
      filename = function() {
        paste0("marker_genes_all_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        req(local_values$markers)
        
        wb <- createWorkbook()
        
        if (input$marker_comparison == "one_vs_all") {
          # Create a sheet for each cluster
          clusters <- unique(local_values$markers$cluster)
          
          for (clust in clusters) {
            clust_markers <- local_values$markers %>%
              filter(cluster == clust) %>%
              arrange(desc(avg_log2FC))
            
            sheet_name <- paste0("Cluster_", clust)
            addWorksheet(wb, sheet_name)
            writeData(wb, sheet_name, clust_markers)
          }
          
          # Add summary sheet
          addWorksheet(wb, "Top_Markers_Summary")
          writeData(wb, "Top_Markers_Summary", local_values$top_markers)
          
        } else {
          # Single sheet for pairwise
          addWorksheet(wb, "Pairwise_Markers")
          writeData(wb, "Pairwise_Markers", local_values$markers)
        }
        
        saveWorkbook(wb, file, overwrite = TRUE)
      }
    )
    
    # Download top markers
    output$downloadTopMarkers <- downloadHandler(
      filename = function() {
        paste0("top_marker_genes_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(local_values$top_markers)
        write.csv(local_values$top_markers, file, row.names = FALSE)
      }
    )
    
    # Return values for export module
    return(list(
      markers = reactive(local_values$markers),
      top_markers = reactive(local_values$top_markers),
      heatmap_genes = reactive(local_values$heatmap_genes)
    ))
    
  })
}
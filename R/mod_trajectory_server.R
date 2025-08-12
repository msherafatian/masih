#' trajectory Server Functions
#'
#' @import Seurat
#' @import slingshot
#' @import SingleCellExperiment
#' @import RColorBrewer
#' @noRd 
mod_trajectory_server <- function(id, seurat_obj, processed, main_values){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
    # Local reactive values
    local_values <- reactiveValues(
      trajectory = NULL,
      trajectory_subset = NULL,
      slingshot_obj = NULL,
      slingshot_sce = NULL
    )
    
    # Update annotation choices when Seurat object changes
    observe({
      req(seurat_obj())
      
      if (!is.null(seurat_obj())) {
        meta_cols <- colnames(seurat_obj()@meta.data)
        
        # Prioritize useful columns for trajectory analysis
        priority_cols <- c("seurat_clusters", "annotate", "cell_type", "cluster")
        available_priority <- priority_cols[priority_cols %in% meta_cols]
        
        # Reorder to put useful columns first
        if(length(available_priority) > 0) {
          other_cols <- meta_cols[!meta_cols %in% priority_cols]
          meta_cols <- c(available_priority, other_cols)
        }
        
        # Update annotation column choices
        updateSelectInput(session, "trajectory_annotation", 
                          choices = meta_cols,
                          selected = if("seurat_clusters" %in% meta_cols) "seurat_clusters" 
                          else if("annotate" %in% meta_cols) "annotate" 
                          else meta_cols[1])
        
        # Update NeoType column choices
        neotype_choices <- c("None" = "none")
        for(col in meta_cols) {
          neotype_choices[col] <- col
        }
        
        # Prioritize common neotype column names
        if("NeoType" %in% meta_cols) {
          updateSelectInput(session, "trajectory_neotype", 
                            choices = neotype_choices, 
                            selected = "NeoType")
        } else if(any(grepl("neo|type", meta_cols, ignore.case = TRUE))) {
          likely_neotype <- meta_cols[grepl("neo|type", meta_cols, ignore.case = TRUE)][1]
          updateSelectInput(session, "trajectory_neotype", 
                            choices = neotype_choices, 
                            selected = likely_neotype)
        } else {
          updateSelectInput(session, "trajectory_neotype", 
                            choices = neotype_choices, 
                            selected = "none")
        }
      }
    })
    
    # Update cell type choices based on annotation column
    observeEvent(input$trajectory_annotation, {
      req(seurat_obj())
      req(input$trajectory_annotation)
      
      if(input$trajectory_annotation %in% colnames(seurat_obj()@meta.data)) {
        # Get unique values from the selected annotation column
        cell_identities <- unique(seurat_obj()@meta.data[[input$trajectory_annotation]])
        cell_identities <- cell_identities[!is.na(cell_identities)]
        
        # For Cell Ranger data, also include cluster identities
        if(input$trajectory_annotation %in% c("orig.ident", "tree.ident")) {
          # Add cluster identities as an option
          cluster_ids <- unique(as.character(Idents(seurat_obj())))
          cluster_ids <- paste0("Cluster_", cluster_ids)
          cell_identities <- c(cell_identities, cluster_ids)
        }
        
        # Sort for better UI
        cell_identities <- sort(cell_identities)
        
        updateSelectInput(session, "trajectory_celltypes", 
                          choices = cell_identities)
        
        # Update start/end cluster choices with the same identities
        start_end_choices <- c("Auto-detect" = "auto")
        for(ct in cell_identities) {
          start_end_choices[ct] <- ct
        }
        updateSelectInput(session, "trajectory_start", choices = start_end_choices)
        updateSelectInput(session, "trajectory_end", choices = start_end_choices)
        
        # Show helpful message for Cell Ranger data
        if(input$trajectory_annotation %in% c("orig.ident", "tree.ident")) {
          showNotification(
            paste("For Cell Ranger data, consider using 'seurat_clusters' as annotation column, or select individual cluster identities for trajectory analysis."),
            type = "message",
            duration = 8
          )
        }
      }
    })
    # Update NeoType choices
    observeEvent(input$trajectory_neotype, {
      req(seurat_obj())
      
      if(input$trajectory_neotype != "none" && 
         input$trajectory_neotype %in% colnames(seurat_obj()@meta.data)) {
        neotypes <- unique(seurat_obj()@meta.data[[input$trajectory_neotype]])
        neotypes <- neotypes[!is.na(neotypes)]
        neotypes <- sort(neotypes)
        
        updateSelectInput(session, "trajectory_neotypes", 
                          choices = neotypes,
                          selected = NULL)
        
        showNotification(
          paste("Found", length(neotypes), "unique values in", 
                input$trajectory_neotype, ":", 
                paste(neotypes, collapse = ", ")),
          type = "message",
          duration = 5
        )
      } else {
        updateSelectInput(session, "trajectory_neotypes", 
                          choices = NULL,
                          selected = NULL)
      }
    })
    
    # Show preview of selected filters
    output$trajectory_filter_preview <- renderText({
      req(seurat_obj())
      
      if(!is.null(input$trajectory_celltypes) && length(input$trajectory_celltypes) > 0) {
        
        # Handle cluster selections differently
        if(input$trajectory_annotation == "seurat_clusters" || 
           any(grepl("Cluster_", input$trajectory_celltypes))) {
          
          # For cluster-based selection
          if(input$trajectory_annotation == "seurat_clusters") {
            cluster_filter <- Idents(seurat_obj()) %in% input$trajectory_celltypes
          } else {
            # Extract cluster numbers from "Cluster_X" format
            cluster_nums <- gsub("Cluster_", "", input$trajectory_celltypes)
            cluster_filter <- as.character(Idents(seurat_obj())) %in% cluster_nums
          }
          
          if(input$trajectory_neotype != "none" && 
             !is.null(input$trajectory_neotypes) && 
             length(input$trajectory_neotypes) > 0) {
            neotype_filter <- seurat_obj()@meta.data[[input$trajectory_neotype]] %in% 
              input$trajectory_neotypes
            final_filter <- cluster_filter & neotype_filter
            
            paste("Selected filters will include", sum(final_filter), 
                  "cells out of", ncol(seurat_obj()), "total cells",
                  "\n- Clusters:", paste(input$trajectory_celltypes, collapse = ", "),
                  "\n- NeoTypes:", paste(input$trajectory_neotypes, collapse = ", "))
          } else {
            paste("Selected filters will include", sum(cluster_filter), 
                  "cells out of", ncol(seurat_obj()), "total cells",
                  "\n- Clusters:", paste(input$trajectory_celltypes, collapse = ", "))
          }
          
        } else {
          # Original logic for other annotation types
          annotation_filter <- seurat_obj()@meta.data[[input$trajectory_annotation]] %in% 
            input$trajectory_celltypes
          
          if(input$trajectory_neotype != "none" && 
             !is.null(input$trajectory_neotypes) && 
             length(input$trajectory_neotypes) > 0) {
            neotype_filter <- seurat_obj()@meta.data[[input$trajectory_neotype]] %in% 
              input$trajectory_neotypes
            final_filter <- annotation_filter & neotype_filter
            
            paste("Selected filters will include", sum(final_filter), 
                  "cells out of", ncol(seurat_obj()), "total cells",
                  "\n- Cell types:", paste(input$trajectory_celltypes, collapse = ", "),
                  "\n- NeoTypes:", paste(input$trajectory_neotypes, collapse = ", "))
          } else {
            paste("Selected filters will include", sum(annotation_filter), 
                  "cells out of", ncol(seurat_obj()), "total cells",
                  "\n- Cell types:", paste(input$trajectory_celltypes, collapse = ", "))
          }
        }
      } else {
        "No cell identities/types selected"
      }
    })
    
    # Trajectory calculation status
    output$trajectoryCalculated <- reactive({
      return(!is.null(local_values$trajectory))
    })
    outputOptions(output, "trajectoryCalculated", suspendWhenHidden = FALSE)
    
    # Run trajectory analysis
    observeEvent(input$run_trajectory, {
      req(processed())
      req(seurat_obj())
      req(input$trajectory_celltypes)
      req(length(input$trajectory_celltypes) >= 2)
      
      withProgress(message = 'Running trajectory analysis...', value = 0, {
        
        tryCatch({
          # Step 1: Filter data
          incProgress(0.1, detail = "Filtering data...")
          
          # Create filter condition
          if(input$trajectory_neotype != "none" && 
             !is.null(input$trajectory_neotypes) && 
             length(input$trajectory_neotypes) > 0) {
            # Filter by both annotation and NeoType
            subset_obj <- seurat_obj()[, 
                                       seurat_obj()@meta.data[[input$trajectory_annotation]] %in% 
                                         input$trajectory_celltypes &
                                         seurat_obj()@meta.data[[input$trajectory_neotype]] %in% 
                                         input$trajectory_neotypes
            ]
          } else {
            # Filter by annotation only
            subset_obj <- seurat_obj()[, 
                                       seurat_obj()@meta.data[[input$trajectory_annotation]] %in% 
                                         input$trajectory_celltypes
            ]
          }
          
          incProgress(0.2, detail = "Processing subset...")
          
          # Process subset
          subset_obj <- process_for_trajectory(subset_obj)
          #subset_obj <- NormalizeData(subset_obj)
          #subset_obj <- FindVariableFeatures(subset_obj,
          #                                   selection.method = "vst",
          #                                   nfeatures = 3000)
          
          # Scale all genes like in app_14.R  
          all.genes <- rownames(subset_obj)
          subset_obj <- ScaleData(subset_obj, features = all.genes)
          
          # Run PCA
          subset_obj <- RunPCA(subset_obj,
                               features = VariableFeatures(subset_obj))
          incProgress(0.4, detail = "Converting to SingleCellExperiment...")
          
          # Convert to SCE
          sce <- as.SingleCellExperiment(subset_obj)
          
          # Get clustering based on annotation
          clustering <- factor(subset_obj@meta.data[[input$trajectory_annotation]])
          
          # Get PCA dimensions
          dimred <- Embeddings(subset_obj, reduction = "pca")[, 1:input$trajectory_dims]
          
          incProgress(0.6, detail = "Running Slingshot...")
          
          # Run slingshot
          set.seed(1)
          
          # Prepare start and end clusters
          start_clus <- if(input$trajectory_start != "auto") input$trajectory_start else NULL
          end_clus <- if(input$trajectory_end != "auto") input$trajectory_end else NULL
          
          # Run lineage identification
          lineages <- getLineages(dimred, clustering, 
                                  start.clus = start_clus, 
                                  end.clus = end_clus)
          
          # Get curves
          crv <- getCurves(lineages)
          
          # Also run full slingshot on SCE
          slingshot_sce <- slingshot(sce, 
                                     clusterLabels = clustering,
                                     reducedDim = dimred,
                                     start.clus = start_clus,
                                     end.clus = end_clus)
          
          # Store results
          local_values$trajectory <- list(
            slingshot = crv,
            slingshot_sce = slingshot_sce,
            subset = subset_obj,
            pca = dimred,
            clustering = clustering,
            lineages = lineages
          )
          
          local_values$trajectory_subset <- subset_obj
          local_values$slingshot_obj <- crv
          local_values$slingshot_sce <- slingshot_sce
          
          # Update main values
          main_values$trajectory <- local_values$trajectory
          
          incProgress(0.8, detail = "Preparing visualization data...")
          
          # Add pseudotime to metadata
          if(!is.null(slingPseudotime(crv))) {
            pseudotime <- slingPseudotime(crv)
            for (i in 1:ncol(pseudotime)) {
              subset_obj[[paste0("slingPseudotime_", i)]] <- pseudotime[, i]
            }
            local_values$trajectory_subset <- subset_obj
          }
          
          incProgress(1, detail = "Complete!")
          
          showNotification("Trajectory analysis completed successfully!", 
                           type = "message", duration = 3)
          
        }, error = function(e) {
          showNotification(paste("Error in trajectory analysis:", e$message), 
                           type = "error", duration = 10)
          local_values$trajectory <- NULL
        })
      })
    })
    
    # Output: Main trajectory plot
    output$trajectoryPlot <- renderPlot({
      req(local_values$trajectory)
      
      safe_plot(
        function() plot_trajectory_main(local_values$trajectory),
        "Error creating trajectory plot"
      )
    })
    
    # Output: Stratified trajectory
    output$trajectorySeuratStyle <- renderPlot({
      req(local_values$trajectory)
      
      plot_trajectory_stratified(
        local_values$trajectory,
        neotype_col = if(input$trajectory_neotype != "none") input$trajectory_neotype else NULL,
        annotation_col = input$trajectory_annotation
      )
    })
    
    # Output: Pseudotime plot
    output$pseudotimePlot <- renderPlot({
      req(local_values$trajectory)
      
      safe_plot(
        function() plot_pseudotime(local_values$trajectory),
        "Error creating pseudotime plot"
      )
    })
    
    # Output: Trajectory statistics
    output$trajectoryStats <- renderPrint({
      req(local_values$trajectory)
      
      print_trajectory_stats(
        local_values$trajectory,
        start_clus = input$trajectory_start,
        end_clus = input$trajectory_end
      )
    })
    
    # Download trajectory data
    output$downloadTrajectory <- downloadHandler(
      filename = function() {
        paste0("trajectory_analysis_", Sys.Date(), ".rds")
      },
      content = function(file) {
        trajectory_export <- list(
          slingshot_object = local_values$trajectory$slingshot,
          subset_seurat = local_values$trajectory$subset,
          pca_embeddings = local_values$trajectory$pca,
          clustering = local_values$trajectory$clustering,
          lineages = local_values$trajectory$lineages,
          parameters = list(
            selected_celltypes = input$trajectory_celltypes,
            annotation_column = input$trajectory_annotation,
            neotype_column = input$trajectory_neotype,
            selected_neotypes = input$trajectory_neotypes,
            n_pcs = input$trajectory_dims,
            start_cluster = input$trajectory_start,
            end_cluster = input$trajectory_end
          )
        )
        saveRDS(trajectory_export, file)
      }
    )
    
  })
}
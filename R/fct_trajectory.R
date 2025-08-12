#' Process Seurat object for trajectory analysis
#'
#' @param seurat_obj Seurat object
#' @param nfeatures Number of variable features
#' @return Processed Seurat object
#' @noRd
process_for_trajectory <- function(seurat_obj, nfeatures = 3000) {
  # EXACT same steps as app_14.R - this ensures identical PCA space
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj,
                                     selection.method = "vst",
                                     nfeatures = nfeatures)
  
  # Scale all genes - EXACT same as app_14.R
  all.genes <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  
  # Run PCA - EXACT same as app_14.R
  seurat_obj <- RunPCA(seurat_obj,
                       features = VariableFeatures(seurat_obj))
  
  return(seurat_obj)
}

#' Plot main trajectory - Enhanced but compatible version
#'
#' @param trajectory Trajectory object from slingshot
#' @return A plot
#' @noRd
plot_trajectory_main <- function(trajectory) {
  # KEEP: Error handling for robustness
  if (is.null(trajectory) || is.null(trajectory$pca) || is.null(trajectory$clustering)) {
    plot.new()
    text(0.5, 0.5, "No trajectory data available", cex = 1.5, col = "red")
    return(invisible(NULL))
  }
  
  # EXACT same data extraction as app_14.R
  dimred <- trajectory$pca[, 1:2]
  clustering <- trajectory$clustering
  
  # ENHANCED: Better color handling for large cluster numbers
  n_clusters <- length(unique(clustering))
  if (n_clusters <= 9) {
    colors <- RColorBrewer::brewer.pal(min(n_clusters, 9), "Set1")
  } else {
    # Fallback for >9 clusters (this was missing in app_14.R)
    base_colors <- RColorBrewer::brewer.pal(9, "Set1")
    colors <- colorRampPalette(base_colors)(n_clusters)
  }
  
  # EXACT same plotting method as app_14.R (this is key for compatibility)
  plot(dimred[, 1:2], 
       col = colors[as.numeric(clustering)], 
       pch = 16,
       xlab = "PC_1",
       ylab = "PC_2",
       main = "Trajectory Analysis")
  
  # EXACT same cluster labels as app_14.R
  for (i in unique(clustering)) {
    text(mean(dimred[clustering == i, 1]), 
         mean(dimred[clustering == i, 2]), 
         labels = i, 
         font = 2)
  }
  
  # ENHANCED: Safe trajectory curve drawing (prevent crashes)
  tryCatch({
    lines(SlingshotDataSet(trajectory$slingshot), lwd = 3, col = "black")
  }, error = function(e) {
    # Silent fallback - don't break the entire plot
    warning("Could not draw trajectory curves: ", e$message)
  })
}

#' Plot stratified trajectory
#'
#' @param trajectory Trajectory object
#' @param neotype_col NeoType column name
#' @param annotation_col Annotation column name
#' @return A ggplot object
#' @noRd
plot_trajectory_stratified <- function(trajectory, neotype_col = NULL, 
                                       annotation_col) {
  
  tryCatch({
    # Check for valid trajectory data
    if (is.null(trajectory) || is.null(trajectory$subset)) {
      p <- ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                 label = "No trajectory data available", 
                 size = 6) +
        theme_void()
      return(p)
    }
    
    # Get curves for overlay
    curves <- slingCurves(trajectory$slingshot_sce, as.df = TRUE)
    
    # Create Seurat DimPlot
    p <- DimPlot(trajectory$subset,
                 label = FALSE,
                 pt.size = 2,
                 shape.by = neotype_col,
                 group.by = annotation_col,
                 dims = c(1, 2),
                 reduction = "pca")
    
    # Add the trajectory curves
    p <- p + 
      geom_path(data = curves %>% arrange(Order),
                aes(x = PC_1, y = PC_2, group = Lineage),
                color = "black",
                size = 2,
                inherit.aes = FALSE) +
      ggtitle("Stratified Trajectory") +
      theme_publication()
    
    return(p)
    
  }, error = function(e) {
    # Return error plot
    p <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = paste("Error:", e$message), 
               size = 4) +
      theme_void()
    return(p)
  })
}

#' Plot pseudotime - Enhanced but compatible version
#'
#' @param trajectory Trajectory object
#' @return A plot
#' @noRd
plot_pseudotime <- function(trajectory) {
  # KEEP: Input validation
  if (is.null(trajectory) || is.null(trajectory$slingshot)) {
    plot.new()
    text(0.5, 0.5, "No trajectory data available", cex = 1.5, col = "red")
    return(invisible(NULL))
  }
  
  pseudotime <- slingPseudotime(trajectory$slingshot)
  if (is.null(pseudotime)) {
    plot.new()
    text(0.5, 0.5, "No pseudotime calculated", cex = 1.5)
    return(invisible(NULL))
  }
  
  # EXACT same data processing as app_14.R
  dimred <- trajectory$pca[, 1:2]
  pt <- pseudotime[, 1]
  
  # ENHANCED: Better validation
  valid_cells <- !is.na(pt)
  if (sum(valid_cells) == 0) {
    plot.new()
    text(0.5, 0.5, "No valid pseudotime values", cex = 1.5, col = "orange")
    return(invisible(NULL))
  }
  
  pt <- pt[valid_cells]
  dimred <- dimred[valid_cells, ]
  
  # EXACT same color generation as app_14.R
  colors <- viridis::viridis(100)
  pt_colors <- colors[cut(pt, breaks = 100)]
  
  # EXACT same plotting method as app_14.R (key for compatibility)
  plot(dimred[, 1:2], 
       col = pt_colors, 
       pch = 16,
       xlab = "PC_1",
       ylab = "PC_2",
       main = "Pseudotime (Lineage 1)")
  
  # ENHANCED: Safer legend positioning with app_14.R method as primary
  tryCatch({
    # Try app_14.R method first
    legend_image <- as.raster(matrix(colors, ncol = 1))
    rasterImage(legend_image, 
                xleft = max(dimred[, 1]) * 0.85,
                ybottom = min(dimred[, 2]),
                xright = max(dimred[, 1]) * 0.9,
                ytop = max(dimred[, 2]))
    text(x = max(dimred[, 1]) * 0.92, 
         y = max(dimred[, 2]), 
         labels = "Late", 
         adj = 0)
    text(x = max(dimred[, 1]) * 0.92, 
         y = min(dimred[, 2]), 
         labels = "Early", 
         adj = 0)
  }, error = function(e) {
    # Fallback: simpler legend if positioning fails
    legend("topright", legend = c("Early", "Late"), 
           fill = c(colors[1], colors[100]), cex = 0.8)
  })
}

#' Print trajectory statistics
#'
#' @param trajectory Trajectory object
#' @param start_clus Starting cluster
#' @param end_clus Ending cluster
#' @noRd
print_trajectory_stats <- function(trajectory, start_clus = NULL, end_clus = NULL) {
  if (is.null(trajectory)) {
    cat("No trajectory analysis has been performed yet.\n")
    return(invisible(NULL))
  }
  
  cat("=== Trajectory Analysis Summary ===\n\n")
  
  if (!is.null(trajectory$subset)) {
    cat("Number of cells:", ncol(trajectory$subset), "\n")
  }
  
  if (!is.null(trajectory$clustering)) {
    cat("Number of cell types:", length(unique(trajectory$clustering)), "\n")
    cat("Cell types included:", paste(sort(unique(trajectory$clustering)), 
                                      collapse = ", "), "\n")
  }
  
  if(start_clus != "auto" && !is.null(start_clus)) {
    cat("Starting cluster:", start_clus, "\n")
  }
  if(end_clus != "auto" && !is.null(end_clus)) {
    cat("Ending cluster:", end_clus, "\n")
  }
  
  # Lineage information
  if (!is.null(trajectory$slingshot)) {
    tryCatch({
      lineages <- slingLineages(trajectory$slingshot)
      cat("\nNumber of lineages:", length(lineages), "\n")
      
      for (i in seq_along(lineages)) {
        cat("\nLineage", i, ":", paste(lineages[[i]], collapse = " -> "), "\n")
      }
      
      # Pseudotime range
      pseudotime <- slingPseudotime(trajectory$slingshot)
      if (!is.null(pseudotime)) {
        cat("\nPseudotime ranges:\n")
        for (i in 1:ncol(pseudotime)) {
          pt <- pseudotime[, i]
          pt <- pt[!is.na(pt)]
          if (length(pt) > 0) {
            cat("  Lineage", i, ": ", 
                round(min(pt), 2), " - ", 
                round(max(pt), 2), "\n")
          }
        }
      }
    }, error = function(e) {
      cat("\nError retrieving trajectory details: ", e$message, "\n")
    })
  }
}
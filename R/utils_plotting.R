#' Create a consistent theme for plots
#'
#' @return A ggplot2 theme object
#' @noRd
theme_publication <- function() {
  theme_minimal() +
    theme(
      text = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}

#' Create a consistent color palette
#'
#' @param n Number of colors needed
#' @param palette Type of palette ("discrete", "continuous", "diverging")
#' @return Vector of colors
#' @noRd
get_color_palette <- function(n, palette = "discrete") {
  if (palette == "discrete") {
    if (n <= 9) {
      RColorBrewer::brewer.pal(n, "Set1")
    } else if (n <= 12) {
      RColorBrewer::brewer.pal(n, "Set3")
    } else {
      colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n)
    }
  } else if (palette == "continuous") {
    viridis::viridis(n)
  } else if (palette == "diverging") {
    colorRampPalette(c("blue", "white", "red"))(n)
  }
}

#' Create a violin plot with consistent styling
#'
#' @param seurat_obj Seurat object
#' @param features Features to plot
#' @param group.by Grouping variable
#' @param pt.size Point size
#' @return A ggplot object
#' @noRd
create_violin_plot <- function(seurat_obj, features, group.by = NULL, pt.size = 0) {
  VlnPlot(seurat_obj, 
          features = features,
          group.by = group.by,
          pt.size = pt.size) + 
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Create a feature plot with consistent styling
#'
#' @param seurat_obj Seurat object
#' @param features Features to plot
#' @param reduction Reduction to use
#' @return A ggplot object
#' @noRd
create_feature_plot <- function(seurat_obj, features, reduction = "umap") {
  FeaturePlot(seurat_obj, 
              features = features,
              reduction = reduction) + 
    scale_color_viridis() +
    theme_publication()
}

#' Create a DimPlot with consistent styling
#'
#' @param seurat_obj Seurat object
#' @param reduction Reduction to use
#' @param group.by Grouping variable
#' @param label Whether to label clusters
#' @return A ggplot object
#' @noRd
create_dim_plot <- function(seurat_obj, reduction = "umap", 
                            group.by = NULL, label = TRUE) {
  DimPlot(seurat_obj, 
          reduction = reduction,
          group.by = group.by,
          label = label,
          label.size = 4) + 
    theme_publication()
}

#' Calculate and format cluster statistics
#'
#' @param seurat_obj Seurat object
#' @return A data frame with cluster statistics
#' @noRd
calculate_cluster_stats <- function(seurat_obj) {
  cluster_counts <- table(Idents(seurat_obj))
  cluster_props <- prop.table(cluster_counts) * 100
  
  data.frame(
    Cluster = names(cluster_counts),
    Cell_Count = as.numeric(cluster_counts),
    Percentage = round(as.numeric(cluster_props), 2),
    stringsAsFactors = FALSE
  )
}

#' Create a heatmap with consistent styling
#'
#' @param data Matrix of data to plot
#' @param scale Whether to scale the data
#' @param color_palette Color palette to use
#' @return A heatmap plot
#' @noRd
create_heatmap <- function(data, scale = TRUE, 
                           color_palette = NULL) {
  if (is.null(color_palette)) {
    color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  }
  
  if (scale) {
    data <- t(scale(t(data)))
    data[is.na(data)] <- 0
  }
  
  heatmap(data, 
          col = color_palette,
          scale = "none",
          margins = c(10, 8))
}

#' Create empty plot with message
#'
#' @param message Message to display
#' @param cex Text size
#' @noRd
create_empty_plot <- function(message = "No data available", cex = 1.2) {
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.5, 0.5, message, cex = cex, col = "gray40")
}

#' Safe plot wrapper
#'
#' @param plot_func Function that creates the plot
#' @param error_message Message to show on error
#' @noRd
safe_plot <- function(plot_func, error_message = "Error creating plot") {
  tryCatch({
    plot_func()
  }, error = function(e) {
    create_empty_plot(paste(error_message, "\n", e$message))
  })
}

#' Calculate pathway averages by cluster (for export functions)
#'
#' @param seurat_obj Seurat object
#' @param score_list Named list of score column names
#' @return Data frame with average scores per cluster
#' @noRd
calculate_pathway_averages <- function(seurat_obj, score_list) {
  clusters <- unique(Idents(seurat_obj))
  
  # Initialize results data frame
  results <- data.frame(Cluster = clusters)
  
  for(pathway_name in names(score_list)) {
    score_col <- score_list[[pathway_name]]
    
    # Calculate average scores per cluster
    avg_scores <- aggregate(
      seurat_obj@meta.data[[score_col]],
      by = list(Idents(seurat_obj)),
      FUN = mean,
      na.rm = TRUE
    )
    
    # Match with results data frame
    results[[pathway_name]] <- avg_scores$x[match(results$Cluster, avg_scores$Group.1)]
  }
  
  return(results)
}
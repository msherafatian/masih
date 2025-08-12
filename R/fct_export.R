#' Create export plot based on type
#'
#' @param plot_type Type of plot to create
#' @param seurat_obj Seurat object
#' @param main_values Main reactive values
#' @param gene Gene name for gene plots
#' @param custom_title Custom title
#' @param remove_axes Remove axes
#' @param white_bg White background
#' @return Plot object or function
#' @noRd
#' @import dplyr
#' @import tidyr
#' @import corrplot
#' @import viridis

create_export_plot <- function(plot_type, seurat_obj, main_values, 
                               gene = NULL, custom_title = "", 
                               remove_axes = FALSE, white_bg = TRUE) {
  
  p <- NULL
  
  # Create plot based on type
  switch(plot_type,
         "cluster" = {
           p <- DimPlot(seurat_obj, 
                        reduction = "umap",
                        label = TRUE,
                        label.size = 4) + 
             theme_minimal() +
             theme(legend.position = "right")
         },
         
         "cellcycle" = {
           if ("Phase" %in% colnames(seurat_obj@meta.data)) {
             p <- DimPlot(seurat_obj, 
                          reduction = "umap",
                          group.by = "Phase") + 
               theme_minimal() +
               scale_color_manual(values = c("G1" = "#E41A1C", 
                                             "G2M" = "#377EB8", 
                                             "S" = "#4DAF4A"))
           }
         },
         
         "cellcycle_distribution" = {
           if ("Phase" %in% colnames(seurat_obj@meta.data)) {
             phase_data <- data.frame(
               Cluster = Idents(seurat_obj),
               Phase = seurat_obj$Phase
             )
             
             phase_props <- phase_data %>%
               group_by(Cluster, Phase) %>%
               summarise(count = n(), .groups = "drop") %>%
               group_by(Cluster) %>%
               mutate(prop = count / sum(count))
             
             p <- ggplot(phase_props, aes(x = Cluster, y = prop, fill = Phase)) +
               geom_bar(stat = "identity") +
               theme_minimal() +
               labs(y = "Proportion", 
                    title = "Cell Cycle Phase Distribution by Cluster") +
               scale_fill_manual(values = c("G1" = "#E41A1C", 
                                            "G2M" = "#377EB8", 
                                            "S" = "#4DAF4A")) +
               theme(axis.text.x = element_text(angle = 45, hjust = 1))
           }
         },
         
         "tree" = {
           # Return function for base R plot
           return(function() PlotClusterTree(seurat_obj))
         },
         
         "gene" = {
           if (!is.null(gene) && gene != "none") {
             p <- FeaturePlot(seurat_obj,
                              features = gene,
                              reduction = "umap") +
               scale_color_viridis() +
               theme_minimal()
           }
         },
         
         "cancersea_feature" = {
           if (length(main_values$cancersea_scores) > 0) {
             # Use first pathway as example
             score_col <- main_values$cancersea_scores[[1]]
             pathway_name <- names(main_values$cancersea_scores)[1]
             
             p <- FeaturePlot(seurat_obj,
                              features = score_col,
                              reduction = "umap") +
               scale_color_viridis() +
               theme_minimal() +
               ggtitle(paste(pathway_name, "Score"))
           }
         },
         
         "cancersea_heatmap" = {
           if (length(main_values$cancersea_scores) >= 2) {
             # Create a matrix of average scores per cluster
             score_matrix <- matrix(
               nrow = length(unique(Idents(seurat_obj))), 
               ncol = length(main_values$cancersea_scores)
             )
             
             for (i in 1:length(main_values$cancersea_scores)) {
               score_col <- main_values$cancersea_scores[[i]]
               avg_scores <- aggregate(
                 seurat_obj@meta.data[[score_col]], 
                 by = list(Idents(seurat_obj)), 
                 FUN = mean
               )
               score_matrix[, i] <- avg_scores$x
             }
             
             rownames(score_matrix) <- paste("Cluster", levels(Idents(seurat_obj)))
             colnames(score_matrix) <- names(main_values$cancersea_scores)
             
             # Scale and plot
             score_matrix_scaled <- t(scale(t(score_matrix)))
             score_matrix_scaled[is.na(score_matrix_scaled)] <- 0
             
             # Return function for base R heatmap
             return(function() {
               heatmap(score_matrix_scaled, 
                       col = viridis::viridis(100),
                       scale = "none",
                       margins = c(10, 8),
                       main = "CancerSEA Pathway Scores by Cluster")
             })
           }
         },
         
         "correlation_matrix" = {
           if (length(main_values$cancersea_scores) >= 2) {
             # Extract scores for correlation
             score_cols <- unlist(main_values$cancersea_scores)
             score_data <- seurat_obj@meta.data[, score_cols, drop = FALSE]
             colnames(score_data) <- names(main_values$cancersea_scores)
             
             # Calculate correlation
             cor_matrix <- cor(score_data, use = "complete.obs")
             
             # Return function for corrplot
             return(function() {
               corrplot::corrplot(cor_matrix, 
                                  method = "color",
                                  type = "upper",
                                  order = "hclust",
                                  tl.cex = 0.8,
                                  tl.col = "black",
                                  col = colorRampPalette(c("blue", "white", "red"))(100),
                                  addCoef.col = "black",
                                  number.cex = 0.7,
                                  title = "Pathway Correlation Matrix")
             })
           }
         },
         
         "pathway_by_cluster" = {
           if (length(main_values$cancersea_scores) >= 2) {
             # Prepare data
             plot_data <- data.frame(Cluster = Idents(seurat_obj))
             
             for (pathway in names(main_values$cancersea_scores)) {
               score_col <- main_values$cancersea_scores[[pathway]]
               plot_data[[pathway]] <- seurat_obj@meta.data[[score_col]]
             }
             
             # Reshape for plotting
             plot_data_long <- pivot_longer(plot_data, 
                                            cols = -Cluster, 
                                            names_to = "Pathway", 
                                            values_to = "Score")
             
             # Calculate mean scores
             plot_summary <- plot_data_long %>%
               group_by(Cluster, Pathway) %>%
               summarise(Mean_Score = mean(Score, na.rm = TRUE),
                         SE = sd(Score, na.rm = TRUE) / sqrt(n()),
                         .groups = "drop")
             
             # Plot
             p <- ggplot(plot_summary, aes(x = Cluster, y = Mean_Score, fill = Pathway)) +
               geom_bar(stat = "identity", position = "dodge") +
               geom_errorbar(aes(ymin = Mean_Score - SE, ymax = Mean_Score + SE),
                             position = position_dodge(0.9), width = 0.25) +
               theme_minimal() +
               labs(y = "Mean Score", title = "Pathway Scores by Cluster") +
               theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
               scale_fill_brewer(palette = "Set2")
           }
         },
         
         "comparison_heatmap" = {
           if (length(main_values$cancersea_scores) >= 2) {
             # Create matrix of average scores per cluster
             score_matrix <- calculate_pathway_averages_export(
               seurat_obj, 
               main_values$cancersea_scores
             )
             
             # Convert to matrix format
             score_mat <- as.matrix(score_matrix[, -1])
             rownames(score_mat) <- score_matrix$Cluster
             
             # Return function for base R heatmap
             return(function() {
               heatmap(score_mat, 
                       col = viridis::viridis(100),
                       scale = "row",
                       margins = c(10, 8),
                       main = "Pathway Comparison Heatmap")
             })
           }
         },
         
         "marker_heatmap" = {
           if (!is.null(main_values$markers)) {
             # Get top markers
             top_markers <- main_values$markers %>%
               group_by(cluster) %>%
               slice_max(avg_log2FC, n = 5) %>%
               pull(gene) %>%
               unique()
             
             if (length(top_markers) > 30) {
               top_markers <- top_markers[1:30]
             }
             
             p <- DoHeatmap(seurat_obj, 
                            features = top_markers,
                            size = 3) + 
               theme(axis.text.y = element_text(size = 8))
           }
         },
         
         "trajectory_main" = {
           if (!is.null(main_values$trajectory)) {
             return(function() {
               tryCatch({
                 plot_trajectory_main(main_values$trajectory)
               }, error = function(e) {
                 plot.new()
                 text(0.5, 0.5, paste("Error creating trajectory plot:", e$message), cex = 1)
               })
             })
           }
         },
         "trajectory_stratified" = {
           if (!is.null(main_values$trajectory)) {
             return(function() {
               tryCatch({
                 plot_trajectory_stratified(
                   main_values$trajectory,
                   neotype_col = NULL,  # You might want to make this configurable
                   annotation_col = "seurat_clusters"  # Or get this from trajectory settings
                 )
               }, error = function(e) {
                 plot.new()
                 text(0.5, 0.5, paste("Error creating stratified trajectory plot:", e$message), cex = 1)
               })
             })
           }
         },
         "pseudotime" = {
           if (!is.null(main_values$trajectory)) {
             return(function() {
               tryCatch({
                 plot_pseudotime(main_values$trajectory)
               }, error = function(e) {
                 plot.new()
                 text(0.5, 0.5, paste("Error creating pseudotime plot:", e$message), cex = 1)
               })
             })
           }
         }
  )
  
  # Apply custom settings if ggplot
  if (!is.null(p) && inherits(p, "ggplot")) {
    if (nchar(custom_title) > 0) {
      p <- p + ggtitle(custom_title)
    }
    
    if (remove_axes) {
      p <- p + theme(axis.text = element_blank(),
                     axis.title = element_blank(),
                     axis.ticks = element_blank())
    }
    
    if (white_bg) {
      p <- p + theme(panel.background = element_rect(fill = "white"),
                     plot.background = element_rect(fill = "white"))
    }
  }
  
  return(p)
}

#' Helper function for export - calculate pathway averages
#'
#' @param seurat_obj Seurat object
#' @param score_list Named list of score column names
#' @return Data frame with average scores per cluster
#' @noRd
calculate_pathway_averages_export <- function(seurat_obj, score_list) {
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

#' Get available plot types based on analyses performed - UPDATED
#'
#' @param main_values Main reactive values
#' @return Vector of plot types
#' @noRd
get_available_plot_types <- function(main_values) {
  plot_types <- c("cluster")
  
  if ("Phase" %in% colnames(main_values$seurat_obj@meta.data)) {
    plot_types <- c(plot_types, "cellcycle", "cellcycle_distribution")
  }
  
  if (length(main_values$cancersea_scores) > 0) {
    plot_types <- c(plot_types, "cancersea_feature")
    
    # Only add comparison plots if we have multiple pathways
    if (length(main_values$cancersea_scores) >= 2) {
      plot_types <- c(plot_types, "cancersea_heatmap", "correlation_matrix", 
                      "pathway_by_cluster", "comparison_heatmap")
    }
  }
  
  if (!is.null(main_values$markers)) {
    plot_types <- c(plot_types, "marker_heatmap")
  }
  
  if (!is.null(main_values$trajectory)) {
    plot_types <- c(plot_types, "trajectory_main", "trajectory_stratified", "pseudotime") 
  }
  
  plot_types <- c(plot_types, "tree")
  
  return(plot_types)
}

#' Export base R plot - UPDATED
#'
#' @param plot_func Function that creates the plot
#' @param file File path
#' @param format Export format
#' @param width Width in inches
#' @param height Height in inches
#' @param dpi Resolution
#' @noRd
export_base_plot <- function(plot_func, file, format, width, height, dpi) {
  # Open device
  if (format == "png") {
    png(file, width = width, height = height, units = "in", res = dpi)
  } else if (format == "pdf") {
    pdf(file, width = width, height = height)
  } else if (format == "svg") {
    svg(file, width = width, height = height)
  } else if (format == "tiff") {
    tiff(file, width = width, height = height, units = "in", res = dpi)
  }
  
  # Try to create plot with error handling
  tryCatch({
    plot_func()
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Error:", e$message), cex = 1)
  })
  
  # Close device
  dev.off()
}

#' Calculate cell cycle proportions
#'
#' @param seurat_obj Seurat object
#' @return Data frame with proportions
#' @noRd
calculate_cellcycle_proportions <- function(seurat_obj) {
  phase_data <- data.frame(
    Cluster = Idents(seurat_obj),
    Phase = seurat_obj$Phase
  )
  
  phase_table <- table(phase_data$Cluster, phase_data$Phase)
  phase_props <- prop.table(phase_table, margin = 1) * 100
  phase_df <- as.data.frame.matrix(phase_props)
  phase_df$Cluster <- rownames(phase_df)
  phase_df <- phase_df[, c("Cluster", "G1", "S", "G2M")]
  
  return(phase_df)
}

#' Summarize trajectory results
#'
#' @param trajectory Trajectory object
#' @return Data frame summary
#' @noRd
summarize_trajectory <- function(trajectory) {
  lineages <- slingLineages(trajectory$slingshot)
  
  summary_list <- list()
  
  for (i in seq_along(lineages)) {
    summary_list[[i]] <- data.frame(
      Lineage = i,
      Path = paste(lineages[[i]], collapse = " -> "),
      n_clusters = length(lineages[[i]]),
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, summary_list)
}

#' Generate methods text
#'
#' @param main_values Main reactive values
#' @return Methods text
#' @noRd
generate_methods_text <- function(main_values) {
  obj <- main_values$seurat_obj
  
  methods <- paste0(
    "Single-cell RNA sequencing data was analyzed using the MASIH pipeline ",
    "(Modular Analysis of Single-cell Insights and Heterogeneity), ",
    "built on Seurat v4 (Hao et al., 2021). ",
    "Data comprised ", ncol(obj), " cells and ", nrow(obj), " genes. ",
    "Cells were normalized using log-normalization and ", 
    length(VariableFeatures(obj)), " highly variable features were identified. ",
    "Principal component analysis was performed, followed by UMAP embedding. ",
    "Clustering was performed using a shared nearest neighbor (SNN) ",
    "modularity optimization-based algorithm, identifying ", 
    length(unique(Idents(obj))), " distinct clusters. "
  )
  
  if ("Phase" %in% colnames(obj@meta.data)) {
    methods <- paste0(methods,
                      "Cell cycle phase was determined using established S and G2/M phase markers ",
                      "(Tirosh et al., 2016). "
    )
  }
  
  if (length(main_values$cancersea_scores) > 0) {
    methods <- paste0(methods,
                      "Cancer cell functional states were assessed using CancerSEA gene signatures ",
                      "(Yuan et al., 2019) for ", length(main_values$cancersea_scores), 
                      " pathways. "
    )
  }
  
  if (!is.null(main_values$markers)) {
    methods <- paste0(methods,
                      "Differential gene expression analysis was performed to identify ",
                      "cluster-specific marker genes. "
    )
  }
  
  if (!is.null(main_values$trajectory)) {
    methods <- paste0(methods,
                      "Trajectory analysis was performed using Slingshot ",
                      "(Street et al., 2018) on selected cell populations. "
    )
  }
  
  methods <- paste0(methods,
                    "All visualizations and analyses were performed using the MASIH application."
  )
  
  return(methods)
}

#' Generate citations
#'
#' @return Citation text
#' @noRd
generate_citations <- function() {
  citations <- paste0(
    "MASIH Pipeline:\n",
    "Custom pipeline for comprehensive single-cell analysis ",
    "(manuscript in preparation).\n\n",
    
    "Seurat v4:\n",
    "Hao, Y. et al. Integrated analysis of multimodal single-cell data. ",
    "Cell 184, 3573-3587 (2021).\n\n",
    
    "CancerSEA:\n",
    "Yuan, H. et al. CancerSEA: a cancer single-cell state atlas. ",
    "Nucleic Acids Res 47, D900-D908 (2019).\n\n",
    
    "Cell Cycle Scoring:\n",
    "Tirosh, I. et al. Dissecting the multicellular ecosystem of ",
    "metastatic melanoma by single-cell RNA-seq. ",
    "Science 352, 189-196 (2016).\n\n",
    
    "Slingshot:\n",
    "Street, K. et al. Slingshot: cell lineage and pseudotime inference ",
    "for single-cell transcriptomics. BMC Genomics 19, 477 (2018).\n\n",
    
    "R packages:\n",
    "R Core Team (2023). R: A language and environment for ",
    "statistical computing.\n",
    "Wickham, H. ggplot2: Elegant Graphics for Data Analysis. ",
    "Springer-Verlag New York (2016)."
  )
  
  return(citations)
}

#' Generate figure legend template
#'
#' @param seurat_obj Seurat object
#' @return Legend text
#' @noRd
generate_figure_legend <- function(seurat_obj) {
  legend <- paste0(
    "UMAP visualization of ", ncol(seurat_obj), " single cells ",
    "colored by [cluster identity/functional state/gene expression]. ",
    "Clustering identified ", length(unique(Idents(seurat_obj))), 
    " distinct populations. ",
    "Color scale represents [description of what colors mean]. ",
    "Each point represents an individual cell. ",
    "Data was processed using MASIH pipeline with log-normalization ",
    "and the first [N] principal components were used for ",
    "dimensionality reduction."
  )
  
  return(legend)
}

#' Print analysis summary
#'
#' @param main_values Main reactive values
#' @noRd
print_analysis_summary <- function(main_values) {
  obj <- main_values$seurat_obj
  
  cat("=== MASIH Analysis Summary ===\n\n")
  cat("Date:", Sys.Date(), "\n")
  cat("Total cells analyzed:", ncol(obj), "\n")
  cat("Total genes detected:", nrow(obj), "\n")
  cat("Number of clusters:", length(unique(Idents(obj))), "\n")
  cat("Variable features:", length(VariableFeatures(obj)), "\n")
  
  reductions <- names(obj@reductions)
  cat("Reductions performed:", paste(reductions, collapse = ", "), "\n")
  
  assays <- names(obj@assays)
  cat("Assays available:", paste(assays, collapse = ", "), "\n")
  
  if (length(main_values$cancersea_scores) > 0) {
    cat("\nCancerSEA pathways analyzed:", length(main_values$cancersea_scores), "\n")
    cat("- ", paste(names(main_values$cancersea_scores), collapse = "\n- "), "\n")
  }
  
  if (!is.null(main_values$markers)) {
    cat("\nMarker gene analysis: Completed\n")
    cat("- Total markers found:", nrow(main_values$markers), "\n")
  }
  
  if (!is.null(main_values$trajectory)) {
    cat("\nTrajectory analysis: Completed\n")
    lineages <- slingLineages(main_values$trajectory$slingshot)
    cat("- Number of lineages:", length(lineages), "\n")
  }
  
  cat("\n=== End of Summary ===\n")
}
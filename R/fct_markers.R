#' Find conserved markers across conditions
#'
#' @param seurat_obj Seurat object
#' @param grouping_var Variable to group by
#' @param ident_var Identity variable
#' @param test_use Statistical test to use
#' @return Data frame of conserved markers
#' @noRd
find_conserved_markers <- function(seurat_obj, grouping_var, 
                                   ident_var = NULL, test_use = "wilcox") {
  
  if (!is.null(ident_var)) {
    Idents(seurat_obj) <- ident_var
  }
  
  markers_list <- list()
  
  for (cluster in levels(Idents(seurat_obj))) {
    tryCatch({
      markers <- FindConservedMarkers(
        seurat_obj,
        ident.1 = cluster,
        grouping.var = grouping_var,
        test.use = test_use,
        verbose = FALSE
      )
      markers$cluster <- cluster
      markers$gene <- rownames(markers)
      markers_list[[cluster]] <- markers
    }, error = function(e) {
      message("Error finding conserved markers for cluster ", 
              cluster, ": ", e$message)
    })
  }
  
  do.call(rbind, markers_list)
}

#' Create a summary of marker statistics
#'
#' @param markers Data frame of markers
#' @return Data frame with summary statistics
#' @noRd
summarize_markers <- function(markers) {
  markers %>%
    group_by(cluster) %>%
    summarise(
      n_markers = n(),
      avg_log2FC = mean(avg_log2FC),
      median_log2FC = median(avg_log2FC),
      max_log2FC = max(avg_log2FC),
      avg_pct1 = mean(pct.1),
      avg_pct2 = mean(pct.2),
      .groups = "drop"
    )
}

#' Filter markers by multiple criteria
#'
#' @param markers Data frame of markers
#' @param min_logfc Minimum log fold change
#' @param min_pct Minimum percent expression
#' @param max_pval Maximum p-value
#' @param only_pos Only positive markers
#' @return Filtered data frame
#' @noRd
filter_markers <- function(markers, min_logfc = 0.25, min_pct = 0.1, 
                           max_pval = 0.05, only_pos = TRUE) {
  
  filtered <- markers
  
  if ("avg_log2FC" %in% colnames(filtered)) {
    if (only_pos) {
      filtered <- filtered %>% filter(avg_log2FC > min_logfc)
    } else {
      filtered <- filtered %>% filter(abs(avg_log2FC) > min_logfc)
    }
  }
  
  if ("pct.1" %in% colnames(filtered)) {
    filtered <- filtered %>% filter(pct.1 > min_pct)
  }
  
  if ("p_val_adj" %in% colnames(filtered)) {
    filtered <- filtered %>% filter(p_val_adj < max_pval)
  }
  
  return(filtered)
}

#' Create marker gene report
#'
#' @param markers Data frame of markers
#' @param seurat_obj Seurat object
#' @param output_dir Output directory
#' @return Path to report file
#' @noRd
create_marker_report <- function(markers, seurat_obj, 
                                 output_dir = tempdir()) {
  
  report_file <- file.path(output_dir, 
                           paste0("marker_report_", Sys.Date(), ".html"))
  
  # This would create an HTML report with marker statistics,
  # plots, and gene lists. Implementation depends on your
  # preferred reporting package (rmarkdown, etc.)
  
  message("Marker report would be created at: ", report_file)
  return(report_file)
}

#' Get marker genes for gene set enrichment
#'
#' @param markers Data frame of markers
#' @param cluster Cluster to get markers for
#' @param n_genes Number of genes to return
#' @param direction Direction of change ("up", "down", "both")
#' @return Vector of gene names
#' @noRd
get_genes_for_enrichment <- function(markers, cluster, 
                                     n_genes = 100, 
                                     direction = "up") {
  
  cluster_markers <- markers %>%
    filter(cluster == !!cluster)
  
  if (direction == "up") {
    genes <- cluster_markers %>%
      filter(avg_log2FC > 0) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = n_genes) %>%
      pull(gene)
  } else if (direction == "down") {
    genes <- cluster_markers %>%
      filter(avg_log2FC < 0) %>%
      arrange(avg_log2FC) %>%
      slice_head(n = n_genes) %>%
      pull(gene)
  } else {
    genes <- cluster_markers %>%
      arrange(desc(abs(avg_log2FC))) %>%
      slice_head(n = n_genes) %>%
      pull(gene)
  }
  
  return(genes)
}
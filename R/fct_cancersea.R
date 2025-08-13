#' Calculate multiple CancerSEA pathways at once
#'
#' @param seurat_obj Seurat object
#' @param pathways Vector of pathway names
#' @param ctrl Number of control features
#' @return Updated Seurat object with pathway scores
#' @noRd
calculate_multiple_pathways <- function(seurat_obj, pathways, ctrl = 20) {
  
  withProgress(message = 'Calculating pathway scores...', 
               value = 0, {
                 
                 all_genes <- rownames(seurat_obj)
                 scores_added <- character()
                 
                 for (i in seq_along(pathways)) {
                   pathway <- pathways[i]
                   
                   incProgress(1/length(pathways), 
                               detail = paste("Processing", pathway))
                   
                   tryCatch({
                     if (!requireNamespace("cancersea", quietly = TRUE)) {
                       message("cancersea package required for pathway: ", pathway_name)
                       return(character())
                     }
                     pathway_data <- get(pathway_name, envir = asNamespace("cancersea"))
                     gene_list <- pathway_data$symbol
                     gene_list_filtered <- gene_list[gene_list %in% all_genes]
                     
                     if (length(gene_list_filtered) > 0) {
                       # FIXED: Use better parameters for AddModuleScore
                       seurat_obj <- AddModuleScore(
                         object = seurat_obj,
                         features = list(gene_list_filtered),
                         name = paste0(pathway, "_"),
                         ctrl = min(100, length(gene_list_filtered) * 5),  # More control genes
                         nbin = 24,  # More bins for better control matching
                         seed = 1,   # Reproducible results
                         assay = "RNA",  # Explicitly use RNA assay
                         search = TRUE   # Search for genes more thoroughly
                       )
                       scores_added <- c(scores_added, pathway)
                     }
                   }, error = function(e) {
                     message("Error calculating ", pathway, ": ", e$message)
                   })
                 }
                 
                 message("Successfully calculated ", length(scores_added), 
                         " pathway scores")
               })
  
  return(list(
    seurat_obj = seurat_obj,
    scores_added = scores_added
  ))
}

#' Get genes for a CancerSEA pathway
#'
#' @param pathway_name Name of the pathway
#' @return Vector of gene symbols
#' @noRd
get_pathway_genes <- function(pathway_name) {
  if (!requireNamespace("cancersea", quietly = TRUE)) {
    message("cancersea package required for pathway: ", pathway_name)
    return(character())
  }
  
  tryCatch({
    pathway_env <- asNamespace("cancersea")
    if (exists(pathway_name, envir = pathway_env)) {
      pathway_data <- get(pathway_name, envir = pathway_env)
      if (is.list(pathway_data) && "symbol" %in% names(pathway_data)) {
        return(pathway_data$symbol)
      }
    }
    character()
  }, error = function(e) {
    message("Error getting genes for ", pathway_name, ": ", e$message)
    character()
  })
}

#' Create pathway score correlation matrix
#'
#' @param seurat_obj Seurat object
#' @param score_columns Named vector of score column names
#' @return Correlation matrix
#' @noRd
create_pathway_correlation <- function(seurat_obj, score_columns) {
  if (length(score_columns) < 2) {
    stop("Need at least 2 pathways for correlation")
  }
  
  score_data <- seurat_obj@meta.data[, score_columns, drop = FALSE]
  colnames(score_data) <- names(score_columns)
  
  cor(score_data, use = "complete.obs")
}

#' Calculate average pathway scores by cluster
#'
#' @param seurat_obj Seurat object
#' @param score_columns Named vector of score column names
#' @return Data frame with average scores
#' @noRd
calculate_pathway_averages <- function(seurat_obj, score_columns) {
  clusters <- levels(Idents(seurat_obj))
  
  avg_scores <- data.frame(
    Cluster = clusters,
    stringsAsFactors = FALSE
  )
  
  for (pathway in names(score_columns)) {
    score_col <- score_columns[[pathway]]
    avg_by_cluster <- aggregate(
      seurat_obj@meta.data[[score_col]], 
      by = list(Idents(seurat_obj)), 
      FUN = mean,
      na.rm = TRUE
    )
    avg_scores[[pathway]] <- avg_by_cluster$x
  }
  
  return(avg_scores)
}
#' Check if Seurat object has required analyses
#'
#' @param seurat_obj A Seurat object
#' @return A list of boolean values indicating which analyses exist
#' @noRd
check_seurat_analyses <- function(seurat_obj) {
  list(
    normalized = !is.null(seurat_obj@assays$RNA@data) && 
      sum(seurat_obj@assays$RNA@data) > 0,
    variable_features = length(VariableFeatures(seurat_obj)) > 0,
    scaled = !is.null(seurat_obj@assays$RNA@scale.data) && 
      nrow(seurat_obj@assays$RNA@scale.data) > 0,
    sct = "SCT" %in% names(seurat_obj@assays),
    pca = "pca" %in% names(seurat_obj@reductions),
    umap = "umap" %in% names(seurat_obj@reductions),
    tsne = "tsne" %in% names(seurat_obj@reductions),
    clusters = !is.null(seurat_obj$seurat_clusters) || 
      length(unique(Idents(seurat_obj))) > 1,
    cell_cycle = all(c("S.Score", "G2M.Score", "Phase") %in% 
                       colnames(seurat_obj@meta.data))
  )
}

#' Get available reductions from Seurat object
#'
#' @param seurat_obj A Seurat object
#' @return Character vector of reduction names
#' @noRd
get_available_reductions <- function(seurat_obj) {
  if (is.null(seurat_obj)) return(character(0))
  names(seurat_obj@reductions)
}

#' Get available metadata columns
#'
#' @param seurat_obj A Seurat object
#' @return Character vector of metadata column names
#' @noRd
get_metadata_columns <- function(seurat_obj) {
  if (is.null(seurat_obj)) return(character(0))
  colnames(seurat_obj@meta.data)
}

#' Create a safe feature plot with error handling
#'
#' @param seurat_obj A Seurat object
#' @param feature Feature to plot
#' @param reduction Reduction to use
#' @return A ggplot object or NULL
#' @noRd
safe_feature_plot <- function(seurat_obj, feature, reduction = "umap") {
  tryCatch({
    FeaturePlot(seurat_obj, 
                features = feature,
                reduction = reduction) + 
      scale_color_viridis() +
      theme_minimal()
  }, error = function(e) {
    message("Error creating feature plot: ", e$message)
    NULL
  })
}
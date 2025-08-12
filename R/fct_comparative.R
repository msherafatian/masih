' Calculate pathway score for a single pathway - FIXED VERSION
#'
#' @param seurat_obj Seurat object
#' @param pathway Pathway name
#' @param ctrl Number of control features
#' @return List with updated Seurat object and score name
#' @noRd
calculate_pathway_score <- function(seurat_obj, pathway, ctrl = 20) {
  all_genes <- rownames(seurat_obj)
  gene_list <- get(pathway)$symbol
  gene_list_filtered <- gene_list[gene_list %in% all_genes]
  
  if (length(gene_list_filtered) > 0) {
    
    # APPLY THE SAME THREE-TIER FALLBACK AS CANCERSEA MODULE
    tryCatch({
      # First attempt: Simple parameters like the old working app
      seurat_obj <- AddModuleScore(
        object = seurat_obj,
        features = list(gene_list_filtered),
        name = paste0(pathway, "_"),
        ctrl = 20  # Same as old app
        # No nbin parameter = uses default (safer)
      )
      
      score_name <- paste0(pathway, "_1")
      
      return(list(
        seurat_obj = seurat_obj,
        score_name = score_name
      ))
      
    }, error = function(e) {
      # If that fails, try with minimal parameters
      if (grepl("Insufficient data", e$message) || grepl("bin", e$message)) {
        
        tryCatch({
          # Second attempt: Minimal parameters
          seurat_obj <- AddModuleScore(
            object = seurat_obj,
            features = list(gene_list_filtered),
            name = paste0(pathway, "_"),
            ctrl = 5,     # Minimal control genes
            nbin = 5      # Minimal bins
          )
          
          score_name <- paste0(pathway, "_1")
          
          return(list(
            seurat_obj = seurat_obj,
            score_name = score_name
          ))
          
        }, error = function(e2) {
          # Third attempt: Disable sophisticated control matching
          tryCatch({
            seurat_obj <- AddModuleScore(
              object = seurat_obj,
              features = list(gene_list_filtered),
              name = paste0(pathway, "_"),
              ctrl = 1,      # Just 1 control gene per feature
              nbin = 1       # No binning
            )
            
            score_name <- paste0(pathway, "_1")
            
            return(list(
              seurat_obj = seurat_obj,
              score_name = score_name
            ))
            
          }, error = function(e3) {
            # If all attempts fail, return NULL
            warning(paste("Failed to calculate", pathway, "score:", e3$message))
            return(NULL)
          })
        })
        
      } else {
        # For non-binning errors, just return NULL
        warning(paste("Error calculating", pathway, "score:", e$message))
        return(NULL)
      }
    })
    
  }
  
  return(NULL)
}
#' Perform pairwise pathway correlation analysis
#'
#' @param seurat_obj Seurat object
#' @param score_columns Named vector of score column names
#' @param method Correlation method
#' @return Correlation results
#' @noRd
pathway_correlation_analysis <- function(seurat_obj, score_columns, 
                                         method = "pearson") {
  
  score_data <- seurat_obj@meta.data[, score_columns, drop = FALSE]
  colnames(score_data) <- names(score_columns)
  
  # Calculate correlations
  cor_matrix <- cor(score_data, use = "complete.obs", method = method)
  
  # Calculate p-values
  n <- nrow(score_data)
  p_matrix <- matrix(NA, nrow = ncol(score_data), ncol = ncol(score_data))
  
  for (i in 1:(ncol(score_data) - 1)) {
    for (j in (i + 1):ncol(score_data)) {
      test_result <- cor.test(score_data[, i], score_data[, j], 
                              method = method)
      p_matrix[i, j] <- p_matrix[j, i] <- test_result$p.value
    }
  }
  
  dimnames(p_matrix) <- dimnames(cor_matrix)
  
  return(list(
    correlation = cor_matrix,
    p_values = p_matrix
  ))
}

#' Create pathway activity profile
#'
#' @param seurat_obj Seurat object
#' @param score_columns Named vector of score column names
#' @param group_by Grouping variable
#' @return Data frame with activity profiles
#' @noRd
create_pathway_profile <- function(seurat_obj, score_columns, 
                                   group_by = "seurat_clusters") {
  
  profile_data <- data.frame(
    Group = seurat_obj@meta.data[[group_by]],
    stringsAsFactors = FALSE
  )
  
  for (pathway in names(score_columns)) {
    score_col <- score_columns[[pathway]]
    profile_data[[pathway]] <- seurat_obj@meta.data[[score_col]]
  }
  
  # Calculate mean and SE for each group
  profile_summary <- profile_data %>%
    pivot_longer(cols = -Group, 
                 names_to = "Pathway", 
                 values_to = "Score") %>%
    group_by(Group, Pathway) %>%
    summarise(
      Mean = mean(Score, na.rm = TRUE),
      SE = sd(Score, na.rm = TRUE) / sqrt(n()),
      Median = median(Score, na.rm = TRUE),
      Q1 = quantile(Score, 0.25, na.rm = TRUE),
      Q3 = quantile(Score, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(profile_summary)
}

#' Identify pathway-specific cell populations
#'
#' @param seurat_obj Seurat object
#' @param score_column Score column name
#' @param threshold Threshold for high expression
#' @param quantile_threshold Use quantile instead of absolute threshold
#' @return Vector of cell names
#' @noRd
identify_pathway_cells <- function(seurat_obj, score_column, 
                                   threshold = NULL, 
                                   quantile_threshold = 0.9) {
  
  scores <- seurat_obj@meta.data[[score_column]]
  
  if (is.null(threshold)) {
    threshold <- quantile(scores, quantile_threshold, na.rm = TRUE)
  }
  
  high_cells <- rownames(seurat_obj@meta.data)[scores > threshold]
  
  return(high_cells)
}

#' Compare pathway distributions between groups
#'
#' @param seurat_obj Seurat object
#' @param score_column Score column name
#' @param group_var Grouping variable
#' @param test_type Type of statistical test
#' @return Test results
#' @noRd
compare_pathway_distributions <- function(seurat_obj, score_column, 
                                          group_var, 
                                          test_type = "wilcox") {
  
  groups <- unique(seurat_obj@meta.data[[group_var]])
  
  if (length(groups) == 2) {
    # Two-group comparison
    group1_scores <- seurat_obj@meta.data[
      seurat_obj@meta.data[[group_var]] == groups[1], score_column
    ]
    group2_scores <- seurat_obj@meta.data[
      seurat_obj@meta.data[[group_var]] == groups[2], score_column
    ]
    
    if (test_type == "wilcox") {
      test_result <- wilcox.test(group1_scores, group2_scores)
    } else if (test_type == "t") {
      test_result <- t.test(group1_scores, group2_scores)
    }
    
    return(list(
      test = test_type,
      statistic = test_result$statistic,
      p_value = test_result$p.value,
      group1 = groups[1],
      group2 = groups[2],
      mean_diff = mean(group1_scores, na.rm = TRUE) - 
        mean(group2_scores, na.rm = TRUE)
    ))
    
  } else {
    # Multiple group comparison
    formula <- as.formula(paste(score_column, "~", group_var))
    
    if (test_type == "kruskal") {
      test_result <- kruskal.test(formula, data = seurat_obj@meta.data)
    } else if (test_type == "anova") {
      test_result <- aov(formula, data = seurat_obj@meta.data)
      test_result <- summary(test_result)
    }
    
    return(test_result)
  }
}
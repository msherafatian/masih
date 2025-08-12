#' compare UI Function
#'
#' @description A shiny Module for comparative pathway analysis.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_compare_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      box(
        title = "Pathway Comparison Setup",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        selectizeInput(ns("compare_pathways"), 
                       "Select pathways to compare:",
                       choices = NULL,
                       multiple = TRUE,
                       options = list(maxItems = 5)),
        actionButton(ns("run_comparison"), 
                     "Run Comparison", 
                     class = "btn-success")
      )
    ),
    
    fluidRow(
      box(
        title = "Correlation Matrix",
        status = "info",
        solidHeader = TRUE,
        width = 6,
        plotOutput(ns("correlationMatrix"), height = "500px")
      ),
      
      box(
        title = "Pathway Expression by Cluster",
        status = "info",
        solidHeader = TRUE,
        width = 6,
        plotOutput(ns("pathwayByCluster"), height = "500px")
      )
    ),
    
    fluidRow(
      box(
        title = "Pathway Comparison Heatmap",
        status = "warning",
        solidHeader = TRUE,
        width = 12,
        plotOutput(ns("comparisonHeatmap"), height = "600px")
      )
    ),
    
    fluidRow(
      box(
        title = "Statistical Summary",
        status = "success",
        solidHeader = TRUE,
        width = 12,
        DTOutput(ns("comparisonStats"))
      )
    )
  )
}

#' compare Server Functions
#'
#' @import corrplot
#' @import tidyr
#' @noRd 
mod_compare_server <- function(id, seurat_obj, processed, cancersea_scores, main_values){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
    # Update pathway choices
    observe({
      data('available_pathways', package = 'cancersea')
      updateSelectizeInput(session, "compare_pathways", 
                           choices = available_pathways)
    })
    
    # Local reactive values
    local_values <- reactiveValues(
      comparison_results = NULL
    )
    
    # Run pathway comparison
    observeEvent(input$run_comparison, {
      req(length(input$compare_pathways) >= 2)
      req(seurat_obj())
      
      withProgress(message = 'Running pathway comparison...', value = 0, {
        
        # Calculate scores for selected pathways if not already done
        incProgress(0.2, detail = "Calculating pathway scores...")
        
        successful_pathways <- c()
        failed_pathways <- c()
        
        for (pathway in input$compare_pathways) {
          if (!(pathway %in% names(main_values$cancersea_scores))) {
            
            incProgress(0.1, detail = paste("Calculating", pathway, "..."))
            
            result <- calculate_pathway_score(seurat_obj(), pathway)
            if (!is.null(result)) {
              main_values$seurat_obj <- result$seurat_obj
              main_values$cancersea_scores[[pathway]] <- result$score_name
              successful_pathways <- c(successful_pathways, pathway)
            } else {
              failed_pathways <- c(failed_pathways, pathway)
            }
          } else {
            successful_pathways <- c(successful_pathways, pathway)
          }
        }
        
        incProgress(0.6, detail = "Analyzing correlations...")
        
        # Only proceed if we have at least 2 successful pathways
        if (length(successful_pathways) >= 2) {
          # Store comparison results using only successful pathways
          local_values$comparison_results <- list(
            pathways = successful_pathways,
            scores = main_values$cancersea_scores[successful_pathways]
          )
          
          incProgress(1, detail = "Complete!")
          
          # Show success message with details
          if (length(failed_pathways) > 0) {
            showNotification(
              paste("Pathway comparison completed for", length(successful_pathways), "pathways.",
                    "Failed pathways:", paste(failed_pathways, collapse = ", ")),
              type = "warning", 
              duration = 5
            )
          } else {
            showNotification("Pathway comparison completed successfully!", 
                             type = "message", duration = 3)
          }
          
        } else {
          # Not enough successful pathways
          showNotification(
            paste("Comparison failed: Only", length(successful_pathways), 
                  "pathways calculated successfully. Need at least 2 pathways.",
                  if(length(failed_pathways) > 0) paste("Failed:", paste(failed_pathways, collapse = ", "))),
            type = "error", 
            duration = 10
          )
          
          local_values$comparison_results <- NULL
        }
      })
    })
    
    # Output: Correlation matrix
    output$correlationMatrix <- renderPlot({
      req(length(input$compare_pathways) >= 2)
      req(seurat_obj())
      req(local_values$comparison_results)
      
      # Extract scores
      score_cols <- unlist(local_values$comparison_results$scores)
      score_data <- seurat_obj()@meta.data[, score_cols, drop = FALSE]
      colnames(score_data) <- names(local_values$comparison_results$scores)
      
      # Calculate correlation
      cor_matrix <- cor(score_data, use = "complete.obs")
      
      # Plot
      corrplot::corrplot(cor_matrix, 
                         method = "color",
                         type = "upper",
                         order = "hclust",
                         tl.cex = 0.8,
                         tl.col = "black",
                         col = colorRampPalette(c("blue", "white", "red"))(100),
                         addCoef.col = "black",
                         number.cex = 0.7)
    })
    
    # Output: Pathway by cluster comparison
    output$pathwayByCluster <- renderPlot({
      req(length(input$compare_pathways) >= 1)
      req(seurat_obj())
      req(local_values$comparison_results)
      
      # Prepare data
      plot_data <- data.frame(Cluster = Idents(seurat_obj()))
      
      for (pathway in names(local_values$comparison_results$scores)) {
        score_col <- local_values$comparison_results$scores[[pathway]]
        plot_data[[pathway]] <- seurat_obj()@meta.data[[score_col]]
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
      ggplot(plot_summary, aes(x = Cluster, y = Mean_Score, fill = Pathway)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_errorbar(aes(ymin = Mean_Score - SE, ymax = Mean_Score + SE),
                      position = position_dodge(0.9), width = 0.25) +
        theme_minimal() +
        labs(y = "Mean Score", title = "Pathway Scores by Cluster") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_brewer(palette = "Set2")
    })
    # Output: Comparison heatmap
    output$comparisonHeatmap <- renderPlot({
      req(local_values$comparison_results)
      req(seurat_obj())
      req(length(local_values$comparison_results$scores) >= 2)
      
      # Create matrix of average scores per cluster
      score_matrix <- calculate_pathway_averages(
        seurat_obj(), 
        local_values$comparison_results$scores
      )
      
      # Convert to matrix format
      score_mat <- as.matrix(score_matrix[, -1])
      rownames(score_mat) <- score_matrix$Cluster
      
      # Scale and plot
      create_heatmap(score_mat, scale = TRUE)
    })
    
    # Output: Comparison statistics
    output$comparisonStats <- renderDT({
      req(local_values$comparison_results)
      req(seurat_obj())
      
      # Calculate statistics for each pathway
      stats_list <- list()
      
      for (pathway in names(local_values$comparison_results$scores)) {
        score_col <- local_values$comparison_results$scores[[pathway]]
        scores <- seurat_obj()@meta.data[[score_col]]
        
        stats_list[[pathway]] <- data.frame(
          Pathway = pathway,
          Mean = mean(scores, na.rm = TRUE),
          Median = median(scores, na.rm = TRUE),
          SD = sd(scores, na.rm = TRUE),
          Min = min(scores, na.rm = TRUE),
          Max = max(scores, na.rm = TRUE),
          CV = sd(scores, na.rm = TRUE) / mean(scores, na.rm = TRUE)
        )
      }
      
      stats_df <- do.call(rbind, stats_list)
      rownames(stats_df) <- NULL
      
      # Format numeric columns
      stats_df[, 2:7] <- round(stats_df[, 2:7], 3)
      
      datatable(stats_df,
                options = list(
                  pageLength = 10,
                  dom = 'Bfrtip',
                  buttons = c('copy', 'csv', 'excel')
                ),
                rownames = FALSE)
    })
    
  })
}
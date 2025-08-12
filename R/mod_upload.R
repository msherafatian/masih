#' upload UI Function
#'
#' @description A shiny Module for data upload and processing.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
#' @import shinydashboard
mod_upload_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      shinydashboard::box(
        title = "Data Upload",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        
        # Input type selection
        radioButtons(ns("input_type"), "Select input type:",
                     choices = c("Seurat Object (.rds)" = "seurat",
                                 "Cell Ranger Output (10X folder)" = "cellranger"),
                     selected = "seurat"),
        
        hr(),
        
        # Conditional panel for Seurat object upload
        conditionalPanel(
          condition = paste0("input['", ns("input_type"), "'] == 'seurat'"),
          fileInput(ns("file"), "Choose RDS file (Seurat Object)",
                    accept = c(".rds"))
        ),
        
        # Conditional panel for Cell Ranger upload
        conditionalPanel(
          condition = paste0("input['", ns("input_type"), "'] == 'cellranger'"),
          h4("Cell Ranger Data Input"),
          p("Select the folder containing the Cell Ranger output files (barcodes.tsv, features.tsv/genes.tsv, matrix.mtx)"),
          
          # Multiple file upload for Cell Ranger files
          fileInput(ns("cellranger_files"), 
                    "Upload Cell Ranger files (select all 3 files):",
                    accept = c(".tsv", ".gz", ".mtx"),
                    multiple = TRUE),
          
          p("OR"),
          
          # Text input for server path
          textInput(ns("cellranger_path"), 
                    "Enter server path to Cell Ranger output folder:",
                    placeholder = "/path/to/filtered_gene_bc_matrices/"),
          
          hr(),
          
          h5("Cell Ranger Object Creation Parameters:"),
          fluidRow(
            column(4,
                   textInput(ns("project_name"), "Project name:", 
                             value = "SingleCell_Project")
            ),
            column(4,
                   numericInput(ns("min_cells"), "Min cells per gene:", 
                                value = 3, min = 1, max = 100)
            ),
            column(4,
                   numericInput(ns("min_features"), "Min genes per cell:", 
                                value = app_config$default_min_features, 
                                min = 50, max = 1000)
            )
          ),
          
          actionButton(ns("load_cellranger"), "Load Cell Ranger Data", 
                       class = "btn-warning")
        ),
        
        conditionalPanel(
          condition = paste0("output['", ns("fileUploaded"), "']"),
          hr(),
          h4("Processing Options"),
          fluidRow(
            column(4,
                   numericInput(ns("nfeatures"), "Number of variable features:", 
                                value = app_config$default_nfeatures, 
                                min = 1000, max = 5000)
            ),
            column(4,
                   numericInput(ns("dims"), "Number of PCs for clustering:", 
                                value = app_config$default_dims, 
                                min = 10, max = 50)
            ),
            column(4,
                   numericInput(ns("resolution"), "Clustering resolution:", 
                                value = app_config$default_resolution, 
                                min = 0.1, max = 2, step = 0.1)
            )
          ),
          fluidRow(
            column(4,
                   checkboxInput(ns("run_sct"), "Run SCTransform", value = TRUE)
            ),
            column(4,
                   checkboxInput(ns("recluster"), "Force re-clustering", value = FALSE)
            ),
            column(4,
                   checkboxInput(ns("run_qc"), "Run QC filtering", value = FALSE)
            )
          ),
          
          # QC parameters (shown only when run_qc is checked)
          conditionalPanel(
            condition = paste0("input['", ns("run_qc"), "'] == true"),
            h5("QC Parameters:"),
            fluidRow(
              column(4,
                     numericInput(ns("max_mt"), "Max mitochondrial %:", 
                                  value = app_config$default_max_mt, 
                                  min = 1, max = 100)
              ),
              column(4,
                     numericInput(ns("min_ncount"), "Min nCount_RNA:", 
                                  value = app_config$default_min_ncount, 
                                  min = 0, max = 10000)
              ),
              column(4,
                     numericInput(ns("max_ncount"), "Max nCount_RNA:", 
                                  value = app_config$default_max_ncount, 
                                  min = 1000, max = 100000)
              )
            )
          ),
          actionButton(ns("process"), "Process Data", class = "btn-success"),
          hr(),
          h4("Data Summary"),
          verbatimTextOutput(ns("dataSummary"))
        )
      )
    ),
    
    fluidRow(
      box(
        title = "Processing Log",
        status = "info",
        solidHeader = TRUE,
        width = 12,
        verbatimTextOutput(ns("processingLog"))
      )
    )
  )
}

#' upload Server Functions
#'
#' @import Seurat
#' @noRd 
mod_upload_server <- function(id, main_values){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
    # Local reactive values for this module
    local_values <- reactiveValues(
      log_text = "",
      existing_analyses = list(
        normalized = FALSE,
        variable_features = FALSE,
        scaled = FALSE,
        sct = FALSE,
        pca = FALSE,
        umap = FALSE,
        tsne = FALSE,
        clusters = FALSE,
        cell_cycle = FALSE
      ),
      is_cellranger = FALSE
    )
    
    # File upload check
    output$fileUploaded <- reactive({
      if (input$input_type == "seurat") {
        return(!is.null(input$file))
      } else {
        return(!is.null(main_values$seurat_obj))
      }
    })
    outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)
    
    # Load Cell Ranger data
    observeEvent(input$load_cellranger, {
      local_values$log_text <- "Loading Cell Ranger data...\n"
      
      tryCatch({
        # Check if files were uploaded or path was provided
        if (!is.null(input$cellranger_files) && nrow(input$cellranger_files) >= 3) {
          # Handle uploaded files
          local_values$log_text <- paste0(local_values$log_text, 
                                          "Processing uploaded Cell Ranger files...\n")
          
          # Create temporary directory for files
          temp_dir <- file.path(tempdir(), "cellranger_data")
          dir.create(temp_dir, showWarnings = FALSE)
          
          # Copy uploaded files to temp directory
          for (i in 1:nrow(input$cellranger_files)) {
            file.copy(input$cellranger_files$datapath[i], 
                      file.path(temp_dir, input$cellranger_files$name[i]))
          }
          
          # Read 10X data
          cellranger_data <- Read10X(data.dir = temp_dir)
          
        } else if (nchar(input$cellranger_path) > 0) {
          # Handle server path
          local_values$log_text <- paste0(local_values$log_text, 
                                          "Loading from path: ", input$cellranger_path, "\n")
          
          if (!dir.exists(input$cellranger_path)) {
            stop("Directory does not exist")
          }
          
          cellranger_data <- Read10X(data.dir = input$cellranger_path)
          
        } else {
          stop("Please either upload files or provide a path")
        }
        
        # Create Seurat object
        local_values$log_text <- paste0(local_values$log_text, "Creating Seurat object...\n")
        main_values$seurat_obj <- CreateSeuratObject(
          counts = cellranger_data, 
          project = input$project_name, 
          min.cells = input$min_cells, 
          min.features = input$min_features
        )
        
        local_values$log_text <- paste0(local_values$log_text, "✓ Seurat object created successfully\n")
        local_values$log_text <- paste0(local_values$log_text, "  - Cells: ", 
                                        ncol(main_values$seurat_obj), "\n")
        local_values$log_text <- paste0(local_values$log_text, "  - Genes: ", 
                                        nrow(main_values$seurat_obj), "\n")
        
        # Mark as Cell Ranger input (needs full processing)
        local_values$is_cellranger <- TRUE
        
        # Reset existing analyses
        local_values$existing_analyses <- list(
          normalized = FALSE,
          variable_features = FALSE,
          scaled = FALSE,
          sct = FALSE,
          pca = FALSE,
          umap = FALSE,
          tsne = FALSE,
          clusters = FALSE,
          cell_cycle = FALSE
        )
        
        local_values$log_text <- paste0(local_values$log_text, 
                                        "\n⚠ Cell Ranger data loaded. Full processing required.\n")
        local_values$log_text <- paste0(local_values$log_text, 
                                        "Click 'Process Data' to run all analyses.\n")
        
      }, error = function(e) {
        local_values$log_text <- paste0(local_values$log_text, 
                                        "✗ Error loading Cell Ranger data: ", e$message, "\n")
      })
    })
    
    # Process uploaded Seurat file
    observeEvent(input$file, {
      req(input$file)
      req(input$input_type == "seurat")
      
      local_values$log_text <- "Loading Seurat object...\n"
      local_values$is_cellranger <- FALSE
      
      tryCatch({
        main_values$seurat_obj <- readRDS(input$file$datapath)
        local_values$log_text <- paste0(local_values$log_text, 
                                        "✓ Seurat object loaded successfully\n")
        
        # Check existing analyses with error catching for each step
        local_values$log_text <- paste0(local_values$log_text, "\nChecking existing analyses:\n")
        
        # Check for normalization
        tryCatch({
          local_values$existing_analyses$normalized <- 
            !is.null(main_values$seurat_obj@assays$RNA@data) && 
            sum(main_values$seurat_obj@assays$RNA@data) > 0
          
          if (local_values$existing_analyses$normalized) {
            local_values$log_text <- paste0(local_values$log_text, "✓ Normalization found\n")
          }
        }, error = function(e) {
          local_values$log_text <- paste0(local_values$log_text, "✗ Error checking normalization: ", e$message, "\n")
        })
        
        # Check for variable features
        tryCatch({
          local_values$existing_analyses$variable_features <- 
            length(VariableFeatures(main_values$seurat_obj)) > 0
          
          if (local_values$existing_analyses$variable_features) {
            local_values$log_text <- paste0(local_values$log_text, "✓ Variable features found (", 
                                            length(VariableFeatures(main_values$seurat_obj)), 
                                            " features)\n")
          }
        }, error = function(e) {
          local_values$log_text <- paste0(local_values$log_text, "✗ Error checking variable features: ", e$message, "\n")
        })
        
        # Check for scaling
        tryCatch({
          local_values$existing_analyses$scaled <- 
            !is.null(main_values$seurat_obj@assays$RNA@scale.data) && 
            nrow(main_values$seurat_obj@assays$RNA@scale.data) > 0
          
          if (local_values$existing_analyses$scaled) {
            local_values$log_text <- paste0(local_values$log_text, "✓ Scaled data found\n")
          }
        }, error = function(e) {
          local_values$log_text <- paste0(local_values$log_text, "✗ Error checking scaled data: ", e$message, "\n")
        })
        
        # Check for SCT
        tryCatch({
          local_values$existing_analyses$sct <- "SCT" %in% names(main_values$seurat_obj@assays)
          if (local_values$existing_analyses$sct) {
            local_values$log_text <- paste0(local_values$log_text, "✓ SCTransform found\n")
          }
        }, error = function(e) {
          local_values$log_text <- paste0(local_values$log_text, "✗ Error checking SCT: ", e$message, "\n")
        })
        
        # Check for reductions
        tryCatch({
          local_values$existing_analyses$pca <- "pca" %in% names(main_values$seurat_obj@reductions)
          if (local_values$existing_analyses$pca) {
            local_values$log_text <- paste0(local_values$log_text, "✓ PCA found\n")
          }
        }, error = function(e) {
          local_values$log_text <- paste0(local_values$log_text, "✗ Error checking PCA: ", e$message, "\n")
        })
        
        tryCatch({
          local_values$existing_analyses$umap <- "umap" %in% names(main_values$seurat_obj@reductions)
          if (local_values$existing_analyses$umap) {
            local_values$log_text <- paste0(local_values$log_text, "✓ UMAP found\n")
          }
        }, error = function(e) {
          local_values$log_text <- paste0(local_values$log_text, "✗ Error checking UMAP: ", e$message, "\n")
        })
        
        tryCatch({
          local_values$existing_analyses$tsne <- "tsne" %in% names(main_values$seurat_obj@reductions)
          if (local_values$existing_analyses$tsne) {
            local_values$log_text <- paste0(local_values$log_text, "✓ t-SNE found\n")
          }
        }, error = function(e) {
          local_values$log_text <- paste0(local_values$log_text, "✗ Error checking t-SNE: ", e$message, "\n")
        })
        
        # Check for clustering
        tryCatch({
          local_values$existing_analyses$clusters <- 
            !is.null(main_values$seurat_obj$seurat_clusters) || 
            length(unique(Idents(main_values$seurat_obj))) > 1
          
          if (local_values$existing_analyses$clusters) {
            local_values$log_text <- paste0(local_values$log_text, "✓ Clustering found (", 
                                            length(unique(Idents(main_values$seurat_obj))), " clusters)\n")
          }
        }, error = function(e) {
          local_values$log_text <- paste0(local_values$log_text, "✗ Error checking clustering: ", e$message, "\n")
        })
        
        # Check for cell cycle scores
        tryCatch({
          meta_names <- colnames(main_values$seurat_obj@meta.data)
          local_values$existing_analyses$cell_cycle <- 
            all(c("S.Score", "G2M.Score", "Phase") %in% meta_names)
          
          if (local_values$existing_analyses$cell_cycle) {
            local_values$log_text <- paste0(local_values$log_text, "✓ Cell cycle scores found\n")
          }
        }, error = function(e) {
          local_values$log_text <- paste0(local_values$log_text, "✗ Error checking cell cycle: ", e$message, "\n")
        })
        
        # Check if fully processed
        if (all(unlist(local_values$existing_analyses[c("normalized", "pca", "umap", 
                                                        "clusters", "cell_cycle")]))) {
          main_values$processed <- TRUE
          local_values$log_text <- paste0(local_values$log_text, 
                                          "\n✓ Object appears to be fully processed!\n")
        } else {
          local_values$log_text <- paste0(local_values$log_text, 
                                          "\n⚠ Some analyses are missing. Click 'Process Data' to run missing steps.\n")
        }
        
      }, error = function(e) {
        local_values$log_text <- paste0(local_values$log_text, 
                                        "✗ Error loading file: ", e$message, "\n")
      })
    })
    
    # Process data - SINGLE observeEvent block
    observeEvent(input$process, {
      req(main_values$seurat_obj)
      
      withProgress(message = 'Processing data...', value = 0, {
        
        tryCatch({
          total_steps <- 10  # Define total steps for progress tracking
          
          # QC filtering (if requested)
          if (input$run_qc || local_values$is_cellranger) {
            incProgress(1/total_steps, detail = "Running QC...")
            local_values$log_text <- paste0(local_values$log_text, "\nRunning QC filtering...\n")
            
            # Calculate mitochondrial percentage
            main_values$seurat_obj[["percent.mt"]] <- PercentageFeatureSet(main_values$seurat_obj, pattern = "^MT-")
            
            # Filter cells
            cells_before <- ncol(main_values$seurat_obj)
            main_values$seurat_obj <- subset(main_values$seurat_obj, 
                                             subset = nFeature_RNA > input$min_features & 
                                               nFeature_RNA < 2500 & 
                                               percent.mt < input$max_mt)
            cells_after <- ncol(main_values$seurat_obj)
            
            local_values$log_text <- paste0(local_values$log_text, 
                                            "✓ QC complete. Cells: ", cells_before, " → ", cells_after, 
                                            " (removed ", cells_before - cells_after, ")\n")
          } else {
            incProgress(1/total_steps, detail = "Skipping QC...")
          }
          
          # Normalize data (if needed)
          if (!local_values$existing_analyses$normalized || local_values$is_cellranger) {
            incProgress(1/total_steps, detail = "Normalizing data...")
            local_values$log_text <- paste0(local_values$log_text, "Normalizing data...\n")
            main_values$seurat_obj <- NormalizeData(main_values$seurat_obj)
            local_values$log_text <- paste0(local_values$log_text, "✓ Normalization complete\n")
          } else {
            local_values$log_text <- paste0(local_values$log_text, "→ Skipping normalization (already done)\n")
            incProgress(1/total_steps, detail = "Skipping normalization...")
          }
          
          # Find variable features (if needed)
          if (!local_values$existing_analyses$variable_features || local_values$is_cellranger) {
            incProgress(1/total_steps, detail = "Finding variable features...")
            local_values$log_text <- paste0(local_values$log_text, "Finding variable features...\n")
            main_values$seurat_obj <- FindVariableFeatures(main_values$seurat_obj, 
                                                           selection.method = "vst", 
                                                           nfeatures = input$nfeatures)
            local_values$log_text <- paste0(local_values$log_text, "✓ Variable features identified\n")
          } else {
            local_values$log_text <- paste0(local_values$log_text, "→ Skipping variable features (already done)\n")
            incProgress(1/total_steps, detail = "Skipping variable features...")
          }
          
          # Scale data (if needed)
          if (!local_values$existing_analyses$scaled || local_values$is_cellranger) {
            incProgress(1/total_steps, detail = "Scaling data...")
            local_values$log_text <- paste0(local_values$log_text, "Scaling data...\n")
            all.genes <- rownames(main_values$seurat_obj)
            main_values$seurat_obj <- ScaleData(main_values$seurat_obj, features = all.genes)
            local_values$log_text <- paste0(local_values$log_text, "✓ Scaling complete\n")
          } else {
            local_values$log_text <- paste0(local_values$log_text, "→ Skipping scaling (already done)\n")
            incProgress(1/total_steps, detail = "Skipping scaling...")
          }
          
          # Run PCA (if needed)
          if (!local_values$existing_analyses$pca || local_values$is_cellranger) {
            incProgress(1/total_steps, detail = "Running PCA...")
            local_values$log_text <- paste0(local_values$log_text, "Running PCA...\n")
            main_values$seurat_obj <- RunPCA(main_values$seurat_obj)
            local_values$log_text <- paste0(local_values$log_text, "✓ PCA complete\n")
          } else {
            local_values$log_text <- paste0(local_values$log_text, "→ Skipping PCA (already done)\n")
            incProgress(1/total_steps, detail = "Skipping PCA...")
          }
          
          # Run UMAP (if needed)
          if (!local_values$existing_analyses$umap || local_values$is_cellranger) {
            incProgress(1/total_steps, detail = "Running UMAP...")
            local_values$log_text <- paste0(local_values$log_text, "Running UMAP...\n")
            main_values$seurat_obj <- RunUMAP(main_values$seurat_obj, dims = 1:input$dims)
            local_values$log_text <- paste0(local_values$log_text, "✓ UMAP complete\n")
          } else {
            local_values$log_text <- paste0(local_values$log_text, "→ Skipping UMAP (already done)\n")
            incProgress(1/total_steps, detail = "Skipping UMAP...")
          }
          
          # Run t-SNE (if needed)
          if (!local_values$existing_analyses$tsne || local_values$is_cellranger) {
            incProgress(1/total_steps, detail = "Running t-SNE...")
            local_values$log_text <- paste0(local_values$log_text, "Running t-SNE...\n")
            main_values$seurat_obj <- RunTSNE(main_values$seurat_obj, dims = 1:input$dims, check_duplicates = FALSE)
            local_values$log_text <- paste0(local_values$log_text, "✓ t-SNE complete\n")
          } else {
            local_values$log_text <- paste0(local_values$log_text, "→ Skipping t-SNE (already done)\n")
            incProgress(1/total_steps, detail = "Skipping t-SNE...")
          }
          
          # Find neighbors and clusters (if needed)
          if (!local_values$existing_analyses$clusters || local_values$is_cellranger) {
            incProgress(1/total_steps, detail = "Finding neighbors and clustering...")
            local_values$log_text <- paste0(local_values$log_text, "Finding neighbors...\n")
            main_values$seurat_obj <- FindNeighbors(main_values$seurat_obj, dims = 1:input$dims)
            main_values$seurat_obj <- FindClusters(main_values$seurat_obj, resolution = input$resolution)
            local_values$log_text <- paste0(local_values$log_text, "✓ Clustering complete\n")
          } else {
            local_values$log_text <- paste0(local_values$log_text, "→ Skipping clustering (already done)\n")
            incProgress(1/total_steps, detail = "Checking clustering...")
            # Option to recluster with different resolution
            if (input$recluster) {
              local_values$log_text <- paste0(local_values$log_text, "Re-clustering with resolution ", input$resolution, "...\n")
              main_values$seurat_obj <- FindClusters(main_values$seurat_obj, resolution = input$resolution)
              local_values$log_text <- paste0(local_values$log_text, "✓ Re-clustering complete\n")
            }
          }
          
          # SCTransform (if needed and user wants it)
          if (!local_values$existing_analyses$sct && input$run_sct) {
            incProgress(1/total_steps, detail = "Running SCTransform...")
            local_values$log_text <- paste0(local_values$log_text, "Running SCTransform...\n")
            
            # IMPROVED: Better SCTransform parameters
            main_values$seurat_obj <- SCTransform(
              main_values$seurat_obj, 
              method = "glmGamPoi",  # Use faster, more robust method
              ncells = 5000,         # Limit cells for faster computation
              n_genes = 3000,        # Limit genes to most variable
              verbose = FALSE,
              return.only.var.genes = FALSE,
              vst.flavor = "v2",     # Use improved version
              conserve.memory = TRUE # Reduce memory usage
            )
            local_values$log_text <- paste0(local_values$log_text, "✓ SCTransform complete\n")
          } else if (local_values$existing_analyses$sct) {
            local_values$log_text <- paste0(local_values$log_text, "→ Skipping SCTransform (already done)\n")
            incProgress(1/total_steps, detail = "Skipping SCTransform...")
          } else {
            local_values$log_text <- paste0(local_values$log_text, "→ Skipping SCTransform (not requested)\n")
            incProgress(1/total_steps, detail = "Skipping SCTransform...")
          }
          
          # Cell cycle scoring (if needed)
          if (!local_values$existing_analyses$cell_cycle || local_values$is_cellranger) {
            incProgress(1/total_steps, detail = "Cell cycle scoring...")
            local_values$log_text <- paste0(local_values$log_text, "Cell cycle scoring...\n")
            s.genes <- cc.genes$s.genes
            g2m.genes <- cc.genes$g2m.genes
            main_values$seurat_obj <- CellCycleScoring(main_values$seurat_obj, 
                                                       s.features = s.genes, 
                                                       g2m.features = g2m.genes, 
                                                       set.ident = FALSE)
            local_values$log_text <- paste0(local_values$log_text, "✓ Cell cycle scoring complete\n")
          } else {
            local_values$log_text <- paste0(local_values$log_text, "→ Skipping cell cycle scoring (already done)\n")
            incProgress(1/total_steps, detail = "Skipping cell cycle...")
          }
          
          # Build cluster tree
          incProgress(1/total_steps, detail = "Building cluster tree...")
          local_values$log_text <- paste0(local_values$log_text, "Building cluster tree...\n")
          main_values$seurat_obj <- BuildClusterTree(main_values$seurat_obj,
                                                     dims = 1:input$dims,
                                                     reorder = TRUE,
                                                     reorder.numeric = TRUE)
          local_values$log_text <- paste0(local_values$log_text, "✓ Cluster tree built\n")
          
          # Final completion
          main_values$processed <- TRUE
          local_values$is_cellranger <- FALSE  # Reset flag after processing
          local_values$log_text <- paste0(local_values$log_text, "\n✓ All processing complete!\n")
          
        }, error = function(e) {
          local_values$log_text <- paste0(local_values$log_text, 
                                          "\n✗ Error during processing: ", e$message, "\n")
        })
      })
    })
    
    # Output: Data summary - Fixed variable names
    output$dataSummary <- renderPrint({
      if (!is.null(main_values$seurat_obj)) {
        cat("Input type:", ifelse(local_values$is_cellranger, "Cell Ranger", "Seurat Object"), "\n")
        cat("Number of cells:", ncol(main_values$seurat_obj), "\n")
        cat("Number of genes:", nrow(main_values$seurat_obj), "\n")
        
        # Show existing analyses status (like in original)
        if (length(local_values$existing_analyses) > 0) {
          cat("\nExisting analyses:\n")
          cat("- Normalized:", ifelse(local_values$existing_analyses$normalized, "Yes", "No"), "\n")
          cat("- Variable features:", ifelse(local_values$existing_analyses$variable_features, "Yes", "No"), "\n")
          cat("- Scaled:", ifelse(local_values$existing_analyses$scaled, "Yes", "No"), "\n")
          cat("- SCTransform:", ifelse(local_values$existing_analyses$sct, "Yes", "No"), "\n")
          cat("- PCA:", ifelse(local_values$existing_analyses$pca, "Yes", "No"), "\n")
          cat("- UMAP:", ifelse(local_values$existing_analyses$umap, "Yes", "No"), "\n")
          cat("- t-SNE:", ifelse(local_values$existing_analyses$tsne, "Yes", "No"), "\n")
          cat("- Clustering:", ifelse(local_values$existing_analyses$clusters, "Yes", "No"), "\n")
          cat("- Cell cycle:", ifelse(local_values$existing_analyses$cell_cycle, "Yes", "No"), "\n")
        }
        
        if (main_values$processed || local_values$existing_analyses$clusters) {
          cat("\nNumber of clusters:", length(unique(Idents(main_values$seurat_obj))), "\n")
          cat("Reductions available:", paste(names(main_values$seurat_obj@reductions), collapse = ", "), "\n")
          cat("Assays available:", paste(names(main_values$seurat_obj@assays), collapse = ", "), "\n")
        }
        
        # Show QC metrics if available
        if ("percent.mt" %in% colnames(main_values$seurat_obj@meta.data)) {
          cat("\nQC Metrics:\n")
          cat("- Mean MT%:", round(mean(main_values$seurat_obj$percent.mt), 2), "\n")
          cat("- Mean nFeature_RNA:", round(mean(main_values$seurat_obj$nFeature_RNA), 0), "\n")
          cat("- Mean nCount_RNA:", round(mean(main_values$seurat_obj$nCount_RNA), 0), "\n")
        }
      }
    })
    
    # Output: Processing log
    output$processingLog <- renderPrint({
      cat(local_values$log_text)
    })
    
  })
}
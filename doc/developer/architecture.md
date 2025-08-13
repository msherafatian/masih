# MASIH Architecture Overview

**Technical documentation for developers and contributors.**

---

## 🏗️ System Architecture

### High-Level Design

MASIH follows a **modular Shiny architecture** built with the golem framework, emphasizing:

- **Modularity**: Independent, reusable components
- **Scalability**: Easy to add new analysis modules
- **Maintainability**: Clear separation of concerns
- **Testability**: Unit and integration testing support
- **Reproducibility**: Consistent analysis workflows

```
┌─────────────────────────────────────────────────────────────────┐
│                         MASIH Application                       │
├─────────────────────────────────────────────────────────────────┤
│  Frontend (Shiny UI)                                           │
│  ├── shinydashboard Layout                                     │
│  ├── Modular UI Components                                     │
│  └── Interactive Visualizations (plotly, DT)                  │
├─────────────────────────────────────────────────────────────────┤
│  Backend (Shiny Server)                                        │
│  ├── Reactive Data Management                                  │
│  ├── Analysis Module Servers                                   │
│  └── State Management (reactiveValues)                        │
├─────────────────────────────────────────────────────────────────┤
│  Analysis Layer                                                │
│  ├── Seurat Integration                                        │
│  ├── CancerSEA Pathway Analysis                               │
│  ├── Slingshot Trajectory Inference                           │
│  └── Statistical Analysis Functions                           │
├─────────────────────────────────────────────────────────────────┤
│  Data Layer                                                    │
│  ├── File I/O (10X, Seurat, CSV)                             │
│  ├── Data Validation                                          │
│  └── Export Functions                                         │
└─────────────────────────────────────────────────────────────────┘
```

### Technology Stack

**Core Framework**:
- **R** (≥4.0.0): Programming language
- **Shiny** (≥1.7.0): Web application framework
- **golem** (≥0.4.0): Shiny application framework
- **shinydashboard**: UI framework

**Analysis Libraries**:
- **Seurat** (≥4.0.0): Single-cell analysis
- **SingleCellExperiment**: Bioconductor single-cell infrastructure
- **slingshot**: Trajectory inference
- **cancersea**: Cancer functional state analysis

**Visualization**:
- **plotly**: Interactive plotting
- **ggplot2**: Static plotting
- **DT**: Interactive data tables
- **viridis**: Color palettes

**Data Handling**:
- **dplyr**: Data manipulation
- **openxlsx**: Excel export
- **Matrix**: Sparse matrix support

---

## 📁 Project Structure

### Directory Organization

```
masih/
├── DESCRIPTION                 # Package metadata
├── NAMESPACE                   # Package namespace
├── NEWS.md                     # Version changelog
├── README.md                   # Project overview
├── LICENSE                     # MIT license
├── .gitignore                  # Git ignore rules
├── .Rbuildignore              # R build ignore rules
├── masih.Rproj                # RStudio project file
│
├── R/                         # Main R code
│   ├── app_config.R           # Application configuration
│   ├── app_server.R           # Main server function
│   ├── app_ui.R               # Main UI function
│   ├── run_app.R              # Application launcher
│   │
│   ├── mod_*.R                # UI modules
│   ├── mod_*_server.R         # Server modules
│   ├── fct_*.R                # Helper functions
│   └── utils_*.R              # Utility functions
│
├── inst/                      # Installed files
│   ├── app/                   # Shiny app assets
│   │   ├── www/               # Web assets (CSS, JS, images)
│   │   └── data/              # Example datasets
│   └── extdata/               # External data files
│
├── man/                       # Documentation (auto-generated)
│   └── *.Rd                   # R documentation files
│
├── tests/                     # Test suite
│   ├── testthat/              # Unit tests
│   │   ├── test-*.R           # Test files
│   │   └── helper-*.R         # Test helpers
│   └── testthat.R             # Test configuration
│
├── vignettes/                 # Package vignettes
│   └── *.Rmd                  # R Markdown vignettes
│
├── docs/                      # Documentation
│   ├── user-guide/            # User documentation
│   ├── tutorials/             # Analysis tutorials
│   ├── developer/             # Developer guides
│   └── api/                   # API documentation
│
└── dev/                       # Development scripts
    ├── 01_start.R             # Package setup
    ├── 02_dev.R               # Development helpers
    └── 03_deploy.R            # Deployment scripts
```

### File Naming Conventions

**R Files**:
- `app_*.R`: Main application files
- `mod_*.R`: Shiny UI modules
- `mod_*_server.R`: Shiny server modules
- `fct_*.R`: Business logic functions
- `utils_*.R`: Utility functions

**Module Naming**:
- Module names match their primary function
- Example: `mod_cluster_analysis.R` & `mod_cluster_analysis_server.R`

---

## 🧩 Module Architecture

### Module Design Principles

**🎯 Single Responsibility**: Each module handles one analysis type
**🔗 Loose Coupling**: Modules communicate through reactive values
**🔄 Reusability**: Modules can be used in different contexts
**📊 Consistent Interface**: All modules follow same pattern

### Module Structure Template

```r
# UI Module (mod_example.R)
#' example UI Function
#'
#' @description A shiny Module for example analysis.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_example_ui <- function(id){
  ns <- NS(id)
  tagList(
    # UI elements here
  )
}

# Server Module (mod_example_server.R)
#' example Server Functions
#'
#' @noRd 
mod_example_server <- function(id, seurat_obj, processed, main_values){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    
    # Local reactive values
    local_values <- reactiveValues()
    
    # Analysis logic
    observeEvent(input$run_analysis, {
      # Perform analysis
      # Update local_values and main_values
    })
    
    # Outputs
    output$plot <- renderPlot({
      # Generate plots
    })
    
    # Return values (optional)
    return(list(
      results = reactive(local_values$results)
    ))
  })
}
```

### Module Communication

**Data Flow Architecture**:

```
Main App
├── Global reactive values (main_values)
│   ├── seurat_obj: Current Seurat object
│   ├── processed: Processing status
│   ├── markers: Marker gene results
│   ├── cancersea_scores: Pathway scores
│   └── trajectory: Trajectory results
│
├── Module A (Data Upload)
│   └── Updates: main_values$seurat_obj
│
├── Module B (Clustering)
│   ├── Reads: main_values$seurat_obj
│   └── Updates: main_values$seurat_obj (with clusters)
│
├── Module C (CancerSEA)
│   ├── Reads: main_values$seurat_obj
│   └── Updates: main_values$cancersea_scores
│
└── Module D (Export)
    └── Reads: All main_values for export
```

**Reactive Value Structure**:
```r
# Main reactive values object
main_values <- reactiveValues(
  # Core data
  seurat_obj = NULL,
  processed = FALSE,
  
  # Analysis results
  markers = NULL,
  cancersea_scores = list(),
  trajectory = NULL,
  
  # UI state
  current_reduction = "umap",
  selected_clusters = NULL
)
```

---

## 📊 Data Management

### Data Flow Pipeline

**1. Data Input**:
```r
# File upload handling
observeEvent(input$file_upload, {
  # Validate file format
  validate_input_file(input$file_upload$datapath)
  
  # Load data based on type
  if (is_10x_data(file_path)) {
    seurat_obj <- load_10x_data(file_path)
  } else if (is_seurat_object(file_path)) {
    seurat_obj <- readRDS(file_path)
  }
  
  # Store in reactive values
  main_values$seurat_obj <- seurat_obj
})
```

**2. Data Processing**:
```r
# Quality control pipeline
run_quality_control <- function(seurat_obj, min_features, max_features, max_mt) {
  # Filter cells and genes
  seurat_obj <- subset(seurat_obj, 
                       subset = nFeature_RNA > min_features & 
                               nFeature_RNA < max_features & 
                               percent.mt < max_mt)
  
  # Normalize data
  seurat_obj <- NormalizeData(seurat_obj)
  
  # Find variable features
  seurat_obj <- FindVariableFeatures(seurat_obj)
  
  return(seurat_obj)
}
```

**3. Data Validation**:
```r
# Validate Seurat object
validate_seurat_object <- function(seurat_obj) {
  # Check object class
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Object is not a Seurat object")
  }
  
  # Check minimum requirements
  if (ncol(seurat_obj) < 10) {
    stop("Too few cells (minimum 10)")
  }
  
  if (nrow(seurat_obj) < 100) {
    stop("Too few genes (minimum 100)")
  }
  
  return(TRUE)
}
```

### Memory Management

**Efficient Data Handling**:
```r
# Use sparse matrices for large datasets
if (ncol(seurat_obj) > 10000) {
  # Convert to sparse format if not already
  if (!inherits(seurat_obj@assays$RNA@counts, "dgCMatrix")) {
    seurat_obj@assays$RNA@counts <- as(seurat_obj@assays$RNA@counts, "dgCMatrix")
  }
}

# Garbage collection for large operations
gc(verbose = FALSE)
```

**Progress Reporting**:
```r
# Progress bars for long operations
withProgress(message = 'Running analysis...', value = 0, {
  incProgress(0.2, detail = "Loading data")
  # Step 1
  
  incProgress(0.4, detail = "Processing")
  # Step 2
  
  incProgress(0.8, detail = "Finalizing")
  # Step 3
  
  incProgress(1, detail = "Complete")
})
```

---

## 🔧 Configuration Management

### Application Configuration

**Configuration File (R/app_config.R)**:
```r
# Application configuration
app_config <- list(
  # Default analysis parameters
  default_min_features = 200,
  default_max_features = 5000,
  default_max_mt = 20,
  default_resolution = 0.5,
  default_marker_logfc = 0.25,
  default_marker_minpct = 0.1,
  default_marker_topn = 10,
  default_marker_display_genes = 50,
  
  # File upload limits
  max_file_size_mb = 2000,
  supported_formats = c(".rds", ".h5", ".mtx", ".csv", ".tsv"),
  
  # UI configuration
  theme = "blue",
  sidebar_width = 300,
  plot_height = "600px",
  
  # Performance settings
  max_cells_for_fast_analysis = 10000,
  default_pca_dims = 10,
  trajectory_max_cells = 20000
)
```

**Environment-Specific Settings**:
```r
# Development vs production settings
if (Sys.getenv("MASIH_ENV") == "production") {
  app_config$max_file_size_mb <- 500  # Smaller limit for production
  app_config$enable_debug <- FALSE
} else {
  app_config$enable_debug <- TRUE
  app_config$verbose_logging <- TRUE
}
```

### Golem Configuration

**Golem Options**:
```r
# In run_app.R
run_app <- function(
  onStart = NULL,
  options = list(), 
  enableBookmarking = NULL,
  ...
) {
  with_golem_options(
    app = shinyApp(
      ui = app_ui,
      server = app_server,
      onStart = onStart,
      options = options, 
      enableBookmarking = enableBookmarking
    ), 
    golem_opts = list(...)
  )
}
```

---

## 🎨 User Interface Architecture

### Dashboard Layout

**Main UI Structure (app_ui.R)**:
```r
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    
    # Main dashboard
    dashboardPage(
      # Dashboard header
      dashboardHeader(title = "MASIH"),
      
      # Sidebar with navigation
      dashboardSidebar(
        sidebarMenu(
          id = "sidebar",
          menuItem("Data Upload", tabName = "data_upload", icon = icon("upload")),
          menuItem("Quality Control", tabName = "qc", icon = icon("filter")),
          menuItem("Cluster Analysis", tabName = "clustering", icon = icon("project-diagram")),
          menuItem("Marker Genes", tabName = "markers", icon = icon("dna")),
          menuItem("Cell Cycle", tabName = "cellcycle", icon = icon("sync-alt")),
          menuItem("CancerSEA", tabName = "cancersea", icon = icon("brain")),
          menuItem("Pathway Comparison", tabName = "compare", icon = icon("balance-scale")),
          menuItem("Trajectory Analysis", tabName = "trajectory", icon = icon("route")),
          menuItem("Export", tabName = "export", icon = icon("download"))
        )
      ),
      
      # Main content area
      dashboardBody(
        # Custom CSS
        tags$head(
          tags$link(rel = "stylesheet", type = "text/css", href = "www/custom.css")
        ),
        
        # Tab items
        tabItems(
          tabItem(tabName = "data_upload", mod_data_upload_ui("data_upload_1")),
          tabItem(tabName = "qc", mod_quality_control_ui("qc_1")),
          tabItem(tabName = "clustering", mod_cluster_analysis_ui("cluster_1")),
          tabItem(tabName = "markers", mod_markers_ui("markers_1")),
          tabItem(tabName = "cellcycle", mod_cell_cycle_ui("cellcycle_1")),
          tabItem(tabName = "cancersea", mod_cancersea_ui("cancersea_1")),
          tabItem(tabName = "compare", mod_compare_ui("compare_1")),
          tabItem(tabName = "trajectory", mod_trajectory_ui("trajectory_1")),
          tabItem(tabName = "export", mod_export_ui("export_1"))
        )
      )
    )
  )
}
```

### Responsive Design

**CSS Framework**:
```css
/* Custom CSS (inst/app/www/custom.css) */
.content-wrapper, .right-side {
  background-color: #f8f9fa;
}

.box {
  border-radius: 8px;
  box-shadow: 0 2px 4px rgba(0,0,0,0.1);
}

.btn-primary {
  background-color: #007bff;
  border-color: #007bff;
}

/* Responsive adjustments */
@media (max-width: 768px) {
  .sidebar-toggle {
    display: block;
  }
  
  .main-sidebar {
    width: 100%;
  }
}
```

---

## ⚡ Performance Optimization

### Reactive Programming Best Practices

**Efficient Reactive Structure**:
```r
# Use reactive expressions for expensive computations
processed_data <- reactive({
  req(main_values$seurat_obj)
  
  # Cache expensive operations
  if (is.null(main_values$processed_cache) || 
      main_values$cache_timestamp < main_values$data_timestamp) {
    
    result <- expensive_processing(main_values$seurat_obj)
    main_values$processed_cache <- result
    main_values$cache_timestamp <- Sys.time()
    
    return(result)
  } else {
    return(main_values$processed_cache)
  }
})

# Use observe() for side effects, reactive() for computations
observe({
  req(input$update_analysis)
  
  # Side effect: update main_values
  main_values$analysis_updated <- TRUE
  
  # Don't return anything
})
```

**Debouncing User Input**:
```r
# Debounce rapid user input
parameter_debounced <- reactive({
  input$parameter_slider
}) %>% debounce(1000)  # Wait 1 second after last change

observeEvent(parameter_debounced(), {
  # Expensive analysis only after user stops changing parameter
  run_analysis(parameter_debounced())
})
```

### Memory Management

**Large Dataset Handling**:
```r
# Subsample for interactive exploration
create_subset_for_ui <- function(seurat_obj, max_cells = 5000) {
  if (ncol(seurat_obj) > max_cells) {
    # Stratified sampling to preserve cluster structure
    if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
      subset_cells <- stratified_sample(seurat_obj, max_cells)
    } else {
      subset_cells <- sample(colnames(seurat_obj), max_cells)
    }
    
    return(seurat_obj[, subset_cells])
  } else {
    return(seurat_obj)
  }
}

# Use subset for UI, full data for analysis
ui_data <- reactive({
  req(main_values$seurat_obj)
  create_subset_for_ui(main_values$seurat_obj)
})

analysis_data <- reactive({
  req(main_values$seurat_obj)
  main_values$seurat_obj  # Full dataset
})
```

**Progress Reporting for Long Operations**:
```r
# Async processing with progress
process_with_progress <- function(data, steps) {
  withProgress(message = 'Processing...', value = 0, {
    
    n_steps <- length(steps)
    
    for (i in seq_along(steps)) {
      incProgress(1/n_steps, detail = steps[[i]]$description)
      
      # Execute step
      data <- steps[[i]]$function(data)
      
      # Allow UI updates
      Sys.sleep(0.1)
    }
    
    return(data)
  })
}
```

---

## 🧪 Testing Architecture

### Test Organization

**Test Structure**:
```
tests/
├── testthat.R                 # Test configuration
└── testthat/
    ├── helper-setup.R          # Test setup functions
    ├── test-data-loading.R     # Data I/O tests
    ├── test-quality-control.R  # QC function tests
    ├── test-clustering.R       # Clustering tests
    ├── test-markers.R          # Marker gene tests
    ├── test-cancersea.R        # CancerSEA tests
    ├── test-trajectory.R       # Trajectory tests
    ├── test-export.R           # Export function tests
    └── test-ui-modules.R       # UI module tests
```

### Unit Testing

**Test Setup (helper-setup.R)**:
```r
# Create test Seurat object
create_test_seurat <- function(n_cells = 100, n_genes = 1000) {
  # Generate synthetic data
  counts <- matrix(rpois(n_cells * n_genes, lambda = 2), 
                   nrow = n_genes, ncol = n_cells)
  
  # Add gene and cell names
  rownames(counts) <- paste0("Gene_", 1:n_genes)
  colnames(counts) <- paste0("Cell_", 1:n_cells)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts)
  
  # Add basic metadata
  seurat_obj$percent.mt <- runif(n_cells, 0, 25)
  
  return(seurat_obj)
}

# Test data paths
get_test_data_path <- function(filename) {
  system.file("testdata", filename, package = "masih")
}
```

**Function Testing Example**:
```r
# Test quality control functions
test_that("quality control filters work correctly", {
  # Setup
  test_obj <- create_test_seurat(n_cells = 200)
  
  # Add extreme values for testing
  test_obj$nFeature_RNA[1:10] <- 50    # Low feature cells
  test_obj$percent.mt[11:20] <- 30     # High MT cells
  
  # Test filtering
  filtered_obj <- filter_cells(test_obj, 
                              min_features = 100, 
                              max_mt = 25)
  
  # Assertions
  expect_s4_class(filtered_obj, "Seurat")
  expect_lt(ncol(filtered_obj), ncol(test_obj))  # Some cells removed
  expect_true(all(filtered_obj$nFeature_RNA >= 100))
  expect_true(all(filtered_obj$percent.mt <= 25))
})
```

### Integration Testing

**Module Integration Tests**:
```r
test_that("full analysis workflow completes successfully", {
  # Load test data
  test_data <- create_test_seurat(n_cells = 500)
  
  # Run complete workflow
  result <- run_complete_workflow(test_data)
  
  # Check workflow completion
  expect_s4_class(result$seurat_obj, "Seurat")
  expect_true(result$qc_completed)
  expect_true(result$clustering_completed)
  expect_true("seurat_clusters" %in% colnames(result$seurat_obj@meta.data))
})
```

### UI Testing

**Shiny Module Testing**:
```r
test_that("clustering module UI generates correctly", {
  # Test UI generation
  ui_output <- mod_cluster_analysis_ui("test")
  
  # Check UI components
  expect_true(inherits(ui_output, "shiny.tag.list"))
  expect_true(any(grepl("reduction_cluster", as.character(ui_output))))
})
```

---

## 📦 Deployment Architecture

### Package Structure

**DESCRIPTION File**:
```r
Package: masih
Title: Modular Analysis Shiny Interface for Heterogeneity
Version: 0.1.0
Authors@R: 
    person("Your", "Name", 
           email = "your.email@institution.edu", 
           role = c("aut", "cre"),
           comment = c(ORCID = "YOUR-ORCID-ID"))
Description: An interactive R Shiny application for comprehensive single-cell 
    RNA sequencing analysis in cancer research.
License: MIT + file LICENSE
Encoding: UTF-8
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.2.3
Depends: 
    R (>= 4.0.0)
Imports:
    shiny (>= 1.7.0),
    shinydashboard,
    golem (>= 0.4.0),
    DT,
    plotly,
    Seurat (>= 4.0.0),
    SingleCellExperiment,
    slingshot,
    dplyr,
    ggplot2,
    viridis,
    RColorBrewer,
    corrplot,
    openxlsx,
    config (>= 0.3.1)
Suggests:
    cancersea,
    testthat (>= 3.0.0),
    knitr,
    rmarkdown
VignetteBuilder: knitr
```

### Deployment Options

**Local Installation**:
```r
# Development installation
devtools::install_github("yourusername/masih")

# Production installation
install.packages("masih")  # When on CRAN
```

**Server Deployment**:
```r
# Shiny Server deployment
# Copy app files to /srv/shiny-server/masih/

# RStudio Connect deployment
rsconnect::deployApp(
  appDir = ".",
  appName = "masih",
  server = "your-connect-server.com"
)

# Docker deployment
# Use provided Dockerfile
```

**Docker Configuration**:
```dockerfile
FROM rocker/shiny-verse:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('BiocManager', 'devtools'))"
RUN R -e "BiocManager::install(c('Seurat', 'SingleCellExperiment', 'slingshot'))"

# Copy app
COPY . /srv/shiny-server/masih/

# Expose port
EXPOSE 3838

# Run app
CMD ["/usr/bin/shiny-server"]
```

---

## 🔒 Security Considerations

### Input Validation

**File Upload Security**:
```r
# Validate file types
validate_file_upload <- function(file_path) {
  # Check file extension
  ext <- tools::file_ext(file_path)
  allowed_extensions <- c("rds", "csv", "tsv", "h5", "mtx")
  
  if (!ext %in% allowed_extensions) {
    stop("File type not allowed: ", ext)
  }
  
  # Check file size
  file_size <- file.info(file_path)$size
  max_size <- 2 * 1024^3  # 2GB
  
  if (file_size > max_size) {
    stop("File too large: ", file_size, " bytes")
  }
  
  # Check file content
  if (ext == "rds") {
    tryCatch({
      obj <- readRDS(file_path)
      if (!inherits(obj, c("Seurat", "SingleCellExperiment"))) {
        stop("Invalid object type in RDS file")
      }
    }, error = function(e) {
      stop("Invalid RDS file: ", e$message)
    })
  }
  
  return(TRUE)
}
```

**Parameter Validation**:
```r
# Validate analysis parameters
validate_parameters <- function(params) {
  # Check numeric ranges
  if (params$min_features < 0 || params$min_features > 10000) {
    stop("min_features must be between 0 and 10000")
  }
  
  if (params$max_mt < 0 || params$max_mt > 100) {
    stop("max_mt must be between 0 and 100")
  }
  
  # Check logical parameters
  if (!is.logical(params$only_pos)) {
    stop("only_pos must be TRUE or FALSE")
  }
  
  return(TRUE)
}
```

### Data Privacy

**Temporary File Management**:
```r
# Clean up temporary files
cleanup_temp_files <- function(session) {
  # Remove uploaded files after processing
  if (exists("temp_files", envir = session$userData)) {
    files_to_remove <- session$userData$temp_files
    file.remove(files_to_remove)
    rm("temp_files", envir = session$userData)
  }
}

# Register cleanup on session end
onSessionEnded(function() {
  cleanup_temp_files(session)
})
```

---

## 📈 Performance Monitoring

### Logging

**Application Logging**:
```r
# Simple logging system
log_message <- function(level, message, details = NULL) {
  timestamp <- Sys.time()
  log_entry <- paste0(
    "[", timestamp, "] ",
    "[", level, "] ",
    message
  )
  
  if (!is.null(details)) {
    log_entry <- paste0(log_entry, " | Details: ", details)
  }
  
  # Log to file in production
  if (app_config$enable_logging) {
    cat(log_entry, "\n", file = "masih.log", append = TRUE)
  }
  
  # Print to console in development
  if (app_config$enable_debug) {
    cat(log_entry, "\n")
  }
}

# Usage
log_message("INFO", "User uploaded file", "cells: 1000, genes: 2000")
log_message("ERROR", "Analysis failed", input$error_message)
```

### Performance Metrics

**Analysis Timing**:
```r
# Time analysis steps
time_analysis <- function(analysis_function, ...) {
  start_time <- Sys.time()
  
  result <- tryCatch({
    analysis_function(...)
  }, error = function(e) {
    log_message("ERROR", "Analysis failed", e$message)
    stop(e)
  })
  
  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  log_message("PERFORMANCE", 
              paste("Analysis completed in", round(duration, 2), "seconds"))
  
  return(result)
}
```

---

## 🔧 Development Tools

### Code Quality Tools

**Linting and Styling**:
```r
# Check code style
lintr::lint_package()

# Auto-format code
styler::style_pkg()

# Check package structure
devtools::check()
```

**Documentation Generation**:
```r
# Generate roxygen documentation
roxygen2::roxygenise()

# Build package website
pkgdown::build_site()

# Create vignettes
usethis::use_vignette("basic-usage")
```

### Continuous Integration

**GitHub Actions Workflow**:
```yaml
# .github/workflows/R-CMD-check.yaml
name: R-CMD-check

on: [push, pull_request]

jobs:
  R-CMD-check:
    runs-on: ubuntu-latest
    
    steps:
    - uses: actions/checkout@v2
    
    - uses: r-lib/actions/setup-r@v2
      with:
        r-version: '4.1.0'
    
    - name: Install dependencies
      run: |
        install.packages(c("remotes", "rcmdcheck"))
        remotes::install_deps(dependencies = TRUE)
      shell: Rscript {0}
    
    - name: Check package
      run: rcmdcheck::rcmdcheck(args = "--no-manual", error_on = "error")
      shell: Rscript {0}
```

---

## 🔮 Future Architecture Considerations

### Scalability Plans

**Microservices Architecture**:
- **Analysis Service**: Heavy computational tasks
- **Data Service**: File management and storage
- **UI Service**: User interface and visualization
- **API Service**: Programmatic access

**Cloud Integration**:
- **AWS/GCP Integration**: Large-scale processing
- **Container Orchestration**: Kubernetes deployment
- **Database Integration**: PostgreSQL for metadata
- **CDN Integration**: Fast static asset delivery

### Technology Evolution

**Planned Upgrades**:
- **WebAssembly**: Client-side analysis for small datasets
- **WebGL**: Advanced 3D visualizations
- **Progressive Web App**: Offline functionality
- **Real-time Collaboration**: Multi-user analysis sessions

---

## 🛠️ Extension Framework

### Plugin Architecture

**Module Registration System**:
```r
# Register new analysis modules
register_analysis_module <- function(module_name, ui_function, server_function) {
  # Add to global module registry
  app_config$modules[[module_name]] <- list(
    ui = ui_function,
    server = server_function,
    dependencies = get_module_dependencies(module_name)
  )
  
  # Update navigation menu
  update_navigation_menu(module_name)
}

# Example plugin registration
register_analysis_module(
  "cell_communication",
  mod_cell_communication_ui,
  mod_cell_communication_server
)
```

**Pathway Database Extensions**:
```r
# Plugin system for new pathway databases
register_pathway_database <- function(db_name, db_loader, db_metadata) {
  app_config$pathway_databases[[db_name]] <- list(
    loader = db_loader,
    metadata = db_metadata,
    version = db_metadata$version
  )
}

# Example: Register Reactome pathways
register_pathway_database(
  "Reactome",
  load_reactome_pathways,
  list(
    description = "Reactome pathway database",
    version = "2023.1",
    n_pathways = 2500
  )
)
```

### API Framework

**RESTful API Structure**:
```r
# API endpoints for programmatic access
#* @apiTitle MASIH API
#* @apiDescription Programmatic access to MASIH functionality

#* Upload and analyze data
#* @param file_path Path to input data
#* @param analysis_params Analysis parameters
#* @post /api/analyze
api_analyze_data <- function(file_path, analysis_params) {
  # Validate inputs
  validate_api_inputs(file_path, analysis_params)
  
  # Run analysis
  results <- run_masih_analysis(file_path, analysis_params)
  
  # Return results
  return(list(
    status = "success",
    results = results,
    timestamp = Sys.time()
  ))
}

#* Get analysis results
#* @param analysis_id Unique analysis identifier
#* @get /api/results/<analysis_id>
api_get_results <- function(analysis_id) {
  results <- load_analysis_results(analysis_id)
  return(results)
}
```

---

## 🔍 Error Handling and Debugging

### Comprehensive Error Management

**Error Handling Strategy**:
```r
# Centralized error handling
handle_analysis_error <- function(error, context) {
  # Log error with context
  log_message("ERROR", paste("Analysis failed:", error$message), 
              paste("Context:", context))
  
  # User-friendly error messages
  user_message <- switch(
    error$type,
    "data_format" = "Data format not supported. Please check your file.",
    "memory" = "Insufficient memory. Try with a smaller dataset.",
    "computation" = "Analysis failed. Please check your parameters.",
    "Unknown error occurred. Please contact support."
  )
  
  # Show notification to user
  showNotification(
    user_message,
    type = "error",
    duration = 10
  )
  
  # Return safe state
  return(get_safe_state())
}

# Wrap analysis functions with error handling
safe_analysis <- function(analysis_function, ..., context = "analysis") {
  tryCatch({
    result <- analysis_function(...)
    return(result)
  }, error = function(e) {
    handle_analysis_error(e, context)
    return(NULL)
  })
}
```

**Debug Mode**:
```r
# Debug information collection
collect_debug_info <- function() {
  debug_info <- list(
    session_info = sessionInfo(),
    system_info = Sys.info(),
    memory_usage = pryr::mem_used(),
    app_config = app_config,
    timestamp = Sys.time()
  )
  
  if (app_config$enable_debug) {
    # Include additional debug information
    debug_info$loaded_packages <- .packages()
    debug_info$search_path <- search()
  }
  
  return(debug_info)
}
```

---

## 🌐 Internationalization Support

### Multi-language Framework

**Translation System**:
```r
# Translation function
translate <- function(key, language = app_config$default_language) {
  translations <- load_translations(language)
  
  if (key %in% names(translations)) {
    return(translations[[key]])
  } else {
    # Fallback to English
    if (language != "en") {
      return(translate(key, "en"))
    } else {
      return(key)  # Return key if no translation found
    }
  }
}

# Usage in UI
h4(translate("cluster_analysis_title"))
actionButton("run_analysis", translate("run_analysis_button"))
```

**Translation Files**:
```json
// inst/translations/en.json
{
  "cluster_analysis_title": "Cluster Analysis",
  "run_analysis_button": "Run Analysis",
  "quality_control_title": "Quality Control",
  "data_upload_title": "Data Upload"
}

// inst/translations/es.json
{
  "cluster_analysis_title": "Análisis de Clusters",
  "run_analysis_button": "Ejecutar Análisis",
  "quality_control_title": "Control de Calidad",
  "data_upload_title": "Cargar Datos"
}
```

---

## 📱 Mobile and Touch Support

### Responsive Design Implementation

**Touch-Friendly Interface**:
```css
/* Touch-friendly buttons and controls */
.btn {
  min-height: 44px;  /* Minimum touch target size */
  min-width: 44px;
  padding: 12px 16px;
}

/* Mobile-optimized plots */
@media (max-width: 768px) {
  .plotly-container {
    height: 300px !important;
  }
  
  .dataTables_wrapper {
    font-size: 14px;
  }
}

/* Swipe gestures for plot navigation */
.plot-container {
  touch-action: pan-x pan-y;
}
```

**Mobile Navigation**:
```r
# Adaptive navigation for mobile
mobile_sidebar <- function() {
  if (input$screen_width < 768) {
    # Collapsed sidebar for mobile
    sidebarMenu(
      menuItem("Upload", tabName = "upload", icon = icon("upload")),
      menuItem("Analysis", tabName = "analysis", icon = icon("chart-bar")),
      menuItem("Results", tabName = "results", icon = icon("download"))
    )
  } else {
    # Full sidebar for desktop
    full_sidebar_menu()
  }
}
```

---

## 🔐 Authentication and Authorization

### User Management System

**Authentication Framework**:
```r
# User authentication module
authenticate_user <- function(username, password) {
  # Hash password for comparison
  password_hash <- digest::digest(password, algo = "sha256")
  
  # Check against user database
  user_record <- get_user_record(username)
  
  if (!is.null(user_record) && user_record$password_hash == password_hash) {
    # Create session token
    session_token <- generate_session_token()
    store_session_token(username, session_token)
    
    return(list(
      success = TRUE,
      token = session_token,
      user_info = user_record
    ))
  } else {
    return(list(success = FALSE))
  }
}

# Authorization levels
check_user_permissions <- function(user, required_permission) {
  user_permissions <- get_user_permissions(user$username)
  return(required_permission %in% user_permissions)
}
```

**Role-Based Access Control**:
```r
# Define user roles and permissions
user_roles <- list(
  "student" = c("upload_data", "basic_analysis", "export_results"),
  "researcher" = c("upload_data", "advanced_analysis", "export_results", "save_projects"),
  "admin" = c("all_permissions", "user_management", "system_config")
)

# Check permissions for specific actions
require_permission <- function(permission) {
  if (!check_user_permissions(current_user(), permission)) {
    showNotification("Insufficient permissions", type = "error")
    return(FALSE)
  }
  return(TRUE)
}
```

---

## 💾 Data Persistence and State Management

### Session State Management

**Persistent Sessions**:
```r
# Save analysis state
save_analysis_state <- function(session_id, state) {
  state_file <- file.path(app_config$state_dir, paste0(session_id, ".rds"))
  saveRDS(state, state_file)
}

# Restore analysis state
restore_analysis_state <- function(session_id) {
  state_file <- file.path(app_config$state_dir, paste0(session_id, ".rds"))
  if (file.exists(state_file)) {
    return(readRDS(state_file))
  } else {
    return(get_default_state())
  }
}

# Auto-save functionality
observe({
  # Auto-save every 5 minutes
  invalidateLater(300000)  # 5 minutes in milliseconds
  
  if (!is.null(main_values$seurat_obj)) {
    save_analysis_state(session$token, main_values)
  }
})
```

**Project Management**:
```r
# Project structure
create_project <- function(project_name, user_id) {
  project_id <- generate_uuid()
  
  project_data <- list(
    id = project_id,
    name = project_name,
    user_id = user_id,
    created_date = Sys.time(),
    last_modified = Sys.time(),
    analyses = list(),
    shared_with = c()
  )
  
  save_project(project_data)
  return(project_id)
}

# Save analysis to project
add_analysis_to_project <- function(project_id, analysis_data) {
  project <- load_project(project_id)
  analysis_id <- generate_uuid()
  
  project$analyses[[analysis_id]] <- analysis_data
  project$last_modified <- Sys.time()
  
  save_project(project)
  return(analysis_id)
}
```

---

## 🎨 Theming and Customization

### Dynamic Theme System

**Theme Management**:
```r
# Theme configuration
app_themes <- list(
  "default" = list(
    primary_color = "#007bff",
    secondary_color = "#6c757d",
    background_color = "#ffffff",
    text_color = "#212529"
  ),
  "dark" = list(
    primary_color = "#0d6efd",
    secondary_color = "#6c757d",
    background_color = "#212529",
    text_color = "#ffffff"
  ),
  "cancer_research" = list(
    primary_color = "#dc3545",
    secondary_color = "#28a745",
    background_color = "#f8f9fa",
    text_color = "#495057"
  )
)

# Apply theme
apply_theme <- function(theme_name) {
  theme <- app_themes[[theme_name]]
  
  css_vars <- paste0(
    ":root {",
    paste0("--", names(theme), ": ", theme, ";", collapse = " "),
    "}"
  )
  
  tags$style(css_vars)
}
```

**Institution Branding**:
```r
# Institutional customization
apply_institution_branding <- function(institution_config) {
  # Custom logo
  if (!is.null(institution_config$logo_url)) {
    update_logo(institution_config$logo_url)
  }
  
  # Custom colors
  if (!is.null(institution_config$primary_color)) {
    update_primary_color(institution_config$primary_color)
  }
  
  # Custom footer
  if (!is.null(institution_config$footer_text)) {
    update_footer(institution_config$footer_text)
  }
}
```

---

## 📊 Analytics and Usage Tracking

### Usage Analytics Framework

**Anonymous Usage Tracking**:
```r
# Track user interactions (anonymized)
track_user_action <- function(action, details = NULL) {
  if (app_config$enable_analytics) {
    event_data <- list(
      timestamp = Sys.time(),
      action = action,
      details = details,
      session_id = digest::digest(session$token),  # Anonymized
      user_agent = session$clientData$user_agent
    )
    
    log_analytics_event(event_data)
  }
}

# Usage examples
track_user_action("data_upload", list(file_size = file.info(file_path)$size))
track_user_action("analysis_complete", list(analysis_type = "clustering"))
track_user_action("export_plot", list(format = "pdf"))
```

**Performance Monitoring**:
```r
# Monitor application performance
monitor_performance <- function() {
  performance_data <- list(
    timestamp = Sys.time(),
    memory_usage = pryr::mem_used(),
    active_sessions = length(session_list),
    cpu_usage = system("ps -o %cpu -p $PPID --no-headers", intern = TRUE),
    response_time = measure_response_time()
  )
  
  log_performance_data(performance_data)
  
  # Alert if performance degrades
  if (performance_data$memory_usage > app_config$memory_alert_threshold) {
    send_performance_alert("High memory usage detected")
  }
}
```

---

## 🔧 Configuration and Environment Management

### Environment-Specific Configuration

**Configuration Hierarchy**:
```r
# Load configuration based on environment
load_app_config <- function() {
  env <- Sys.getenv("MASIH_ENV", "development")
  
  # Base configuration
  base_config <- load_base_config()
  
  # Environment-specific overrides
  env_config <- switch(env,
    "development" = load_dev_config(),
    "testing" = load_test_config(),
    "production" = load_prod_config(),
    list()
  )
  
  # Merge configurations
  final_config <- merge_configs(base_config, env_config)
  
  # Validate configuration
  validate_config(final_config)
  
  return(final_config)
}

# Environment-specific settings
load_prod_config <- function() {
  list(
    max_file_size_mb = 1000,
    enable_debug = FALSE,
    enable_analytics = TRUE,
    database_url = Sys.getenv("DATABASE_URL"),
    redis_url = Sys.getenv("REDIS_URL")
  )
}
```

---

## 🔄 Version Management and Migrations

### Data Migration Framework

**Version Migration System**:
```r
# Handle data format migrations between versions
migrate_data <- function(data, from_version, to_version) {
  migration_path <- get_migration_path(from_version, to_version)
  
  for (migration in migration_path) {
    data <- migration$function(data)
    log_message("INFO", paste("Applied migration:", migration$name))
  }
  
  return(data)
}

# Example migration for Seurat object updates
migration_seurat_v4_to_v5 <- function(seurat_obj) {
  if (inherits(seurat_obj, "Seurat") && seurat_obj@version < "5.0.0") {
    # Update object structure for Seurat v5
    seurat_obj <- UpdateSeuratObject(seurat_obj)
    log_message("INFO", "Updated Seurat object to v5 format")
  }
  return(seurat_obj)
}
```

**Backward Compatibility**:
```r
# Maintain compatibility with older analysis results
load_legacy_analysis <- function(file_path) {
  analysis_data <- readRDS(file_path)
  
  # Detect version
  version <- detect_analysis_version(analysis_data)
  
  # Apply necessary migrations
  if (version < current_version()) {
    analysis_data <- migrate_analysis_data(analysis_data, version)
  }
  
  return(analysis_data)
}
```

---

**This architecture supports MASIH's mission** 🎯 **of making advanced single-cell analysis accessible while maintaining scientific rigor and extensibility.**
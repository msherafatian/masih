# Contributing to MASIH

**Welcome contributors! Help us make single-cell cancer analysis accessible to all researchers.**

---

## 🎯 Ways to Contribute

We welcome contributions of all kinds:

### 🐛 Bug Reports and Feature Requests
- Report bugs through [GitHub Issues](https://github.com/yourusername/masih/issues)
- Request new features or enhancements
- Suggest improvements to documentation

### 📖 Documentation
- Improve existing documentation
- Add new tutorials or examples
- Translate documentation
- Create video tutorials

### 🧪 Code Contributions
- Fix bugs and add features
- Improve performance
- Add new analysis modules
- Enhance user interface

### 🧬 Scientific Contributions
- Add new pathway databases
- Validate methods on new datasets
- Contribute example datasets
- Improve biological interpretations

### 🎨 Design and UX
- Improve user interface design
- Enhance user experience
- Create logos and graphics
- Improve accessibility

---

## 🚀 Getting Started

### Prerequisites

**Development Environment**:
- **R** ≥ 4.0.0
- **RStudio** (recommended)
- **Git** for version control
- **GitHub account**

**Required R Packages**:
```r
# Development tools
install.packages(c("devtools", "roxygen2", "testthat", "usethis"))

# MASIH dependencies
install.packages(c("shiny", "shinydashboard", "golem"))
BiocManager::install(c("Seurat", "SingleCellExperiment", "slingshot"))
```

### Development Setup

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/yourusername/masih.git
   cd masih
   ```

3. **Set up development environment**:
   ```r
   # In RStudio, open masih.Rproj
   # Install development dependencies
   devtools::install_deps(dependencies = TRUE)
   
   # Load package for development
   devtools::load_all()
   ```

4. **Create a new branch** for your feature:
   ```bash
   git checkout -b feature/your-feature-name
   ```

---

## 📋 Development Workflow

### Code Style Guidelines

**R Code Style**:
```r
# Follow tidyverse style guide
# Use descriptive variable names
calculate_pathway_score <- function(seurat_obj, pathway_name) {
  # Good: descriptive function name
  # Good: snake_case for variables
  
  gene_list <- get_pathway_genes(pathway_name)
  filtered_genes <- gene_list[gene_list %in% rownames(seurat_obj)]
  
  # Good: clear variable names and logical flow
  return(filtered_genes)
}
```

**Shiny Module Structure**:
```r
# UI Function
mod_analysis_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # UI elements with consistent naming
    selectInput(ns("parameter"), "Parameter:", choices = NULL),
    actionButton(ns("run_analysis"), "Run Analysis"),
    plotOutput(ns("result_plot"))
  )
}

# Server Function  
mod_analysis_server <- function(id, shared_data) {
  moduleServer(id, function(input, output, session) {
    # Server logic with clear separation of concerns
  })
}
```

**Documentation Standards**:
```r
#' Calculate CancerSEA pathway scores
#'
#' @description Calculate functional state scores for cancer cells using
#'   CancerSEA gene signatures.
#'
#' @param seurat_obj Seurat object containing single-cell data
#' @param pathway_name Character string specifying CancerSEA pathway
#' @param ctrl Number of control features for scoring (default: 100)
#'
#' @return Seurat object with pathway scores added to metadata
#'
#' @examples
#' \dontrun{
#' seurat_obj <- calculate_cancersea_score(seurat_obj, "Stemness")
#' }
#'
#' @export
calculate_cancersea_score <- function(seurat_obj, pathway_name, ctrl = 100) {
  # Function implementation
}
```

### Testing Guidelines

**Unit Tests**:
```r
# tests/testthat/test-pathway-scoring.R
library(testthat)
library(masih)

test_that("pathway scoring works correctly", {
  # Create test data
  test_data <- create_test_seurat_object()
  
  # Test pathway scoring
  result <- calculate_cancersea_score(test_data, "Stemness")
  
  # Assertions
  expect_s4_class(result, "Seurat")
  expect_true("Stemness_1" %in% colnames(result@meta.data))
  expect_true(all(!is.na(result$Stemness_1)))
})
```

**Integration Tests**:
```r
test_that("full analysis workflow works", {
  # Test complete workflow
  app_data <- load_example_data()
  processed_data <- run_quality_control(app_data)
  clustered_data <- run_clustering(processed_data)
  
  expect_s4_class(clustered_data, "Seurat")
  expect_true("seurat_clusters" %in% colnames(clustered_data@meta.data))
})
```

---

## 🧪 Adding New Features

### Adding a New Analysis Module

**1. Create Module Files**:
```
R/
├── mod_new_analysis.R          # UI function
├── mod_new_analysis_server.R   # Server function  
└── fct_new_analysis.R          # Helper functions
```

**2. UI Module Template**:
```r
#' new_analysis UI Function
#'
#' @description A shiny Module for new analysis functionality.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_new_analysis_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      box(
        title = "New Analysis",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        # Your UI elements here
      )
    )
  )
}
```

**3. Server Module Template**:
```r
#' new_analysis Server Functions
#'
#' @noRd 
mod_new_analysis_server <- function(id, seurat_obj, processed, main_values){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    
    # Your server logic here
    
  })
}
```

**4. Integration Steps**:
```r
# In app_ui.R - Add menu item
menuItem("New Analysis", tabName = "new_analysis", icon = icon("chart-line"))

# In app_ui.R - Add tab content  
tabItem(tabName = "new_analysis",
        mod_new_analysis_ui("new_analysis_1"))

# In app_server.R - Call server module
mod_new_analysis_server("new_analysis_1", 
                       reactive(values$seurat_obj),
                       reactive(values$processed),
                       values)
```

### Adding New Pathway Databases

**1. Create Database Interface**:
```r
#' Load custom pathway database
#'
#' @param database_name Name of the pathway database
#' @return List of pathways with gene symbols
#'
#' @export
load_pathway_database <- function(database_name) {
  switch(database_name,
    "CancerSEA" = load_cancersea_pathways(),
    "Hallmark" = load_hallmark_pathways(),
    "Custom" = load_custom_pathways(),
    stop("Unknown database: ", database_name)
  )
}
```

**2. Add to UI Options**:
```r
# Update pathway selection dropdown
selectInput("pathway_database", "Pathway Database:",
            choices = c("CancerSEA", "Hallmark", "Custom"))
```

---

## 🐛 Bug Reports

### Good Bug Report Template

```markdown
**Bug Description**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected Behavior**
A clear description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Environment Information:**
 - OS: [e.g. Windows 10, macOS 11.0, Ubuntu 20.04]
 - R Version: [e.g. 4.1.0]
 - MASIH Version: [e.g. 0.1.0]
 - Browser: [e.g. Chrome 91, Firefox 89]

**Data Information:**
 - Data type: [e.g. 10X, Seurat object]
 - Number of cells: [approximate]
 - File size: [approximate]

**Additional Context**
Add any other context about the problem here.

**Session Info**
```r
sessionInfo()
```

**Error Messages**
Paste any error messages here.
```

### Bug Triage Labels

We use these labels for bug triage:
- **bug**: Confirmed bug
- **enhancement**: New feature request
- **documentation**: Documentation improvements
- **good first issue**: Good for newcomers
- **help wanted**: Extra attention needed
- **priority-high**: Critical bugs
- **priority-low**: Minor issues

---

## 💡 Feature Requests

### Feature Request Template

```markdown
**Is your feature request related to a problem?**
A clear description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear description of what you want to happen.

**Describe alternatives you've considered**
A clear description of alternative solutions you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.

**Biological Rationale**
Explain why this feature would be valuable for cancer research.

**Implementation Ideas**
If you have ideas about how to implement this, please share.
```

---

## 📚 Documentation Contributions

### Documentation Standards

**Writing Style**:
- **Clear and concise**: Use simple language
- **Step-by-step**: Provide clear instructions
- **Biological context**: Explain why steps matter
- **Examples**: Include practical examples
- **Screenshots**: Add visual aids where helpful

**Documentation Types Needed**:
1. **User guides**: How to use MASIH features
2. **Tutorials**: Complete analysis walkthroughs
3. **API documentation**: Function references
4. **Developer guides**: Technical implementation details
5. **Troubleshooting**: Common problems and solutions

### Adding New Documentation

**1. Create markdown file** in appropriate directory:
```
docs/
├── user-guide/
├── tutorials/
├── developer/
└── api/
```

**2. Follow existing templates** for consistency

**3. Test all code examples** to ensure they work

**4. Add to navigation** in README.md and other relevant files

---

## 🧬 Scientific Contributions

### Validating Methods

**Dataset Requirements**:
- **Public datasets**: Use published, accessible data
- **Documentation**: Clearly document data source and processing
- **Expected results**: Define what results should look like
- **Validation**: Compare with published findings

**Validation Process**:
1. **Select dataset**: Choose appropriate published data
2. **Reproduce analysis**: Follow published methods
3. **Compare results**: MASIH vs published findings
4. **Document findings**: Create validation report
5. **Submit contribution**: Via pull request with documentation

### Adding Example Datasets

**Dataset Criteria**:
- **Size**: <50MB for inclusion in package
- **Quality**: High-quality, well-processed data
- **Documentation**: Clear metadata and expected results
- **Licensing**: Appropriate license for redistribution

**Dataset Format**:
```r
# Save example dataset
example_data <- list(
  seurat_object = processed_seurat_obj,
  metadata = data.frame(...),
  expected_clusters = 8,
  description = "Melanoma progression dataset...",
  source = "Published study XYZ",
  citation = "Author et al. Journal 2023"
)

usethis::use_data(example_data, overwrite = TRUE)
```

---

## 🔄 Pull Request Process

### Before Submitting

**Checklist**:
- [ ] **Code follows style guidelines**
- [ ] **Tests pass**: `devtools::test()`
- [ ] **Documentation updated**: New functions documented
- [ ] **News updated**: Add entry to NEWS.md
- [ ] **Examples work**: All code examples run successfully
- [ ] **Branch is up to date**: Merge latest main branch

### Pull Request Template

```markdown
## Description
Brief description of changes made.

## Type of Change
- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)
- [ ] Documentation update

## Testing
- [ ] Unit tests pass
- [ ] Integration tests pass
- [ ] Manual testing completed
- [ ] Examples in documentation work

## Screenshots (if applicable)
Add screenshots of UI changes.

## Checklist
- [ ] My code follows the style guidelines
- [ ] I have performed a self-review of my code
- [ ] I have commented my code, particularly in hard-to-understand areas
- [ ] I have made corresponding changes to the documentation
- [ ] My changes generate no new warnings
- [ ] I have added tests that prove my fix is effective or that my feature works
- [ ] New and existing unit tests pass locally with my changes
```

### Review Process

**Review Criteria**:
1. **Code quality**: Follows style guidelines, well-documented
2. **Testing**: Adequate test coverage, tests pass
3. **Documentation**: Updated and accurate
4. **Biological validity**: Makes scientific sense
5. **User experience**: Improves or maintains usability

**Review Timeline**:
- **Initial response**: Within 1 week
- **Full review**: Within 2 weeks
- **Revision cycles**: As needed
- **Merge**: After approval from maintainers

---

## 🏗️ Development Environment

### Required Tools

**Essential Software**:
```bash
# Git for version control
git --version

# R development environment
R --version
# Should be ≥ 4.0.0

# RStudio (recommended)
# Download from rstudio.com
```

**R Development Packages**:
```r
# Install development tools
install.packages(c(
  "devtools",      # Development tools
  "roxygen2",      # Documentation generation
  "testthat",      # Unit testing
  "usethis",       # Package development helpers
  "pkgdown",       # Documentation websites
  "lintr",         # Code linting
  "styler"         # Code formatting
))
```

### Useful Development Commands

```r
# Load package for development
devtools::load_all()

# Run tests
devtools::test()

# Check package
devtools::check()

# Generate documentation
roxygen2::roxygenise()

# Build package website
pkgdown::build_site()

# Check code style
lintr::lint_package()
styler::style_pkg()
```

---

## 🎯 Maintainer Guidelines

### Release Process

**Version Numbering**: Follow semantic versioning (MAJOR.MINOR.PATCH)
- **MAJOR**: Breaking changes
- **MINOR**: New features, backwards compatible
- **PATCH**: Bug fixes, backwards compatible

**Release Checklist**:
1. **Update version** in DESCRIPTION
2. **Update NEWS.md** with changes
3. **Run full test suite**
4. **Update documentation**
5. **Create release on GitHub**
6. **Announce on social media/forums**

### Community Management

**Communication Channels**:
- **GitHub Issues**: Bug reports and feature requests
- **GitHub Discussions**: General questions and ideas
- **Email**: Direct support for complex issues
- **Social Media**: Announcements and community building

**Response Guidelines**:
- **Be welcoming**: Encourage all contributors
- **Be patient**: Help newcomers learn
- **Be constructive**: Focus on solutions
- **Be responsive**: Aim for timely responses

---

## 📞 Getting Help

### For Contributors

**Development Questions**:
- **GitHub Discussions**: Ask development questions
- **Email**: your.email@institution.edu for complex issues
- **Documentation**: Check existing developer docs

**Resources**:
- **R Packages book**: [r-pkgs.org](https://r-pkgs.org/)
- **Shiny development**: [mastering-shiny.org](https://mastering-shiny.org/)
- **Golem framework**: [thinkr-open.github.io/golem](https://thinkr-open.github.io/golem/)

### For Maintainers

**Maintainer Resources**:
- **Code review guidelines**: Focus on quality and usability
- **Community management**: Foster inclusive environment
- **Release management**: Coordinate releases and testing

---

## 🙏 Recognition

### Contributors

All contributors will be recognized in:
- **README.md**: Contributor list
- **Package documentation**: Author/contributor credits
- **Release notes**: Acknowledgment of contributions
- **Website**: Contributor gallery

### Types of Recognition

**Code Contributors**:
- Listed as contributors in DESCRIPTION file
- Mentioned in release announcements
- Invited to contribute to papers using MASIH

**Documentation Contributors**:
- Credited in documentation sections
- Listed in acknowledgments
- Invited to co-author tutorials or reviews

**Scientific Contributors**:
- Co-authorship opportunities on methods papers
- Collaboration on validation studies
- Conference presentation opportunities

---

**Thank you for contributing to MASIH!** 🎉 Together, we're making single-cell cancer analysis accessible to researchers worldwide.
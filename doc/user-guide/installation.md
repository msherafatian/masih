# MASIH Installation Guide

**Complete step-by-step instructions for installing MASIH and all dependencies.**

---

## 📋 System Requirements

Before installing MASIH, ensure your system meets these requirements:

### Minimum Requirements
- **R**: Version 4.0.0 or higher ([Download R](https://cran.r-project.org/))
- **RStudio**: Latest version recommended ([Download RStudio](https://rstudio.com/products/rstudio/download/))
- **Operating System**: Windows 10+, macOS 10.14+, or Linux
- **Memory**: 8GB RAM minimum
- **Storage**: 2GB free space for installation and dependencies

### Recommended Requirements
- **Memory**: 16GB RAM (for large datasets >10,000 cells)
- **Storage**: 5GB free space (for example data and exports)
- **Internet**: Stable connection for package downloads

---

## 🚀 Quick Installation (5 minutes)

### Option 1: One-Line Install (Recommended)

```r
# Install all dependencies and MASIH in one go
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

devtools::install_github("yourusername/masih", dependencies = TRUE)
```

### Test Your Installation

```r
library(masih)
run_app()
```

**✅ Success**: If a web browser opens showing the MASIH interface, you're all set!

**❌ Problems**: Continue to the detailed installation below.

---

## 🔧 Detailed Installation (Step-by-Step)

### Step 1: Install Base R Packages

Open RStudio and run each code block:

```r
# Essential CRAN packages
required_packages <- c(
    "shiny", "shinydashboard", "DT", "plotly",
    "dplyr", "ggplot2", "viridis", "RColorBrewer",
    "corrplot", "openxlsx", "devtools", "config",
    "golem", "tidyr", "magrittr"
)

# Install packages that aren't already installed
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
```

### Step 2: Install Bioconductor Packages

```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Set Bioconductor version (use latest)
BiocManager::install(version = "3.18")

# Install required Bioconductor packages
bioc_packages <- c(
    "SingleCellExperiment", 
    "slingshot"
)

BiocManager::install(bioc_packages)
```

### Step 3: Install Seurat

```r
# Install Seurat (large package, may take several minutes)
install.packages("Seurat")

# Verify Seurat installation
library(Seurat)
packageVersion("Seurat")  # Should be >= 4.0.0
```

### Step 4: Install CancerSEA

```r
# Install CancerSEA for functional analysis
if (!requireNamespace("cancersea", quietly = TRUE)) {
    # Try from CRAN first
    install.packages("cancersea")
    
    # If that fails, try from GitHub
    if (!requireNamespace("cancersea", quietly = TRUE)) {
        devtools::install_github("Moonerss/cancersea")
    }
}
```

### Step 5: Install MASIH

```r
# Install MASIH from GitHub
devtools::install_github("yourusername/masih")

# Load and test
library(masih)
```

### Step 6: Final Verification

```r
# Test all major components
library(masih)
library(Seurat)
library(SingleCellExperiment)
library(slingshot)

# Launch MASIH
run_app()
```

---

## 🛠️ Troubleshooting Common Issues

### Issue 1: R Version Too Old

**Error**: `R version 4.0.0 or higher is required`

**Solution**:
1. Download latest R from [CRAN](https://cran.r-project.org/)
2. Install new R version
3. Restart RStudio
4. Reinstall packages

### Issue 2: Package Installation Failures

**Error**: `installation of package 'XXX' had non-zero exit status`

**Solutions**:

#### For Windows Users:
```r
# Install Rtools if missing
install.packages("installr")
installr::install.Rtools()
```

#### For macOS Users:
```bash
# Install Xcode command line tools in Terminal
xcode-select --install
```

#### For Linux Users:
```bash
# Ubuntu/Debian - install system dependencies
sudo apt-get update
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev

# CentOS/RedHat
sudo yum install openssl-devel libcurl-devel libxml2-devel fontconfig-devel
```

#### Alternative Package Installation:
```r
# Try different repository
options(repos = c(CRAN = "https://cloud.r-project.org/"))
install.packages("package_name")

# Or try binary packages (faster)
install.packages("package_name", type = "binary")
```

### Issue 3: Bioconductor Installation Problems

**Error**: `Bioconductor version X.X is not supported`

**Solution**:
```r
# Update BiocManager to latest version
remove.packages("BiocManager")
install.packages("BiocManager")

# Install latest Bioconductor
BiocManager::install(version = "devel")  # For latest features
# OR
BiocManager::install()  # For stable release
```

### Issue 4: Memory Issues During Installation

**Error**: `cannot allocate vector of size X`

**Solutions**:
1. **Close other applications** to free memory
2. **Restart R session**: Session → Restart R
3. **Install packages one by one**:
   ```r
   install.packages("Seurat")  # Wait for completion
   install.packages("SingleCellExperiment")  # Then this one
   ```
4. **Increase memory limit** (Windows only):
   ```r
   memory.limit(size = 8000)  # 8GB limit
   ```

### Issue 5: Seurat Installation Problems

**Error**: `Seurat installation failed`

**Solutions**:
```r
# Try installing dependencies first
install.packages(c("Matrix", "Rcpp", "RcppAnnoy", "RcppHNSW"))

# Then install Seurat
install.packages("Seurat")

# Alternative: Install from GitHub (development version)
devtools::install_github("satijalab/seurat", ref = "release/4.3.0")
```

### Issue 6: CancerSEA Not Found

**Error**: `package 'cancersea' is not available`

**Solutions**:
```r
# Method 1: Install from Bioconductor
BiocManager::install("cancersea")

# Method 2: Install from GitHub
devtools::install_github("Moonerss/cancersea")

# Method 3: Skip for now (MASIH will work without it)
# CancerSEA features will be disabled
```

### Issue 7: MASIH Won't Launch

**Error**: `Shiny application failed to start`

**Diagnostic Steps**:
```r
# Check package loading
library(masih)  # Look for error messages

# Check dependencies
library(shiny)
library(shinydashboard)
library(Seurat)

# Try launching with error details
options(shiny.error = browser)
run_app()
```

**Common Solutions**:
```r
# Try different port
run_app(options = list(port = 3839))

# Disable browser auto-launch
run_app(options = list(launch.browser = FALSE))
# Then manually open: http://127.0.0.1:3838

# Clear package cache
remove.packages("masih")
devtools::install_github("yourusername/masih", force = TRUE)
```

---

## 🔍 Installation Verification Checklist

Run this checklist to ensure everything is working:

```r
# ✅ Check R version
R.version.string

# ✅ Check key packages
required_libs <- c("masih", "Seurat", "shiny", "SingleCellExperiment", "slingshot")
missing_libs <- required_libs[!sapply(required_libs, requireNamespace, quietly = TRUE)]

if(length(missing_libs) == 0) {
    cat("✅ All packages installed successfully!\n")
} else {
    cat("❌ Missing packages:", paste(missing_libs, collapse = ", "), "\n")
}

# ✅ Test MASIH launch
cat("Testing MASIH launch...\n")
tryCatch({
    library(masih)
    cat("✅ MASIH loaded successfully!\n")
    cat("Ready to run: run_app()\n")
}, error = function(e) {
    cat("❌ MASIH loading failed:", e$message, "\n")
})
```

---

## 📱 Platform-Specific Notes

### Windows
- **Rtools**: Required for compiling packages from source
- **Antivirus**: May block package downloads (temporarily disable)
- **Firewall**: Allow R/RStudio through firewall for Shiny apps

### macOS
- **Xcode**: Command line tools required for compilation
- **Gatekeeper**: May prevent running downloaded software
- **Rosetta**: M1 Macs may need Rosetta for some packages

### Linux
- **System libraries**: Many R packages need system dependencies
- **Permissions**: Ensure write permissions in R library directory
- **Package managers**: Use system package manager for dependencies

---

## 🚀 Next Steps

Once installation is complete:

1. **📖 Read the Getting Started Guide**: [getting-started.md](getting-started.md)
2. **🧪 Try a Tutorial**: [Basic Analysis Tutorial](../tutorials/basic-analysis.md)
3. **📊 Prepare Your Data**: [Data Preparation Guide](data-preparation.md)
4. **❓ Get Help**: [Troubleshooting Guide](troubleshooting.md)

---

## 💬 Getting Help

**Still having trouble?**

1. **Check our FAQ**: [troubleshooting.md](troubleshooting.md)
2. **Search existing issues**: [GitHub Issues](https://github.com/yourusername/masih/issues)
3. **Create new issue**: Include your sessionInfo() output
4. **Email support**: your.email@institution.edu

### When Reporting Issues:

```r
# Include this information in bug reports
sessionInfo()
```

---

**Installation successful? Great!** 🎉 You're ready to explore single-cell data with MASIH!
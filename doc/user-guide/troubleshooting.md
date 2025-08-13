# MASIH Troubleshooting Guide

**Solutions to common problems and how to get help when you're stuck.**

---

## 🚨 Quick Fixes (Try These First)

### App Not Responding
1. **Refresh browser tab** (F5 or Ctrl+R)
2. **Restart MASIH**: Close browser, run `run_app()` again
3. **Restart R session**: In RStudio → Session → Restart R
4. **Check memory usage**: Close other applications

### Error Messages
1. **Read the error carefully** - often tells you exactly what's wrong
2. **Check file formats** - ensure data is properly formatted
3. **Verify file sizes** - keep under 2GB limit
4. **Clear browser cache** - Ctrl+Shift+Delete

### Slow Performance
1. **Use smaller datasets** for initial exploration
2. **Close unused browser tabs**
3. **Restart the app** periodically
4. **Check available RAM** - need 8GB minimum

---

## 📂 Data Upload Issues

### Problem: "File upload failed"

**Possible Causes**:
- File too large (>2GB)
- Network connection interrupted
- Unsupported file format
- Corrupted file

**Solutions**:
```r
# Check file size in R
file.info("your_file.rds")$size / 1024^3  # Size in GB

# Test file integrity
test_obj <- readRDS("your_file.rds")
class(test_obj)  # Should be "Seurat"

# Reduce file size if needed
seurat_subset <- seurat_obj[, sample(colnames(seurat_obj), 5000)]
saveRDS(seurat_subset, "smaller_file.rds")
```

### Problem: "Data format not recognized"

**10X Data Issues**:
```r
# Verify 10X file structure
list.files("path/to/10x/data")
# Should show: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz

# Test loading in R first
library(Seurat)
test_data <- Read10X("path/to/10x/data")
dim(test_data)  # Should show genes × cells
```

**Seurat Object Issues**:
```r
# Check Seurat object validity
seurat_obj <- readRDS("your_file.rds")
class(seurat_obj)                    # Should be "Seurat"
length(seurat_obj@assays)           # Should have at least "RNA"
head(rownames(seurat_obj))          # Should be gene symbols
```

### Problem: "Gene names not recognized"

**Ensembl ID Issue**:
```r
# Convert Ensembl IDs to gene symbols
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Example conversion (simplified)
head(rownames(seurat_obj))  # Check current format
# If they look like "ENSG00000000003", convert them
```

**Species Mismatch**:
```r
# Check for mixed human/mouse genes
human_pattern <- "^[A-Z][A-Z0-9-]+$"     # GAPDH, TP53
mouse_pattern <- "^[A-Z][a-z0-9-]+$"     # Gapdh, Tp53

human_genes <- grep(human_pattern, rownames(seurat_obj), value = TRUE)
mouse_genes <- grep(mouse_pattern, rownames(seurat_obj), value = TRUE)

length(human_genes)  # Should be majority for human data
length(mouse_genes)  # Should be majority for mouse data
```

---

## 🔍 Quality Control Problems

### Problem: "Too many/few cells after filtering"

**Too Few Cells Remaining**:
```r
# Check your QC thresholds - they might be too strict
summary(seurat_obj$nFeature_RNA)  # Gene count distribution
summary(seurat_obj$percent.mt)    # Mitochondrial percentage

# Relax filters
# Instead of: min_features = 500, try min_features = 200
# Instead of: max_mt = 15, try max_mt = 25
```

**Too Many Cells (Suspected Doublets)**:
```r
# Check for high gene count cells (possible doublets)
hist(seurat_obj$nFeature_RNA, breaks = 50)
# Look for cells with >5000 genes - likely doublets

# Consider using DoubletFinder before MASIH
library(DoubletFinder)
# (Run DoubletFinder workflow)
```

### Problem: "Strange QC plots"

**Bimodal Distributions**:
- **Cause**: Likely mixing cell types or qualities
- **Solution**: Check experimental setup, consider batch effects

**No Clear Cutoffs**:
- **Cause**: Dataset may need platform-specific thresholds
- **Solution**: Compare to published data from same platform

---

## 🎯 Clustering Issues

### Problem: "Clustering results don't make sense"

**Too Many Clusters**:
```r
# Try lower resolution
# In MASIH: Change resolution from 0.8 to 0.3 or 0.5
```

**Too Few Clusters**:
```r
# Try higher resolution  
# In MASIH: Change resolution from 0.5 to 0.8 or 1.0
```

**Artificial Clusters**:
- **Check QC**: Poor quality cells create artificial clusters
- **Check batch effects**: Different samples clustering separately
- **Check cell cycle**: G2/M cells might cluster together

### Problem: "UMAP plot looks strange"

**All cells in one blob**:
- **Increase PCA dimensions**: Try 20-30 instead of 10
- **Check feature selection**: Ensure highly variable genes calculated
- **Verify normalization**: Data should be log-normalized

**Cells scattered everywhere**:
- **Check for doublets**: High gene count cells
- **Reduce PCA dimensions**: Try 10-15 instead of 30
- **Check data quality**: Might have technical artifacts

---

## 🧬 Marker Gene Problems

### Problem: "No significant marker genes found"

**Possible Causes**:
- Clusters too similar (over-clustering)
- Statistical thresholds too strict
- Insufficient cells per cluster

**Solutions**:
```r
# In MASIH marker gene module:
# - Reduce logfc.threshold from 0.25 to 0.1
# - Increase max p-value from 0.05 to 0.1
# - Reduce min.pct from 0.1 to 0.05
```

### Problem: "Marker genes don't make biological sense"

**Check for Technical Artifacts**:
- **Mitochondrial genes**: High in dying cells
- **Ribosomal genes**: Technical variation
- **Cell cycle genes**: Proliferating cells

**Validation Steps**:
```r
# Look up unknown genes
# Use databases like GeneCards, UniProt
# Check if genes are known cell type markers
```

---

## 🧠 CancerSEA Issues

### Problem: "CancerSEA pathway calculation failed"

**Error Messages**:
- `"Insufficient data for binning"`: Dataset too small
- `"No genes found"`: Gene name mismatch
- `"Score calculation failed"`: Technical error

**Solutions**:
```r
# Check gene overlap with pathway
data('available_pathways', package = 'cancersea')
pathway_genes <- get("Stemness")$symbol
overlap <- intersect(pathway_genes, rownames(seurat_obj))
length(overlap)  # Should be >10 for reliable scoring

# If low overlap, check gene naming
head(rownames(seurat_obj))  # Should be like "GAPDH", "TP53"
```

### Problem: "CancerSEA scores seem wrong"

**All Scores Similar**:
- **Normal**: Scores are relative within dataset
- **Check**: Compare scores between clusters, not absolute values

**Unexpected High/Low Scores**:
- **Verify cell types**: Make sure you have cancer cells
- **Check normalization**: Use log-normalized data
- **Compare literature**: Expected scores for your cancer type

---

## 📈 Trajectory Analysis Problems

### Problem: "Trajectory analysis failed to run"

**Common Issues**:
- **Too few cells selected**: Need >100 cells per cell type
- **No clear trajectory**: Cells might not be on developmental path
- **Parameter issues**: PCA dimensions or cell type selection

**Solutions**:
```r
# Check cell type selection
table(seurat_obj$cell_type)  # Ensure enough cells per type

# Verify biological relationship
# Are selected cell types expected to be on same trajectory?

# Try different parameters
# - Increase PCA dimensions: 10 → 20
# - Include more cell types
# - Check start/end cluster specification
```

### Problem: "Trajectory doesn't match biology"

**Unexpected Paths**:
- **Check cell type annotations**: Might be mislabeled
- **Consider batch effects**: Different samples on different paths
- **Review literature**: Expected developmental relationships

**Multiple Trajectories**:
- **Normal**: Complex systems often have multiple paths
- **Focus**: Analyze one trajectory at a time

---

## 💾 Export Problems

### Problem: "Plot export failed"

**Large File Issues**:
- **Reduce plot size**: 8×6 inches instead of 12×10
- **Lower DPI**: 150 instead of 300 for testing
- **Try different format**: PNG instead of PDF

**Memory Issues**:
```r
# Restart R session before export
# In RStudio: Session → Restart R
library(masih)
run_app()
# Then try export again
```

### Problem: "Data export incomplete"

**Excel File Issues**:
- **Check file size**: Large datasets might timeout
- **Export sections separately**: One module at a time
- **Use CSV format**: For very large data tables

---

## 🖥️ Technical Issues

### Problem: "App crashes or freezes"

**Browser Issues**:
1. **Try different browser**: Chrome, Firefox, Safari
2. **Clear browser cache**: Ctrl+Shift+Delete
3. **Disable browser extensions**: AdBlockers might interfere
4. **Check JavaScript**: Ensure JavaScript is enabled

**Memory Issues**:
```r
# Check memory usage
gc()  # Garbage collection
memory.limit()  # Windows only

# Increase memory (Windows)
memory.limit(size = 16000)  # 16GB

# Monitor memory usage
pryr::mem_used()
```

**R Session Issues**:
```r
# Check R session info
sessionInfo()

# Update packages if needed
update.packages()

# Reinstall MASIH if corrupted
remove.packages("masih")
devtools::install_github("yourusername/masih")
```

### Problem: "Port already in use"

**Error**: `listen tcp 127.0.0.1:3838: bind: address already in use`

**Solutions**:
```r
# Try different port
run_app(options = list(port = 3839))

# Or kill existing process (advanced)
# Windows: netstat -ano | findstr :3838
# Mac/Linux: lsof -ti:3838 | xargs kill
```

---

## 🔬 Platform-Specific Issues

### Windows Specific

**Path Issues**:
```r
# Use forward slashes or double backslashes
"C:/Users/name/data.rds"          # ✅ Good
"C:\\Users\\name\\data.rds"       # ✅ Good  
"C:\Users\name\data.rds"          # ❌ Bad
```

**Encoding Issues**:
```r
# If you see strange characters
Sys.setlocale("LC_ALL", "English")
```

### macOS Specific

**Gatekeeper Issues**:
- **Problem**: "App can't be opened because it's from unidentified developer"
- **Solution**: System Preferences → Security → Allow apps downloaded from anywhere

**M1 Mac Issues**:
```r
# Some packages need Rosetta
# Install Rosetta 2 if prompted

# Check architecture
Sys.info()["machine"]  # Should show architecture
```

### Linux Specific

**Missing System Libraries**:
```bash
# Ubuntu/Debian
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev

# CentOS/RedHat  
sudo yum install openssl-devel libcurl-devel libxml2-devel
```

---

## 📊 Performance Optimization

### For Large Datasets (>20,000 cells)

```r
# Subsample for initial exploration
set.seed(123)
subset_cells <- sample(colnames(seurat_obj), 10000)
seurat_subset <- seurat_obj[, subset_cells]

# Use fewer PCA dimensions initially
# Start with 10 PCs instead of 30

# Calculate fewer pathways initially
# Start with 2-3 CancerSEA pathways instead of all
```

### For Slow Analysis

```r
# Reduce clustering resolution for speed
# Use 0.3-0.5 instead of 0.8-1.0

# Use fewer highly variable genes
# 2000 instead of 3000

# Restart R session periodically
# Prevents memory accumulation
```

---

## 🆘 Getting Help

### Before Asking for Help

**Gather Information**:
```r
# System information
sessionInfo()

# MASIH version
packageVersion("masih")

# Data information
ncol(seurat_obj)  # Number of cells
nrow(seurat_obj)  # Number of genes
```

**Document the Problem**:
1. **What were you trying to do?**
2. **What happened instead?**
3. **Any error messages?** (copy exact text)
4. **What data are you using?** (size, type, source)

### Where to Get Help

**1. Check Documentation First**:
- [Installation Guide](installation.md)
- [Getting Started](getting-started.md)
- [Data Preparation](data-preparation.md)
- [Tutorials](../tutorials/)

**2. Search Existing Issues**:
- [GitHub Issues](https://github.com/yourusername/masih/issues)
- Use search to find similar problems

**3. Community Help**:
- [GitHub Discussions](https://github.com/yourusername/masih/discussions)
- Biostars forum
- Stack Overflow (use `masih` and `single-cell` tags)

**4. Direct Support**:
- **Email**: your.email@institution.edu
- **Include**: sessionInfo(), error messages, data description

### Creating Good Bug Reports

**Template**:
```
**Problem Description**
Brief description of what went wrong

**Steps to Reproduce**
1. Load data: [describe data]
2. Go to [module name]
3. Click [button name]
4. Error occurs

**Expected Behavior**
What you expected to happen

**Actual Behavior**  
What actually happened (include error messages)

**System Information**
- Operating System: [Windows 10/macOS/Linux]
- R Version: [from R.version.string]
- MASIH Version: [from packageVersion("masih")]
- Browser: [Chrome/Firefox/Safari + version]

**Data Information**
- Data type: [10X/Seurat/Matrix]
- Number of cells: [approximate]
- File size: [approximate]

**Additional Context**
Any other information that might help
```

---

## ✅ Diagnostic Checklist

**When something goes wrong, check:**

### Basic System
- [ ] **R version** ≥ 4.0.0
- [ ] **MASIH installed** and up to date
- [ ] **All dependencies** installed
- [ ] **Sufficient memory** available (8GB+)

### Data Issues
- [ ] **File format** supported
- [ ] **File size** under limits
- [ ] **Gene names** are symbols
- [ ] **Data quality** looks reasonable

### App Function
- [ ] **Browser compatibility** (Chrome/Firefox recommended)
- [ ] **JavaScript enabled**
- [ ] **No browser extensions** interfering
- [ ] **Port not blocked** by firewall

### Analysis Parameters
- [ ] **Reasonable QC thresholds**
- [ ] **Appropriate clustering resolution**
- [ ] **Sufficient cells** for analysis
- [ ] **Expected biological relationships**

---

## 🚀 Prevention Tips

**Avoid Problems**:

1. **Start with example data** to learn the interface
2. **Use recommended file formats** (10X or Seurat)
3. **Check data quality** before upload
4. **Save work frequently** by downloading results
5. **Document parameters** that work well
6. **Keep backups** of important data
7. **Update software** regularly
8. **Monitor system resources** during analysis

**Best Practices**:
- **Read error messages carefully** - they often contain the solution
- **Try simple solutions first** - restart, refresh, clear cache
- **Test with smaller data** if having problems
- **Save intermediate results** to avoid losing work
- **Ask specific questions** when seeking help

---

**Remember**: Most problems have simple solutions! 🎯 Don't hesitate to ask for help if you're stuck.
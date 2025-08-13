# Data Preparation Guide

**Optimize your single-cell data for the best MASIH analysis results.**

---

## 📋 Overview

Proper data preparation is crucial for meaningful single-cell analysis. This guide covers:

- **Supported data formats** and requirements
- **Pre-processing recommendations** before MASIH
- **Quality considerations** for optimal results
- **Common data issues** and solutions
- **Best practices** for different data types

**Time to read**: ~10 minutes  
**Prerequisites**: Basic understanding of single-cell data

---

## 📊 Supported Data Formats

### 1. 10X Genomics Data (Recommended)

**Cell Ranger Output Directory**:
```
your_sample/
├── filtered_feature_bc_matrix/
│   ├── barcodes.tsv.gz          # Cell barcodes
│   ├── features.tsv.gz          # Gene information
│   └── matrix.mtx.gz            # Expression matrix
└── raw_feature_bc_matrix/       # Optional: unfiltered data
```

**Requirements**:
- ✅ **Cell Ranger v3+** output format
- ✅ **Gzipped files** (.gz extension)
- ✅ **Standard file names** (barcodes, features, matrix)
- ✅ **Gene symbols** in features file

**File Size Limits**:
- **Maximum**: 2GB total upload size
- **Recommended**: <500MB for smooth performance
- **Large datasets**: Consider subsampling initially

### 2. Seurat Objects (.rds)

**Pre-processed Seurat Objects**:
```r
# Example of properly formatted Seurat object
seurat_obj <- CreateSeuratObject(
    counts = expression_matrix,
    min.cells = 3,
    min.features = 200,
    project = "MyProject"
)

# Save for MASIH
saveRDS(seurat_obj, "my_data.rds")
```

**Requirements**:
- ✅ **Seurat v4+** compatible
- ✅ **Gene symbols** as row names
- ✅ **Metadata included** (if available)
- ✅ **Raw counts** preserved

### 3. Gene Expression Matrices

**CSV/TSV Format**:
```
Gene        Cell1    Cell2    Cell3    ...
GAPDH       150      89       234      ...
ACTB        89       156      67       ...
TP53        12       45       23       ...
```

**Requirements**:
- ✅ **Genes as rows**, cells as columns
- ✅ **Gene symbols** in first column
- ✅ **Raw or normalized counts** (specify which)
- ✅ **No special characters** in gene names

### 4. SingleCellExperiment Objects

**Bioconductor SCE Objects**:
```r
# Convert to Seurat for MASIH compatibility
library(Seurat)
seurat_obj <- as.Seurat(sce_object)
saveRDS(seurat_obj, "converted_data.rds")
```

---

## 🔬 Pre-Processing Recommendations

### For Raw 10X Data (Recommended Workflow)

If you have raw Cell Ranger output, minimal pre-processing is best:

```r
library(Seurat)

# Load data
data <- Read10X("path/to/cellranger/output/")
seurat_obj <- CreateSeuratObject(
    counts = data,
    min.cells = 3,      # Filter genes in <3 cells
    min.features = 200, # Filter cells with <200 genes
    project = "MyProject"
)

# Add mitochondrial gene percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
    seurat_obj, 
    pattern = "^MT-"    # Human: ^MT-, Mouse: ^mt-
)

# Save for MASIH - let MASIH do the rest!
saveRDS(seurat_obj, "for_masih.rds")
```

**✅ What MASIH will handle**:
- Quality control filtering
- Normalization and scaling
- Feature selection
- Dimensionality reduction
- Clustering

### For Pre-Processed Data

If your data is already processed, ensure these steps were done properly:

```r
# Check your pre-processed object
seurat_obj

# Verify data integrity
head(rownames(seurat_obj))           # Should be gene symbols
head(colnames(seurat_obj))           # Should be cell barcodes
max(seurat_obj[["RNA"]]@counts)      # Should be reasonable counts

# Check metadata
head(seurat_obj@meta.data)

# Save if everything looks good
saveRDS(seurat_obj, "preprocessed_for_masih.rds")
```

---

## 🎯 Data Quality Guidelines

### Minimum Dataset Requirements

**Cell Numbers**:
- **Minimum**: 100 cells (for testing)
- **Recommended**: 1,000-10,000 cells
- **Large studies**: 10,000+ cells (may need subsampling)

**Gene Coverage**:
- **Minimum**: 500 genes per cell on average
- **Typical**: 1,000-3,000 genes per cell
- **High quality**: 2,000+ genes per cell

**Sequencing Depth**:
- **Minimum**: 1,000 UMIs per cell
- **Recommended**: 2,000-10,000 UMIs per cell
- **Deep sequencing**: 10,000+ UMIs per cell

### Quality Metrics to Check

Before uploading to MASIH, verify:

```r
# Basic statistics
ncol(seurat_obj)  # Number of cells
nrow(seurat_obj)  # Number of genes

# Quality metrics per cell
summary(seurat_obj$nFeature_RNA)  # Genes per cell
summary(seurat_obj$nCount_RNA)    # UMIs per cell
summary(seurat_obj$percent.mt)    # Mitochondrial %

# Look for outliers
hist(seurat_obj$nFeature_RNA, breaks = 50)
hist(seurat_obj$percent.mt, breaks = 50)
```

**Warning Signs**:
- **Very few genes per cell** (<200): Poor cell capture
- **Too many genes per cell** (>6000): Possible doublets
- **High mitochondrial %** (>25%): Dying cells
- **Low total UMIs** (<500): Poor sequencing

---

## 🐛 Common Data Issues and Solutions

### Issue 1: Gene Name Problems

**Problem**: Ensembl IDs instead of gene symbols
```
ENSG00000000003    # ❌ Not ideal
ENSG00000000005    # ❌ Not ideal
```

**Solution**: Convert to gene symbols
```r
# Using biomaRt to convert Ensembl to symbols
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get gene symbols
gene_map <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = rownames(seurat_obj),
    mart = mart
)

# Update gene names (simplified example)
# Note: Handle duplicates and missing symbols appropriately
```

### Issue 2: Mixed Species Data

**Problem**: Human and mouse genes mixed
```
GAPDH      # Human
Gapdh      # Mouse
```

**Solution**: Filter to single species
```r
# For human data, keep only human genes
human_genes <- rownames(seurat_obj)[grepl("^[A-Z]", rownames(seurat_obj))]
seurat_obj <- seurat_obj[human_genes, ]

# For mouse data, keep only mouse genes  
mouse_genes <- rownames(seurat_obj)[grepl("^[A-Z][a-z]", rownames(seurat_obj))]
seurat_obj <- seurat_obj[mouse_genes, ]
```

### Issue 3: Batch Effects

**Problem**: Multiple samples with technical differences

**Solution**: Add batch information
```r
# Add batch/sample information to metadata
seurat_obj$batch <- c(rep("Batch1", 1000), rep("Batch2", 1000))
seurat_obj$sample <- c(rep("Sample_A", 500), rep("Sample_B", 500), 
                       rep("Sample_C", 500), rep("Sample_D", 500))

# MASIH can use this information for visualization
```

### Issue 4: Extremely Large Datasets

**Problem**: >50,000 cells causing memory issues

**Solution**: Intelligent subsampling
```r
# Option 1: Random subsampling
set.seed(123)
subset_cells <- sample(colnames(seurat_obj), 10000)
seurat_subset <- seurat_obj[, subset_cells]

# Option 2: Stratified sampling (if you have cell type annotations)
library(dplyr)
metadata <- seurat_obj@meta.data
metadata$cell_id <- rownames(metadata)

# Sample 1000 cells per cell type
balanced_sample <- metadata %>%
    group_by(cell_type) %>%
    sample_n(min(1000, n())) %>%
    pull(cell_id)

seurat_subset <- seurat_obj[, balanced_sample]
```

### Issue 5: Missing Metadata

**Problem**: No cell type annotations or other metadata

**Solution**: Add available information
```r
# Add basic metadata
seurat_obj$tissue <- "tumor"
seurat_obj$condition <- "untreated"
seurat_obj$patient_id <- "Patient_001"

# If you have experimental conditions
seurat_obj$treatment <- ifelse(
    grepl("treated", colnames(seurat_obj)), 
    "treated", 
    "control"
)
```

---

## 📁 File Organization Best Practices

### Directory Structure

```
my_scrnaseq_project/
├── raw_data/
│   ├── sample1_cellranger_output/
│   ├── sample2_cellranger_output/
│   └── README.txt                    # Document data source
├── processed_data/
│   ├── sample1_for_masih.rds
│   ├── sample2_for_masih.rds
│   └── combined_for_masih.rds
├── masih_results/
│   ├── plots/
│   ├── data_exports/
│   └── analysis_notes.txt
└── documentation/
    ├── sample_metadata.csv
    ├── experimental_design.txt
    └── analysis_log.txt
```

### File Naming Conventions

**Good naming examples**:
- `patient001_tumor_treated_masih.rds`
- `batch1_control_preprocessed.rds`
- `melanoma_sample_A_cellranger.zip`

**Avoid**:
- `data.rds` (not descriptive)
- `final_final_v2.rds` (version confusion)
- `my data with spaces.rds` (spaces cause issues)

---

## 🔧 Platform-Specific Tips

### 10X Genomics Cell Ranger

**Optimal Cell Ranger Settings**:
```bash
# Recommended Cell Ranger command
cellranger count \
    --id=sample_name \
    --transcriptome=/path/to/refdata \
    --fastqs=/path/to/fastqs \
    --sample=sample_name \
    --localcores=8 \
    --localmem=64
```

**Important notes**:
- Use **latest reference genome** (GRCh38 for human)
- **Force cells parameter** only if you know expected cell number
- **Check Cell Ranger web summary** before analysis

### Drop-seq/inDrop Data

**Convert to 10X format**:
```r
# If you have a digital expression matrix
library(DropletUtils)
write10xCounts("output_directory", expression_matrix)

# Then load as 10X data
data <- Read10X("output_directory")
seurat_obj <- CreateSeuratObject(data)
```

### Smart-seq2/Smart-seq3 Data

**Full-length transcript data**:
```r
# Usually comes as FPKM/TPM - convert if possible to counts
# Or proceed with normalized data but specify in metadata
seurat_obj$data_type <- "TPM"
seurat_obj$protocol <- "Smart-seq2"
```

---

## ✅ Pre-Upload Checklist

Before uploading to MASIH, verify:

### Data Format
- [ ] **File format** is supported (10X, Seurat .rds, or matrix)
- [ ] **File size** is under 2GB
- [ ] **Gene names** are symbols (not Ensembl IDs)
- [ ] **No special characters** in gene names

### Data Quality  
- [ ] **Cell count** >100 (preferably >1000)
- [ ] **Gene count** per cell >200 on average
- [ ] **Mitochondrial genes** calculated (if applicable)
- [ ] **Obvious outliers** identified

### Metadata
- [ ] **Sample information** included
- [ ] **Experimental conditions** documented
- [ ] **Batch effects** considered
- [ ] **Cell types** annotated (if known)

### Documentation
- [ ] **Data source** documented
- [ ] **Processing steps** recorded
- [ ] **Expected results** considered
- [ ] **Backup files** saved

---

## 🎯 Optimization Tips

### For Better Clustering
- **Remove ambient RNA** using tools like SoupX
- **Filter doublets** using DoubletFinder or Scrublet
- **Include cell cycle genes** for downstream analysis

### For Trajectory Analysis
- **Keep developmental intermediates** (don't over-filter)
- **Include known markers** for your system
- **Consider time-series data** if available

### For CancerSEA Analysis
- **Verify gene symbols** match CancerSEA database
- **Include control samples** if available
- **Document tumor type** and stage

---

## 📞 Getting Help with Data Prep

**Data preparation issues?**

1. **Check our troubleshooting**: [troubleshooting.md](troubleshooting.md)
2. **Review tutorials**: Detailed examples in [tutorials](../tutorials/)
3. **GitHub discussions**: Ask the community
4. **Email support**: your.email@institution.edu with data details

**When asking for help, include**:
- Data type and source (10X, Smart-seq, etc.)
- File sizes and cell/gene counts
- Any error messages
- Preprocessing steps already performed

---

## 🚀 Next Steps

**Data ready for MASIH?**

1. **Upload your data**: Follow the [Getting Started Guide](getting-started.md)
2. **Try basic analysis**: [Basic Analysis Tutorial](../tutorials/basic-analysis.md)
3. **Explore advanced features**: [Advanced tutorials](../tutorials/)

**Need more data prep help?**

1. **Platform-specific guides**: Check 10X, Smart-seq documentation
2. **Seurat tutorials**: [satijalab.org](https://satijalab.org/seurat/)
3. **Community resources**: Biostars, Stack Overflow

---

**Proper preparation leads to better insights!** 🎯 Take time to prepare your data well for the best MASIH experience.
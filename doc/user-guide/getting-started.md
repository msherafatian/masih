# Getting Started with MASIH

**Your first single-cell analysis in 15 minutes - from data upload to biological insights.**

---

## 🎯 What You'll Learn

By the end of this guide, you'll know how to:
- Launch MASIH and navigate the interface
- Upload your single-cell data
- Perform quality control and clustering
- Calculate CancerSEA pathway scores
- Export your results

**Time required**: ~15 minutes  
**Prerequisites**: MASIH installed ([Installation Guide](installation.md))

---

## 🚀 Step 1: Launch MASIH

### Start the Application

```r
# Load MASIH
library(masih)

# Launch the application
run_app()
```

**What happens**: A web browser window opens showing the MASIH interface.

### Navigate the Interface

MASIH has **8 main tabs** in the left sidebar:

- **📁 Data Upload** - Load your data
- **🔍 Quality Control** - Filter cells and genes  
- **🎯 Cluster Analysis** - Identify cell populations
- **🧬 Marker Genes** - Find cluster-specific genes
- **🔄 Cell Cycle** - Analyze cell cycle phases
- **🧠 CancerSEA** - Functional state analysis
- **📊 Pathway Comparison** - Compare multiple pathways
- **📈 Trajectory Analysis** - Pseudotime inference
- **📤 Export** - Download results

---

## 📂 Step 2: Upload Your Data

### Option A: Use Example Data (Quickest)

1. **Click "Data Upload" tab**
2. **Select "Load Example Data"** button
3. **Choose** one of the provided cancer datasets
4. **Click "Load Data"**

### Option B: Upload Your Own Data

#### For 10X Genomics Data:
1. **Click "Data Upload" tab**
2. **Select "10X Genomics Data"**
3. **Upload files**:
   - `matrix.mtx.gz` (gene expression matrix)
   - `barcodes.tsv.gz` (cell identifiers)  
   - `features.tsv.gz` (gene information)
4. **Click "Process Data"**

#### For Seurat Objects:
1. **Select "Seurat Object (.rds)"**
2. **Browse and upload** your `.rds` file
3. **Click "Load Object"**

### Data Upload Tips:
- **File size limit**: 2GB maximum
- **Supported formats**: 10X, Seurat objects, CSV matrices
- **Example datasets**: Start with these to learn the interface

---

## 🔍 Step 3: Quality Control

### Automatic QC Metrics

Once data is loaded, MASIH automatically calculates:
- **nFeature_RNA**: Number of genes per cell
- **nCount_RNA**: Total UMI counts per cell  
- **percent.mt**: Mitochondrial gene percentage

### Set QC Filters

1. **Review QC plots** showing cell distributions
2. **Adjust filter parameters**:
   - **Min genes per cell**: 200 (default)
   - **Max genes per cell**: 5000 (default)
   - **Max mitochondrial %**: 20% (default)
3. **Click "Apply Filters"**

### What to Look For:
- **Low gene counts** (< 200): Likely empty droplets or dead cells
- **High gene counts** (> 5000): Possible doublets
- **High mitochondrial %** (> 20%): Stressed or dying cells

### QC Results:
- See **before/after cell counts**
- Review **filtered data statistics**
- **Proceed to clustering** when satisfied

---

## 🎯 Step 4: Perform Clustering

### Run Standard Analysis

1. **Go to "Cluster Analysis" tab**
2. **Click "Run Standard Analysis"** (uses default parameters)
3. **Wait for processing** (1-3 minutes for most datasets)

### What Happens:
- **Normalization**: Log-normalization of gene expression
- **Feature selection**: Identify highly variable genes
- **Scaling**: Z-score scaling of expression values
- **PCA**: Principal component analysis
- **UMAP**: Non-linear dimensionality reduction
- **Clustering**: Graph-based clustering

### View Results:
- **UMAP plot**: Cells colored by cluster
- **Cluster tree**: Hierarchical relationships
- **Cluster statistics**: Cell counts per cluster

### Adjust Parameters (Optional):
- **Resolution**: Higher = more clusters (0.5 default)
- **Dimensions**: Number of PCs to use (10 default)
- **k-parameter**: Nearest neighbors (20 default)

---

## 🧠 Step 5: CancerSEA Functional Analysis

### Calculate Pathway Scores

1. **Go to "CancerSEA" tab**
2. **Select a pathway** from dropdown:
   - Start with **"Stemness"** or **"Proliferation"**
3. **Click "Calculate/Update"**
4. **Wait for calculation** (30 seconds - 2 minutes)

### Available Pathways:
- **Angiogenesis** - Blood vessel formation
- **Apoptosis** - Programmed cell death
- **Cell cycle** - Cell division processes
- **Differentiation** - Cell maturation
- **EMT** - Epithelial-mesenchymal transition
- **Hypoxia** - Low oxygen response
- **Inflammation** - Immune response
- **Invasion** - Tissue infiltration
- **Metastasis** - Distant spread
- **Proliferation** - Cell growth
- **Quiescence** - Cell dormancy
- **Stemness** - Stem cell properties

### Interpret Results:
- **Feature plot**: Pathway scores on UMAP
- **Violin plot**: Score distribution by cluster
- **Heatmap**: Compare multiple pathways

---

## 📊 Step 6: Explore Your Results

### Identify Interesting Clusters

1. **Look at cluster sizes** in statistics table
2. **Check pathway scores** - which clusters are high/low?
3. **Compare different pathways** for the same clusters

### Find Marker Genes

1. **Go to "Marker Genes" tab**
2. **Click "Find Marker Genes"**
3. **Review marker table** - top genes per cluster
4. **Visualize specific genes** in the plot section

### Example Interpretation:
- **Cluster 0**: High stemness, low differentiation → Stem-like cells
- **Cluster 3**: High proliferation, high cell cycle → Actively dividing
- **Cluster 5**: High EMT, high invasion → Mesenchymal/invasive cells

---

## 📤 Step 7: Export Your Results

### Quick Export

1. **Go to "Export" tab**
2. **Select plot type**: e.g., "Current Cluster Plot"
3. **Choose format**: PNG (for presentations), PDF (for papers)
4. **Click "Download Plot"**

### Comprehensive Export

1. **Select data options**:
   - ✅ Cell metadata
   - ✅ Cluster statistics  
   - ✅ CancerSEA scores
2. **Click "Download Excel File"**

### For Publications:
- **Copy methods text** from the Citation section
- **Download high-resolution plots** (300 DPI, PDF format)
- **Save Seurat object** for future analysis

---

## 🎉 Congratulations!

You've completed your first MASIH analysis! You now have:

- **Clustered cells** into distinct populations
- **Characterized functional states** using CancerSEA
- **Identified marker genes** for each cluster
- **Export-ready results** for further analysis

---

## 🔄 Next Steps

### Immediate Next Steps:
1. **Try more pathways**: Calculate additional CancerSEA scores
2. **Explore markers**: Find genes specific to interesting clusters  
3. **Adjust parameters**: Fine-tune clustering resolution

### Advanced Analysis:
1. **Trajectory Analysis**: Trace cell development paths
2. **Pathway Comparison**: Correlate multiple functional states
3. **Cell Cycle Analysis**: Understand proliferation patterns

### Learn More:
- **[Data Preparation Guide](data-preparation.md)** - Optimize your input data
- **[Basic Analysis Tutorial](../tutorials/basic-analysis.md)** - Detailed walkthrough
- **[Advanced Trajectory Tutorial](../tutorials/advanced-trajectory.md)** - Pseudotime analysis

---

## 💡 Tips for Success

### Data Quality:
- **Start with good data** - proper 10X processing is crucial
- **Check QC metrics carefully** - don't skip quality control
- **Remove doublets** if you have tools like DoubletFinder

### Analysis Strategy:
- **Use example data first** - learn the interface before your data
- **Start simple** - basic clustering before advanced features
- **Document parameters** - record what settings work well

### Biological Interpretation:
- **Know your system** - understand the expected cell types
- **Validate findings** - check marker genes make biological sense
- **Use pathway scores wisely** - combine with marker gene analysis

### Technical Tips:
- **Save frequently** - download intermediate results
- **Try different parameters** - clustering resolution affects results
- **Check memory usage** - restart if the app becomes slow

---

## ❓ Common Questions

### "My clusters don't make biological sense"
- **Try different resolution**: 0.1-0.3 for fewer clusters, 0.8-1.2 for more
- **Check quality control**: Poor QC can create artificial clusters
- **Validate with markers**: Do cluster markers match expected cell types?

### "CancerSEA scores seem too high/low"
- **Normal variation**: Scores are relative within your dataset
- **Compare across clusters**: Focus on relative differences
- **Check gene overlap**: Some pathways may have few genes in your data

### "The app is running slowly"
- **Large dataset**: Consider subsampling for exploration
- **Restart browser**: Close and reopen the browser tab
- **Restart R**: In RStudio: Session → Restart R, then relaunch MASIH

### "I can't reproduce my results"
- **Set random seed**: Some functions use randomization
- **Document parameters**: Note all settings used
- **Save Seurat object**: Export and reload to continue later

---

## 🆘 Getting Help

**Stuck on something?**

1. **Check troubleshooting**: [troubleshooting.md](troubleshooting.md)
2. **Review tutorials**: [Basic tutorial](../tutorials/basic-analysis.md) has more details
3. **GitHub issues**: [Report problems](https://github.com/yourusername/masih/issues)
4. **Email support**: your.email@institution.edu

**When asking for help, include:**
- What step you're on
- Any error messages
- Your data type and size
- Screenshots if helpful

---

**Ready to dive deeper?** 🚀 Check out our [detailed tutorials](../tutorials/) or start analyzing your own data!
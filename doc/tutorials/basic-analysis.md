# Basic Analysis Tutorial

**Complete walkthrough: From raw data to biological insights in 30 minutes.**

---

## 🎯 Tutorial Overview

This comprehensive tutorial guides you through a complete MASIH analysis using example cancer data. You'll learn:

- **Data upload and quality control** with biological interpretation
- **Clustering and visualization** of cell populations
- **Marker gene analysis** and biological validation
- **CancerSEA functional characterization** 
- **Results interpretation** and biological insights
- **Export for publication** and further analysis

**⏱️ Time required**: 30-45 minutes  
**Prerequisites**: MASIH installed ([Installation Guide](../user-guide/installation.md))  
**Dataset**: Example melanoma single-cell data (included with MASIH)

---

## 📋 Learning Objectives

By completing this tutorial, you will be able to:

1. **Perform quality control** and interpret QC metrics biologically
2. **Identify cell clusters** and validate them with known markers
3. **Characterize functional states** using CancerSEA pathways
4. **Find cluster-specific genes** and assess their biological relevance
5. **Export publication-ready results** with proper documentation
6. **Interpret results** in the context of cancer biology

---

## 🚀 Step 1: Launch MASIH and Load Data

### Start Your Analysis

```r
# Load MASIH
library(masih)

# Launch the application
run_app()
```

**What you'll see**: The MASIH interface opens in your browser with the sidebar showing analysis modules.

### Load Example Dataset

1. **Navigate to "Data Upload" tab** (first tab in sidebar)
2. **Click "Load Example Data"** button
3. **Select "Melanoma Dataset"** from the dropdown
4. **Click "Load Data"**

**📊 Dataset Information**:
- **Cell type**: Primary melanoma tumor cells
- **Sample size**: ~2,500 cells after QC
- **Source**: Simulated based on published melanoma scRNA-seq data
- **Expected cell types**: Melanoma cells, immune cells, stromal cells

### Initial Data Overview

Once loaded, you'll see:
- **Total cells**: 2,847 cells detected
- **Total genes**: 18,617 genes detected
- **Metadata columns**: Sample, Cell_type, Treatment

**🔍 What to notice**:
- Data is already minimally processed (Cell Ranger output)
- Some metadata is available (sample origin, treatment status)
- Ready for quality control analysis

---

## 🔍 Step 2: Quality Control Analysis

### Navigate to Quality Control

1. **Click "Quality Control" tab** in the sidebar
2. **Review the QC plots** that appear automatically

### Interpret QC Metrics

**📊 Plot 1: Number of Features (Genes) per Cell**
- **X-axis**: Cells (ordered by gene count)
- **Y-axis**: Number of genes detected
- **Expected range**: 500-4,000 genes per cell
- **Red flags**: 
  - Very low (<200): Dead/empty cells
  - Very high (>6,000): Potential doublets

**📊 Plot 2: UMI Counts per Cell**
- **X-axis**: Cells
- **Y-axis**: Total UMI counts
- **Expected range**: 1,000-15,000 UMIs per cell
- **Interpretation**: Higher counts = more RNA captured

**📊 Plot 3: Mitochondrial Gene Percentage**
- **X-axis**: Cells
- **Y-axis**: % mitochondrial genes
- **Expected range**: 5-20% for healthy cells
- **Red flags**: >25% indicates dying/stressed cells

### Set QC Filters

**Recommended filters for this dataset**:

1. **Minimum features per cell**: 200
2. **Maximum features per cell**: 5,000  
3. **Maximum mitochondrial %**: 20%

**Why these cutoffs?**:
- **Min 200 genes**: Ensures sufficient data per cell
- **Max 5,000 genes**: Removes likely doublets
- **Max 20% MT**: Removes dying cells

### Apply Filters

1. **Set the filter values** in the input boxes
2. **Click "Apply Filters"**
3. **Review filtering results**:
   - **Before filtering**: 2,847 cells
   - **After filtering**: ~2,156 cells (76% retained)
   - **Genes remaining**: ~16,500 genes

**🧬 Biological Interpretation**:
- Lost ~691 cells (24%) - normal for cancer samples
- Retained cells have good RNA quality
- Ready for downstream analysis

---

## 🎯 Step 3: Clustering Analysis

### Run Standard Clustering

1. **Navigate to "Cluster Analysis" tab**
2. **Click "Run Standard Analysis"** (uses optimal defaults)
3. **Wait for processing** (~2-3 minutes)

**What happens during processing**:
1. **Normalization**: Log-normalize gene expression
2. **Feature selection**: Find 2,000 highly variable genes
3. **Scaling**: Z-score transformation
4. **PCA**: Reduce to top 10 principal components
5. **Neighbor graph**: Build k-nearest neighbor graph
6. **Clustering**: Leiden algorithm with resolution 0.5
7. **UMAP**: Generate 2D visualization

### Interpret Clustering Results

**📊 UMAP Plot Interpretation**:
- **Number of clusters**: 8 clusters identified (0-7)
- **Cluster sizes**: Range from 89 to 584 cells
- **Spatial organization**: Distinct, well-separated clusters
- **Resolution**: Good balance - not over/under-clustered

**🔍 What each cluster might represent**:
- **Large clusters (0, 1, 2)**: Likely different melanoma cell states
- **Medium clusters (3, 4, 5)**: Possibly immune cells or stromal cells
- **Small clusters (6, 7)**: Rare cell types or transitional states

### Validate Clustering Quality

**📊 Cluster Tree**: Shows hierarchical relationships
- **Close clusters**: Similar expression profiles
- **Distant clusters**: Very different cell types
- **Branch lengths**: Degree of similarity

**📊 Cluster Statistics Table**:
```
Cluster    Cell_Count    Percentage
0          584           27.1%
1          421           19.5%
2          386           17.9%
3          289           13.4%
4          198           9.2%
5          156           7.2%
6          89            4.1%
7          33            1.5%
```

**🧬 Biological Assessment**:
- **Cluster sizes look reasonable** - no tiny or huge clusters
- **Distribution is expected** - few rare types, many common types
- **Good separation in UMAP** - distinct cell populations

---

## 🧬 Step 4: Marker Gene Analysis

### Find Cluster Markers

1. **Navigate to "Marker Genes" tab**
2. **Keep default parameters**:
   - Statistical test: Wilcoxon
   - Min log fold change: 0.25
   - Min percent expressed: 0.1
   - Only positive markers: Yes
3. **Click "Find Marker Genes"**
4. **Wait for analysis** (~1-2 minutes)

### Interpret Top Markers

**📊 Marker Gene Table** - Top 5 genes per cluster:

**Cluster 0** (584 cells, 27.1%):
- **MITF** (p < 0.001, logFC = 1.2) - Melanocyte master regulator
- **TYRP1** (p < 0.001, logFC = 1.1) - Melanin synthesis enzyme  
- **DCT** (p < 0.001, logFC = 1.0) - Melanocyte differentiation
- **MLANA** (p < 0.001, logFC = 0.9) - Melanoma antigen
- **SOX10** (p < 0.001, logFC = 0.8) - Neural crest transcription factor

**🧬 Interpretation**: **Differentiated Melanoma Cells**
- Classic melanocyte/melanoma markers
- Maintained differentiation program
- Likely primary tumor cells

**Cluster 1** (421 cells, 19.5%):
- **MKI67** (p < 0.001, logFC = 2.1) - Proliferation marker
- **TOP2A** (p < 0.001, logFC = 1.8) - DNA replication
- **BIRC5** (p < 0.001, logFC = 1.6) - Anti-apoptotic factor
- **CCNB1** (p < 0.001, logFC = 1.4) - Cell cycle regulation
- **AURKA** (p < 0.001, logFC = 1.3) - Mitotic kinase

**🧬 Interpretation**: **Proliferating Melanoma Cells**
- High proliferation markers
- Active cell cycle
- Aggressive tumor cells

**Cluster 2** (386 cells, 17.9%):
- **AXL** (p < 0.001, logFC = 1.5) - Receptor tyrosine kinase
- **WNT5A** (p < 0.001, logFC = 1.4) - Wnt signaling
- **TWIST1** (p < 0.001, logFC = 1.2) - EMT transcription factor
- **VIM** (p < 0.001, logFC = 1.1) - Mesenchymal marker
- **CDH2** (p < 0.001, logFC = 1.0) - N-cadherin

**🧬 Interpretation**: **Mesenchymal-like Melanoma Cells**
- EMT-associated genes
- Loss of epithelial features
- Potential for invasion/metastasis

**Cluster 3** (289 cells, 13.4%):
- **CD3D** (p < 0.001, logFC = 2.0) - T cell marker
- **CD3E** (p < 0.001, logFC = 1.9) - T cell marker
- **CD2** (p < 0.001, logFC = 1.7) - T cell activation
- **IL7R** (p < 0.001, logFC = 1.5) - T cell survival
- **LCK** (p < 0.001, logFC = 1.4) - T cell signaling

**🧬 Interpretation**: **T Cells**
- Clear T cell identity
- Tumor-infiltrating lymphocytes
- Important for immune response

### Validate with Known Markers

**Use the gene visualization tool**:
1. **Select "MITF"** from marker dropdown
2. **Choose "Feature Plot"** 
3. **Click "Update Plot"**

**Expected result**: MITF expression should be highest in Cluster 0, confirming melanoma cell identity.

**Try additional validations**:
- **MKI67**: Should be high in Cluster 1 (proliferating cells)
- **CD3D**: Should be high in Cluster 3 (T cells)
- **AXL**: Should be high in Cluster 2 (mesenchymal cells)

---

## 🧠 Step 5: CancerSEA Functional Analysis

### Calculate Stemness Scores

1. **Navigate to "CancerSEA" tab**
2. **Select "Stemness"** from pathway dropdown
3. **Click "Calculate/Update"**
4. **Wait for calculation** (~30 seconds)

### Interpret Stemness Results

**📊 Feature Plot**: 
- **Color scale**: Blue (low) to Red (high) stemness
- **Expected pattern**: Some clusters higher than others
- **Biological meaning**: Red cells have stem-cell-like properties

**📊 Violin Plot by Cluster**:
- **Y-axis**: Stemness score
- **X-axis**: Cluster identity
- **Interpretation**: Compare median scores between clusters

**Expected Results**:
- **Cluster 0**: Moderate stemness (differentiated but plastic)
- **Cluster 1**: Lower stemness (committed to proliferation)
- **Cluster 2**: Higher stemness (EMT-associated stemness)
- **Cluster 3**: Very low stemness (immune cells)

### Calculate Additional Pathways

**Add Proliferation scores**:
1. **Select "Proliferation"** from dropdown
2. **Click "Calculate/Update"**
3. **Compare with clustering results**

**Expected Results**:
- **Cluster 1**: Highest proliferation scores
- **Cluster 0, 2**: Moderate proliferation
- **Cluster 3**: Low proliferation (immune cells)

**Add EMT scores**:
1. **Select "EMT"** from dropdown  
2. **Click "Calculate/Update"**
3. **Look for EMT-high clusters**

**Expected Results**:
- **Cluster 2**: Highest EMT scores
- **Cluster 0, 1**: Lower EMT scores
- **Cluster 3**: Very low EMT (immune cells)

### Compare Multiple Pathways

**📊 Pathway Heatmap** (appears after calculating ≥2 pathways):
- **Rows**: Clusters
- **Columns**: CancerSEA pathways
- **Colors**: Red (high activity) to Blue (low activity)
- **Patterns**: Each cluster has distinct functional profile

**🧬 Biological Interpretation**:
- **Cluster 0**: Balanced profile (differentiated melanoma)
- **Cluster 1**: Proliferation-dominant (aggressive melanoma)
- **Cluster 2**: EMT/Stemness-high (invasive melanoma)
- **Cluster 3**: All pathways low (immune infiltrate)

---

## 📊 Step 6: Pathway Comparison Analysis

### Run Comparative Analysis

1. **Navigate to "Pathway Comparison" tab**
2. **Select multiple pathways**: Stemness, Proliferation, EMT
3. **Click "Run Comparison"**
4. **Wait for correlation analysis**

### Interpret Correlation Results

**📊 Correlation Matrix**:
- **Stemness vs EMT**: Positive correlation (r = 0.65)
  - *Biological meaning*: EMT promotes stemness
- **Proliferation vs Stemness**: Negative correlation (r = -0.42)
  - *Biological meaning*: Proliferating cells less stem-like
- **Proliferation vs EMT**: Weak correlation (r = 0.1)
  - *Biological meaning*: Independent processes

**📊 Pathway by Cluster Bar Plot**:
- **X-axis**: Clusters
- **Y-axis**: Mean pathway scores
- **Colors**: Different pathways
- **Error bars**: Standard error

**🧬 Key Insights**:
1. **Cluster 2**: High stemness + EMT = invasive potential
2. **Cluster 1**: High proliferation + low stemness = aggressive but less plastic
3. **Cluster 0**: Balanced profile = stable differentiated state

---

## 🔄 Step 7: Cell Cycle Analysis

### Analyze Cell Cycle States

1. **Navigate to "Cell Cycle" tab**
2. **Review cell cycle distribution plot**
3. **Examine phase distribution by cluster**

### Interpret Cell Cycle Results

**📊 Cell Cycle Phase Distribution**:
- **G1 phase**: 65% of cells (resting/gap phase)
- **S phase**: 20% of cells (DNA synthesis)
- **G2/M phase**: 15% of cells (mitosis preparation)

**📊 Phase Distribution by Cluster**:
- **Cluster 1**: Enriched in S and G2/M phases (proliferating)
- **Cluster 0, 2**: Mostly G1 phase (quiescent)
- **Cluster 3**: G1 dominant (immune cells typically quiescent)

**🧬 Biological Validation**:
- **Cluster 1** high proliferation scores + high S/G2M = consistent
- **Immune clusters** in G1 = expected (not actively dividing)
- **EMT cluster** (2) mostly G1 = focus on invasion not proliferation

---

## 📈 Step 8: Results Summary and Interpretation

### Biological Summary

**🧬 Cell Population Identified**:

1. **Cluster 0 - Differentiated Melanoma Cells (27.1%)**:
   - **Markers**: MITF, TYRP1, MLANA
   - **Function**: Maintained melanocyte program
   - **State**: Differentiated, moderate stemness
   - **Clinical relevance**: Primary tumor bulk

2. **Cluster 1 - Proliferating Melanoma Cells (19.5%)**:
   - **Markers**: MKI67, TOP2A, CCNB1
   - **Function**: Active cell division
   - **State**: High proliferation, low stemness
   - **Clinical relevance**: Aggressive growth, therapy target

3. **Cluster 2 - Mesenchymal-like Melanoma Cells (17.9%)**:
   - **Markers**: AXL, WNT5A, TWIST1
   - **Function**: EMT program active
   - **State**: High stemness, invasive potential
   - **Clinical relevance**: Metastasis risk, therapy resistance

4. **Cluster 3 - Tumor-infiltrating T Cells (13.4%)**:
   - **Markers**: CD3D, CD3E, IL7R
   - **Function**: Immune surveillance
   - **State**: Quiescent, low cancer pathways
   - **Clinical relevance**: Immunotherapy targets

### Clinical Implications

**🎯 Therapeutic Targets**:
- **Cluster 1**: Target proliferation (CDK4/6 inhibitors)
- **Cluster 2**: Target EMT/stemness (AXL inhibitors)
- **Cluster 3**: Enhance T cell function (checkpoint inhibitors)

**📊 Prognostic Markers**:
- **High Cluster 2 proportion**: Poor prognosis (invasive potential)
- **High Cluster 3 proportion**: Better prognosis (immune active)
- **Cluster 1 dominance**: Aggressive but potentially targetable

---

## 📤 Step 9: Export Results

### Export Publication Plots

1. **Navigate to "Export" tab**
2. **Select plot type**: "Current Cluster Plot"
3. **Set parameters**:
   - Format: PDF
   - Resolution: 300 DPI
   - Size: 8×6 inches
   - White background: Yes
4. **Click "Download Plot"**

### Export Data Tables

1. **Select data export options**:
   - ✅ Cell metadata
   - ✅ Cluster statistics
   - ✅ CancerSEA scores
   - ✅ Marker genes
2. **Click "Download Excel File"**

### Generate Methods Text

**Copy the auto-generated methods text**:
```
Single-cell RNA sequencing data were analyzed using MASIH (Modular Analysis 
Shiny Interface for Heterogeneity). Quality control was performed to filter 
cells with <200 or >5000 detected genes and >20% mitochondrial gene expression. 
Data were log-normalized and scaled, followed by principal component analysis. 
Clustering was performed using the Leiden algorithm with resolution 0.5. 
Functional states were characterized using CancerSEA pathway scoring...
```

---

## 🎯 Key Learning Points

### Technical Skills Gained

1. **Quality Control**: Set biologically meaningful filters
2. **Clustering**: Interpret UMAP plots and validate clusters
3. **Marker Analysis**: Find and validate cluster-specific genes
4. **Functional Analysis**: Use CancerSEA for biological insights
5. **Export**: Generate publication-ready outputs

### Biological Insights

1. **Tumor Heterogeneity**: Identified 4 distinct cell populations
2. **Functional States**: Linked molecular profiles to biological functions
3. **Clinical Relevance**: Connected findings to therapeutic opportunities
4. **Immune Context**: Characterized tumor-immune interactions

### Best Practices Learned

1. **Start with QC**: Good quality control is crucial
2. **Validate findings**: Use known markers to confirm clusters
3. **Think biologically**: Interpret results in biological context
4. **Multiple approaches**: Combine clustering, markers, and pathways
5. **Document everything**: Export parameters and results

---

## 🚀 Next Steps

### Immediate Next Steps

1. **Try your own data**: Apply these skills to your research
2. **Explore more pathways**: Calculate additional CancerSEA scores
3. **Advanced analysis**: Try trajectory analysis tutorial
4. **Parameter optimization**: Fine-tune clustering resolution

### Advanced Tutorials

1. **[Advanced Trajectory Analysis](advanced-trajectory.md)**: Pseudotime inference
2. **[Comparative Pathway Analysis](comparative-pathways.md)**: Multi-sample comparison
3. **Cell-cell Communication**: Coming in future versions

### Further Learning

1. **Seurat tutorials**: [satijalab.org](https://satijalab.org/seurat/)
2. **CancerSEA database**: [biocc.hrbmu.edu.cn/CancerSEA](http://biocc.hrbmu.edu.cn/CancerSEA/)
3. **Single-cell best practices**: Luecken & Theis, 2019

---

## ❓ Tutorial Q&A

### "My results look different"

**Possible reasons**:
- **Different random seeds**: Some variation is normal
- **Different parameters**: Check QC filters and clustering resolution
- **Different data**: Results will vary with different datasets

### "I don't see the expected cell types"

**Troubleshooting**:
- **Check marker genes**: Do they match known cell type markers?
- **Adjust resolution**: Try higher (more clusters) or lower (fewer clusters)
- **Review QC**: Poor quality cells can create artificial clusters

### "CancerSEA scores seem wrong"

**Remember**:
- **Scores are relative**: Compare between clusters, not absolute values
- **Dataset specific**: Scores depend on your specific data
- **Biological context**: Consider your cancer type and experimental conditions

---

**Congratulations! 🎉** You've completed your first comprehensive MASIH analysis. You're now ready to tackle your own single-cell cancer data with confidence!
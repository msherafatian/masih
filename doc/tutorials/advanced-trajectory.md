# Advanced Trajectory Analysis Tutorial

**Master pseudotime inference to trace cancer cell evolution and developmental paths.**

---

## 🎯 Tutorial Overview

Trajectory analysis reveals how cells transition between states over time, making it crucial for understanding cancer progression, treatment resistance, and cellular development. This advanced tutorial covers:

- **Trajectory inference theory** and biological applications
- **Cell type selection** for meaningful trajectories
- **Slingshot algorithm** parameter optimization
- **Pseudotime interpretation** and validation
- **Gene expression dynamics** along trajectories
- **Clinical implications** of trajectory patterns

**⏱️ Time required**: 45-60 minutes  
**Prerequisites**: Completed [Basic Analysis Tutorial](basic-analysis.md)  
**Dataset**: Melanoma progression dataset (differentiation → EMT trajectory)

---

## 📚 Trajectory Analysis Background

### What is Trajectory Analysis?

**Trajectory inference** reconstructs developmental or progression paths from single-cell snapshots. Unlike clustering (which identifies discrete states), trajectory analysis:

- **Connects cell states** along continuous paths
- **Orders cells** in pseudotime (inferred temporal order)
- **Identifies transition genes** that change along trajectories
- **Reveals branching points** where cells make fate decisions

### Biological Applications in Cancer

**🧬 Cancer Progression Trajectories**:
1. **Tumor evolution**: Normal → Pre-malignant → Malignant
2. **EMT progression**: Epithelial → Partial EMT → Mesenchymal
3. **Drug resistance**: Sensitive → Intermediate → Resistant
4. **Metastatic cascade**: Primary → Circulating → Metastatic

**🎯 Clinical Relevance**:
- **Identify early events** in cancer progression
- **Find intervention points** to block progression
- **Predict resistance mechanisms** before they emerge
- **Personalize therapy** based on progression state

### When to Use Trajectory Analysis

**✅ Good candidates**:
- Cell populations with **continuous variation**
- **Known biological relationships** between cell types
- **Sufficient cell numbers** (>100 per major state)
- **Clear developmental/progression question**

**❌ Poor candidates**:
- **Discrete, unrelated** cell types
- **Technical batch effects** dominating variation
- **Too few cells** for stable inference
- **No biological expectation** of continuity

---

## 🚀 Step 1: Data Preparation and Setup

### Load Processed Data

```r
# Launch MASIH with processed data from basic tutorial
library(masih)
run_app()

# Or load fresh example data
# Navigate to "Data Upload" → "Load Example Data" → "Melanoma Progression Dataset"
```

**📊 Dataset characteristics**:
- **3,200 cells** from melanoma progression model
- **Expected trajectory**: Differentiated → EMT-like → Invasive
- **Key cell types**: Melanoma cells at different progression stages
- **Controls**: Immune infiltrate (should not be on main trajectory)

### Verify Clustering Results

**Navigate to "Cluster Analysis" tab**:
- Should show **6-8 clusters** from previous analysis
- **Expected pattern**: Gradual transition between melanoma clusters
- **Validation**: Check that immune clusters are distinct

**🔍 Key clusters for trajectory**:
- **Cluster 0**: Differentiated melanoma (MITF+, TYRP1+)
- **Cluster 2**: Transitional state (moderate EMT markers)
- **Cluster 4**: EMT-like melanoma (AXL+, TWIST1+)
- **Cluster 6**: Invasive melanoma (high stemness + EMT)

---

## 📈 Step 2: Configure Trajectory Analysis

### Navigate to Trajectory Module

1. **Click "Trajectory Analysis" tab**
2. **Review the parameter panel**
3. **Check cell type selection options**

### Select Cell Types for Trajectory

**🎯 Strategic cell type selection**:

1. **Annotation column**: Select "seurat_clusters"
2. **Select cell identities**: Choose clusters 0, 2, 4, 6
   - **Rationale**: These represent the melanoma progression path
   - **Exclude**: Immune clusters (1, 3, 5) - not part of progression

**Why exclude immune cells?**:
- **Different biology**: Immune cells follow different developmental paths
- **Confounding signal**: May create artificial trajectories
- **Focus**: Keep analysis on cancer cell progression

### Set Trajectory Parameters

**📊 Key parameters to configure**:

1. **Number of PCs**: 15 (default: 10)
   - **Rationale**: Capture more variation in complex progression
   - **Rule**: Use more PCs for complex trajectories

2. **Starting cluster**: "Cluster 0" 
   - **Biological rationale**: Most differentiated state
   - **Markers**: High MITF, TYRP1 (classic melanocyte program)

3. **Ending cluster**: "Auto-detect"
   - **Let Slingshot determine**: Most progressed state
   - **Expected**: Should identify Cluster 6 (invasive)

**💡 Parameter tips**:
- **Start conservative**: Begin with defaults, then optimize
- **Biological knowledge**: Use known biology to guide choices
- **Iterative process**: May need multiple attempts

---

## 🔬 Step 3: Run Trajectory Analysis

### Execute Slingshot Analysis

1. **Review filter preview**: Should show ~1,800 selected cells
2. **Click "Run Trajectory Analysis"**
3. **Monitor progress**: Should take 2-3 minutes
4. **Wait for completion**: "Trajectory analysis completed successfully!"

### Initial Results Overview

**📊 What you'll see**:
- **Trajectory visualization**: Curves overlaid on UMAP
- **Pseudotime coloring**: Cells colored by inferred time
- **Lineage paths**: Multiple potential trajectories
- **Statistics summary**: Trajectory characteristics

**🔍 Expected trajectory pattern**:
- **Start**: Cluster 0 (differentiated melanoma)
- **Middle**: Cluster 2 (transitional)
- **Branches**: 
  - Path 1: → Cluster 4 (EMT-focused)
  - Path 2: → Cluster 6 (invasive/stemness-focused)

---

## 📊 Step 4: Interpret Trajectory Results

### Main Trajectory Visualization

**📈 Trajectory Plot Interpretation**:

1. **Smooth curves**: Well-fitted trajectories (good!)
2. **Starting point**: Should begin in Cluster 0
3. **Branching pattern**: 
   - **Single trunk**: Early progression shared
   - **Branch point**: Around Cluster 2 (transition state)
   - **Two endpoints**: Different final states

**🧬 Biological interpretation**:
- **Trunk (Cluster 0→2)**: Early EMT initiation
- **Branch 1 (→Cluster 4)**: EMT-dominant path
- **Branch 2 (→Cluster 6)**: Stemness + invasion path

### Pseudotime Analysis

**📊 Pseudotime Plot Features**:
- **Color gradient**: Blue (early) → Red (late)
- **Smooth transition**: Gradual color change is good
- **No jumps**: Avoid sudden color changes within clusters

**🔍 Validation checks**:
1. **Cluster 0**: Should be mostly blue (early pseudotime)
2. **Cluster 6**: Should be mostly red (late pseudotime)
3. **Intermediate clusters**: Should show gradient colors

### Stratified Trajectory View

**📊 Panel-by-cluster view**:
- **Each panel**: One cluster's cells on trajectory
- **Trajectory curves**: Consistent across panels
- **Cell positions**: Validate cluster assignments make sense

**Key insights**:
- **Cluster 0**: Cells at trajectory start
- **Cluster 2**: Cells at branching point
- **Clusters 4 & 6**: Cells at trajectory endpoints

---

## 🧬 Step 5: Gene Expression Dynamics

### Identify Trajectory-Associated Genes

**🔍 What to look for**:
- **Early genes**: High at trajectory start, decrease over time
- **Late genes**: Low initially, increase toward endpoints
- **Branch-specific genes**: High in one branch, low in another

**📊 Expected patterns**:

**Early progression genes (decreasing)**:
- **MITF**: Melanocyte master regulator (lost during EMT)
- **TYRP1**: Melanin synthesis (differentiation marker)
- **MLANA**: Melanoma antigen (differentiated state)

**Late progression genes (increasing)**:
- **AXL**: RTK associated with progression
- **TWIST1**: EMT transcription factor
- **VIM**: Mesenchymal marker

**Branch-specific patterns**:
- **Branch 1 (EMT)**: High SNAI1, ZEB1, CDH2
- **Branch 2 (Invasion)**: High MMP2, TIMP1, stemness markers

### Validate Key Transition Genes

**Use marker gene visualization**:

1. **Navigate to "Marker Genes" tab**
2. **Select gene**: "MITF"
3. **Plot type**: "Feature Plot"
4. **Expected**: High in Cluster 0, decreasing along trajectory

**Try additional genes**:
- **AXL**: Should increase along trajectory
- **TWIST1**: Should be high in EMT branch
- **SOX2**: Should be high in stemness branch

---

## 🎯 Step 6: Functional Analysis Along Trajectory

### CancerSEA Pathway Dynamics

**Navigate back to trajectory results**:
1. **Overlay CancerSEA scores** on trajectory plot
2. **Expected pathway changes**:

**Stemness pathway**:
- **Early (Cluster 0)**: Moderate stemness
- **Middle (Cluster 2)**: Increased stemness
- **Late (Cluster 6)**: Highest stemness

**EMT pathway**:
- **Early (Cluster 0)**: Low EMT
- **Transition (Cluster 2)**: Moderate EMT
- **Branch 1 (Cluster 4)**: High EMT
- **Branch 2 (Cluster 6)**: Very high EMT

**Differentiation pathway**:
- **Pattern**: Should decrease along trajectory
- **Early high**: Cluster 0 (differentiated state)
- **Late low**: Clusters 4 & 6 (dedifferentiated)

### Trajectory Statistics Summary

**📊 Key metrics to review**:

```
Trajectory Statistics:
- Number of lineages: 2
- Total trajectory length: 12.3 units
- Branching point: Pseudotime 0.34
- Cells on trajectory: 1,847 (85.7%)
- Average curve length: 6.8 units
```

**🔍 Interpretation**:
- **2 lineages**: Confirms branching biology
- **Branching at 0.34**: Early in progression (good)
- **85.7% on trajectory**: Most selected cells fit trajectory (excellent)

---

## 🔬 Step 7: Advanced Trajectory Validation

### Cross-Reference with Known Biology

**✅ Biological validation checklist**:

1. **Trajectory direction**: Differentiated → EMT ✓
2. **Branching makes sense**: EMT vs stemness paths ✓
3. **Marker gene patterns**: MITF decreases, AXL increases ✓
4. **Pathway dynamics**: Stemness and EMT increase ✓
5. **Cell cycle patterns**: Check if progression affects proliferation

### Statistical Validation

**📊 Trajectory confidence measures**:
- **Curve smoothness**: Avoid jagged trajectories
- **Cell density**: Avoid gaps in trajectory
- **Branch confidence**: Strong separation at branches

**🔍 Red flags to watch for**:
- **Artificial loops**: May indicate parameter issues
- **Disconnected segments**: Suggests missing intermediate states
- **Reversed biology**: Check gene expression makes sense

### Alternative Parameter Testing

**If results seem suboptimal, try**:

1. **Different PC numbers**:
   - More PCs (20-25): For complex datasets
   - Fewer PCs (8-12): For simpler trajectories

2. **Different starting points**:
   - Try different starting clusters
   - Compare trajectory shapes

3. **Cell selection refinement**:
   - Remove transition doublets
   - Include additional intermediate clusters

---

## 🎯 Step 8: Clinical Interpretation

### Progression Timing Analysis

**🕐 Pseudotime interpretation**:
- **Pseudotime 0-0.3**: Early progression (targetable)
- **Pseudotime 0.3-0.6**: Transition phase (intervention window)
- **Pseudotime 0.6-1.0**: Advanced progression (poor prognosis)

**🎯 Clinical implications**:
- **Early detection**: Focus on pseudotime 0-0.3 markers
- **Intervention timing**: Target transition phase
- **Prognosis**: Late pseudotime = worse outcomes

### Therapeutic Target Identification

**📊 Trajectory-based target discovery**:

**Early intervention targets** (Pseudotime 0-0.3):
- **Maintain differentiation**: MITF pathway activation
- **Prevent EMT**: Target early EMT signals
- **Block progression**: Interfere with AXL signaling

**Late intervention targets** (Pseudotime 0.6-1.0):
- **Target stemness**: CSC-specific therapies
- **Block invasion**: MMP inhibitors
- **Overcome resistance**: Combination therapies

### Prognostic Signatures

**🧬 Trajectory-derived biomarkers**:
1. **Early progression score**: Genes changing early in trajectory
2. **EMT branch score**: Branch 1-specific genes
3. **Stemness branch score**: Branch 2-specific genes
4. **Transition risk score**: Genes at branching point

---

## 📤 Step 9: Export Trajectory Results

### Export Trajectory Plots

1. **Navigate to "Export" tab**
2. **Select "Main Trajectory"** plot type
3. **Set high quality parameters**:
   - Format: PDF
   - Resolution: 600 DPI
   - Size: 10×8 inches
4. **Download plot**

### Export Trajectory Data

**Comprehensive trajectory export**:
1. **Download trajectory data** (includes pseudotime scores)
2. **Export gene expression** along trajectories
3. **Save analysis parameters** for reproducibility

### Generate Trajectory Methods

**Sample methods text**:
```
Trajectory analysis was performed using Slingshot algorithm implemented 
in MASIH. Melanoma cell clusters (0, 2, 4, 6) were selected based on 
expression of melanocyte markers and EMT signatures. Principal component 
analysis using 15 dimensions was performed prior to trajectory inference. 
Cluster 0 was specified as the starting point based on high expression 
of differentiation markers (MITF, TYRP1). Two trajectory branches were 
identified, representing EMT-focused and stemness-focused progression paths...
```

---

## 💡 Advanced Tips and Best Practices

### Trajectory Design Principles

**🎯 Planning your trajectory analysis**:

1. **Clear biological hypothesis**: What progression do you expect?
2. **Sufficient cell sampling**: >100 cells per major state
3. **Quality over quantity**: Better to analyze fewer, well-defined cell types
4. **Validation strategy**: How will you confirm biological relevance?

### Common Pitfalls and Solutions

**❌ Problem**: Trajectory connects unrelated cell types
**✅ Solution**: More careful cell type selection

**❌ Problem**: Artificial branching in homogeneous population  
**✅ Solution**: Reduce PCA dimensions or change resolution

**❌ Problem**: Trajectory reverses expected biology
**✅ Solution**: Check starting cluster selection and gene expression validation

**❌ Problem**: Gaps in trajectory
**✅ Solution**: Include more intermediate cell states or adjust parameters

### Parameter Optimization Strategy

**🔄 Iterative approach**:

1. **Start simple**: Default parameters, obvious cell types
2. **Validate biology**: Check known marker genes
3. **Refine selection**: Add/remove cell types as needed
4. **Optimize parameters**: Adjust PCs and resolution
5. **Cross-validate**: Compare multiple parameter sets

### Integration with Other Analyses

**🔗 Combine trajectory with**:
- **Differential expression**: Find trajectory-associated genes
- **Gene set enrichment**: Identify pathway changes
- **Cell-cell communication**: How signaling changes along trajectory
- **Spatial analysis**: If spatial coordinates available

---

## 🚀 Next Steps and Advanced Applications

### Immediate Extensions

1. **Multiple trajectories**: Analyze tumor-immune interaction trajectories
2. **Time-course data**: If available, compare with real time points
3. **Drug response**: How do treatments affect trajectory paths?
4. **Patient stratification**: Use trajectory positions for prognosis

### Advanced Trajectory Methods

**🔬 Beyond Slingshot**:
- **RNA velocity**: Directional information from splicing
- **Palantir**: Multifurcating trajectories
- **CytoTRACE**: Differentiation potential scoring
- **scVelo**: RNA velocity-informed trajectories

### Research Applications

**📚 Trajectory analysis for**:
- **Development biology**: Normal developmental paths
- **Cancer evolution**: Tumor progression modeling
- **Drug resistance**: Resistance emergence paths
- **Regenerative medicine**: Reprogramming trajectories

---

## ❓ Troubleshooting Trajectory Analysis

### "Trajectory looks wrong"

**Diagnostic steps**:
1. **Check cell selection**: Are all selected cells biologically related?
2. **Validate markers**: Do known genes show expected patterns?
3. **Review clustering**: Is underlying clustering reasonable?
4. **Try different parameters**: Adjust PCs, starting points

### "No clear trajectory structure"

**Possible causes**:
- **Discrete cell types**: Some biology is truly discrete, not continuous
- **Technical variation**: Batch effects masking biological signal
- **Insufficient resolution**: Need more cells or better clustering
- **Wrong question**: Maybe trajectory analysis isn't appropriate

### "Multiple conflicting trajectories"

**Interpretation**:
- **May be real**: Complex biology often has multiple paths
- **Focus on question**: Which trajectory answers your biological question?
- **Validate separately**: Analyze each trajectory independently

---

## 🎓 Learning Summary

### Key Concepts Mastered

1. **Trajectory theory**: Understanding pseudotime and developmental inference
2. **Parameter optimization**: Systematic approach to trajectory analysis
3. **Biological validation**: Ensuring results match known biology
4. **Clinical interpretation**: Translating trajectories to therapeutic insights
5. **Quality control**: Recognizing good vs problematic trajectories

### Biological Insights Gained

1. **Cancer progression**: EMT and stemness as parallel paths
2. **Intervention windows**: Early transition states as therapeutic targets
3. **Heterogeneity**: Multiple progression routes in same tumor
4. **Biomarker development**: Trajectory-based prognostic signatures

### Technical Skills Developed

1. **Cell type curation**: Strategic selection for meaningful trajectories
2. **Parameter tuning**: Optimizing analysis for biological questions
3. **Validation workflows**: Confirming trajectory biological relevance
4. **Integration analysis**: Combining trajectories with functional data

---

**Congratulations! 🎉** You've mastered advanced trajectory analysis in MASIH. You can now trace cellular evolution, identify intervention points, and discover progression biomarkers in your own cancer research!
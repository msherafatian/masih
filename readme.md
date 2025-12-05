# MASIH: Modular Analysis Suite for Interactive Heterogeneity in single cells<img src="https://raw.githubusercontent.com/msherafatian/masih/main/man/figures/logo.PNG" align="right" width="120" />

[![DOI](https://zenodo.org/badge/1036947264.svg)](https://doi.org/10.5281/zenodo.17824082)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E=4.0-blue.svg)](https://cran.r-project.org/)
[![GitHub issues](https://img.shields.io/github/issues/msherafatian/masih.svg)](https://github.com/msherafatian/masih/issues)
**An interactive R Shiny application for comprehensive single-cell RNA sequencing analysis to focus on the multidimensional biological portrait of each cell.**

> 🐍 **Python version available!** Check out [MASIH-Python](https://github.com/msherafatian/masih-python) for a Dash/Scanpy implementation with Docker support.
---

## 📋 Table of Contents

- [Overview](#-overview)
- [Key Features](#-key-features)
- [Why MASIH?](#-why-masih)
- [Quick Start](#-quick-start)
- [Documentation](#-documentation)
- [Analysis Modules](#-analysis-modules)
- [Supported Data Formats](#-supported-data-formats)
- [Use Cases](#-use-cases)
- [System Requirements](#-system-requirements)
- [Contributing](#-contributing)
- [Citation](#-citation)
- [Roadmap](#-roadmap)
- [Related Projects](#-related-projects)
- [Acknowledgments](#-acknowledgments)
- [License](#-license)

---

## 📖 Overview

MASIH provides a user-friendly web interface for analyzing single-cell RNA sequencing (scRNA-seq) data to focus on the **multidimensional biological portrait of each cancer cell**.  

It offers a complete workflow from **raw and previously analysed single-cell gene expression** to **many cell biological insights**.

> 💡 **Looking for Python?** A Python/Dash implementation is available at [MASIH-Python](https://github.com/msherafatian/masih-python), maintaining architectural parity with this R version.

---

## 🌟 Key Features

### 📂 Flexible Input & Preprocessing
- Accepts both raw 10x Genomics outputs and processed Seurat objects  
- Tailored workflows for each input type, preserving prior results and metadata  
- Automatically detects and completes missing analysis steps  

### 🔍 Cellular Landscape Mapping
- **High-Resolution Clustering** — Graph-based clustering with statistical validation for resolving rare subpopulations  
- **Dimensionality Reduction** — PCA, t-SNE, UMAP for revealing data structure  
- **Marker Gene Profiling** — Differential expression analysis with integrated marker references  

### 🧠 Functional & Dynamic State Inference
- **Functional State Characterization** — CancerSEA-based scoring of cancer-related pathways  
- **Trajectory Mapping** — Pseudotime analysis to trace developmental or oncogenic progression  
- **Cell Cycle Deconvolution** — Phase scoring and integration into downstream analyses  

### 📤 Insight Sharing & Reporting
- Export publication-quality plots and structured analysis results  
- Auto-generate dataset-specific methods text for manuscripts  

---

## 🎯 Why MASIH?

- **Cancer-Focused**: Specialized cancer research tools with CancerSEA integration  
- **No Coding Required**: Accessible to all researchers  
- **Comprehensive Workflow**: From raw data to publication-ready figures  
- **Reproducible**: Exports parameters and generates methods text  
- **Modular Design**: Easily extend with new analysis modules  
- **🐍 Python Version**: Also available as [MASIH-Python](https://github.com/msherafatian/masih-python) for Python users


---

## 🚀 Quick Start

### Installation

```r
# Install from GitHub
devtools::install_github("camlab-bioml/cancersea")
devtools::install_github("msherafatian/masih")
```

### Launch MASIH

```r
library(cancersea)
library(masih)
run_app()
```

### Load Example Data

```r
# Load provided example dataset
data("example_cancer_data")

# Or launch with pre-loaded data
run_app(seurat_object = example_cancer_data)
```

---

## 📖 Documentation

### User Guides
- **[Installation Guide](docs/user-guide/installation.md)** - Detailed installation instructions
- **[Getting Started](docs/user-guide/getting-started.md)** - Your first analysis with MASIH
- **[Data Preparation](docs/user-guide/data-preparation.md)** - How to prepare your data
- **[Troubleshooting](docs/user-guide/troubleshooting.md)** - Common issues and solutions

### Tutorials
- **[Basic Analysis](docs/tutorials/basic-analysis.md)** - Complete walkthrough (30 min)
- **[Advanced Trajectory Analysis](docs/tutorials/advanced-trajectory.md)** - Pseudotime analysis
- **[Comparative Pathway Analysis](docs/tutorials/comparative-pathways.md)** - CancerSEA workflows

### For Developers
- **[Contributing Guidelines](docs/developer/contributing.md)** - How to contribute
- **[Architecture Overview](docs/developer/architecture.md)** - Technical details

---

## 🔬 Analysis Modules

### Core Analysis
- **Data Upload**: Multiple format support (10X, Seurat, matrices)
- **Quality Control**: Interactive filtering and validation
- **Clustering**: Graph-based clustering with resolution optimization
- **Marker Genes**: Statistical testing with multiple methods

### Cancer-Specific Features
- **CancerSEA Integration**: 14 functional state pathways
- **Pathway Comparison**: Correlation and comparative analysis
- **Trajectory Analysis**: Pseudotime inference for cancer progression
- **Cell Cycle Scoring**: G1/S/G2M phase identification

### Visualization & Export
- **Interactive Plots**: Plotly-powered visualizations
- **High-Quality Export**: Publication-ready figures
- **Comprehensive Data Export**: Excel, CSV, Seurat objects
- **Methods Generation**: Automatic methods text for papers

---

## 📊 Supported Data Formats

- **10X Genomics**: Cell Ranger outputs (h5, mtx, barcodes, features)
- **Seurat Objects**: .rds files containing processed Seurat objects
- **SingleCellExperiment**: SCE objects from Bioconductor
- **Gene Expression Matrices**: CSV/TSV format

---

## 🧪 Example Workflows

### Basic Cancer Analysis (30 minutes)
1. Upload 10X data → 2. Quality control → 3. Clustering → 4. CancerSEA analysis → 5. Export results

### Advanced Trajectory Analysis (45 minutes)
1. Basic workflow → 2. Cell type selection → 3. Trajectory inference → 4. Pseudotime analysis → 5. Publication figures

### Comparative Study (60 minutes)
1. Load multiple samples → 2. Batch correction → 3. Comparative clustering → 4. Pathway comparison → 5. Statistical analysis

---

## 🏥 Use Cases

MASIH is designed for cancer researchers studying:

- **Tumor Heterogeneity**: Identify and characterize cancer cell subpopulations
- **Treatment Response**: Analyze single-cell responses to therapy
- **Cancer Progression**: Trace developmental trajectories and metastasis
- **Functional States**: Characterize stemness, invasion, drug resistance
- **Microenvironment**: Analyze tumor-immune interactions

---

## 📋 System Requirements

- **R**: Version 4.0.0 or higher
- **Operating System**: Windows 10+, macOS 10.14+, or Linux
- **Memory**: 8GB RAM minimum (16GB recommended for large datasets)
- **Storage**: 2GB free space for installation

---

## 🤝 Contributing

We welcome contributions! Please read our [Contributing Guidelines](docs/developer/contributing.md) for details.

### Ways to Contribute
- 🐛 Report bugs and request features
- 📖 Improve documentation
- 🧪 Add new analysis modules
- 🎨 Enhance user interface
- 🧬 Add new pathway databases

---

## 📄 Citation

If you use MASIH, please cite the Zenodo DOI:

[![DOI](https://zenodo.org/badge/1036947264.svg)](https://doi.org/10.5281/zenodo.17824082)

**Concept DOI (latest release):**  
https://doi.org/10.5281/zenodo.17824082

**Version-specific DOI (for reproducibility, e.g., manuscript-linked release):**  
v1.0.0 → https://doi.org/10.5281/zenodo.17824081

### BibTeX
```bibtex
@software{masih_2024,
  author       = {Sherafatian, Masih},
  title        = {MASIH: Modular Analysis Suite for Interactive Heterogeneity in Single Cells},
  year         = 2024,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.17824082},
  url          = {https://doi.org/10.5281/zenodo.17824082}
}

---

## 📞 Support

- **📖 Documentation**: Check our comprehensive guides
- **🐛 Issues**: Report bugs on [GitHub Issues](https://github.com/msherafatian/masih/issues)
- **✉️ Email**: your.email@institution.edu
- **💬 Discussions**: Join our [GitHub Discussions](https://github.com/msherafatian/masih/discussions)

---

## 📈 Roadmap

### Version 1.1 (Coming Soon)
- Cell-cell communication analysis
- Spatial transcriptomics support
- Enhanced pathway databases
- Cloud deployment options

### Future Features
- Multi-sample integration tools
- Machine learning cell type prediction
- Real-time collaboration features
- API for programmatic access

---

## 🙏 Acknowledgments

MASIH is built on the shoulders of giants:

- **[Seurat](https://satijalab.org/seurat/)** - Single-cell analysis framework
- **[Slingshot](https://bioconductor.org/packages/slingshot/)** - Trajectory inference
- **[CancerSEA](http://biocc.hrbmu.edu.cn/CancerSEA/)** - Cancer functional states
- **[golem](https://github.com/ThinkR-open/golem)** - Shiny app framework
- **[Shiny](https://shiny.rstudio.com/)** - Web application framework

---

## 🔗 Related Projects

- **[MASIH-Python](https://github.com/msherafatian/masih-python)** — Python/Dash implementation with Scanpy backend and Docker support
- **[Seurat](https://satijalab.org/seurat/)** — Single-cell analysis framework (R)
- **[Scanpy](https://scanpy.readthedocs.io/)** — Single-cell analysis toolkit (Python)
- **[CancerSEA](http://biocc.hrbmu.edu.cn/CancerSEA/)** — Cancer functional state database

---


## 📜 License

This project is licensed under the GPL-3.0 license - see the [LICENSE](LICENSE) file for details.

---

**MASIH**: Making single-cell cancer analysis accessible to all researchers.

[![Made with ❤️ for Cancer Research](https://img.shields.io/badge/Made%20with%20%E2%9D%A4%EF%B8%8F%20for-Cancer%20Research-red.svg)](https://github.com/msherafatian/masih)
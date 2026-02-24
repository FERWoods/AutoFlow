# AutoFlow

## Table of Contents
- [Introduction](#introduction)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Getting Started](#getting-started)
- [Automatic Data Processing](#automatic-data-processing)
- [Channel Identification and Renaming](#channel-identification-and-renaming)
- [Seurat Processing](#seurat-processing)
- [Supervised Cell Identification](#supervised-cell-identification)
- [Model Bundles (Supervised)](#model-bundles-supervised)
- [Feature Matching & Standardization](#feature-matching--standardization)
- [Output and Results](#output-and-results)
- [Troubleshooting](#troubleshooting)
- [Support and Contact Information](#support-and-contact-information)

## Introduction

**AutoFlow** is an R Shiny application for automated analysis of flow cytometry data.  
It supports **unsupervised clustering** and **supervised classification** workflows, automating preprocessing, visualization, and summary reporting.

---

## System Requirements

- **R ≥ 4.5.0**  
- **Windows users:** Rtools for building packages from source  
- Suggested R packages:  
  `Seurat`, `flowCore`, `PeacoQC`, `mclust`, `data.table`, `plotly`, `DT`  
- **Optional for Leiden clustering:**  
  - Python (3.8+) available to `reticulate`  
  - Python package `leidenalg`  

---

## Installation

```r
# Requires devtools or remotes
install.packages("devtools")
devtools::install_github("FERWoods/AutoFlow")

# OR
install.packages("remotes")
remotes::install_github("FERWoods/AutoFlow")

# Launch the app
library(AutoFlow)
run_app()
```

## Getting Started
(With either some example data (available on Zenodo DOI: 10.5281/zenodo.18235796), or your own flow cytometry data.)
Selecting FCS Files
Click "Select FCS Files".
Choose one or more FCS files.
The application extracts metadata from FCS files, using file paths for reference.

## Automatic Data Processing (Optional)
*Experimental: this pipeline is provided for convenience and may not be robust across all panels or instruments. Users are encouraged to validate QC and gating behaviour on their own datasets.*
- Single Cell Identification: RemoveDoublets() via PeacoQC.
- Compensation: Uses spillover matrix for compensation (if not already applied).
- Logicle Transformation: Manages negative values and dynamic range.
- Quality Control with PeacoQC: Removes low-quality events.
- Debris Removal: Uses Gaussian mixture model (GMM via Mclust), with k-means fallback.
- Viability Identification: Threshold generated using GMM on Live/Dead marker, with user option to tune.
- Channel Identification and Renaming
- Renames channels by extracting descriptions from flow cytometry files.

**Selecting FCS Files**
1. Click **Select FCS Files**.
2. Choose one or more `.FCS` files.  
   AutoFlow extracts basic metadata (file name/path) and carries it forward for summaries.


If **Pre-process files = No**, AutoFlow still standardises channel names using the FCS `desc` where available (falling back to `name`), while leaving raw values intact.

## Channel Identification and Renaming

AutoFlow reads per-channel metadata (`flowCore::parameters(ff)`) and:
- Uses **`desc`** where available (fallback to `name`).
- Prefers **area channels** (`*-A`) and **compensated channels** (`Comp-*`).
- User selects which columns are used for downstream unsupervised learning.

Duplicate `desc` values are disambiguated by appending the original `name`.

## Seurat Processing
Dimensionality Reduction
Uses UMAP for visualisation.

Cluster Analysis
Employs unsupervised clustering to group cells.
**Dimensionality Reduction**  
UMAP (3D) on variable features.

Differential Expression Analysis
Identifies marker genes with a Wilcox ranked sum test.
**Clustering**  
Nearest-neighbors graph → Louvain/Leiden clustering with adjustable resolution.

**Markers**  
Differential expression (Wilcoxon rank sum) to auto-annotate clusters with “marker strings” (e.g., `CD34+CD38-`).

## Supervised Cell Identification
Uses a pretrained Random Forest model for cell type annotations.


## Model Bundles (Supervised)

1. **Install Rtools for your version of R**  
   - Visit the official Rtools page: [https://cran.r-project.org/bin/windows/Rtools/](https://cran.r-project.org/bin/windows/Rtools/)  
   - Download the version that matches your installed version of R (e.g., Rtools44 for R 4.4.x).  
   - During installation, ensure you check the box to **"Add Rtools to system PATH."**
A model bundle is a simple R list saved with `saveRDS()` that contains everything AutoFlow needs:

2. **Ensure the correct folder is in your PATH**  
   R needs access to the compiler tools inside this folder:
C:\rtools44\usr\bin
```r
bundle <- list(
  model    = <trained model object>,     # e.g., ranger::ranger(probability=TRUE)
  features = <character vector>,         # feature names expected by the model (DESC names)
  scaling  = list(means = <named numeric>, sds = <named numeric>),  # z-score params (TRAIN only)
  levels   = <levels vector>,            # factor levels order used during training
  meta     = list(
    dataset           = "BM-MPS Controls (multiclass)",
    model             = "ranger",
    splitrule         = "gini",
    mtry              = 8,
    min.node.size     = 5,
    num.trees         = 1000,
    sample.fraction   = 0.8,
    class.weights     = c(HSCs=..., ...),
    zscale            = TRUE,
    features_are_desc = TRUE,
    created           = Sys.time()
  )
)
saveRDS(bundle, "my_model_bundle.rds")
```

**Example bundling scripts** (installed with the package):
- `inst/extdata/Supervised_model_build_helper_scripts/Bmmps_model_build.R`
- `inst/extdata/Supervised_model_build_helper_scripts/Mosmann_rare_model_build.R`
- `inst/extdata/Supervised_model_build_helper_scripts/Nilsson_rare_model_build.R`


Each example:
- **Builds** a model (caret + ranger),
- **Stores** feature names, class levels, **TRAIN-only** `means`/`sds` for z-scaling, and model parameters used.
- **Saves** the bundle for direct use in AutoFlow’s Supervised mode.

## Feature Matching & Standardization

1. **Feature Matching**
   - AutoFlow first tries **exact matches** between the model’s `features` and your uploaded column names.
   - If some are missing, it attempts **synonym/alias** matches (e.g., small formatting differences).
   - If anything is still unresolved, a **feature mapping UI** lets you manually map your file’s columns to the model’s required features.  
     - Already-matched columns are pre-selected so you *only* map the ambiguous ones.

2. **Standardisation**
   - AutoFlow applies **z-scaling** using the bundle’s `scaling$means` and `scaling$sds` (computed on the training set that produced the model).
   - If the bundle sets `meta$zscale = TRUE`, scaling is always applied; otherwise raw values are used.
   - Optional **min–max** normalization can be enabled in the app if your model expects it (rare; most models just need z-scaling).

3. **Prediction**
   - Predictions are made with the bundled model (e.g., `ranger` with `probability=TRUE` for multiclass).
   - Predicted labels are written into the Seurat metadata and used in downstream counts/plots.

## Output and Results

**Visualisations**
The application provides the following visual outputs:
- **Unsupervised analysis**  
  3D UMAP embedding coloured by cluster assignment.
- **Marker quality control (QC)**  
  Density plot for a selected marker. By default, this is used to identify a viability marker for thresholding, but any marker can be selected when preprocessing is enabled (`preprocess = TRUE`)
- **Unsupervised marker distributions**  
  Marker density plots stratified by cluster ID, enabling inspection of marker expression patterns across unsupervised clusters.
- **Supervised analysis**  
  Bar plot showing predicted cell type counts for the uploaded dataset.

## Tables
The following tabular outputs are generated by the application:
- **Supervised analysis**  
  Per-file cell counts by predicted cell type, based on model predictions (downloadable as CSV).
- **Unsupervised analysis (counts)**  
  Per-file cell counts by unsupervised cluster assignment, with options to export in **long** or **wide** format.

**Files**
- **Processed `.fcs`** (optional) with:
  - Channels renamed to **DESC** (e.g., `CD34`, `CD13`), using only **`Comp-* -A`** channels for phenotyping,
  - Extra numeric channels: `label` (class code), `timepoint`, `replicate`,
  - A `label_map` keyword that maps integer `label` codes to class names.
- **Seurat object** (`.rds`) with UMAP, clusters, metadata.

## Troubleshooting

### Supervised model will not predict or requests feature mapping
- Ensure the uploaded data represent the **same biological system** and use **compatible preprocessing** to the data used to train the model.
- Use the **feature mapping UI** to align any non-matching columns.
- Confirm the model bundle contains `scaling$means` and `scaling$sds`, and that feature names in `scaling` match those in `bundle$features`.

### “Could not find tools necessary to compile a package”
- On Windows, install **Rtools** and ensure it is available on your system `PATH`.
- Refer to **System Requirements** for the correct Rtools version for your installed version of R.
- You can verify installation with:
  ```r
  pkgbuild::check_build_tools(debug = TRUE)
  ```
  
## Support and Contact Information
For issues, visit GitHub Issues. 
Please open an issue on the project’s GitHub **Issues** page with a minimal reproducible example (FCS header, a few rows of metadata, and any relevant console output).

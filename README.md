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

Welcome to AutoFlow. This application processes flow cytometry data, performs various analyses, and extracts insights. It is designed for flexibility and automation.
Welcome to AutoFlow. This application processes flow cytometry data, performs unsupervised and supervised analyses, and outputs results suitable for downstream visualization and reporting.

## System Requirements

- R version > 4.1.0
- R ≥ 4.1.0
- (Recommended) Rtools on Windows for building from source packages

## Installation

Install the package from github with:
(Requires install.packages("devtools"))
Install the package from GitHub:

```r
# Requires install.packages("devtools") or "remotes"
devtools::install_github("FERWoods/AutoFlow")
or
(Requires install.packages("remotes"))
# or
remotes::install_github("FERWoods/AutoFlow")

library(AutoFlow) 

library(AutoFlow)
run_app()
```

## Getting Started
Selecting FCS Files
Click "Select FCS Files".
Choose one or more FCS files.
Metadata Extraction
The application extracts metadata from FCS files, using file paths for reference.

## Automatic Data Processing (Optional)
- Single Cell Identification: Fits a linear model to FSC.A and FSC.H.
- Compensation: Uses "SPILLOVER" matrix for compensation.
- Logicle Transformation: Manages negative values and dynamic range.
- Quality Control with PeacoQC: Removes low-quality events.
- Debris Removal: Uses Gaussian mixture model.
- Viability Identification: Threshold on "Viability" variable (<2).
- Channel Identification and Renaming
- Renames channels by extracting descriptions from flow cytometry files.

**Selecting FCS Files**
1. Click **Select FCS Files**.
2. Choose one or more `.FCS` files.  
   AutoFlow extracts basic metadata (file name/path) and carries it forward for summaries.

**Metadata (optional)**  
Upload a `.csv`/`.xlsx` with sample annotations to enable treatment plots and grouped summaries.

## Automatic Data Processing

If **Pre-process files = Yes**, AutoFlow performs:

- **Compensation** using the file’s `SPILL`/`SPILLOVER` keyword.
- **Logicle transformation** (biexponential) on compensated channels.
- **QC with PeacoQC** (margin removal and quality filtering).
- **Debris removal** via a 2-component GMM on FSC.A.
- **Singlet clean-up** using the FSC.H ~ FSC.A relationship.
- **Viability filter** (if a “Viability” channel is present).
- **Channel selection** (see below) and renaming to human-readable labels.

If **Pre-process files = No**, AutoFlow still standardizes channel names using the FCS `desc` where available (falling back to `name`), while leaving raw values intact.

## Channel Identification and Renaming

AutoFlow reads per-channel metadata (`flowCore::parameters(ff)`) and:
- Uses **`desc`** where available (fallback to `name`).
- Prefers **area channels** (`*-A`) and **compensated channels** (`Comp-*`).
- **Excludes Live/Dead/Viability** from feature sets used for phenotyping by default.
- Preserves **EdU**/**proliferation** markers for downstream analysis.

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

## Output and Results
Visualisation includes UMAP representation and cell type counts.
Outputs include summary tables, processed .fcs files, and Seurat objects.
AutoFlow supports loading **pre-trained models** to classify cells, requiring a pre-labelled dataset to label new data from:

## Troubleshooting
## Common Issues
## FAQ
### Q: I got an error saying “Could not find tools necessary to compile a package.” What should I do?
- Upload a **Model Bundle** (`.rds`) via **Model Type → Supervised**.
- AutoFlow will:
  1. **Match features** between the uploaded data and the model’s required feature set.
  2. **Standardize** the data using the **means/SDs stored in the bundle**.
  3. **Predict** cell types (e.g., multiclass Random Forest/Ranger).
- If any model features are missing from the uploaded data, AutoFlow:
  - **Automatically maps** exact name matches and obvious synonyms (when available).
  - If anything remains unresolved, it opens an **interactive feature mapper** so you can point-and-click to map your columns to the model’s features.

This error usually means R cannot find the system tools required to build packages from source — particularly on Windows.
> Tip: Example scripts for **creating model bundles** (including **BM-MPS leave-one-day-out** training and classic binary models) are provided in `inst/examples/` (see below: [Model Bundles (Supervised)](#model-bundles-supervised)).

#### Solution: Install and Configure Rtools
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
- `vignettes/Bmmps_model_build.R`
- `vignettes/Mosmann_rare_model_build.R`
- `vignettes/Nilsson_rare_model_build.R`


Each example:
- **Builds** a model (caret + ranger),
- **Stores** feature names, class levels, and **TRAIN-only** `means`/`sds`,
- **Saves** the bundle for direct use in AutoFlow’s Supervised mode.

## Feature Matching & Standardization

If this folder is not already in your system `PATH`, add it manually:
- Open **Environment Variables** in Windows
- Under **System variables**, find `Path` → click **Edit**
- Click **New** and add:
When you upload a bundle in **Supervised** mode:

  ```
  C:\rtools44\usr\bin
  ```
1. **Feature Matching**
   - AutoFlow first tries **exact matches** between the model’s `features` and your uploaded column names.
   - If some are missing, it attempts **synonym/alias** matches (e.g., small formatting differences).
   - If anything is still unresolved, a **feature mapping UI** lets you manually map your file’s columns to the model’s required features.  
     - Already-matched columns are pre-selected so you *only* map the ambiguous ones.

2. **Standardization**
   - AutoFlow applies **z-scaling** using the bundle’s `scaling$means` and `scaling$sds` (computed on the training set that produced the model).
   - If the bundle sets `meta$zscale = TRUE`, scaling is always applied; otherwise raw values are used.
   - Optional **min–max** normalization can be enabled in the app if your model expects it (rare; most models just need z-scaling).

3. **Prediction**
   - Predictions are made with the bundled model (e.g., `ranger` with `probability=TRUE` for multiclass).
   - Predicted labels are written into the Seurat metadata and used in downstream counts/plots.

## Output and Results

3. **Restart R or RStudio**, then run the following command to verify the tools are available:
**Visualizations**
- UMAP (3D) colored by cluster or predicted label.
- Treatment plot (when metadata is provided).

**Tables**
- Per-file **cell counts** by predicted/assigned class (downloadable CSV).
- DeltaFlow-format count tables (wide format).

**Files**
- **Processed `.fcs`** (optional) with:
  - Channels renamed to **DESC** (e.g., `CD34`, `CD13`), using only **`Comp-* -A`** channels for phenotyping,
  - Extra numeric channels: `label` (class code), `timepoint`, `replicate`,
  - A `label_map` keyword that maps integer `label` codes to class names.
- **Seurat object** (`.rds`) with UMAP, clusters, and predictions.
- **Model bundle** examples (see *Model Bundles* section).

## Troubleshooting

### Supervised model won’t predict / asks for feature mapping
- Ensure your uploaded data has the **same biology** and **compatible preprocessing** as the model’s training data.
- Use the **feature mapping UI** to align any non-matching columns.
- Confirm the bundle contains `scaling$means` and `scaling$sds` and that feature names in `scaling` match `bundle$features`.

### “Could not find tools necessary to compile a package.”
- Install Rtools and ensure it’s on your PATH (Windows). See the detailed FAQ below.

## FAQ
**Q: I got an error saying “Could not find tools necessary to compile a package.”**  
A: See **System Requirements** and the Rtools checklist; verify with:
```r
pkgbuild::check_build_tools(debug = TRUE)
```

## Support and Contact Information
For issues, visit GitHub Issues. 
Please open an issue on the project’s GitHub **Issues** page with a minimal reproducible example (FCS header, a few rows of metadata, and any relevant console output).

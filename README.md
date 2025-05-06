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
- [Output and Results](#output-and-results)
- [Troubleshooting](#troubleshooting)
- [Support and Contact Information](#support-and-contact-information)

## Introduction

Welcome to AutoFlow. This application processes flow cytometry data, performs various analyses, and extracts insights. It is designed for flexibility and automation.

## System Requirements

- R version > 4.1.0

## Installation

Install the package from github with:
(Requires install.packages("devtools"))
devtools::install_github("FERWoods/AutoFlow")
or
(Requires install.packages("remotes"))
remotes::install_github("FERWoods/AutoFlow")

library(AutoFlow) 

run_app()

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

## Seurat Processing
Dimensionality Reduction
Uses UMAP for visualisation.

Cluster Analysis
Employs unsupervised clustering to group cells.

Differential Expression Analysis
Identifies marker genes with a Wilcox ranked sum test.

## Supervised Cell Identification
Uses a pretrained Random Forest model for cell type annotations.

## Output and Results
Visualisation includes UMAP representation and cell type counts.
Outputs include summary tables, processed .fcs files, and Seurat objects.

## Troubleshooting
## Common Issues
## FAQ
### Q: I got an error saying “Could not find tools necessary to compile a package.” What should I do?

This error usually means R cannot find the system tools required to build packages from source — particularly on Windows.

#### Solution: Install and Configure Rtools

1. **Install Rtools for your version of R**  
   - Visit the official Rtools page: [https://cran.r-project.org/bin/windows/Rtools/](https://cran.r-project.org/bin/windows/Rtools/)  
   - Download the version that matches your installed version of R (e.g., Rtools44 for R 4.4.x).  
   - During installation, ensure you check the box to **"Add Rtools to system PATH."**

2. **Ensure the correct folder is in your PATH**  
   R needs access to the compiler tools inside this folder:
C:\rtools44\usr\bin


If this folder is not already in your system `PATH`, add it manually:
- Open **Environment Variables** in Windows
- Under **System variables**, find `Path` → click **Edit**
- Click **New** and add:

  ```
  C:\rtools44\usr\bin
  ```

3. **Restart R or RStudio**, then run the following command to verify the tools are available:

```r
pkgbuild::check_build_tools(debug = TRUE)
```

## Support and Contact Information
For issues, visit GitHub Issues. 

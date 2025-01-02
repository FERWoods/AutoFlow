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

Install the package locally in R with:
r
devtools::install_github("AutoFlowApp") library(AutoFlowApp) run_app() ```

Getting Started
Selecting FCS Files
Click "Select FCS Files".
Choose one or more FCS files.
Metadata Extraction
The application extracts metadata from FCS files, using file paths for reference.

Automatic Data Processing (Optional)
Single Cell Identification: Fits a linear model to FSC.A and FSC.H.
Compensation: Uses "SPILLOVER" matrix for compensation.
Logicle Transformation: Manages negative values and dynamic range.
Quality Control with PeacoQC: Removes low-quality events.
Debris Removal: Uses Gaussian mixture model.
Viability Identification: Threshold on "Viability" variable (<2).
Channel Identification and Renaming
Renames channels by extracting descriptions from flow cytometry files.

Seurat Processing
Dimensionality Reduction
Uses UMAP for visualisation.

Cluster Analysis
Employs unsupervised clustering to group cells.

Differential Expression Analysis
Identifies marker genes with a Wilcox ranked sum test.

Supervised Cell Identification
Uses a pretrained Random Forest model for cell type annotations.

Output and Results
Visualisation includes UMAP representation and cell type counts.
Outputs include summary tables, processed .fcs files, and Seurat objects.
Troubleshooting
Common Issues
FAQ
Support and Contact Information
For issues, visit GitHub Issues. 

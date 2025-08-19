
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `{autoflow}`

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test
coverage](https://codecov.io/gh/FERWoods/AutoFlow/graph/badge.svg)](https://app.codecov.io/gh/FERWoods/AutoFlow)
[![R-CMD-check](https://github.com/FERWoods/AutoFlow/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/FERWoods/AutoFlow/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Installation

You can install the development version of `{autoflow}` like so:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

## Run

You can launch the application by running:

``` r
autoflow::run_app()
```

## About

You are reading the doc about version : 3.0.0

This README has been compiled on the

``` r
Sys.time()
#> [1] "2025-08-19 13:46:26 BST"
```

Here are the tests results and package coverage:

``` r
devtools::check(quiet = TRUE)
#> ℹ Loading autoflow
#> Warning: replacing previous import 'flowViz::contour' by 'graphics::contour'
#> when loading 'flowStats'
#> Warning: replacing previous import 'dplyr::filter' by 'flowCore::filter' when
#> loading 'autoflow'
#> ── R CMD check results ───────────────────────────────────── autoflow 3.0.0 ────
#> Duration: 1m 3.1s
#> 
#> ❯ checking for portable file names ... WARNING
#>   Found the following files with non-portable file names:
#>     outputs_bmmps/Samples 1_C5_C05_016_TEST.fcs
#>     outputs_bmmps/Samples 1_C5_C05_016_TRAIN.fcs
#>     outputs_bmmps/Samples 2_E10_E10_045_TEST.fcs
#>     outputs_bmmps/Samples 2_E10_E10_045_TRAIN.fcs
#>     outputs_bmmps/Samples 2_E7_E07_042_TEST.fcs
#>     outputs_bmmps/Samples 2_E7_E07_042_TRAIN.fcs
#>     outputs_bmmps/Samples 2_F9_F09_056_TEST.fcs
#>     outputs_bmmps/Samples 2_F9_F09_056_TRAIN.fcs
#>     outputs_bmmps/Samples 3_G5_G05_064_TEST.fcs
#>     outputs_bmmps/Samples 3_G5_G05_064_TRAIN.fcs
#>   These are not fully portable file names.
#>   See section ‘Package structure’ in the ‘Writing R Extensions’ manual.
#> 
#> ❯ checking whether package ‘autoflow’ can be installed ... WARNING
#>   See below...
#> 
#> ❯ checking package subdirectories ... WARNING
#>   Subdirectory ‘data’ contains no data sets.
#> 
#> ❯ checking code files for non-ASCII characters ... WARNING
#>   Found the following file with non-ASCII characters:
#>     R/app_server.R
#>   Portable packages must use only ASCII characters in their R code and
#>   NAMESPACE directives, except perhaps in comments.
#>   Use \uxxxx escapes for other characters.
#>   Function ‘tools::showNonASCIIfile’ can help in finding non-ASCII
#>   characters in files.
#> 
#> ❯ checking dependencies in R code ... WARNING
#>   'library' or 'require' calls not declared from:
#>     ‘Seurat’ ‘data.table’ ‘dplyr’ ‘readr’ ‘readxl’ ‘stringr’ ‘tidyr’
#>   'library' or 'require' calls in package code:
#>     ‘Seurat’ ‘data.table’ ‘dplyr’ ‘readr’ ‘readxl’ ‘stringr’ ‘tidyr’
#>     Please use :: or requireNamespace() instead.
#>     See section 'Suggested packages' in the 'Writing R Extensions' manual.
#>   Namespaces in Imports field not imported from:
#>     ‘biocthis’ ‘e1071’ ‘flowWorkspace’ ‘ggplot2’ ‘markdown’
#>     ‘randomForest’ ‘ranger’ ‘readr’ ‘readxl’ ‘shinycssloaders’ ‘shinyjs’
#>     ‘stringr’
#>     All declared Imports should be used.
#> 
#> ❯ checking for missing documentation entries ... WARNING
#>   Undocumented code objects:
#>     ‘run_unsupervised_func’
#>   All user-level objects in a package should have documentation entries.
#>   See chapter ‘Writing R documentation files’ in the ‘Writing R
#>   Extensions’ manual.
#> 
#> ❯ checking contents of ‘data’ directory ... WARNING
#>   Files not of a type allowed in a ‘data’ directory:
#>     ‘BM-MPS’ ‘benchmarking’
#>   Please use e.g. ‘inst/extdata’ for non-R data files
#> 
#> ❯ checking files in ‘vignettes’ ... WARNING
#>   Files in the 'vignettes' directory but no files in 'inst/doc':
#>     ‘Bmmps_model_build.R’ ‘Mosmann_rare_model_build.R’
#>     ‘Nilsson_rare_model_build.R’ ‘Supervised_model_building.Rmd’
#>   Files named as vignettes but with no recognized vignette engine:
#>      ‘vignettes/Supervised_model_building.Rmd’
#>   (Is a VignetteBuilder field missing?)
#> 
#> ❯ checking for hidden files and directories ... NOTE
#>   Found the following hidden files and directories:
#>     .idea
#>   These were most likely included in error. See section ‘Package
#>   structure’ in the ‘Writing R Extensions’ manual.
#> 
#> ❯ checking top-level files ... NOTE
#>   Non-standard files/directories found at top level:
#>     ‘PeacoQC_results’ ‘main.R’ ‘outputs’ ‘outputs_bmmps’ ‘test_file’
#> 
#> ❯ checking R code for possible problems ... NOTE
#>   app_server : auto_threshold_gmm: no visible global function definition
#>     for ‘median’
#>   app_server : auto_threshold_gmm : f1: no visible global function
#>     definition for ‘dnorm’
#>   app_server : auto_threshold_gmm : f2: no visible global function
#>     definition for ‘dnorm’
#>   app_server : auto_threshold_gmm: no visible global function definition
#>     for ‘quantile’
#>   app_server : auto_threshold_gmm: no visible global function definition
#>     for ‘uniroot’
#>   app_server: no visible global function definition for ‘setNames’
#>   app_server: no visible global function definition for ‘density’
#>   app_server: no visible binding for global variable ‘assignment’
#>   app_server: no visible global function definition for
#>     ‘auto_map_features’
#>   app_server: no visible global function definition for ‘predict’
#>   app_server: no visible binding for global variable ‘.’
#>   app_server: no visible binding for global variable ‘Sample’
#>   app_server: no visible binding for global variable ‘proliferation’
#>   app_server: no visible binding for '<<-' assignment to ‘summary_tab’
#>   app_server: no visible binding for global variable ‘summary_tab’
#>   app_server: no visible binding for global variable ‘mtcars’
#>   app_server : <anonymous>: no visible binding for global variable
#>     ‘summary_tab’
#>   app_server : <anonymous>: no visible binding for global variable
#>     ‘assignment’
#>   buffer: no visible global function definition for ‘as2’
#>   cell_labelling_bm: no visible global function definition for
#>     ‘str_detect’
#>   channel_select: no visible global function definition for ‘%like%’
#>   read_file: no visible global function definition for ‘read_excel’
#>   read_file: no visible global function definition for ‘read_csv’
#>   remove_outliers: no visible global function definition for ‘quantile’
#>   workspace_cell_labels: no visible binding for global variable
#>     ‘cell_types_fixed’
#>   workspace_cell_labels: no visible global function definition for
#>     ‘fj_ws_get_sample_groups’
#>   workspace_cell_labels: no visible global function definition for
#>     ‘open_flowjo_xml’
#>   workspace_cell_labels: no visible global function definition for
#>     ‘GetFlowJoLabels’
#>   Undefined global functions or variables:
#>     %like% . GetFlowJoLabels Sample as2 assignment auto_map_features
#>     cell_types_fixed density dnorm fj_ws_get_sample_groups median mtcars
#>     open_flowjo_xml predict proliferation quantile read_csv read_excel
#>     setNames str_detect summary_tab uniroot
#>   Consider adding
#>     importFrom("datasets", "mtcars")
#>     importFrom("stats", "density", "dnorm", "median", "predict",
#>                "quantile", "setNames", "uniroot")
#>   to your NAMESPACE file.
#> 
#> ❯ checking package vignettes ... NOTE
#>   Package has ‘vignettes’ subdirectory but apparently no vignettes.
#>   Perhaps the ‘VignetteBuilder’ information is missing from the
#>   DESCRIPTION file?
#> 
#> 0 errors ✔ | 8 warnings ✖ | 4 notes ✖
#> Error: R CMD check found WARNINGs
```

``` r
covr::package_coverage()
#> Error in loadNamespace(x): there is no package called 'covr'
```

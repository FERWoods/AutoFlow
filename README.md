
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `{autoflow}`

<!-- badges: start -->
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

You are reading the doc about version : 2.0.0

This README has been compiled on the

``` r
Sys.time()
#> [1] "2025-01-27 00:10:00 GMT"
```

Here are the tests results and package coverage:

``` r
devtools::check(quiet = TRUE)
#> ℹ Loading autoflow
#> Warning: replacing previous import 'dplyr::filter' by 'flowCore::filter' when
#> loading 'autoflow'
#> ── R CMD check results ───────────────────────────────────── autoflow 2.0.0 ────
#> Duration: 2.4s
#> 
#> ❯ checking package dependencies ... ERROR
#>   Namespace dependencies missing from DESCRIPTION Imports/Depends entries:
#>     'Seurat', 'dplyr', 'flowCore'
#>   
#>   See section 'The DESCRIPTION file' in the 'Writing R Extensions'
#>   manual.
#> 
#> 1 error ✖ | 0 warnings ✔ | 0 notes ✔
#> Error: R CMD check found ERRORs
```

``` r
covr::package_coverage()
#> Error in loadNamespace(x): there is no package called 'covr'
```


<!-- README.md is generated from README.Rmd. Please edit that file -->

# `{autoflow}`

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
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
#> [1] "2025-01-27 00:50:54 GMT"
```

Here are the tests results and package coverage:

``` r
devtools::check(quiet = TRUE)
#> ℹ Loading autoflow
#> Warning: replacing previous import 'flowViz::contour' by 'graphics::contour'
#> when loading 'flowStats'
#> Warning: replacing previous import 'dplyr::filter' by 'flowCore::filter' when
#> loading 'autoflow'
#> ── R CMD check results ───────────────────────────────────── autoflow 2.0.0 ────
#> Duration: 3m 0.4s
#> 
#> ❯ checking tests ...
#>   See below...
#> 
#> ❯ checking whether package 'autoflow' can be installed ... [28s] WARNING
#>   See below...
#> 
#> ❯ checking dependencies in R code ... WARNING
#>   'library' or 'require' calls not declared from:
#>     'Seurat' 'data.table' 'dplyr' 'readr' 'readxl' 'stringr' 'tidyr'
#>   'library' or 'require' calls in package code:
#>     'Seurat' 'data.table' 'dplyr' 'readr' 'readxl' 'stringr' 'tidyr'
#>     Please use :: or requireNamespace() instead.
#>     See section 'Suggested packages' in the 'Writing R Extensions' manual.
#>   Namespaces in Imports field not imported from:
#>     'data.table' 'pkgload' 'readr' 'readxl' 'stringr' 'tidyr'
#>     All declared Imports should be used.
#> 
#> ❯ checking for missing documentation entries ... WARNING
#>   Undocumented code objects:
#>     'app_server' 'run_unsupervised_func'
#>   All user-level objects in a package should have documentation entries.
#>   See chapter 'Writing R documentation files' in the 'Writing R
#>   Extensions' manual.
#> 
#> ❯ checking package dependencies ... NOTE
#>   Imports includes 21 non-default packages.
#>   Importing from so many packages makes the package vulnerable to any of
#>   them becoming unavailable.  Move as many as possible to Suggests and
#>   use conditionally.
#> 
#> ❯ checking for future file timestamps ... NOTE
#>   unable to verify current time
#> 
#> ❯ checking R code for possible problems ... [32s] NOTE
#>   app_server : <anonymous>: no visible global function definition for
#>     'lm'
#>   app_server : <anonymous>: no visible global function definition for
#>     'coef'
#>   app_server: no visible binding for global variable 'processFCS'
#>   app_server: no visible binding for global variable 'scaling_parameters'
#>   app_server: no visible global function definition for 'predict'
#>   app_server: no visible binding for '<<-' assignment to 'seurat_to_df'
#>   app_server: no visible binding for global variable 'seurat_to_df'
#>   app_server: no visible binding for global variable '.'
#>   app_server: no visible binding for global variable 'Sample'
#>   app_server: no visible binding for global variable 'assignment'
#>   app_server: no visible binding for global variable 'proliferation'
#>   app_server: no visible global function definition for 'setNames'
#>   app_server: no visible binding for '<<-' assignment to 'summary_tab'
#>   app_server: no visible binding for global variable 'nCount_RNA'
#>   app_server: no visible binding for global variable 'summary_tab'
#>   app_server : <anonymous>: no visible binding for global variable
#>     'summary_tab'
#>   app_server : <anonymous>: no visible global function definition for
#>     'write.csv'
#>   app_server : <anonymous>: no visible global function definition for
#>     'data.table'
#>   app_server : <anonymous>: no visible binding for global variable
#>     'assignment'
#>   app_server : <anonymous>: no visible global function definition for
#>     'pivot_wider'
#>   app_server: no visible global function definition for 'plot_ly'
#>   app_server: no visible binding for global variable 'mtcars'
#>   buffer: no visible global function definition for 'as2'
#>   cell_labelling_bm: no visible global function definition for
#>     'str_detect'
#>   channel_select: no visible global function definition for '%like%'
#>   read_file: no visible global function definition for 'read_excel'
#>   read_file: no visible global function definition for 'read_csv'
#>   remove_outliers: no visible global function definition for 'quantile'
#>   workspace_cell_labels: no visible binding for global variable
#>     'cell_types_fixed'
#>   workspace_cell_labels: no visible global function definition for
#>     'fj_ws_get_sample_groups'
#>   workspace_cell_labels: no visible global function definition for
#>     'open_flowjo_xml'
#>   workspace_cell_labels: no visible global function definition for
#>     'GetFlowJoLabels'
#>   Undefined global functions or variables:
#>     %like% . GetFlowJoLabels Sample as2 assignment cell_types_fixed coef
#>     data.table fj_ws_get_sample_groups lm mtcars nCount_RNA
#>     open_flowjo_xml pivot_wider plot_ly predict processFCS proliferation
#>     quantile read_csv read_excel scaling_parameters setNames seurat_to_df
#>     str_detect summary_tab write.csv
#>   Consider adding
#>     importFrom("datasets", "mtcars")
#>     importFrom("stats", "coef", "lm", "predict", "quantile", "setNames")
#>     importFrom("utils", "write.csv")
#>   to your NAMESPACE file.
#> 
#> ── Test failures ───────────────────────────────────────────────── testthat ────
#> 
#> > # This file is part of the standard setup for testthat.
#> > # It is recommended that you do not modify it.
#> > #
#> > # Where should you do additional test configuration?
#> > # Learn more about the roles of various files in:
#> > # * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
#> > # * https://testthat.r-lib.org/articles/special-files.html
#> > 
#> > library(testthat)
#> Warning message:
#> package 'testthat' was built under R version 4.3.3 
#> > library(autoflow)
#> Warning message:
#> replacing previous import 'dplyr::filter' by 'flowCore::filter' when loading 'autoflow' 
#> > 
#> > test_check("autoflow")
#> Loading required package: shiny
#> [ FAIL 2 | WARN 2 | SKIP 0 | PASS 86 ]
#> 
#> ══ Warnings ════════════════════════════════════════════════════════════════════
#> ── Warning ('test-golem-recommended.R:55:1'): (code run outside of `test_that()`) ──
#> package 'shiny' was built under R version 4.3.3
#> Backtrace:
#>     ▆
#>  1. └─shiny::testServer(...) at test-golem-recommended.R:55:1
#>  2.   └─base::require(shiny)
#>  3.     ├─base::tryCatch(...)
#>  4.     │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#>  5.     │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#>  6.     │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#>  7.     └─base::library(...)
#>  8.       └─base (local) testRversion(pkgInfo, package, pkgpath)
#> ── Warning ('test-golem-recommended.R:55:1'): (code run outside of `test_that()`) ──
#> cannot open file 'R/helper_functions.R': No such file or directory
#> Backtrace:
#>      ▆
#>   1. └─shiny::testServer(...) at test-golem-recommended.R:55:1
#>   2.   ├─shiny:::withMockContext(...)
#>   3.   │ ├─shiny::isolate(...)
#>   4.   │ │ ├─shiny::..stacktraceoff..(...)
#>   5.   │ │ └─ctx$run(...)
#>   6.   │ │   ├─promises::with_promise_domain(...)
#>   7.   │ │   │ └─domain$wrapSync(expr)
#>   8.   │ │   ├─shiny::withReactiveDomain(...)
#>   9.   │ │   │ └─promises::with_promise_domain(...)
#>  10.   │ │   │   └─domain$wrapSync(expr)
#>  11.   │ │   │     └─base::force(expr)
#>  12.   │ │   ├─shiny::captureStackTraces(...)
#>  13.   │ │   │ └─promises::with_promise_domain(...)
#>  14.   │ │   │   └─domain$wrapSync(expr)
#>  15.   │ │   │     └─base::withCallingHandlers(expr, error = doCaptureStack)
#>  16.   │ │   └─env$runWith(self, func)
#>  17.   │ │     └─shiny (local) contextFunc()
#>  18.   │ │       └─shiny::..stacktraceon..(expr)
#>  19.   │ ├─shiny::withReactiveDomain(...)
#>  20.   │ │ └─promises::with_promise_domain(...)
#>  21.   │ │   └─domain$wrapSync(expr)
#>  22.   │ │     └─base::force(expr)
#>  23.   │ └─withr::with_options(...)
#>  24.   │   └─base::force(code)
#>  25.   └─autoflow (local) server(input = session$input, output = session$output, session = session)
#>  26.     └─base::source("R/helper_functions.R")
#>  27.       └─base::file(filename, "r", encoding = encoding)
#> 
#> ══ Failed tests ════════════════════════════════════════════════════════════════
#> ── Error ('test-golem-recommended.R:55:1'): (code run outside of `test_that()`) ──
#> Error in `file(filename, "r", encoding = encoding)`: cannot open the connection
#> Backtrace:
#>      ▆
#>   1. ├─shiny::testServer(...) at test-golem-recommended.R:55:1
#>   2. │ ├─shiny:::withMockContext(...)
#>   3. │ │ ├─shiny::isolate(...)
#>   4. │ │ │ ├─shiny::..stacktraceoff..(...)
#>   5. │ │ │ └─ctx$run(...)
#>   6. │ │ │   ├─promises::with_promise_domain(...)
#>   7. │ │ │   │ └─domain$wrapSync(expr)
#>   8. │ │ │   ├─shiny::withReactiveDomain(...)
#>   9. │ │ │   │ └─promises::with_promise_domain(...)
#>  10. │ │ │   │   └─domain$wrapSync(expr)
#>  11. │ │ │   │     └─base::force(expr)
#>  12. │ │ │   ├─shiny::captureStackTraces(...)
#>  13. │ │ │   │ └─promises::with_promise_domain(...)
#>  14. │ │ │   │   └─domain$wrapSync(expr)
#>  15. │ │ │   │     └─base::withCallingHandlers(expr, error = doCaptureStack)
#>  16. │ │ │   └─env$runWith(self, func)
#>  17. │ │ │     └─shiny (local) contextFunc()
#>  18. │ │ │       └─shiny::..stacktraceon..(expr)
#>  19. │ │ ├─shiny::withReactiveDomain(...)
#>  20. │ │ │ └─promises::with_promise_domain(...)
#>  21. │ │ │   └─domain$wrapSync(expr)
#>  22. │ │ │     └─base::force(expr)
#>  23. │ │ └─withr::with_options(...)
#>  24. │ │   └─base::force(code)
#>  25. │ └─autoflow (local) server(input = session$input, output = session$output, session = session)
#>  26. │   └─base::source("R/helper_functions.R")
#>  27. │     └─base::file(filename, "r", encoding = encoding)
#>  28. └─base::.handleSimpleError(...)
#>  29.   └─shiny (local) h(simpleError(msg, call))
#> ── Error ('test-mod_helper_functions.R:1:1'): (code run outside of `test_that()`) ──
#> Error in `eval(code, test_env)`: object 'mod_helper_functions_server' not found
#> Backtrace:
#>     ▆
#>  1. └─shiny::testServer(...) at test-mod_helper_functions.R:1:1
#>  2.   └─shiny:::isModuleServer(app)
#> 
#> [ FAIL 2 | WARN 2 | SKIP 0 | PASS 86 ]
#> Error: Test failures
#> Execution halted
#> 
#> 1 error ✖ | 3 warnings ✖ | 3 notes ✖
#> Error: R CMD check found ERRORs
```

``` r
covr::package_coverage()
#> Error in loadNamespace(x): there is no package called 'covr'
```

testthat::test_that("preprocess_flowframes returns stable empty structures on NULL/empty input", {
  out1 <- preprocess_flowframes(NULL)
  testthat::expect_true(is.list(out1))
  testthat::expect_true(all(c("frames","flags","counts") %in% names(out1)))
  testthat::expect_equal(out1$frames, list())
  testthat::expect_equal(out1$flags,  list())
  testthat::expect_s3_class(out1$counts, "data.frame")
  testthat::expect_equal(nrow(out1$counts), 0)

  out2 <- preprocess_flowframes(list())
  testthat::expect_equal(out2$frames, list())
  testthat::expect_equal(out2$flags,  list())
  testthat::expect_equal(nrow(out2$counts), 0)
})

testthat::test_that("preprocess_flowframes preserves list length; non-flowFrames become NULL slots", {
  out <- preprocess_flowframes(list("nope", 123))
  testthat::expect_length(out$frames, 2)
  testthat::expect_length(out$flags,  2)
  testthat::expect_true(is.null(out$frames[[1]]) && is.null(out$flags[[1]]))
  testthat::expect_true(is.null(out$frames[[2]]) && is.null(out$flags[[2]]))
  testthat::expect_equal(nrow(out$counts), 2)
})

testthat::test_that("preprocess_flowframes returns transformed frames + per-event flags aligned to rows", {
  testthat::skip_if_not_installed("flowCore")
  testthat::skip_if_not_installed("Biobase")
  testthat::skip_if_not_installed("mclust")   # debris GMM

  make_ff <- function(expr, desc = NULL) {
    stopifnot(is.matrix(expr), !is.null(colnames(expr)))
    cn <- colnames(expr)
    if (is.null(desc)) desc <- cn

    rng  <- apply(expr, 2, function(v) range(v, na.rm = TRUE))
    minR <- as.numeric(rng[1, ])
    maxR <- as.numeric(rng[2, ])
    span <- pmax(1, ceiling(maxR - minR))

    params_df <- data.frame(
      name     = cn,
      desc     = as.character(desc),
      range    = span,
      minRange = minR,
      maxRange = maxR,
      stringsAsFactors = FALSE
    )
    params <- Biobase::AnnotatedDataFrame(params_df)
    flowCore::flowFrame(exprs = expr, parameters = params)
  }

  set.seed(1)
  n_low  <- 80
  n_high <- 220
  n <- n_low + n_high

  # Make sure we have >=50 finite for debris detection.
  expr <- cbind(
    `FSC.A` = c(rnorm(n_low,  50,   10),  rnorm(n_high, 5000, 500)),
    `FSC.H` = c(rnorm(n_low,  40,   10),  rnorm(n_high, 4500, 500)),
    `SSC.A` = rnorm(n, 2000, 200),
    `CD3`   = rnorm(n, 2, 0.5),
    `Time`  = seq_len(n)
  )
  ff <- make_ff(expr)

  out <- preprocess_flowframes(list(ff), min_events_peacoqc = 999999) # skip badqc branch deterministically

  testthat::expect_length(out$frames, 1)
  testthat::expect_length(out$flags,  1)
  testthat::expect_equal(nrow(out$counts), 1)

  ff2 <- out$frames[[1]]
  fl  <- out$flags[[1]]

  testthat::expect_s4_class(ff2, "flowFrame")
  testthat::expect_s3_class(fl, "data.frame")
  testthat::expect_true(all(c("AutoFlow_debris","AutoFlow_doublet","AutoFlow_badqc") %in% colnames(fl)))

  # alignment: flags rows = events
  testthat::expect_equal(nrow(fl), nrow(flowCore::exprs(ff2)))

  # flags are 0/1 integers
  for (nm in colnames(fl)) {
    testthat::expect_type(fl[[nm]], "integer")
    testthat::expect_true(all(fl[[nm]] %in% c(0L, 1L)))
  }

  # counts consistent with flags
  testthat::expect_equal(out$counts$n_before[1], nrow(flowCore::exprs(ff)))
  testthat::expect_equal(out$counts$n_after_transform[1], nrow(flowCore::exprs(ff2)))
  testthat::expect_equal(out$counts$n_debris[1],  sum(fl$AutoFlow_debris  == 1L))
  testthat::expect_equal(out$counts$n_doublet[1], sum(fl$AutoFlow_doublet == 1L))
  testthat::expect_equal(out$counts$n_badqc[1],   sum(fl$AutoFlow_badqc   == 1L))

  # in this synthetic bimodal FSC, we should flag *some* debris but not all
  prop <- mean(fl$AutoFlow_debris == 1L)
  testthat::expect_gt(prop, 0.01)
  testthat::expect_lt(prop, 0.95)

  # we forced badqc skipping -> should be all zeros
  testthat::expect_true(all(fl$AutoFlow_badqc == 0L))
})

testthat::test_that("preprocess_flowframes handles a mix of valid/invalid frames without crashing", {
  testthat::skip_if_not_installed("flowCore")
  testthat::skip_if_not_installed("Biobase")
  testthat::skip_if_not_installed("mclust")

  make_ff <- function(expr) {
    cn <- colnames(expr)
    rng  <- apply(expr, 2, function(v) range(v, na.rm = TRUE))
    minR <- as.numeric(rng[1, ])
    maxR <- as.numeric(rng[2, ])
    span <- pmax(1, ceiling(maxR - minR))
    params <- Biobase::AnnotatedDataFrame(
      data.frame(
        name     = cn,
        desc     = cn,
        range    = span,
        minRange = minR,
        maxRange = maxR,
        stringsAsFactors = FALSE
      )
    )
    flowCore::flowFrame(exprs = expr, parameters = params)
  }

  set.seed(2)
  n <- 120
  ff_ok <- make_ff(cbind(`FSC.A` = rnorm(n, 5000, 500),
                         `FSC.H` = rnorm(n, 4500, 500)))

  out <- preprocess_flowframes(list(ff_ok, "bad", NULL), min_events_peacoqc = 999999)

  testthat::expect_length(out$frames, 3)
  testthat::expect_length(out$flags,  3)
  testthat::expect_equal(nrow(out$counts), 3)

  testthat::expect_s4_class(out$frames[[1]], "flowFrame")
  testthat::expect_s3_class(out$flags[[1]], "data.frame")

  testthat::expect_true(is.null(out$frames[[2]]) && is.null(out$flags[[2]]))
  testthat::expect_true(is.null(out$frames[[3]]) && is.null(out$flags[[3]]))
})

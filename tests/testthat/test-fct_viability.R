test_that("get_marker_vector returns numeric(0) on NULL/empty inputs", {
  expect_equal(get_marker_vector(NULL, "CD3"), numeric(0))
  expect_equal(get_marker_vector(list(), "CD3"), numeric(0))
  expect_equal(get_marker_vector(list("not a flowFrame"), NULL), numeric(0))
})

test_that("get_marker_vector concatenates marker values across frames and skips missing channels", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("Biobase")

  make_ff <- function(expr, desc = NULL) {
    stopifnot(is.matrix(expr), !is.null(colnames(expr)))
    cn <- colnames(expr)
    if (is.null(desc)) desc <- cn

    rng <- apply(expr, 2, function(v) range(v, na.rm = TRUE))
    minR <- as.numeric(rng[1, ])
    maxR <- as.numeric(rng[2, ])
    span <- pmax(1, ceiling(maxR - minR))

    params <- Biobase::AnnotatedDataFrame(
      data.frame(
        name     = cn,
        desc     = as.character(desc),
        range    = span,
        minRange = minR,
        maxRange = maxR,
        stringsAsFactors = FALSE
      )
    )
    flowCore::flowFrame(exprs = expr, parameters = params)
  }

  ff1 <- make_ff(cbind(CD3 = c(1, 2, 3), FSC.A = c(10, 11, 12)))
  ff2 <- make_ff(cbind(CD3 = c(4, 5),     SSC.A = c(20, 21)))
  ff3 <- make_ff(cbind(CD4 = c(9, 9, 9)))  # missing CD3

  v <- get_marker_vector(list(ff1, ff2, ff3), "CD3")

  expect_type(v, "double")
  expect_equal(v, c(1, 2, 3, 4, 5))
})

test_that("norm_minmax scales to [0,1] and handles NAs", {
  x <- c(2, 4, 6, NA)
  y <- norm_minmax(x)

  expect_length(y, length(x))
  expect_true(all(is.na(y) == is.na(x)))
  expect_equal(y[1:3], c(0, 0.5, 1))
})

test_that("norm_minmax returns zeros when range is zero or non-finite", {
  expect_equal(norm_minmax(c(5, 5, 5)), c(0, 0, 0))
  expect_equal(norm_minmax(c(NA, NA)), c(0, 0))
  expect_equal(norm_minmax(c(Inf, 1, 2)), c(0, 0, 0)) # range non-finite -> zeros
})

test_that("auto_threshold_gmm returns NA for <200 finite values", {
  skip_if_not_installed("mclust")
  set.seed(1)
  x <- rnorm(199)
  expect_true(is.na(auto_threshold_gmm(x)))
})

test_that("auto_threshold_gmm returns a finite threshold for clear bimodal data", {
  skip_if_not_installed("mclust")
  set.seed(1)

  # Two well-separated modes
  x <- c(rnorm(500, mean = 0, sd = 0.5),
         rnorm(500, mean = 5, sd = 0.5))

  thr <- auto_threshold_gmm(x)

  expect_type(thr, "double")
  expect_true(is.finite(thr))

  # Should lie between the two component means in a clear bimodal setup
  expect_gt(thr, 1)
  expect_lt(thr, 4)
})

test_that("auto_threshold_gmm ignores non-finite values", {
  skip_if_not_installed("mclust")
  set.seed(2)

  x <- c(rnorm(500, 0, 0.5),
         rnorm(500, 5, 0.5),
         NA, Inf, -Inf)

  thr <- auto_threshold_gmm(x)

  expect_true(is.finite(thr))
})


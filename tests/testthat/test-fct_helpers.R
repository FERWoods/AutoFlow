test_that("%||% returns left-hand side when not NULL", {
  expect_equal(1 %||% 2, 1)
  expect_equal("a" %||% "b", "a")
  expect_equal(FALSE %||% TRUE, FALSE)
})

test_that("%||% returns right-hand side when left is NULL", {
  expect_equal(NULL %||% 2, 2)
  expect_equal(NULL %||% "b", "b")
})

test_that("canon lowercases and strips non-alphanumeric", {
  expect_equal(canon("Live/Dead-A"), "livedeada")
  expect_equal(canon("  CD34 (PE) "), "cd34pe")
  expect_equal(canon(c("A_B", "A-B")), c("ab", "ab"))
})

test_that("dens_or_hist adds a density-like trace when x has variance", {
  skip_if_not_installed("plotly")

  p0 <- plotly::plot_ly()
  x  <- c(1, 2, 3, 4, 5)

  p1 <- dens_or_hist(p0, x, "mytrace")
  b1 <- plotly::plotly_build(p1)

  tr <- b1$x$data %||% list()
  expect_true(length(tr) >= 1)

  has_density_like <- any(vapply(tr, function(tt) {
    is.list(tt) &&
      # density line should have numeric x/y vectors of length > 1
      is.numeric(tt$x) && length(tt$x) > 1 &&
      is.numeric(tt$y) && length(tt$y) > 1 &&
      # and it's typically a scatter/lines trace (but don't over-specify)
      (is.null(tt$type) || tt$type %in% c("scatter", "scattergl")) &&
      (!is.null(tt$mode) && grepl("lines", tt$mode))
  }, logical(1)))

  expect_true(has_density_like)
})

test_that("dens_or_hist drops non-finite values before deciding", {
  skip_if_not_installed("plotly")

  p0 <- plotly::plot_ly()
  x  <- c(1, NA, Inf, 2, NaN)

  p1 <- dens_or_hist(p0, x, "finite_only")
  b1 <- plotly::plotly_build(p1)

  tr <- b1$x$data %||% list()
  expect_true(length(tr) >= 1)

  has_density_like <- any(vapply(tr, function(tt) {
    is.list(tt) &&
      is.numeric(tt$x) && length(tt$x) > 1 &&
      is.numeric(tt$y) && length(tt$y) > 1 &&
      (is.null(tt$type) || tt$type %in% c("scatter", "scattergl")) &&
      (!is.null(tt$mode) && grepl("lines", tt$mode))
  }, logical(1)))

  expect_true(has_density_like)
})

test_that("dens_or_hist adds a 1-bin histogram when exactly one finite value remains", {
  skip_if_not_installed("plotly")

  p0 <- plotly::plot_ly()
  x  <- c(NA, Inf, 7)

  p1 <- dens_or_hist(p0, x, "single")
  b1 <- plotly::plotly_build(p1)

  tr <- b1$x$data %||% list()
  expect_true(length(tr) >= 1)

  has_hist_1bin <- any(vapply(tr, function(tt) {
    is.list(tt) &&
      # plotly histogram traces usually have type == "histogram"
      (
        identical(tt$type, "histogram") ||
          # but sometimes nbinsx is the giveaway
          (!is.null(tt$nbinsx) && as.integer(tt$nbinsx) == 1L)
      )
  }, logical(1)))

  expect_true(has_hist_1bin)
})

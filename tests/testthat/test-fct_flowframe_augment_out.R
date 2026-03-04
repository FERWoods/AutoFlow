test_that("ff_add_cols_safe produces a structurally valid flowFrame", {

  skip_if_not_installed("flowCore")
  skip_if_not_installed("Biobase")

  n <- 100

  expr <- cbind(
    FSC.A = rnorm(n),
    SSC.A = rnorm(n),
    FL1.A = rnorm(n)
  )

  pd <- data.frame(
    name     = colnames(expr),
    desc     = colnames(expr),
    range    = rep(1024, ncol(expr)),
    minRange = rep(0,    ncol(expr)),
    maxRange = rep(1023, ncol(expr)),
    stringsAsFactors = FALSE
  )
  rownames(pd) <- paste0("$P", seq_len(nrow(pd)))

  params <- Biobase::AnnotatedDataFrame(pd)

  ff <- flowCore::flowFrame(
    exprs       = expr,
    parameters  = params,
    description = list(
      "$PAR" = as.character(ncol(expr)),
      "$TOT" = as.character(nrow(expr))
    )
  )

  add <- data.frame(
    AutoFlow_debris    = sample(0:1, n, TRUE),
    AutoFlow_doublet   = sample(0:1, n, TRUE),
    AutoFlow_badqc     = sample(0:1, n, TRUE),
    AutoFlow_viable    = sample(0:1, n, TRUE),
    AutoFlow_singlecell= sample(0:1, n, TRUE),
    check.names = FALSE
  )

  ff2 <- autoflow:::ff_add_cols_safe(ff, add)

  expr2 <- flowCore::exprs(ff2)
  pd2   <- flowCore::pData(flowCore::parameters(ff2))
  kw2   <- flowCore::keyword(ff2)

  # Core invariants that prevent "subscript out of bounds"
  expect_equal(ncol(expr2), nrow(pd2))
  expect_identical(colnames(expr2), as.character(pd2$name))

  # Keyword $PAR must match number of params/channels
  expect_equal(as.integer(kw2[["$PAR"]]), ncol(expr2))
  expect_equal(nrow(expr2), n)
})

test_that("augmented flowFrame writes and reads correctly", {

  skip_if_not_installed("flowCore")
  skip_if_not_installed("Biobase")

  n <- 200
  p <- 5

  expr <- matrix(rnorm(n * p), ncol = p)
  colnames(expr) <- paste0("Ch", seq_len(p))

  pd <- data.frame(
    name     = colnames(expr),
    desc     = colnames(expr),
    range    = rep(1024, p),
    minRange = rep(0,    p),
    maxRange = rep(1023, p),
    stringsAsFactors = FALSE
  )
  rownames(pd) <- paste0("$P", seq_len(p))

  ff <- flowCore::flowFrame(
    exprs       = expr,
    parameters  = Biobase::AnnotatedDataFrame(pd),
    description = list("$PAR" = as.character(p), "$TOT" = as.character(n))
  )

  add <- data.frame(
    AutoFlow_debris  = sample(0:1, n, TRUE),
    AutoFlow_doublet = sample(0:1, n, TRUE),
    check.names = FALSE
  )

  ff2 <- autoflow:::ff_add_cols_safe(ff, add)

  tf <- tempfile(fileext = ".fcs")

  expect_error(flowCore::write.FCS(ff2, tf), NA)

  ff3 <- flowCore::read.FCS(tf, transformation = FALSE)

  expect_equal(nrow(flowCore::exprs(ff3)), n)
  expect_equal(ncol(flowCore::exprs(ff3)), ncol(flowCore::exprs(ff2)))
})

test_that("ff_add_cols_safe survives stale $Pn keywords (real-world FCS)", {

  skip_if_not_installed("flowCore")
  skip_if_not_installed("Biobase")

  n <- 50
  expr <- cbind(A = rnorm(n), B = rnorm(n))
  p <- ncol(expr)

  pd <- data.frame(
    name     = colnames(expr),
    desc     = colnames(expr),
    range    = rep(1024, p),
    minRange = rep(0,    p),
    maxRange = rep(1023, p),
    stringsAsFactors = FALSE
  )
  rownames(pd) <- paste0("$P", seq_len(p))

  # Deliberately include stale $P* keys as if from an old FCS header
  dirty_kw <- list(
    "$PAR" = as.character(p),
    "$TOT" = as.character(n),
    "$P1N" = "A",
    "$P2N" = "B",
    "$P9N" = "NONEXISTENT"  # stale / wrong
  )

  ff <- flowCore::flowFrame(
    exprs       = expr,
    parameters  = Biobase::AnnotatedDataFrame(pd),
    description = dirty_kw
  )

  add <- data.frame(Flag1 = sample(0:1, n, TRUE), check.names = FALSE)

  ff2 <- autoflow:::ff_add_cols_safe(ff, add)

  tf <- tempfile(fileext = ".fcs")
  expect_error(flowCore::write.FCS(ff2, tf), NA)
})

test_that("param_table returns name/desc and replaces NA desc with empty string", {
  make_ff <- function(names, descs = names) {
    testthat::skip_if_not_installed("flowCore")
    testthat::skip_if_not_installed("Biobase")

    n <- length(names)
    stopifnot(length(descs) == n)

    expr <- matrix(runif(n * 10), nrow = 10, ncol = n)
    colnames(expr) <- names

    pd <- data.frame(
      name     = names,
      desc     = descs,
      range    = rep(262144, n),
      minRange = rep(0, n),
      maxRange = rep(262144, n),
      stringsAsFactors = FALSE
    )

    params <- Biobase::AnnotatedDataFrame(pd)
    flowCore::flowFrame(exprs = expr, parameters = params)
  }
  ff <- make_ff(c("FSC-A", "SSC-A"), c(NA, "SSC-A"))
  tb <- param_table(ff)

  expect_s3_class(tb, "data.frame")
  expect_equal(tb$name, c("FSC-A", "SSC-A"))
  expect_equal(tb$desc[1], "")
  expect_equal(tb$desc[2], "SSC-A")
})

test_that("build_name_label_maps uses union of names and majority (case-insensitive) desc", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("Biobase")

  make_ff <- function(names, descs = names) {
    testthat::skip_if_not_installed("flowCore")
    testthat::skip_if_not_installed("Biobase")

    n <- length(names)
    stopifnot(length(descs) == n)

    expr <- matrix(runif(n * 10), nrow = 10, ncol = n)
    colnames(expr) <- names

    pd <- data.frame(
      name     = names,
      desc     = descs,
      range    = rep(262144, n),
      minRange = rep(0, n),
      maxRange = rep(262144, n),
      stringsAsFactors = FALSE
    )

    params <- Biobase::AnnotatedDataFrame(pd)
    flowCore::flowFrame(exprs = expr, parameters = params)
  }

  # Same NAME appears with desc variants in different frames
  ff1 <- make_ff(
    names = c("CD3", "CD4"),
    descs = c("CD3", "CD4")
  )
  ff2 <- make_ff(
    names = c("CD3", "CD8"),
    descs = c("cd3", NA)          # lower-case variant; CD8 desc missing
  )
  ff3 <- make_ff(
    names = c("CD3", "CD8"),
    descs = c("CD3", "CD8 marker") # CD3 again
  )

  maps <- build_name_label_maps(list(ff1, ff2, ff3))

  # union includes CD3, CD4, CD8
  expect_true(all(c("CD3", "CD4", "CD8") %in% maps$union_names))

  # CD3: majority label should resolve to "CD3" (case-insensitive)
  expect_equal(maps$name_to_label[["CD3"]], "CD3")

  # CD8: desc missing in one frame but present in another -> should use present desc
  expect_equal(maps$name_to_label[["CD8"]], "CD8 marker")

  # label_to_name should be inverse mapping for the chosen labels
  expect_equal(maps$label_to_name[["CD3"]], "CD3")
  expect_equal(maps$label_to_name[["CD8 marker"]], "CD8")
})

test_that("build_name_label_maps disambiguates duplicate display labels by appending [NAME]", {
  make_ff <- function(names, descs = names) {
    testthat::skip_if_not_installed("flowCore")
    testthat::skip_if_not_installed("Biobase")

    n <- length(names)
    stopifnot(length(descs) == n)

    expr <- matrix(runif(n * 10), nrow = 10, ncol = n)
    colnames(expr) <- names

    pd <- data.frame(
      name     = names,
      desc     = descs,
      range    = rep(262144, n),
      minRange = rep(0, n),
      maxRange = rep(262144, n),
      stringsAsFactors = FALSE
    )

    params <- Biobase::AnnotatedDataFrame(pd)
    flowCore::flowFrame(exprs = expr, parameters = params)
  }

  # Two different channels share the same DESC
  ff <- make_ff(
    names = c("FL1.A", "FL2.A"),
    descs = c("CD3", "CD3")  # duplicate display label
  )

  maps <- build_name_label_maps(list(ff))

  # both should be disambiguated
  expect_equal(maps$name_to_label[["FL1.A"]], "CD3 [FL1.A]")
  expect_equal(maps$name_to_label[["FL2.A"]], "CD3 [FL2.A]")

  # and label_to_name should map back
  expect_equal(maps$label_to_name[["CD3 [FL1.A]"]], "FL1.A")
  expect_equal(maps$label_to_name[["CD3 [FL2.A]"]], "FL2.A")
})

test_that("pick_fsc finds FSC-A and FSC-H across common naming variants", {
  cn <- c("Time", "FSC.A", "FSC.H", "SSC.A")
  expect_equal(pick_fsc(cn, "A"), "FSC.A")
  expect_equal(pick_fsc(cn, "H"), "FSC.H")

  cn2 <- c("FSC-A", "FSC-H")
  expect_equal(pick_fsc(cn2, "A"), "FSC-A")
  expect_equal(pick_fsc(cn2, "H"), "FSC-H")

  cn3 <- c("FSC_Area", "FSC_Height")
  expect_equal(pick_fsc(cn3, "A"), "FSC_Area")
  expect_equal(pick_fsc(cn3, "H"), "FSC_Height")

  cn4 <- c("SSC.A", "Time")
  expect_true(is.na(pick_fsc(cn4, "A")))
  expect_true(is.na(pick_fsc(cn4, "H")))
})

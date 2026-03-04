testthat::test_that("scale_with_bundle returns matrix with columns ordered like means", {
  df <- data.frame(
    B = c(4, 5, 6),
    A = c(1, 2, 3),
    C = c(10, 11, 12),
    check.names = FALSE
  )

  means <- c(A = 2, B = 5)
  sds   <- c(A = 1, B = 2)

  X <- scale_with_bundle(df, means = means, sds = sds)

  testthat::expect_true(is.matrix(X))
  testthat::expect_equal(dim(X), c(3, 2))
  testthat::expect_equal(colnames(X), names(means))

  # A: (1,2,3) - 2 / 1 = (-1,0,1)
  testthat::expect_equal(as.numeric(X[, "A"]), c(-1, 0, 1))

  # B: (4,5,6) - 5 / 2 = (-0.5,0,0.5)
  testthat::expect_equal(as.numeric(X[, "B"]), c(-0.5, 0, 0.5))
})

testthat::test_that("scale_with_bundle replaces NA/0 SDs with 1 (no division by zero)", {
  df <- data.frame(A = c(1, 2, 3), B = c(4, 5, 6))

  means <- c(A = 2, B = 5)
  sds   <- c(A = 0, B = NA_real_)  # should become 1,1

  X <- scale_with_bundle(df, means = means, sds = sds)

  # With sd=1:
  # A -> (-1,0,1), B -> (-1,0,1)
  testthat::expect_equal(as.numeric(X[, "A"]), c(-1, 0, 1))
  testthat::expect_equal(as.numeric(X[, "B"]), c(-1, 0, 1))
})

testthat::test_that("scale_with_bundle errors if required columns are missing", {
  df <- data.frame(A = 1:3, B = 4:6)
  means <- c(A = 2, D = 0)  # D not present
  sds   <- c(A = 1, D = 1)

  testthat::expect_error(
    scale_with_bundle(df, means = means, sds = sds),
    "names\\(means\\).*colnames\\(df_mat\\)|names\\(sds\\).*colnames\\(df_mat\\)"
  )
})

testthat::test_that("auto_map_features returns NA mapping when inputs are empty", {
  out1 <- auto_map_features(character(0), c("A", "B"))
  testthat::expect_true(is.character(out1))
  testthat::expect_length(out1, 0)

  out2 <- auto_map_features(c("A", "B"), character(0))
  testthat::expect_equal(out2, c(A = NA_character_, B = NA_character_))

  out3 <- auto_map_features(character(0), character(0))
  testthat::expect_length(out3, 0)
})

testthat::test_that("auto_map_features matches case-insensitive exact names", {
  model_feats <- c("CD3", "CD19")
  dataset_cols <- c("cd3", "Cd19", "Other")

  out <- auto_map_features(model_feats, dataset_cols)

  testthat::expect_equal(out[["CD3"]], "cd3")
  testthat::expect_equal(out[["CD19"]], "Cd19")
})

testthat::test_that("auto_map_features matches via make.names heuristic", {
  model_feats <- c("Live/Dead", "HLA-DR")
  dataset_cols <- c("Live.Dead", "HLA.DR", "X")

  out <- auto_map_features(model_feats, dataset_cols)

  testthat::expect_equal(out[["Live/Dead"]], "Live.Dead")
  testthat::expect_equal(out[["HLA-DR"]], "HLA.DR")
})

testthat::test_that("auto_map_features matches via canonicalised strings (strip punctuation)", {
  model_feats <- c("CD45RA", "CD45-RA")
  dataset_cols <- c("CD45_RA", "foo")

  out <- auto_map_features(model_feats, dataset_cols)

  # Both canon -> cd45ra, should match first hit "CD45_RA"
  testthat::expect_equal(out[["CD45RA"]], "CD45_RA")
  testthat::expect_equal(out[["CD45-RA"]], "CD45_RA")
})

testthat::test_that("auto_map_features soft-matches startsWith when no exact match", {
  model_feats <- c("FSC-A", "SSC-A")
  dataset_cols <- c("FSC.A", "SSC.A.Width")  # SSC-A should soft match SSC.A.Width

  out <- auto_map_features(model_feats, dataset_cols)

  testthat::expect_equal(out[["FSC-A"]], "FSC.A")
  testthat::expect_equal(out[["SSC-A"]], "SSC.A.Width")
})

testthat::test_that("auto_map_features returns NA when no match found", {
  model_feats <- c("Nope1", "Nope2")
  dataset_cols <- c("A", "B", "C")

  out <- auto_map_features(model_feats, dataset_cols)

  testthat::expect_true(all(is.na(out)))
  testthat::expect_equal(names(out), model_feats)
})

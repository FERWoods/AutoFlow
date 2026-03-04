testthat::test_that("assign_by_knn returns correct lengths and scores in [0,1]", {
  testthat::skip_if_not_installed("RANN")

  set.seed(1)
  pca_train <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
  labels_train <- rep(c("A", "B"), each = 50)

  pca_query <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)

  out <- assign_by_knn(pca_train, labels_train, pca_query, k = 7)

  testthat::expect_type(out, "list")
  testthat::expect_true(all(c("label", "score") %in% names(out)))
  testthat::expect_length(out$label, nrow(pca_query))
  testthat::expect_length(out$score, nrow(pca_query))
  testthat::expect_true(all(is.finite(out$score)))
  testthat::expect_true(all(out$score >= 0 & out$score <= 1))
})

testthat::test_that("assign_by_knn errors on invalid k", {
  pca_train <- matrix(0, nrow = 5, ncol = 2)
  labels_train <- letters[1:5]
  pca_query <- matrix(0, nrow = 2, ncol = 2)

  testthat::expect_error(
    assign_by_knn(pca_train, labels_train, pca_query, k = 0),
    "k must be"
  )
})

testthat::test_that("run_unsupervised_func adds Timepoint and assignment columns (no subsample branch)", {
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("SeuratObject")

  set.seed(1)
  counts <- matrix(abs(rpois(30 * 80, lambda = 5)), nrow = 30, ncol = 80)
  rownames(counts) <- paste0("M", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell_", seq_len(ncol(counts)))

  meta <- data.frame(dummy = rep("x", ncol(counts)), row.names = colnames(counts))
  obj <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta)

  # Force no-subsample branch (quiet)
  out <- suppressWarnings(suppressMessages(
    run_unsupervised_func(
      obj,
      res = 0.2,
      logfold = 0.25,
      percentage_cells = 0.1,
      max_umap_train = NULL
    )
  ))

  testthat::expect_s4_class(out, "Seurat")
  testthat::expect_true("Timepoint" %in% colnames(out@meta.data))
  testthat::expect_true("assignment" %in% colnames(out@meta.data))
  testthat::expect_true("assignment_conf" %in% colnames(out@meta.data))

  testthat::expect_equal(length(out@meta.data$assignment), ncol(out))
  testthat::expect_true(all(!is.na(out@meta.data$assignment)))
  testthat::expect_true(all(out@meta.data$assignment_conf == 1))
})

testthat::test_that("run_unsupervised_func subsample branch returns assignment for all cells and confidence in [0,1]", {
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("RANN")

  set.seed(2)
  counts <- matrix(abs(rpois(30 * 120, lambda = 5)), nrow = 30, ncol = 120)
  rownames(counts) <- paste0("M", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell_", seq_len(ncol(counts)))

  meta <- data.frame(dummy = rep("x", ncol(counts)), row.names = colnames(counts))
  obj <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta)

  max_train <- 50L

  # Subsample branch (quiet)
  out <- suppressWarnings(suppressMessages(
    run_unsupervised_func(
      obj,
      res = 0.2,
      logfold = 0.25,
      percentage_cells = 0.1,
      max_umap_train = max_train,
      seed = 123,
      k_transfer = 7,
      min_transfer_conf = 0.6,
      low_conf_label = "Uncertain"
    )
  ))

  testthat::expect_s4_class(out, "Seurat")
  testthat::expect_true("assignment" %in% colnames(out@meta.data))
  testthat::expect_true("assignment_conf" %in% colnames(out@meta.data))
  testthat::expect_true("seurat_clusters" %in% colnames(out@meta.data))

  testthat::expect_equal(length(out@meta.data$assignment), ncol(out))
  testthat::expect_equal(length(out@meta.data$assignment_conf), ncol(out))
  testthat::expect_equal(length(out@meta.data$seurat_clusters), ncol(out))

  conf <- out@meta.data$assignment_conf
  testthat::expect_true(all(is.finite(conf)))
  testthat::expect_true(all(conf >= 0 & conf <= 1))

  low_idx <- which(conf < 0.6)
  if (length(low_idx)) {
    testthat::expect_true(all(out@meta.data$assignment[low_idx] == "Uncertain"))
  }
})

testthat::test_that("run_unsupervised_func treats non-sensical max_umap_train as 'no subsample'", {
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("SeuratObject")

  set.seed(3)
  counts <- matrix(abs(rpois(20 * 60, lambda = 5)), nrow = 20, ncol = 60)
  rownames(counts) <- paste0("M", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell_", seq_len(ncol(counts)))

  meta <- data.frame(dummy = rep("x", ncol(counts)), row.names = colnames(counts))
  obj <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta)

  out <- suppressWarnings(suppressMessages(
    run_unsupervised_func(
      obj,
      max_umap_train = -1
    )
  ))

  testthat::expect_true(all(out@meta.data$assignment_conf == 1))
})

test_that("assign_by_knn recovers labels on an easy separable problem", {
  skip_if_not_installed("RANN")

  set.seed(1)

  # Two well-separated clusters in 2D
  n_train <- 200
  n_query <- 100

  A_train <- matrix(rnorm(n_train * 2, mean = 0,  sd = 0.3), ncol = 2)
  B_train <- matrix(rnorm(n_train * 2, mean = 5,  sd = 0.3), ncol = 2)
  pca_train <- rbind(A_train, B_train)
  labels_train <- c(rep("A", n_train), rep("B", n_train))

  A_query <- matrix(rnorm(n_query * 2, mean = 0,  sd = 0.3), ncol = 2)
  B_query <- matrix(rnorm(n_query * 2, mean = 5,  sd = 0.3), ncol = 2)
  pca_query <- rbind(A_query, B_query)
  true_labels <- c(rep("A", n_query), rep("B", n_query))

  out <- assign_by_knn(pca_train, labels_train, pca_query, k = 15)

  expect_type(out$label, "character")
  expect_type(out$score, "double")
  expect_length(out$label, nrow(pca_query))
  expect_length(out$score, nrow(pca_query))
  expect_true(all(out$score >= 0 & out$score <= 1))

  acc <- mean(out$label == true_labels)
  expect_gt(acc, 0.98)
  expect_gt(mean(out$score), 0.9)
})

test_that("assign_by_knn respects k bounds and errors on bad inputs", {
  skip_if_not_installed("RANN")

  pca_train <- matrix(rnorm(20), ncol = 2)
  labels_train <- rep("A", nrow(pca_train))
  pca_query <- matrix(rnorm(10), ncol = 2)

  expect_error(assign_by_knn(pca_train, labels_train[-1], pca_query, k = 3))
  expect_error(assign_by_knn(pca_train, labels_train, pca_query, k = 0))
})

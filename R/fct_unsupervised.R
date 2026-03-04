#' Unsupervised clustering + UMAP with optional fast subsampling and kNN transfer
#'
#' Runs an unsupervised Seurat workflow (Normalize → HVGs → Scale → PCA → UMAP → clustering)
#' and assigns human-readable cluster labels (`assignment`) using marker strings when possible.
#'
#' If `max_umap_train` is set and the dataset is larger than this value, a fast mode is used:
#' UMAP and clustering are trained on a random subset of cells, the learned UMAP model is
#' projected to all remaining cells, and `seurat_clusters` / `assignment` for the remaining
#' cells are transferred by k-nearest neighbours (kNN) in PCA space.
#'
#' @param flow_data A Seurat object (cells as columns) containing per-cell expression/features.
#' @param res Numeric. Clustering resolution passed to `Seurat::FindClusters()`.
#' @param logfold Numeric. Log-fold-change threshold passed to `Seurat::FindAllMarkers()`.
#' @param percentage_cells Numeric in [0,1]. Minimum fraction of cells expressing a feature
#'   used by `Seurat::FindAllMarkers()` (`min.pct`).
#' @param batch_correct Logical. Reserved for future use (currently not applied).
#' @param max_umap_train Integer or NULL. If NULL, run the full workflow on all cells.
#'   If set (e.g. `100000`) and `ncol(flow_data) > max_umap_train`, UMAP + clustering are trained
#'   on a random subset of size `max_umap_train`. The UMAP model is then projected to all cells
#'   via `Seurat::ProjectUMAP()`, and labels for non-training cells are assigned using kNN in
#'   PCA space (see `k_transfer`).
#' @param seed Integer. RNG seed used for subsampling when `max_umap_train` triggers fast mode.
#' @param k_transfer Integer. Number of nearest neighbours used for kNN label transfer from
#'   training cells to query cells in PCA space.
#' @param min_transfer_conf Numeric in [0,1] or NULL. If set, query cells with
#'   `assignment_conf < min_transfer_conf` are labelled as `low_conf_label`.
#' @param low_conf_label Character. Label used for low-confidence transferred assignments.
#'
#' @details
#' The returned object contains:
#' \itemize{
#'   \item A UMAP reduction named `"umap"` (3 components are requested; if fewer are available,
#'         the reduction will contain fewer dimensions).
#'   \item `meta.data$seurat_clusters` (cluster IDs as character).
#'   \item `meta.data$assignment` (marker-string label if marker assignment is successful,
#'         otherwise falls back to cluster IDs).
#'   \item `meta.data$assignment_conf` (confidence score in [0,1] for query cells in fast mode,
#'         defined as the fraction of kNN votes supporting the assigned label; training cells
#'         are set to 1).
#' }
#'
#' Marker-based assignment requires at least two clusters; if only one cluster is detected,
#' marker discovery and marker-string assignment are skipped with a warning.
#'
#' @return A Seurat object with UMAP, clustering, and per-cell labels stored in `meta.data`.
#' @export
run_unsupervised_func <- function(flow_data,
                                  res = 0.5,
                                  logfold = 0.25,
                                  percentage_cells = 0.25,
                                  batch_correct = FALSE,
                                  max_umap_train = NULL,
                                  seed = 1L,
                                  k_transfer = 15L,
                                  min_transfer_conf = NULL,
                                  low_conf_label = "Uncertain") {

  if (!("Timepoint" %in% colnames(flow_data@meta.data))) {
    flow_data@meta.data$Timepoint <- NA
  }

  ref <- flow_data

  # ----------------------------
  # Preprocess on ALL cells
  # ----------------------------
  ref <- Seurat::NormalizeData(ref)
  ref <- Seurat::FindVariableFeatures(ref)
  ref <- Seurat::ScaleData(ref)

  # Choose a conservative npcs for tiny test objects
  nvar <- length(Seurat::VariableFeatures(ref))
  npcs <- max(2L, min(30L, nvar, ncol(ref) - 1L))
  ref <- Seurat::RunPCA(ref, features = Seurat::VariableFeatures(ref), npcs = npcs)

  # Safe dims for downstream
  ndims <- min(npcs, ncol(ref[["pca"]]@cell.embeddings))
  dims_use <- seq_len(ndims)

  n_cells   <- ncol(ref)
  all_cells <- colnames(ref)

  do_subsample <- is.numeric(max_umap_train) &&
    is.finite(max_umap_train) &&
    max_umap_train > 0 &&
    n_cells > max_umap_train

  # Helper: can we safely do markers?
  .can_do_markers <- function(obj, stage_label = "") {
    clu <- Seurat::Idents(obj)
    n_clu <- length(unique(as.character(clu)))
    if (n_clu < 2) {
      warning(sprintf(
        "run_unsupervised_func: only %s cluster detected%s; skipping FindAllMarkers + marker-based assignment.",
        n_clu,
        if (nzchar(stage_label)) paste0(" (", stage_label, ")") else ""
      ), call. = FALSE)
      return(FALSE)
    }
    TRUE
  }

  .safe_markers_join <- function(markers.pos_neg, markers.pos) {
    if (is.null(markers.pos_neg) || is.null(markers.pos)) return(NULL)
    if (!is.data.frame(markers.pos_neg) || !is.data.frame(markers.pos)) return(NULL)
    if (nrow(markers.pos_neg) == 0 || nrow(markers.pos) == 0) return(NULL)

    need <- c("cluster", "gene")
    if (!all(need %in% names(markers.pos_neg))) return(NULL)
    if (!all(need %in% names(markers.pos))) return(NULL)

    dplyr::full_join(markers.pos_neg, markers.pos, by = c("cluster", "gene"))
  }

  # ----------------------------
  # Marker assignment WITHOUT dplyr joins into meta.data
  # ----------------------------
  .apply_marker_assignment <- function(obj, all_markers) {

    # Fallback helper: assignments = clusters
    .fallback <- function() {
      obj@meta.data$assignment <- as.character(Seurat::Idents(obj))
      obj@meta.data$assignment_conf <- 1
      obj
    }

    if (is.null(all_markers) || !nrow(all_markers)) return(.fallback())

    all_markers$gene <- sub("__dup.*$", "", all_markers$gene)

    # unify FC column naming
    if ("avg_log2FC.x" %in% colnames(all_markers)) {
      fc <- all_markers$avg_log2FC.x
    } else if ("avg_log2FC" %in% colnames(all_markers)) {
      fc <- all_markers$avg_log2FC
    } else if ("avg_logFC" %in% colnames(all_markers)) {
      fc <- all_markers$avg_logFC
    } else {
      fc <- rep(NA_real_, nrow(all_markers))
    }
    all_markers$avg_log2FC_join <- fc

    if (all(is.na(all_markers$avg_log2FC_join))) {
      warning("run_unsupervised_func: marker table has no logFC column; using clusters as assignment.", call. = FALSE)
      return(.fallback())
    }

    neg <- all_markers[all_markers$avg_log2FC_join < 0, , drop = FALSE]
    pos <- all_markers[all_markers$avg_log2FC_join > 0, , drop = FALSE]
    if (!nrow(neg) || !nrow(pos)) return(.fallback())

    all_table_neg <- table(neg$cluster, neg$gene)
    all_table_pos <- table(pos$cluster, pos$gene)

    clusters <- intersect(rownames(all_table_neg), rownames(all_table_pos))
    if (!length(clusters)) return(.fallback())

    all_labels <- lapply(clusters, function(cl) {
      labs_neg <- paste0(colnames(all_table_neg)[all_table_neg[cl, ] == 1], "-")
      labs_pos <- paste0(colnames(all_table_pos)[all_table_pos[cl, ] == 1], "+")
      data.frame(
        cluster = as.character(cl),
        assignment = paste(paste(labs_pos, collapse = ""), paste(labs_neg, collapse = ""), collapse = ""),
        stringsAsFactors = FALSE
      )
    })

    manual <- do.call(rbind, all_labels)

    # --- SAFE join: match seurat_clusters -> manual$cluster ---
    sc <- as.character(obj@meta.data$seurat_clusters)
    m  <- match(sc, manual$cluster)

    obj@meta.data$assignment <- manual$assignment[m]

    if (!"assignment" %in% colnames(obj@meta.data) || all(is.na(obj@meta.data$assignment))) {
      obj@meta.data$assignment <- as.character(Seurat::Idents(obj))
    }
    obj@meta.data$assignment_conf <- 1
    obj
  }

  # ----------------------------
  # Branch 1: Standard (no subsample)
  # ----------------------------
  if (!do_subsample) {

    ref <- Seurat::RunUMAP(ref, dims = dims_use, n.components = 3L, return.model = TRUE)
    ref <- Seurat::FindNeighbors(ref, dims = dims_use, reduction = "pca")

    algo_to_use <- if ("algorithm" %in% names(formals(Seurat::FindClusters))) 4 else 1
    message(sprintf("Running clustering with algorithm = %s (%s)",
                    algo_to_use, ifelse(algo_to_use == 4, "Leiden", "Louvain")))
    ref <- Seurat::FindClusters(ref, resolution = res, algorithm = algo_to_use)

    if (.can_do_markers(ref, "full")) {
      suppressWarnings({
        markers.pos     <- Seurat::FindAllMarkers(ref, only.pos = TRUE,  min.pct = percentage_cells, logfc.threshold = logfold)
        markers.pos_neg <- Seurat::FindAllMarkers(ref, only.pos = FALSE, min.pct = percentage_cells, logfc.threshold = logfold)
      })
      all_markers <- .safe_markers_join(markers.pos_neg, markers.pos)
      ref <- .apply_marker_assignment(ref, all_markers)
    } else {
      ref@meta.data$assignment <- as.character(Seurat::Idents(ref))
      ref@meta.data$assignment_conf <- 1
    }

    rownames(ref@meta.data) <- colnames(ref)
    if (!"assignment" %in% colnames(ref@meta.data)) {
      ref@meta.data$assignment <- as.character(Seurat::Idents(ref))
    }
    if (!"assignment_conf" %in% colnames(ref@meta.data)) {
      ref@meta.data$assignment_conf <- 1
    }

    return(ref)
  }

  # ----------------------------
  # Branch 2: Subsample + project + kNN transfer
  # ----------------------------
  set.seed(as.integer(seed))
  train_cells <- sample(all_cells, size = as.integer(max_umap_train))
  query_cells <- setdiff(all_cells, train_cells)

  # Seurat::subset() not exported
  ref_train <- subset(ref, cells = train_cells)

  ref_train <- Seurat::RunUMAP(
    ref_train,
    dims = dims_use,
    n.components = 3L,
    return.model = TRUE
  )

  ref <- Seurat::ProjectUMAP(
    reference = ref_train,
    query = ref,
    reference.reduction = "pca",
    query.reduction = "pca",
    reduction.model = "umap",
    reduction.name = "umap",     # <- make the name explicit
    reduction.key  = "UMAP_"     # <- gives UMAP_1, UMAP_2, ...
  )

  # --- Ensure plotting can find the embeddings (many apps look in meta.data) ---
  um <- Seurat::Embeddings(ref, "umap")
  # defensive: ensure rownames match cells
  um <- um[colnames(ref), , drop = FALSE]

  ref@meta.data$UMAP_1 <- um[, 1]
  ref@meta.data$UMAP_2 <- um[, 2]
  if (ncol(um) >= 3) ref@meta.data$UMAP_3 <- um[, 3]
  ref_train <- Seurat::FindNeighbors(ref_train, dims = dims_use, reduction = "pca")

  algo_to_use <- if ("algorithm" %in% names(formals(Seurat::FindClusters))) 4 else 1
  message(sprintf("Running clustering with algorithm = %s (%s) on training subset (n=%s)",
                  algo_to_use, ifelse(algo_to_use == 4, "Leiden", "Louvain"),
                  format(length(train_cells))))

  ref_train <- Seurat::FindClusters(ref_train, resolution = res, algorithm = algo_to_use)

  if (.can_do_markers(ref_train, "train")) {
    suppressWarnings({
      markers.pos     <- Seurat::FindAllMarkers(ref_train, only.pos = TRUE,  min.pct = percentage_cells, logfc.threshold = logfold)
      markers.pos_neg <- Seurat::FindAllMarkers(ref_train, only.pos = FALSE, min.pct = percentage_cells, logfc.threshold = logfold)
    })
    all_markers <- .safe_markers_join(markers.pos_neg, markers.pos)
    ref_train <- .apply_marker_assignment(ref_train, all_markers)
  } else {
    ref_train@meta.data$assignment <- as.character(Seurat::Idents(ref_train))
    ref_train@meta.data$assignment_conf <- 1
  }

  rownames(ref_train@meta.data) <- colnames(ref_train)
  if (!"assignment" %in% colnames(ref_train@meta.data)) {
    ref_train@meta.data$assignment <- as.character(Seurat::Idents(ref_train))
  }

  # --- kNN transfer in PCA space ---
  pca_all   <- ref[["pca"]]@cell.embeddings
  pca_train <- pca_all[train_cells, , drop = FALSE]
  pca_query <- pca_all[query_cells, , drop = FALSE]

  knn_asg <- assign_by_knn(
    pca_train = pca_train,
    labels_train = ref_train@meta.data[train_cells, "assignment", drop = TRUE],
    pca_query = pca_query,
    k = k_transfer
  )

  knn_clu <- assign_by_knn(
    pca_train = pca_train,
    labels_train = ref_train@meta.data[train_cells, "seurat_clusters", drop = TRUE],
    pca_query = pca_query,
    k = k_transfer
  )

  ref@meta.data$assignment <- NA_character_
  ref@meta.data$assignment_conf <- NA_real_
  ref@meta.data$seurat_clusters <- NA

  ref@meta.data[train_cells, "assignment"] <- ref_train@meta.data[train_cells, "assignment", drop = TRUE]
  ref@meta.data[train_cells, "assignment_conf"] <- 1
  ref@meta.data[train_cells, "seurat_clusters"] <- as.character(ref_train@meta.data[train_cells, "seurat_clusters", drop = TRUE])

  ref@meta.data[query_cells, "assignment"] <- knn_asg$label
  ref@meta.data[query_cells, "assignment_conf"] <- knn_asg$score
  ref@meta.data[query_cells, "seurat_clusters"] <- as.character(knn_clu$label)

  if (!is.null(min_transfer_conf)) {
    min_transfer_conf <- as.numeric(min_transfer_conf)
    low <- which(!is.na(ref@meta.data$assignment_conf) & ref@meta.data$assignment_conf < min_transfer_conf)
    if (length(low)) ref@meta.data$assignment[low] <- low_conf_label
  }

  rownames(ref@meta.data) <- colnames(ref)
  ref
}
#' Assign labels by k-nearest neighbours in PCA space
#'
#' @param pca_train Matrix (n_train x d) PCA embeddings for training cells.
#' @param labels_train Character/factor length n_train; labels for training cells.
#' @param pca_query Matrix (n_query x d) PCA embeddings for query cells.
#' @param k Integer; number of neighbours.
#' @return List with `label` (character vector length n_query) and `score`
#'   (fraction of neighbours voting for assigned label, in [0,1]).
#' @noRd
assign_by_knn <- function(pca_train, labels_train, pca_query, k = 15L) {
  stopifnot(is.matrix(pca_train), is.matrix(pca_query))
  stopifnot(nrow(pca_train) == length(labels_train))
  k <- as.integer(k)
  if (k < 1L) stop("k must be >= 1")

  # RANN gives fast exact NN; if not installed, you could fall back to FNN.
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Package 'RANN' is required for kNN assignment.")
  }

  nn <- RANN::nn2(data = pca_train, query = pca_query, k = min(k, nrow(pca_train)))
  idx <- nn$nn.idx

  labs <- as.character(labels_train)
  out_label <- character(nrow(pca_query))
  out_score <- numeric(nrow(pca_query))

  for (i in seq_len(nrow(idx))) {
    neigh <- labs[idx[i, ]]
    tab <- sort(table(neigh), decreasing = TRUE)
    out_label[i] <- names(tab)[1]
    out_score[i] <- as.numeric(tab[1]) / length(neigh)
  }

  list(label = out_label, score = out_score)
}

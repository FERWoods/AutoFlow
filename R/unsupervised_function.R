#' @export

run_unsupervised <- function(flow_data, res = 0.5, logfold = 0.25, percentage_cells = 0.25, batch_correct = FALSE) {
  library(Seurat)
  library(dplyr)
  if (!("Timepoint" %in% colnames(flow_data))) {
    cat("Creating 'Timepoint' column with NULL values\n")
    flow_data$Timepoint <- rep(NA, nrow(flow_data))
  }

  ref_seurat <- flow_data

  if (batch_correct && length(unique(ref_seurat$Timepoint)) > 1) {
    cat("Running batch correction based on timepoint\n")
    day.list <- SplitObject(ref_seurat, split.by = "Timepoint")
    features <- SelectIntegrationFeatures(object.list = day.list, nfeatures = 3000)
    anchors <- FindIntegrationAnchors(object.list = day.list, anchor.features = features,
                                      normalization.method = "LogNormalize", reduction = "cca",
                                      dims = 1:(nrow(ref_seurat) - 1), verbose = TRUE, k.filter = 7)
    day.integrated <- IntegrateData(anchorset = anchors,
                                    normalization.method = "LogNormalize",
                                    dims = 1:(nrow(ref_seurat) - 1), features = features,
                                    features.to.integrate = features,
                                    k.weight = 5)
    day.integrated <- ScaleData(day.integrated)
    day.integrated <- RunPCA(day.integrated, verbose = TRUE, features = features)
    day.integrated <- RunUMAP(day.integrated, dims = 1:ncol(ref_seurat$pca), n.components = 3L)
    day.integrated <- FindNeighbors(day.integrated, dims = 1:ncol(ref_seurat$pca), reduction = "pca", verbose = TRUE)
    day.integrated <- FindClusters(day.integrated, resolution = res, verbose = TRUE)
    ref_seurat <- day.integrated
  } else if (!batch_correct || (batch_correct && length(unique(ref_seurat$Timepoint)) == 1)) {
    if (batch_correct && length(unique(ref_seurat$Timepoint)) == 1) {
      cat("Cannot run batch correction as only 1 study day found")
    }
    ref_seurat <- NormalizeData(ref_seurat)
    ref_seurat <- FindVariableFeatures(ref_seurat)#, loess.span = 10)
    ref_seurat <- ScaleData(ref_seurat)
    ref_seurat <- RunPCA(ref_seurat, "RNA")
    ref_seurat <- RunUMAP(ref_seurat,dims=1:ncol(ref_seurat$pca), n.components = 3L, return.model = TRUE)
    ref_seurat <- FindNeighbors(ref_seurat, dims = 1:ncol(ref_seurat$pca), reduction = "pca", verbose = TRUE)
    ref_seurat <- FindClusters(ref_seurat, resolution = res, verbose = TRUE)
  }
  markers.pos <- FindAllMarkers(ref_seurat, only.pos = TRUE, min.pct = percentage_cells, logfc.threshold = logfold)
  markers.pos_neg <- FindAllMarkers(ref_seurat, only.pos = FALSE, min.pct = percentage_cells, logfc.threshold = logfold)

  all_markers <- merge(markers.pos_neg, markers.pos, by = c("cluster", "gene"), all = TRUE)
  all_table_neg <- table(all_markers[all_markers$avg_log2FC.x < 0, ]$cluster, all_markers[all_markers$avg_log2FC.x < 0, ]$gene)
  all_table_pos <- table(all_markers[all_markers$avg_log2FC.x > 0, ]$cluster, all_markers[all_markers$avg_log2FC.x > 0, ]$gene)

  all_labels <- lapply(1:nrow(all_table_neg), function(i) {
    labels_tmp_neg <- paste0(colnames(all_table_neg)[all_table_neg[i, ] == 1], "-")
    labels_tmp_pos <- paste0(colnames(all_table_pos)[all_table_pos[i, ] == 1], "+")

    c("cluster" = rownames(all_table_pos)[i],
      "assignment" = paste(paste(labels_tmp_pos, collapse = ""), paste(labels_tmp_neg, collapse = ""), collapse = "")
    )
  })

  manual_labelling <- data.frame(do.call("rbind", all_labels))
  manual_labelling_sorted <- manual_labelling[order(as.numeric(manual_labelling$cluster)), ]

  new.cluster.ids <- manual_labelling_sorted$assignment
  names(new.cluster.ids) <- levels(ref_seurat)

  ref_seurat@meta.data <- ref_seurat@meta.data %>%
    left_join(manual_labelling_sorted, by = c("seurat_clusters" = "cluster"))

  rownames(ref_seurat@meta.data) <- colnames(ref_seurat)

  return(ref_seurat)
}

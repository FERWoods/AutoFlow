# functions

# pre-process functions
preprocess_2 <- function(ff){
  ff <- flowCore::compensate(ff, flowCore::keyword(ff)$SPILL)
  ff <- flowCore::transform(ff,
                            flowCore::estimateLogicle(ff,
                                                      colnames(flowCore::keyword(ff)$SPILL)))
  #scale(ff)
}
preprocess_3 <- function(ff){
  ff <- flowCore::compensate(ff, flowCore::keyword(ff)$SPILL)
  ff <- flowCore::transform(ff,
                            flowCore::estimateLogicle(ff,
                                                      colnames(flowCore::keyword(ff)$SPILL)))
  ff <- flowStats::gaussNorm(flowSet(ff), channel.names = colnames(flowCore::keyword(ff)$SPILL))
  #scale(ff)
}
channel_select <- function(ff){
  library(data.table)
  ff = ff[,c((flowCore::colnames(ff) %like% ".A") | flowCore::colnames(ff) %like% "Time") | flowCore::colnames(ff) %like% "FSC.H"]

  ff
}
# Functions for .fsc data importing and mapping to metadata

# Below sections are to allow for cbind.fill -- package origin had depreciated
cbind_fill<-function(...,fill=NULL)
{
  inputs<-list(...)
  inputs<-lapply(inputs,vert)
  maxlength<-max(unlist(lapply(inputs,len)))
  bufferedInputs<-lapply(inputs,buffer,length.out=maxlength,fill,preserveClass=FALSE)
  return(Reduce(cbind.data.frame,bufferedInputs))
}

vert<-function(object)
{
  #result<-as.data.frame(cbind(as.matrix(object)))
  if(is.list(object))
    object<-cbind(object)
  object<-data.frame(object)

  return(object)
}

len <- function(data)
{
  result<-ifelse(is.null(nrow(data)),length(data),nrow(data))
  return(result)
}

buffer<-function(x,length.out=len(x),fill=NULL,preserveClass=TRUE)
{
  xclass<-class(x)
  input<-lapply(vert(x),unlist)
  results<-as.data.frame(lapply(input,rep,length.out=length.out))
  if(length.out>len(x) && !is.null(fill))
  {
    results<-t(results)
    results[(length(unlist(x))+1):length(unlist(results))]<-fill
    results<-t(results)
  }
  if(preserveClass)
    results<-as2(results,xclass)
  return(results)
}

# Now process each with the preprocess function which applies the set spillover file and normalises with biexponential logicle function
preprocess <- function(ff) {
  ff <- compensate(ff, ff@description$SPILL)
  ff <- transform(ff, flowCore::transformList(colnames(ff@description$SPILL), flowCore::logicleTransform()))
  ff@exprs <- scale(ff@exprs)
  return(ff)
}

# function to convert wellIDs to match those attached to flowjo objects
checkWellNames = function(wellNames) {
  #' Check and convert well names to the appropriate format
  #' E.g. A1 -> A01
  #'
  #' @param wellNames String vector (e.g. c("A1", "A2", "B3", "B4"))
  #'
  #' @return String vector
  #'
  o = wellNames
  rows = substr(wellNames, 1, 1)
  stopifnot(all(rows %in% toupper(letters)[1:8]))
  columns = as.integer(substr(wellNames, 2, 10))
  stopifnot(all(columns >= 1 & columns <= 12))
  columns = as.character(columns)
  columns = sapply(columns, function(x) {
    if (nchar(x) == 1) {
      return(paste0("0", x))
    } else {
      return(x)
    }
  })
  return(paste0(rows, columns))
}

# Function to extract labels from fcs file paired with workspace --
workspace_cell_labels <- function(flowset, workspace, cell_types=cell_types_fixed, ws_group){
  groups = fj_ws_get_sample_groups(open_flowjo_xml(workspace))#, execute = FALSE)
  group_pos = groups[match(ws_group, groups$groupName), 2] + 1 # returns the group ID for the matched group +1 as the numbering begins at 0
  gating_manual = GetFlowJoLabels(flowset,
                                  workspace,
                                  group = group_pos,
                                  cellTypes = cell_types,
                                  get_data = TRUE)
  manual_labels = do.call("c",lapply(gating_manual,function(x){x$manual}))
  rm(gating_manual)
  return(manual_labels)
}


# wiht this the condition is a on points that deviate from FSC.H ~ FSC.A linearly, plus removing debrit (v. small vals in both)

# Function to remove data 5% from the line
remove_outliers <- function(x, y, p = 0.05, slope, intercept) {
  # Get the predicted values
  predicted <- slope * x + intercept

  # Calculate the residuals
  residuals <- y - predicted

  # Calculate the threshold for outliers
  threshold <- quantile(abs(residuals), 1 - p)

  # Identify the outliers
  outliers <- which(abs(residuals) > threshold)

  # Return the new data without outliers
  return(outliers)
}

log_transform_with_shift <- function(data) {
  min_value <- min(data)

  # Calculate the shift based on the minimum value
  shift <- ifelse(min_value < 0, abs(min_value) + 1, 0)

  # Log-transform the adjusted data
  log_transformed_data <- log(data + shift)

  return(log_transformed_data)
}


cell_labelling_bm <- function(seur) {
  library(stringr)
  seur$cell_type <- ifelse(str_detect(seur$assignment, "CD34\\+"), "HSC",
                           ifelse(str_detect(seur$assignment, "CD38\\+"), "LDPs",
                                  ifelse(str_detect(seur$assignment, "CD13\\+") & str_detect(seur$assignment, "CD16\\+"), "Late Granulocytes",
                                         ifelse(str_detect(seur$assignment, "CD13\\+") & str_detect(seur$assignment, "CD36\\+"), "Late Monocytes",
                                                ifelse(str_detect(seur$assignment, "CD13\\+"), "Early Granulocytes",
                                                       ifelse(str_detect(seur$assignment, "CD36\\+") & str_detect(seur$assignment, "CD71\\+"), "Early Erythroid",
                                                              ifelse(str_detect(seur$assignment, "CD41\\+") | str_detect(seur$assignment, "CD42b\\+"), "Megakaryocytes",
                                                                     ifelse(str_detect(seur$assignment, "CD235a\\+"), "Late Erythroid", seur$assignment))))))))

  return(seur)
}

norm_minmax <- function(x){
  tmp = (x - min(x)) / (max(x) - min(x))
  return(tmp)
}

#' @export

run_unsupervised_func <- function(flow_data, res = 0.5, logfold = 0.25, percentage_cells = 0.25, batch_correct = FALSE) {
  library(Seurat)
  library(dplyr)
  if (!("Timepoint" %in% colnames(flow_data@meta.data))) {
    cat("Creating 'Timepoint' column with NULL values\n")
    flow_data@meta.data$Timepoint <- NA
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
    ref_seurat <- RunPCA(ref_seurat, "RNA", features = VariableFeatures(ref_seurat))
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

# Function to read either .xlsx or .csv based on the file extension
read_file <- function(file_path) {
  library(readxl)  # For reading Excel files
  library(readr)   # For reading CSV files
  # Detect the file extension
  file_ext <- tools::file_ext(file_path)

  # Read file based on the extension
  if (file_ext == "xlsx") {
    data <- read_excel(file_path)  # Read Excel file
  } else if (file_ext == "csv") {
    data <- read_csv(file_path)    # Read CSV file
  } else {
    stop("Unsupported file format. Please provide a .xlsx or .csv file.")
  }

  return(data)
}

find_common_columns <- function(df1, df2) {
  common_cols <- list()  # Initialize an empty list to store matched columns
  matched_cols <- c()    # Keep track of already matched columns in df2

  for (col1 in colnames(df1)) {
    for (col2 in setdiff(colnames(df2), matched_cols)) {  # Exclude already matched columns
      # Check for common values between columns
      if (length(intersect(unique(df1[[col1]]), unique(df2[[col2]]))) > 0) {
        common_cols[[col1]] <- col2  # Map the match
        matched_cols <- c(matched_cols, col2)  # Mark the column in df2 as matched
        break  # Stop looking for matches for this col1 once a match is found
      }
    }
  }

  return(common_cols)
}



## Setup
suppressPackageStartupMessages({
  library(flowCore)
  library(Biobase)     # for AnnotatedDataFrame
  library(dplyr)
  library(caret)
  library(ranger)
  library(stringr)
  library(gtools)
})
suppressWarnings({ ok <- requireNamespace("data.table", quietly = TRUE) })

`%||%` <- function(a,b) if (!is.null(a)) a else b
set.seed(123)

##Paths
root_dir <- "~/OneDrive - Swansea University/Documents/AutoFlow/data/BM-MPS/"
outdir   <- "~/OneDrive - Swansea University/Documents/AutoFlow/outputs_bmmps"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Hold out this day for testing (character, e.g. "28")
test_day <- "28"

# Expected layout per day:
#   1_FCS_combined_D{DAY}.FCS ... 5_FCS_combined_D{DAY}.FCS
#   Gating_Matrix_D{DAY}.csv
get_day <- function(x) { sub("^.*_D(\\d+).*$", "\\1", x, perl = TRUE) }
get_rep <- function(x) { as.integer(sub("^([0-9]+)_FCS_combined_.*$", "\\1", basename(x))) }

all_fcs <- list.files(root_dir, pattern = "_FCS_combined_D\\d+\\.FCS$", ignore.case = TRUE, full.names = TRUE)
all_csv <- list.files(root_dir, pattern = "^Gating_Matrix_D\\d+\\.csv$", ignore.case = TRUE, full.names = TRUE)

fcs_by_day <- split(all_fcs, sapply(all_fcs, get_day))
csv_by_day <- split(all_csv, sapply(all_csv, get_day))
days <- mixedsort(intersect(names(fcs_by_day), names(csv_by_day)))
stopifnot(length(days) > 0, test_day %in% days)

##Robust reader for gating matrices
read_gating_file <- function(path) {
  if (ok) {
    df <- try(data.table::fread(path, data.table = FALSE, showProgress = FALSE), silent = TRUE)
    if (!inherits(df, "try-error") && ncol(df) > 1) {
      message(sprintf("Read (fread): %s  (cols=%d, rows=%d)", basename(path), ncol(df), nrow(df)))
      return(df)
    }
  }
  raw4 <- try(readBin(path, what = "raw", n = 4), silent = TRUE)
  enc  <- "UTF-8"
  if (!inherits(raw4, "try-error") && length(raw4) >= 2) {
    if (raw4[1] == as.raw(0xFF) && raw4[2] == as.raw(0xFE)) enc <- "UTF-16LE"
    if (raw4[1] == as.raw(0xFE) && raw4[2] == as.raw(0xFF)) enc <- "UTF-16BE"
  }
  peek <- try(readLines(path, n = 5, encoding = if (enc == "UTF-8") "UTF-8" else enc), silent = TRUE)
  peek <- if (inherits(peek, "try-error")) character(0) else peek
  first <- if (length(peek)) peek[which(nchar(trimws(peek)) > 0)[1]] else ""

  delim <- if (grepl("\t", first, fixed = TRUE)) "\t" else if (grepl(";", first, fixed = TRUE)) ";" else ","
  rd <- function(sep, fe="") try(read.table(path, header=TRUE, sep=sep,
                                            fileEncoding=fe, check.names=FALSE,
                                            quote="\"", fill=TRUE,
                                            comment.char="", stringsAsFactors=FALSE),
                                 silent = TRUE)
  fe_use <- if (enc %in% c("UTF-16LE","UTF-16BE")) enc else ""
  df <- rd(delim, fe_use); if (inherits(df, "try-error") || ncol(df) == 1) df <- rd("\t", fe_use)
  if (inherits(df, "try-error") || ncol(df) == 1) df <- rd(";", fe_use)
  if (inherits(df, "try-error") || ncol(df) == 1) df <- rd(",")
  message(sprintf("Read (base): %s  (cols=%d, rows=%d)", basename(path), ncol(df), nrow(df)))
  if (ncol(df) == 1) warning("Gating file still parsed as 1 column. Check delimiter/encoding.")
  df
}

## Gate resolution + coercion helpers
.canon <- function(x) tolower(gsub("[^a-z0-9]+", "", x))

resolve_col <- function(gdf, key) {
  # exact case-insensitive first
  hit_exact <- which(tolower(colnames(gdf)) == tolower(key))
  if (length(hit_exact) == 1) return(colnames(gdf)[hit_exact])

  # canonical exact
  can_key  <- .canon(key); can_cols <- .canon(colnames(gdf))
  hit <- match(can_key, can_cols); if (!is.na(hit)) return(colnames(gdf)[hit])

  # canonical contains (unique)
  cand <- grep(can_key, can_cols, fixed = TRUE); if (length(cand) == 1) return(colnames(gdf)[cand])

  NA_character_
}

to_logi <- function(v) {
  if (is.logical(v)) return(replace(v, is.na(v), FALSE))
  if (is.numeric(v)) return(!is.na(v) & v != 0)
  v <- trimws(tolower(as.character(v)))
  true_vals  <- c("true","t","1","yes","y","pos","+")
  false_vals <- c("false","f","0","no","n","neg","-","")
  out <- ifelse(v %in% true_vals, TRUE, ifelse(v %in% false_vals, FALSE, NA))
  replace(out, is.na(out), FALSE)
}

gget <- function(gdf, key) {
  col <- resolve_col(gdf, key)
  if (is.na(col)) return(rep(FALSE, nrow(gdf)))
  to_logi(gdf[[col]])
}

pick_viability <- function(gdf) {
  col_v  <- resolve_col(gdf, "Viable Cells")
  col_nv <- resolve_col(gdf, "Non Viable Cells")
  if (!is.na(col_v))  return(to_logi(gdf[[col_v]]))
  if (!is.na(col_nv)) return(!to_logi(gdf[[col_nv]]))
  # heuristic fallback
  cols <- colnames(gdf); can <- .canon(cols)
  pats <- c("viable","viability","live","dead","7aad","amine","dapi","dump")
  hits <- unique(unlist(lapply(pats, function(p) which(grepl(p, can)))))
  if (!length(hits)) return(rep(TRUE, nrow(gdf)))
  best <- NULL; best_score <- -Inf
  for (i in hits) {
    v <- to_logi(gdf[[cols[i]]]); if (grepl("dead|dump", can[i])) v <- !v
    prop <- mean(v); score <- -abs(prop - 0.7)
    if (is.finite(score) && score > best_score) { best <- v; best_score <- score }
  }
  if (is.null(best)) best <- rep(TRUE, nrow(gdf))
  best
}

pick_singlets <- function(gdf) {
  col_sc <- resolve_col(gdf, "Single Cells")
  if (!is.na(col_sc)) return(to_logi(gdf[[col_sc]]))
  cols <- colnames(gdf); can <- .canon(cols)
  pats <- c("singlecells","singlets","singlet","single")
  hits <- unique(unlist(lapply(pats, function(p) which(grepl(p, can)))))
  if (!length(hits)) return(rep(TRUE, nrow(gdf)))
  best <- to_logi(gdf[[cols[hits[1]]]])
  if (mean(best) < 0.1) best <- rep(TRUE, nrow(gdf))
  best
}

## EXACT gating logic + precedence -- aligns with experimental defs
mk_labels_one_day <- function(gdf) {
  sc  <- ifelse(gget(gdf, "Single Cells"), 1L, 0L)
  vc  <- ifelse(gget(gdf, "Viable Cells"), 1L, 0L)

  hsc <- ifelse(gget(gdf, "CD34+ CD38- HSCs V3") &
                  gget(gdf, "CD41- CD16-") &
                  gget(gdf, "CD13-"), 1L, 0L)

  mgk <- ifelse(gget(gdf, "CD36+- Megakaryocytes V3") &
                  gget(gdf, "CD13-") &
                  gget(gdf, "CD13-/CD41+ CD16-"), 1L, 0L)

  ee  <- ifelse(gget(gdf, "CD71+ CD235+- Early Erythroid V3") &
                  gget(gdf, "CD13-") &
                  gget(gdf, "CD13-/CD36+"), 1L, 0L)

  le  <- ifelse(gget(gdf, "CD71+- CD235+ Later Erythroid V3") &
                  gget(gdf, "CD13-") &
                  gget(gdf, "CD13-/CD36-"), 1L, 0L)

  eg  <- ifelse(gget(gdf, "CD41- CD16- Early Granulocytes") &
                  gget(gdf, "CD13+ CD235-"), 1L, 0L)

  lg  <- ifelse(gget(gdf, "CD41- CD16+ Late Granulocytes") &
                  gget(gdf, "CD13+ CD235-"), 1L, 0L)

  lm  <- ifelse(gget(gdf, "CD36+ CD16+- Late Monocytes") &
                  gget(gdf, "CD13+ CD235-") &
                  gget(gdf, "CD41-"), 1L, 0L)

  ldp <- ifelse(gget(gdf, "CD34+ CD38+ Lineage Differentiated Progenitors"), 1L, 0L)

  # precedence
  out <- ifelse(vc == 0L, "non-viable",
                ifelse(sc == 0L, "debris/doublet",
                       ifelse(lm == 1L, "Late Monocytes",
                              ifelse(lg == 1L, "Late Granulocytes",
                                     ifelse(eg == 1L, "Early Granulocytes",
                                            ifelse(ee == 1L, "Early Erythroid",
                                                   ifelse(le == 1L, "Later Erythroid",
                                                          ifelse(mgk== 1L, "Megakaryocytes",
                                                                 ifelse(hsc== 1L, "HSCs",
                                                                        ifelse(ldp== 1L, "Progenitors", "Unknown"))))))))))

  # factor with meaningful level order (precedence on top, then others)
  lvl_order <- c("non-viable","debris/doublet",
                 "Late Monocytes","Late Granulocytes","Early Granulocytes",
                 "Early Erythroid","Later Erythroid","Megakaryocytes","HSCs","Progenitors","Unknown")
  factor(out, levels = lvl_order)
}

## Read FCS and build per-day data
read_day_flowframes <- list()
day_exprs <- list()
row_ranges <- list()
rep_vec_by_day <- list()

for (d in days) {
  files_d <- fcs_by_day[[d]]
  reps <- sapply(files_d, get_rep); ord <- order(reps)
  files_d <- files_d[ord]; reps <- reps[ord]

  read_day_flowframes[[d]] <- lapply(files_d, function(fp) read.FCS(fp, transformation = NULL))
  mats <- lapply(read_day_flowframes[[d]], function(ff) exprs(ff))

  starts <- c(1L, cumsum(sapply(mats, nrow))[-length(mats)] + 1L)
  ends   <- cumsum(sapply(mats, nrow))
  row_ranges[[d]] <- data.frame(rep = reps, start = starts, end = ends)

  day_exprs[[d]]  <- do.call(rbind, mats)

  rv <- integer(nrow(day_exprs[[d]]))
  for (k in seq_len(nrow(row_ranges[[d]]))) {
    rv[row_ranges[[d]]$start[k]:row_ranges[[d]]$end[k]] <- row_ranges[[d]]$rep[k]
  }
  rep_vec_by_day[[d]] <- rv
}

#Read gating CSVs
gating_by_day <- lapply(days, function(d) read_gating_file(csv_by_day[[d]][1]))
names(gating_by_day) <- days

# Align row counts
for (d in days) {
  n_fcs <- nrow(day_exprs[[d]]); n_csv <- nrow(gating_by_day[[d]])
  if (n_fcs != n_csv) {
    warning(sprintf("Row mismatch day %s: FCS=%s CSV=%s; truncating to min", d, n_fcs, n_csv))
    n_min <- min(n_fcs, n_csv)
    day_exprs[[d]]     <- day_exprs[[d]][seq_len(n_min), , drop = FALSE]
    gating_by_day[[d]] <- gating_by_day[[d]][seq_len(n_min), , drop = FALSE]
    rr <- row_ranges[[d]]
    rr$start <- pmin(rr$start, n_min + 1L); rr$end <- pmin(rr$end, n_min)
    row_ranges[[d]] <- rr[rr$start <= rr$end, , drop = FALSE]
    rep_vec_by_day[[d]] <- rep_vec_by_day[[d]][seq_len(n_min)]
  }
}

#Decide channels (desc; Comp-* -A; no Live/Dead
template_ff <- read_day_flowframes[[days[1]]][[1]]
par_tbl <- pData(parameters(template_ff))

is_area   <- grepl("-A$", par_tbl$name)
is_comp   <- grepl("^Comp", par_tbl$name)
desc_txt  <- ifelse(is.na(par_tbl$desc), par_tbl$name, par_tbl$desc)
can_desc  <- .canon(desc_txt)

# Exclude only viability / live-dead; KEEP EdU
is_viab   <- can_desc %in% c("viability","livedead","live","livedeaddye", "VIABILITY", "EdU")

keep_idx  <- which(is_area & is_comp & !is_viab)
stopifnot(length(keep_idx) > 0)

selected_names <- par_tbl$name[keep_idx]                       # raw matrix column names
selected_desc  <- make.names(desc_txt[keep_idx], unique = TRUE) # desc labels as new names
names(selected_desc) <- selected_names

#Combine days with desired channels
dfs_by_day <- lapply(days, function(d) {
  X <- day_exprs[[d]][, selected_names, drop = FALSE]
  X <- as.data.frame(X)
  # rename to DESC names for modeling
  colnames(X) <- unname(selected_desc[colnames(X)])

  X$label <- mk_labels_one_day(gating_by_day[[d]])
  X$.day  <- d
  X$.rep  <- rep_vec_by_day[[d]]
  X$.row  <- seq_len(nrow(X))
  X
})
df_all <- bind_rows(dfs_by_day)

cat("Overall label % BEFORE filtering:\n")
print(round(100 * prop.table(table(df_all$label, useNA = "ifany")), 2))

# Optional filter (keep if meaningful; else revert)
drop_labels <- c("debris/doublet", "non-viable", "Unknown")
df_filt <- df_all[!(df_all$label %in% drop_labels), , drop = FALSE]
if (nrow(df_filt) < 100) df_filt <- df_all
df_filt$label <- droplevels(df_filt$label)

#Leave-one-day-out split
test_df  <- df_filt[df_filt$.day == test_day, , drop = FALSE]
train_df <- df_filt[df_filt$.day != test_day, , drop = FALSE]
stopifnot(nrow(test_df) > 0, nrow(train_df) > 0)

feature_names <- setdiff(colnames(train_df), c("label",".day",".rep",".row"))

#Scaling on training set
mu  <- sapply(train_df[, feature_names, drop = FALSE], mean, na.rm = TRUE)
sds <- sapply(train_df[, feature_names, drop = FALSE], sd,   na.rm = TRUE)
sds[is.na(sds) | sds == 0] <- 1

scale_df <- function(X, mu, sds) {
  X <- as.data.frame(X)
  X[, names(mu)]  <- sweep(X[, names(mu),  drop = FALSE], 2, mu,  "-")
  X[, names(sds)] <- sweep(X[, names(sds), drop = FALSE], 2, sds, "/")
  X
}
train_scaled <- train_df; train_scaled[, feature_names] <- scale_df(train_df[, feature_names, drop = FALSE], mu, sds)
test_scaled  <- test_df;  test_scaled[,  feature_names] <- scale_df(test_df[,  feature_names, drop = FALSE], mu, sds)

#Train & Evaluate
ctrl <- caret::trainControl(
  method = "repeatedcv", number = 3, repeats = 1,
  classProbs = FALSE, summaryFunction = defaultSummary,
  verboseIter = TRUE, allowParallel = TRUE
)

set.seed(123)
fit <- caret::train(
  x = train_scaled[, feature_names, drop = FALSE],
  y = droplevels(train_scaled$label),
  method = "ranger",
  trControl = ctrl,
  tuneLength = 8,
  metric = "Accuracy",
  num.trees = 800,
  importance = "impurity"
)
print(fit$bestTune)

pred_test <- predict(fit, newdata = test_scaled[, feature_names, drop = FALSE])
print(caret::confusionMatrix(pred_test, droplevels(test_scaled$label)))

#Final refit + model bundle
best <- fit$bestTune
tbl <- table(train_scaled$label)
wts <- as.numeric(1 / tbl); names(wts) <- names(tbl)

set.seed(123)
final_fit <- ranger::ranger(
  formula         = label ~ .,
  data            = train_scaled[, c(feature_names, "label")],
  num.trees       = 1000,
  mtry            = if (!is.null(best$mtry)) best$mtry else floor(sqrt(length(feature_names))),
  min.node.size   = if (!is.null(best$min.node.size)) best$min.node.size else 5,
  splitrule       = if (!is.null(best$splitrule)) best$splitrule else "gini",
  probability     = TRUE,
  importance      = "impurity",
  sample.fraction = 0.8,
  class.weights   = wts
)

bundle <- list(
  model    = final_fit,
  features = feature_names,                 # DESC-named features used by the model
  scaling  = list(means = mu, sds = sds),   # z-score params from TRAIN only
  levels   = levels(train_scaled$label),    # class level order
  meta     = list(
    dataset           = "BM-MPS Controls (multiclass)",
    test_day          = test_day,
    model             = "ranger",
    splitrule         = if (!is.null(best$splitrule)) best$splitrule else "gini",
    mtry              = if (!is.null(best$mtry)) best$mtry else NA_integer_,
    min.node.size     = if (!is.null(best$min.node.size)) best$min.node.size else NA_integer_,
    num.trees         = 1000,
    sample.fraction   = 0.8,
    class.weights     = wts,
    zscale            = TRUE,
    features_are_desc = TRUE,               # channels renamed to 'desc' before modeling
    created           = Sys.time()
  )
)

saveRDS(bundle, file.path(outdir, sprintf("bmmps_ranger_bundle_LODO_D%s.rds", test_day)))
cat("Saved model bundle to: ",
    file.path(outdir, sprintf("bmmps_ranger_bundle_LODO_D%s.rds", test_day)), "\n", sep = "")

# Write combined TRAIN/TEST FCS
# Build parameter data with DESC names and append numeric meta channels.

# Build selected parameter set (renamed to DESC)
template_par <- pData(parameters(template_ff))
pdat <- template_par[keep_idx, , drop = FALSE]
pdat$name <- unname(selected_desc)                 # channel names = DESC
rownames(pdat) <- paste0("$P", seq_len(nrow(pdat)))
par_selected <- AnnotatedDataFrame(pdat)           # requires Biobase

make_combined_ff <- function(df, par_selected) {
  feat_cols <- pData(par_selected)$name
  stopifnot(all(feat_cols %in% colnames(df)))

  classes <- levels(df$label)
  label_code <- as.integer(df$label)
  tp_num  <- as.numeric(df$.day)
  rep_num <- as.numeric(df$.rep)

  expr <- as.matrix(df[, feat_cols, drop = FALSE])
  expr <- cbind(expr,
                label     = as.numeric(label_code),
                timepoint = as.numeric(tp_num),
                replicate = as.numeric(rep_num))

  add_row <- function(name, desc, vals) {
    rng <- range(vals, na.rm = TRUE)
    data.frame(name = name, desc = desc,
               range   = max(1, ceiling(rng[2] - rng[1])),
               minRange= rng[1],
               maxRange= rng[2],
               stringsAsFactors = FALSE)
  }
  p_extra <- rbind(
    add_row("label",     "class label code", expr[, "label"]),
    add_row("timepoint", "day number",       expr[, "timepoint"]),
    add_row("replicate", "replicate id",     expr[, "replicate"])
  )
  par2 <- AnnotatedDataFrame(rbind(pData(par_selected), p_extra))

  desc_kw <- list()
  lab_map <- paste(sprintf("%d=%s", seq_along(classes), classes), collapse = ";")
  desc_kw[["label_map"]] <- lab_map
  desc_kw[["$TOT"]] <- nrow(expr)

  flowFrame(exprs = expr, parameters = par2, description = desc_kw)
}

ord_cols <- pData(par_selected)$name
train_common <- train_df[, c(ord_cols, "label",".day",".rep",".row")]
test_common  <- test_df[,  c(ord_cols, "label",".day",".rep",".row")]

ff_train <- make_combined_ff(train_common, par_selected)
ff_test  <- make_combined_ff(test_common,  par_selected)

write.FCS(ff_train, file.path(outdir, sprintf("BM_MPS_TRAIN_Dnot%s.fcs", test_day)))
write.FCS(ff_test,  file.path(outdir, sprintf("BM_MPS_TEST_D%s.fcs",    test_day)))

cat("Wrote combined FCS (channels = DESC names + meta):\n",
    file.path(outdir, sprintf("BM_MPS_TRAIN_Dnot%s.fcs", test_day)), "\n",
    file.path(outdir, sprintf("BM_MPS_TEST_D%s.fcs",    test_day)), "\n", sep = "")

# Optional sanity summaries
for (d in days) {
  cat(sprintf("Day %s label distribution (%%):\n", d))
  print(round(100 * prop.table(table(mk_labels_one_day(gating_by_day[[d]]))), 2))
}
cat("Overall (after filtering) label distribution (%%):\n")
print(round(100 * prop.table(table(df_filt$label)), 2))

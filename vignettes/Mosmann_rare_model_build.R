## Setup
# Packages
library(flowCore)
library(dplyr)
library(caret)
library(randomForest)
library(ranger)
library(pROC)

set.seed(123)

# Load + prep data
# Paths
dir    <- "~/OneDrive - Swansea University/Documents/AutoFlow/data/benchmarking/"
outdir <- "~/OneDrive - Swansea University/Documents/AutoFlow/outputs"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load FCS (no transformation)
mosmann <- read.FCS(file.path(dir, "Mosmann_rare.fcs"), transformation = NULL)
X_all   <- exprs(mosmann)                     # numeric matrix [events x channels]
all_cols <- colnames(X_all)

# Build data.frame for modeling
df <- as.data.frame(X_all)

# Keep label, drop purely technical channels (edit to taste)
drop_cols <- c("FSC-A","FSC-W","FSC-H","SSC-A","SSC-W","SSC-H","Live_Dead","Time")
keep <- setdiff(colnames(df), drop_cols)
df <- df[, keep, drop = FALSE]

# Ensure 'label' exists and is factor
stopifnot("label" %in% colnames(df))
df$label <- as.factor(df$label)

# Preserve original names (for mapping if needed later)
orig_names <- colnames(df)

# Make syntactically valid names for modeling APIs
colnames(df) <- make.names(colnames(df), unique = TRUE)

# Stratified split (no leakage)
train_idx <- caret::createDataPartition(df$label, p = 0.70, list = FALSE)
train <- df[train_idx, , drop = FALSE]
test  <- df[-train_idx, , drop = FALSE]

# Upsample *training set only* (if imbalanced)
train_up <- caret::upSample(
  x = subset(train, select = -label),
  y = train$label,
  yname = "label"
)

# Scale using TRAIN stats
num_cols <- setdiff(colnames(train_up), "label")

mu  <- sapply(train_up[, num_cols, drop = FALSE], mean, na.rm = TRUE)
sds <- sapply(train_up[, num_cols, drop = FALSE], sd,   na.rm = TRUE)
sds[is.na(sds) | sds == 0] <- 1  # guard

scale_df <- function(X, mu, sds) {
  X <- as.data.frame(X)
  X[, names(mu)]  <- sweep(X[, names(mu),  drop = FALSE], 2, mu,  "-")
  X[, names(sds)] <- sweep(X[, names(sds), drop = FALSE], 2, sds, "/")
  X
}

train_scaled <- train_up
train_scaled[, num_cols] <- scale_df(train_up[, num_cols, drop = FALSE], mu, sds)

test_scaled <- test
test_scaled[, num_cols] <- scale_df(test[, num_cols, drop = FALSE], mu, sds)

# Label handling (pos = X1 first)
# If original labels are "1" / "0", this puts positive first.
train_scaled$label <- factor(make.names(train_scaled$label), levels = c("X1","X0"))
test_scaled$label  <- factor(make.names(test_scaled$label),  levels = c("X1","X0"))

# caret tuning (ranger)
feature_names <- setdiff(colnames(train_scaled), "label")
p <- length(feature_names)

# Class weights (handle imbalance)
tbl <- table(train_scaled$label)                  # X1 (pos), X0 (neg)
wts <- as.numeric(1 / tbl); names(wts) <- names(tbl)

ctrl <- caret::trainControl(
  method = "repeatedcv",
  number = 3, repeats = 1,
  classProbs = TRUE,
  summaryFunction = caret::twoClassSummary,
  sampling = "up",         # upsample within CV folds too
  verboseIter = TRUE,
  allowParallel = TRUE
)

set.seed(123)
fit_tuned <- caret::train(
  x = train_scaled[, feature_names, drop = FALSE],
  y = train_scaled$label,
  method = "ranger",
  trControl = ctrl,
  tuneLength = 6,
  metric = "Spec",    # pick your target metric (Spec/Sens/ROC)
  num.trees = 1000,
  importance = "impurity",
  class.weights = wts
)

print(fit_tuned$bestTune)
best <- fit_tuned$bestTune

## Final fit on full TRAIN (scaled)
set.seed(123)
final_fit <- ranger::ranger(
  formula = label ~ .,
  data = train_scaled[, c(feature_names, "label")],
  num.trees = 1000,
  mtry = best$mtry,
  min.node.size = best$min.node.size,
  splitrule = best$splitrule,
  class.weights = wts,
  importance = "impurity",
  probability = TRUE,
  sample.fraction = 0.8,
  num.threads = max(1, parallel::detectCores() - 1)
)

## Evaluation
pr <- predict(final_fit, data = test_scaled)$predictions   # matrix: cols X1/X0
pred_class <- factor(colnames(pr)[max.col(pr, ties.method = "first")],
                     levels = levels(test_scaled$label))

cm <- caret::confusionMatrix(pred_class, test_scaled$label)
print(cm)

# ROC AUC (positive = X1)
roc_obj <- pROC::roc(
  response  = test_scaled$label,
  predictor = pr[, "X1"],
  levels    = rev(levels(test_scaled$label))
)
cat("ROC AUC:", pROC::auc(roc_obj), "\n")

# Save model bundle
bundle <- list(
  model    = final_fit,
  features = feature_names,
  scaling  = list(means = mu, sds = sds),     # z-scaling params from TRAIN
  levels   = levels(train_scaled$label),      # c("X1","X0")
  meta     = list(
    source          = "Mosmann_rare.fcs",
    model           = "ranger",
    splitrule       = best$splitrule,
    mtry            = best$mtry,
    min.node.size   = best$min.node.size,
    num.trees       = 1000,
    sample.fraction = 0.8,
    class.weights   = wts,
    zscale          = TRUE,
    created         = Sys.time()
  )
)
saveRDS(bundle, file.path(outdir, "autoflow_mosmann_ranger_bundle.rds"))

# Save raw TRAIN/TEST subsets as FCS files
# We’ll write unscaled data with ALL original channels preserved.
# Use the train/test row indices created earlier (train_idx).

# Helperto build a flowFrame with inherited parameters/keywords and updated $TOT
make_subset_ff <- function(ff, rows) {
  mat <- exprs(ff)[rows, , drop = FALSE]
  par <- parameters(ff)       # carry original parameter metadata
  desc <- description(ff)     # carry original keywords
  desc[["$TOT"]] <- nrow(mat) # update event count
  flowFrame(exprs = mat, parameters = par, description = desc)
}

ff_train <- make_subset_ff(mosmann, train_idx)
ff_test  <- make_subset_ff(mosmann, setdiff(seq_len(nrow(X_all)), train_idx))

write.FCS(ff_train, file.path(outdir, "Mosmann_rare_TRAIN.fcs"))
write.FCS(ff_test,  file.path(outdir, "Mosmann_rare_TEST.fcs"))

cat("Wrote:\n",
    file.path(outdir, "Mosmann_rare_TRAIN.fcs"), "\n",
    file.path(outdir, "Mosmann_rare_TEST.fcs"),  "\n")

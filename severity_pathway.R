# RStudio工作目录设为本文件夹
# 先运行ev-pathway-score.R
library(data.table)
library(glmnet)


dm <- "inputs"
out_dir <- "generated/clinical_context"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_repeats <- 10L
outer_k <- 5L
inner_k <- 4L
base_seed <- 20260817L
alpha_grid <- c(0, 0.25, 0.5, 0.75, 1)

meta <- fread(file.path(
  dm, "refine-logs", "clinical_rna_meta_20260815",
  "analysis_ready_clinical_rna_meta.tsv"
))
score <- fread(cmd = paste(
  "gzip -cd",
  shQuote(file.path("generated", "pathway_scores", "Hallmark_single_sample_rank_scores_wide.tsv.gz"))
))

hallmark <- grep("^HALLMARK_", names(score), value = TRUE)
keep_meta <- c(
  "sample_id", "sample_library_id", "group", "age", "sex",
  "disease_duration_years", "FVC_pct", "MMT8", "steroid_daily_mg",
  "DLCO_pct", "DLCO_severity", "DLCO_severity_code"
)
dat <- merge(
  meta[, ..keep_meta],
  score[, c(
    "sample_id", "library_id", "batch", "log10_library_size_z", "detected_genes_z",
    hallmark
  ), with = FALSE],
  by = "sample_id",
  all = FALSE
)
dat <- dat[group %in% c("MDA5", "ARS") & !is.na(DLCO_pct)]
stopifnot(nrow(dat) == 87L, !anyDuplicated(dat$sample_id))
stopifnot(all(dat$sample_library_id == dat$library_id))
stopifnot(!anyNA(dat[, ..hallmark]))
dat[, outcome := as.integer(DLCO_severity_code >= 2)]
dat[, group_mda5 := as.integer(group == "MDA5")]
dat[, sex_male := as.integer(sex %in% c("男", "male", "Male", "M"))]
dat[, batch := factor(batch)]

clinical_vars <- c(
  "age", "sex_male", "disease_duration_years", "FVC_pct", "MMT8", "steroid_daily_mg"
)

make_stratified_folds <- function(group, y, k, seed) {
  set.seed(seed)
  strata <- interaction(group, y, drop = TRUE)
  fold <- integer(length(y))
  for (lev in levels(strata)) {
    idx <- sample(which(strata == lev))
    fold[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold
}

clinical_matrix <- function(train, test, vars) {
  x_train <- matrix(nrow = nrow(train), ncol = 0)
  x_test <- matrix(nrow = nrow(test), ncol = 0)
  for (v in vars) {
    tr <- as.numeric(train[[v]])
    te <- as.numeric(test[[v]])
    med <- median(tr, na.rm = TRUE)
    if (!is.finite(med)) med <- 0
    miss_tr <- as.integer(is.na(tr))
    miss_te <- as.integer(is.na(te))
    tr[is.na(tr)] <- med
    te[is.na(te)] <- med
    x_train <- cbind(x_train, tr)
    x_test <- cbind(x_test, te)
    colnames(x_train)[ncol(x_train)] <- v
    colnames(x_test)[ncol(x_test)] <- v
    if (any(miss_tr) || any(miss_te)) {
      x_train <- cbind(x_train, miss_tr)
      x_test <- cbind(x_test, miss_te)
      colnames(x_train)[ncol(x_train)] <- paste0(v, "_missing")
      colnames(x_test)[ncol(x_test)] <- paste0(v, "_missing")
    }
  }
  list(train = x_train, test = x_test)
}

residualize_hallmark <- function(train, test, pathways) {
  batch_levels <- levels(dat$batch)
  train$batch <- factor(train$batch, levels = batch_levels)
  test$batch <- factor(test$batch, levels = batch_levels)
  z_train <- model.matrix(~ batch + log10_library_size_z + detected_genes_z, data = train)
  z_test <- model.matrix(~ batch + log10_library_size_z + detected_genes_z, data = test)
  x_train <- as.matrix(train[, ..pathways])
  x_test <- as.matrix(test[, ..pathways])
  coef <- qr.solve(crossprod(z_train) + diag(1e-8, ncol(z_train)), crossprod(z_train, x_train))
  list(
    train = x_train - z_train %*% coef,
    test = x_test - z_test %*% coef
  )
}

auc_rank <- function(y, p) {
  pos <- p[y == 1]
  neg <- p[y == 0]
  if (!length(pos) || !length(neg)) return(NA_real_)
  z <- outer(pos, neg, "-")
  mean(z > 0) + 0.5 * mean(z == 0)
}

lambda_score <- function(fit, s) {
  target <- if (s == "lambda.1se") fit$lambda.1se else fit$lambda.min
  idx <- which.min(abs(log(fit$lambda) - log(target)))
  fit$cvm[idx]
}

fit_glmnet <- function(x_train, y_train, x_test, inner_fold) {
  weights <- ifelse(
    y_train == 1,
    length(y_train) / (2 * sum(y_train == 1)),
    length(y_train) / (2 * sum(y_train == 0))
  )
  fits <- vector("list", length(alpha_grid))
  score_min <- score_1se <- rep(Inf, length(alpha_grid))
  for (i in seq_along(alpha_grid)) {
    fit <- try(cv.glmnet(
      x_train, y_train,
      family = "binomial",
      alpha = alpha_grid[i],
      foldid = inner_fold,
      weights = weights,
      type.measure = "deviance",
      standardize = TRUE,
      nlambda = 80
    ), silent = TRUE)
    if (inherits(fit, "try-error")) next
    fits[[i]] <- fit
    score_min[i] <- lambda_score(fit, "lambda.min")
    score_1se[i] <- lambda_score(fit, "lambda.1se")
  }
  if (!any(is.finite(score_min))) stop("All inner glmnet fits failed")
  i_min <- which.min(score_min)
  i_1se <- which.min(score_1se)
  fit_min <- fits[[i_min]]
  fit_1se <- fits[[i_1se]]
  beta_min <- as.matrix(coef(fit_min, s = "lambda.min"))
  beta_1se <- as.matrix(coef(fit_1se, s = "lambda.1se"))
  selected <- function(beta, alpha) {
    if (alpha == 0) return(character())
    setdiff(rownames(beta)[beta[, 1] != 0], "(Intercept)")
  }
  list(
    probability_min = as.numeric(predict(fit_min, x_test, s = "lambda.min", type = "response")),
    probability_1se = as.numeric(predict(fit_1se, x_test, s = "lambda.1se", type = "response")),
    alpha_min = alpha_grid[i_min],
    alpha_1se = alpha_grid[i_1se],
    lambda_min = fit_min$lambda.min,
    lambda_1se = fit_1se$lambda.1se,
    selected_min = selected(beta_min, alpha_grid[i_min]),
    selected_1se = selected(beta_1se, alpha_grid[i_1se])
  )
}

fit_autoantibody <- function(x_train, y_train, x_test) {
  d <- data.frame(y = y_train, group_MDA5 = as.numeric(x_train[, 1]))
  fit <- glm(y ~ group_MDA5, data = d, family = binomial())
  p <- as.numeric(predict(
    fit,
    newdata = data.frame(group_MDA5 = as.numeric(x_test[, 1])),
    type = "response"
  ))
  list(
    probability_min = p,
    probability_1se = p,
    alpha_min = NA_real_,
    alpha_1se = NA_real_,
    lambda_min = NA_real_,
    lambda_1se = NA_real_,
    selected_min = "group_MDA5",
    selected_1se = "group_MDA5"
  )
}

model_matrix <- function(model, train, test) {
  if (model == "autoantibody_only") {
    return(list(
      train = matrix(train$group_mda5, ncol = 1, dimnames = list(NULL, "group_MDA5")),
      test = matrix(test$group_mda5, ncol = 1, dimnames = list(NULL, "group_MDA5"))
    ))
  }
  clinical <- clinical_matrix(train, test, clinical_vars)
  rna <- residualize_hallmark(train, test, hallmark)
  colnames(rna$train) <- hallmark
  colnames(rna$test) <- hallmark
  if (model == "clinical_only") return(clinical)
  if (model == "EV_RNA_pathway_only") return(rna)
  if (model == "clinical_autoantibody") {
    return(list(
      train = cbind(group_MDA5 = train$group_mda5, clinical$train),
      test = cbind(group_MDA5 = test$group_mda5, clinical$test)
    ))
  }
  if (model == "combined") {
    return(list(
      train = cbind(group_MDA5 = train$group_mda5, clinical$train, rna$train),
      test = cbind(group_MDA5 = test$group_mda5, clinical$test, rna$test)
    ))
  }
  stop("Unknown model")
}

models <- c(
  "autoantibody_only", "clinical_only", "clinical_autoantibody",
  "EV_RNA_pathway_only", "combined"
)
pred_rows <- list()
fold_rows <- list()
selected_min <- list()
selected_1se <- list()
ii <- 0L
inner_allocations <- list()

for (r in seq_len(n_repeats)) {
  outer_fold <- make_stratified_folds(dat$group, dat$outcome, outer_k, base_seed + r * 1000L)
  for (f in seq_len(outer_k)) {
    te <- which(outer_fold == f)
    tr <- setdiff(seq_len(nrow(dat)), te)
    inner_fold <- make_stratified_folds(
      dat$group[tr], dat$outcome[tr], inner_k,
      base_seed + r * 1000L + f * 10L
    )
    inner_allocations[[length(inner_allocations)+1L]] <- data.table(
      repeat_id=r, outer_fold=f, sample_id=dat$sample_id[tr], inner_fold=inner_fold)
    for (model in models) {
      ii <- ii + 1L
      x <- model_matrix(model, dat[tr], dat[te])
      fit <- if (model == "autoantibody_only") {
        fit_autoantibody(x$train, dat$outcome[tr], x$test)
      } else {
        fit_glmnet(x$train, dat$outcome[tr], x$test, inner_fold)
      }
      pred_rows[[ii]] <- data.table(
        repeat_id = r,
        fold_id = f,
        model = model,
        sample_id = dat$sample_id[te],
        library_id = dat$library_id[te],
        group = dat$group[te],
        truth = dat$outcome[te],
        probability_min = fit$probability_min,
        probability_1se = fit$probability_1se
      )
      fold_rows[[ii]] <- data.table(
        repeat_id = r,
        fold_id = f,
        model = model,
        train_n = length(tr),
        test_n = length(te),
        input_features = ncol(x$train),
        alpha_min = fit$alpha_min,
        alpha_1se = fit$alpha_1se,
        lambda_min = fit$lambda_min,
        lambda_1se = fit$lambda_1se,
        selected_min_n = length(fit$selected_min),
        selected_1se_n = length(fit$selected_1se),
        train_sample_ids = paste(dat$sample_id[tr], collapse = ";"),
        test_sample_ids = paste(dat$sample_id[te], collapse = ";")
      )
      selected_min[[ii]] <- data.table(model = model, feature = fit$selected_min)
      selected_1se[[ii]] <- data.table(model = model, feature = fit$selected_1se)
    }
  }
}

pred <- rbindlist(pred_rows)
fwrite(rbindlist(inner_allocations), file.path(out_dir,"inner_fold_assignments.tsv"), sep="\t")
folds <- rbindlist(fold_rows)
stopifnot(nrow(pred) == nrow(dat) * n_repeats * length(models))

long_pred <- melt(
  pred,
  id.vars = c("repeat_id", "fold_id", "model", "sample_id", "library_id", "group", "truth"),
  measure.vars = c("probability_min", "probability_1se"),
  variable.name = "penalty_rule",
  value.name = "probability"
)
long_pred[, penalty_rule := sub("probability_", "lambda_", penalty_rule)]

repeat_auc <- long_pred[, .(auc = auc_rank(truth, probability)), by = .(model, penalty_rule, repeat_id)]
participant <- long_pred[, .(
  library_id = unique(library_id),
  group = unique(group),
  truth = unique(truth),
  mean_oof_score = mean(probability),
  sd_oof_score = sd(probability),
  prediction_repeats = .N
), by = .(model, penalty_rule, sample_id)]

performance <- repeat_auc[, .(
  primary_mean_repeated_cv_auc = mean(auc),
  repeat_auc_sd = sd(auc),
  repeat_auc_min = min(auc),
  repeat_auc_max = max(auc)
), by = .(model, penalty_rule)]
ensemble_auc <- participant[, .(
  secondary_participant_mean_oof_auc = auc_rank(truth, mean_oof_score)
), by = .(model, penalty_rule)]
performance <- merge(performance, ensemble_auc, by = c("model", "penalty_rule"))

confusion <- participant[, {
  pred_class <- as.integer(mean_oof_score >= 0.5)
  .(
    tp = sum(pred_class == 1 & truth == 1),
    tn = sum(pred_class == 0 & truth == 0),
    fp = sum(pred_class == 1 & truth == 0),
    fn = sum(pred_class == 0 & truth == 1),
    sensitivity = sum(pred_class == 1 & truth == 1) / sum(truth == 1),
    specificity = sum(pred_class == 0 & truth == 0) / sum(truth == 0),
    accuracy = mean(pred_class == truth),
    brier = mean((mean_oof_score - truth)^2)
  )
}, by = .(model, penalty_rule)]

stability_table <- function(x, rule) {
  z <- rbindlist(x, fill = TRUE)
  if (!nrow(z)) return(data.table())
  z[, .(selected_outer_models = .N), by = .(model, feature)][
    , `:=`(penalty_rule = rule, total_outer_models = n_repeats * outer_k)
  ][, selection_frequency := selected_outer_models / total_outer_models]
}
stability <- rbind(
  stability_table(selected_min, "lambda_min"),
  stability_table(selected_1se, "lambda_1se"),
  fill = TRUE
)

manifest <- dat[, .(
  sample_id, library_id, group, outcome, DLCO_pct, DLCO_severity,
  age, sex, disease_duration_years, FVC_pct, MMT8, steroid_daily_mg,
  batch, log10_library_size_z, detected_genes_z
)]
missingness <- data.table(
  variable = c(clinical_vars, hallmark),
  nonmissing_n = sapply(c(clinical_vars, hallmark), function(v) sum(!is.na(dat[[v]]))),
  missing_n = sapply(c(clinical_vars, hallmark), function(v) sum(is.na(dat[[v]])))
)

fwrite(manifest, file.path(out_dir, "analysis_manifest_87.tsv"), sep = "\t")
fwrite(missingness, file.path(out_dir, "input_missingness.tsv"), sep = "\t")
fwrite(long_pred, file.path(out_dir, "all_outer_fold_predictions.tsv"), sep = "\t")
fwrite(folds, file.path(out_dir, "outer_fold_details.tsv"), sep = "\t")
fwrite(repeat_auc, file.path(out_dir, "auc_by_repeat.tsv"), sep = "\t")
fwrite(participant, file.path(out_dir, "participant_mean_oof_scores.tsv"), sep = "\t")
fwrite(performance, file.path(out_dir, "performance_summary.tsv"), sep = "\t")
fwrite(confusion, file.path(out_dir, "fixed_0.5_confusion_summary.tsv"), sep = "\t")
fwrite(stability, file.path(out_dir, "feature_stability.tsv"), sep = "\t")

writeLines(c(
  "# Pooled binary DLCO-severity nested-CV analysis",
  "",
  "- Complete-case outcome cohort: 87 patients; 42 MDA5 and 45 ARS.",
  "- Outcome: mild versus moderate/severe DLCO impairment.",
  "- Outer validation: 10 repeats of stratified 5-fold CV, stratified jointly by subtype and outcome.",
  "- Inner tuning: stratified 4-fold CV in each outer-training set.",
  "- Penalized logistic candidates: ridge, elastic net, and lasso (alpha 0, 0.25, 0.5, 0.75, 1).",
  "- lambda.min is the performance-oriented rule; lambda.1se is the more strongly regularized sensitivity rule.",
  "- All alpha/lambda selection uses only outer-training data.",
  "- Clinical missing values are median-imputed using the outer-training set and accompanied by missingness indicators.",
  "- Hallmark scores are residualized for sequencing-date batch, log10 library size, and detected-gene count using coefficients estimated in the outer-training set.",
  "- Clinical-only inputs: age, sex, disease duration, concurrent FVC%, MMT8, and daily steroid dose.",
  "- Clinical-plus-autoantibody inputs: the clinical variables plus anti-MDA5/anti-ARS subtype.",
  "- Combined inputs: autoantibody subtype, clinical variables, and 50 QC/batch-residualized Hallmark scores.",
  "- Because FVC% and DLCO were measured concurrently, this is cross-sectional severity classification, not prospective prediction.",
  "- This is internal validation from one cohort; it is not external validation."
), file.path(out_dir, "METHODS_AND_LIMITATIONS.md"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
saveRDS(list(
  manifest = manifest,
  predictions = long_pred,
  folds = folds,
  performance = performance,
  confusion = confusion,
  stability = stability
), file.path(out_dir, "compact_results.rds"))

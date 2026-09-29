# RStudio工作目录设为本文件夹
library(data.table)
library(edgeR)
library(e1071)


options(stringsAsFactors = FALSE)
set.seed(20260826)

dm <- "inputs"
root <- "generated/transcript"
quick_test <- FALSE
input_mode <- "raw_count"
run_dir <- root
result_dir <- file.path(run_dir, "results")
log_dir <- file.path(run_dir, "logs")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

n_repeat <- 10L
n_outer <- 5L
n_inner <- 4L
n_boot <- 2000L
n_cores <- 4L
n_grid <- 1:8
cost_grid <- c(0.1, 1, 10, 100)
gamma_mult_grid <- c(0.5, 1, 2)

anno <- read.delim(file.path(dm, "annotation.txt"), check.names = FALSE)
colnames(anno)[2] <- "library_id"
anno <- anno[anno$mistake == 0, , drop = FALSE]
anno$batch <- substr(anno$library_id, 1, 6)
stopifnot(nrow(anno) == 155L)
stopifnot(!anyDuplicated(anno$sample_id), !anyDuplicated(anno$library_id))

if (input_mode == "combat_global") {
  input_path <- file.path(dm, "IM EV", "gencode_rmbatch.all.txt")
  raw <- fread(input_path, data.table = FALSE, check.names = FALSE)
  feature_id <- raw[[1]]
  raw <- raw[, -1, drop = FALSE]
  rownames(raw) <- feature_id
} else if (input_mode == "raw_count") {
  input_path <- file.path(dm, "1204", "gencode.txt")
  raw <- fread(input_path, data.table = FALSE, check.names = FALSE)
  feature_id <- raw[[1]]
  raw <- raw[, -1, drop = FALSE]
  raw <- raw[, -c(40:44, 50:54), drop = FALSE]
  colnames(raw)[35:44] <- sub("-GP", "", colnames(raw)[35:44], fixed = TRUE)
  rownames(raw) <- feature_id
} else {
  stop("DM_INPUT_MODE must be combat_global or raw_count")
}

stopifnot(setequal(colnames(raw), anno$library_id))
raw <- raw[, anno$library_id, drop = FALSE]

part <- tstrsplit(rownames(raw), "|", fixed = TRUE, fill = "")
feature <- data.frame(
  feature_id = rownames(raw),
  gene_symbol = part[[3]],
  biotype = part[[4]],
  stringsAsFactors = FALSE
)

rna_type <- c(
  "protein_coding", "lncRNA", "Mt_tRNA", "scaRNA", "sRNA",
  "IG_C_gene", "IG_D_gene", "IG_V_gene", "IG_J_gene",
  "TR_C_gene", "TR_D_gene", "TR_V_gene", "TR_J_gene",
  "miRNA", "snRNA", "snoRNA"
)
keep <- feature$biotype %in% rna_type & feature$gene_symbol != "Y_RNA"
count_all <- as.matrix(raw[keep, , drop = FALSE])
storage.mode(count_all) <- "double"
feature <- feature[keep, , drop = FALSE]
rownames(feature) <- feature$feature_id
rm(raw)
gc()

auc_value <- function(y, score) {
  ok <- is.finite(score) & !is.na(y)
  y <- y[ok]
  score <- score[ok]
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

average_precision <- function(y, score) {
  ok <- is.finite(score) & !is.na(y)
  y <- y[ok]
  score <- score[ok]
  if (!sum(y == 1)) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  yy <- y[ord]
  precision <- cumsum(yy == 1) / seq_along(yy)
  mean(precision[yy == 1])
}

make_foldid <- function(strata, k, seed) {
  set.seed(seed)
  fold <- integer(length(strata))
  for (s in unique(strata)) {
    idx <- sample(which(strata == s))
    fold[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold
}

log_cpm <- function(count) {
  lib <- colSums(count)
  stopifnot(all(lib > 0))
  log2(sweep(count, 2, lib / 1e6, "/") + 1)
}

filter_train_features <- function(count, y) {
  yy <- factor(y, levels = c(0, 1))
  dge <- DGEList(counts = count, group = yy)
  keep <- filterByExpr(dge, group = yy, min.count = 2, min.prop = 0.2)
  rownames(count)[keep]
}

rank_train_features <- function(x, y) {
  x1 <- x[, y == 1, drop = FALSE]
  x0 <- x[, y == 0, drop = FALSE]
  d <- rowMeans(x1) - rowMeans(x0)
  se <- sqrt(apply(x1, 1, var) / ncol(x1) + apply(x0, 1, var) / ncol(x0))
  z <- abs(d) / (se + 1e-8)
  z[!is.finite(z)] <- 0
  names(sort(z, decreasing = TRUE))
}

scale_train_test <- function(x_train, x_test) {
  center <- colMeans(x_train)
  scale <- apply(x_train, 2, sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  list(
    train = sweep(sweep(x_train, 2, center, "-"), 2, scale, "/"),
    test = sweep(sweep(x_test, 2, center, "-"), 2, scale, "/")
  )
}

class_weight <- function(y) {
  n <- length(y)
  c("0" = n / (2 * sum(y == 0)), "1" = n / (2 * sum(y == 1)))
}

svm_probability <- function(x_train, y_train, x_test, kernel, cost, gamma, seed) {
  set.seed(seed)
  fit <- svm(
    x = x_train,
    y = factor(y_train, levels = c(0, 1)),
    type = "C-classification",
    kernel = kernel,
    cost = cost,
    gamma = gamma,
    class.weights = class_weight(y_train),
    probability = TRUE,
    scale = FALSE
  )
  pr <- predict(fit, x_test, probability = TRUE)
  prob <- attr(pr, "probabilities")
  if (is.null(prob) || !"1" %in% colnames(prob)) {
    stop("SVM probability output is incomplete")
  }
  as.numeric(prob[, "1"])
}

tune_svm <- function(count_train, meta_train, seed) {
  inner_fold <- make_foldid(
    meta_train$fold_stratum,
    n_inner, seed
  )
  fold_output <- file.path("generated", "transcript", "folds")
  dir.create(fold_output,recursive=TRUE,showWarnings=FALSE)
  fwrite(data.table(sample_id=meta_train$sample_id, seed=seed, inner_fold=inner_fold),
    file.path(fold_output,paste0("inner_",seed,".tsv")),sep="\t")
  radial_config <- expand.grid(
    n_feature = n_grid,
    kernel = "radial",
    cost = cost_grid,
    gamma_mult = gamma_mult_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  linear_config <- expand.grid(
    n_feature = n_grid,
    kernel = "linear",
    cost = cost_grid,
    gamma_mult = 1,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  config <- rbind(radial_config, linear_config)
  pred <- matrix(NA_real_, nrow(meta_train), nrow(config))

  for (fold in seq_len(n_inner)) {
    ii_train <- which(inner_fold != fold)
    ii_test <- which(inner_fold == fold)
    y_train <- meta_train$y[ii_train]
    eligible <- filter_train_features(count_train[, ii_train, drop = FALSE], y_train)
    x_train_all <- log_cpm(count_train[eligible, ii_train, drop = FALSE])
    x_test_all <- log_cpm(count_train[eligible, ii_test, drop = FALSE])
    ranked <- rank_train_features(x_train_all, y_train)

    for (j in seq_len(nrow(config))) {
      n_use <- min(config$n_feature[j], length(ranked))
      selected <- ranked[seq_len(n_use)]
      xtr <- t(x_train_all[selected, , drop = FALSE])
      xte <- t(x_test_all[selected, , drop = FALSE])
      ss <- scale_train_test(xtr, xte)
      gamma <- config$gamma_mult[j] / ncol(ss$train)
      ans <- try(
        svm_probability(
          ss$train, y_train, ss$test,
          config$kernel[j], config$cost[j], gamma,
          seed + fold * 10000 + j
        ),
        silent = TRUE
      )
      if (!inherits(ans, "try-error")) pred[ii_test, j] <- ans
    }
  }

  score <- vapply(seq_len(nrow(config)), function(j) {
    auc_value(meta_train$y, pred[, j])
  }, numeric(1))
  brier <- vapply(seq_len(nrow(config)), function(j) {
    mean((pred[, j] - meta_train$y)^2, na.rm = TRUE)
  }, numeric(1))
  config$inner_auc <- score
  config$inner_brier <- brier
  config <- config[order(-config$inner_auc, config$inner_brier,
                         config$n_feature, config$kernel,
                         config$cost, config$gamma_mult), ]
  config[1, , drop = FALSE]
}

fit_outer <- function(count, meta, train_idx, test_idx, seed) {
  count_train <- count[, train_idx, drop = FALSE]
  count_test <- count[, test_idx, drop = FALSE]
  meta_train <- meta[train_idx, , drop = FALSE]
  best <- tune_svm(count_train, meta_train, seed)

  eligible <- filter_train_features(count_train, meta_train$y)
  x_train_all <- log_cpm(count_train[eligible, , drop = FALSE])
  x_test_all <- log_cpm(count_test[eligible, , drop = FALSE])
  ranked <- rank_train_features(x_train_all, meta_train$y)
  n_use <- min(best$n_feature, length(ranked))
  selected <- ranked[seq_len(n_use)]
  xtr <- t(x_train_all[selected, , drop = FALSE])
  xte <- t(x_test_all[selected, , drop = FALSE])
  ss <- scale_train_test(xtr, xte)
  gamma <- best$gamma_mult / ncol(ss$train)
  pred <- svm_probability(
    ss$train, meta_train$y, ss$test,
    best$kernel, best$cost, gamma, seed + 900000
  )

  list(pred = pred, selected = selected, best = best, gamma = gamma,
       eligible_n = length(eligible))
}

bootstrap_auc <- function(y, score, B, seed) {
  set.seed(seed)
  pos <- which(y == 1)
  neg <- which(y == 0)
  val <- replicate(B, {
    idx <- c(sample(pos, length(pos), replace = TRUE),
             sample(neg, length(neg), replace = TRUE))
    auc_value(y[idx], score[idx])
  })
  unname(quantile(val, c(0.025, 0.975), na.rm = TRUE))
}

summarize_patient_prediction <- function(pred, comparison, validation, seed) {
  agg <- aggregate(
    score ~ sample_id + library_id + group + batch + truth,
    pred, mean
  )
  auc <- auc_value(agg$truth, agg$score)
  ci <- bootstrap_auc(agg$truth, agg$score, n_boot, seed)
  cls <- as.integer(agg$score >= 0.5)
  perf <- data.frame(
    comparison = comparison,
    validation = validation,
    n = nrow(agg),
    AUC = auc,
    AUC_low = ci[1],
    AUC_high = ci[2],
    PR_AUC = average_precision(agg$truth, agg$score),
    Brier = mean((agg$score - agg$truth)^2),
    accuracy_0.5 = mean(cls == agg$truth),
    sensitivity_0.5 = mean(cls[agg$truth == 1] == 1),
    specificity_0.5 = mean(cls[agg$truth == 0] == 0),
    stringsAsFactors = FALSE
  )
  agg$comparison <- comparison
  agg$validation <- validation
  list(performance = perf, patient = agg)
}

run_comparison <- function(comparison) {
  if (comparison == "Pooled_severity") {
    severity_path <- file.path(
      dm, "pooled_binary_severity_20260817", "output", "analysis_manifest_87.tsv"
    )
    severity <- fread(severity_path, data.table = FALSE)
    meta <- severity[, c("sample_id", "library_id", "group", "batch", "outcome")]
    meta$y <- as.integer(meta$outcome)
    meta$fold_stratum <- paste(meta$group, meta$y, sep = "_")
    stopifnot(nrow(meta) == 87L, !anyDuplicated(meta$sample_id))
  } else {
    group_map <- list(
      MDA5_vs_HC = c("MDA5", "HC"),
      ARS_vs_HC = c("ARS", "HC"),
      MDA5_vs_ARS = c("MDA5", "ARS")
    )
    groups <- group_map[[comparison]]
    if (is.null(groups)) stop("Unknown comparison: ", comparison)
    meta <- anno[anno$group %in% groups,
                 c("sample_id", "library_id", "group", "batch"), drop = FALSE]
    meta$y <- as.integer(meta$group == groups[1])
    meta$fold_stratum <- paste(meta$group, meta$batch, sep = "_")
  }
  meta <- meta[order(meta$group, meta$batch, meta$sample_id), , drop = FALSE]
  stopifnot(all(meta$library_id %in% colnames(count_all)))
  count <- count_all[, meta$library_id, drop = FALSE]

  task <- list()
  for (repeat_id in seq_len(n_repeat)) {
    outer <- make_foldid(
      meta$fold_stratum,
      n_outer, 100000 + repeat_id
    )
    for (fold in seq_len(n_outer)) {
      task[[length(task) + 1L]] <- list(
        repeat_id = repeat_id,
        fold = fold,
        train_idx = which(outer != fold),
        test_idx = which(outer == fold),
        seed = 200000 + repeat_id * 100 + fold
      )
    }
  }

  run_task <- function(z) {
    ans <- fit_outer(count, meta, z$train_idx, z$test_idx, z$seed)
    pred <- data.frame(
      comparison = comparison,
      validation = "repeated_nested_5fold",
      repeat_id = z$repeat_id,
      fold = z$fold,
      sample_id = meta$sample_id[z$test_idx],
      library_id = meta$library_id[z$test_idx],
      group = meta$group[z$test_idx],
      batch = meta$batch[z$test_idx],
      truth = meta$y[z$test_idx],
      score = ans$pred,
      stringsAsFactors = FALSE
    )
    detail <- data.frame(
      comparison = comparison,
      repeat_id = z$repeat_id,
      fold = z$fold,
      n_train = length(z$train_idx),
      n_test = length(z$test_idx),
      n_eligible = ans$eligible_n,
      n_feature = ans$best$n_feature,
      kernel = ans$best$kernel,
      cost = ans$best$cost,
      gamma_mult = ans$best$gamma_mult,
      gamma = ans$gamma,
      inner_auc = ans$best$inner_auc,
      inner_brier = ans$best$inner_brier,
      outer_auc = auc_value(meta$y[z$test_idx], ans$pred),
      stringsAsFactors = FALSE
    )
    selected <- data.frame(
      comparison = comparison,
      repeat_id = z$repeat_id,
      fold = z$fold,
      feature_id = ans$selected,
      gene_symbol = feature[ans$selected, "gene_symbol"],
      biotype = feature[ans$selected, "biotype"],
      rank = seq_along(ans$selected),
      stringsAsFactors = FALSE
    )
    list(pred = pred, detail = detail, selected = selected)
  }

  ans <- parallel::mclapply(task, run_task, mc.cores = n_cores,
                            mc.preschedule = FALSE)
  pred <- rbindlist(lapply(ans, `[[`, "pred"))
  detail <- rbindlist(lapply(ans, `[[`, "detail"))
  selected <- rbindlist(lapply(ans, `[[`, "selected"))

  repeat_auc <- pred[, .(AUC = auc_value(truth, score)),
                     by = .(comparison, repeat_id)]
  s1 <- summarize_patient_prediction(
    pred, comparison, "repeated_nested_5fold",
    700000 + match(comparison, c("MDA5_vs_HC", "ARS_vs_HC", "MDA5_vs_ARS", "Pooled_severity"))
  )

  if (grepl("_vs_HC$", comparison)) {
    patient_group <- sub("_vs_HC$", "", comparison)
    patient_batch <- sort(unique(meta$batch[meta$group == patient_group]))
    hc_batch <- sort(unique(meta$batch[meta$group == "HC"]))
    pair <- expand.grid(patient_batch = patient_batch, hc_batch = hc_batch,
                        stringsAsFactors = FALSE)
    batch_task <- lapply(seq_len(nrow(pair)), function(i) {
      test_idx <- which(
        (meta$group == patient_group & meta$batch == pair$patient_batch[i]) |
        (meta$group == "HC" & meta$batch == pair$hc_batch[i])
      )
      list(
        pair_id = i,
        held_out = paste(pair$patient_batch[i], pair$hc_batch[i], sep = "__"),
        train_idx = setdiff(seq_len(nrow(meta)), test_idx),
        test_idx = test_idx,
        seed = 800000 + i
      )
    })
    batch_validation <- "held_out_batch_pair"
  } else {
    hold_batch <- sort(unique(meta$batch))
    batch_task <- lapply(seq_along(hold_batch), function(i) {
      test_idx <- which(meta$batch == hold_batch[i])
      list(
        pair_id = i,
        held_out = hold_batch[i],
        train_idx = setdiff(seq_len(nrow(meta)), test_idx),
        test_idx = test_idx,
        seed = 800000 + i
      )
    })
    batch_validation <- "held_out_batch"
  }

  run_batch_task <- function(z) {
    ans <- fit_outer(count, meta, z$train_idx, z$test_idx, z$seed)
    pred <- data.frame(
      comparison = comparison,
      validation = batch_validation,
      repeat_id = 1,
      fold = z$pair_id,
      held_out = z$held_out,
      sample_id = meta$sample_id[z$test_idx],
      library_id = meta$library_id[z$test_idx],
      group = meta$group[z$test_idx],
      batch = meta$batch[z$test_idx],
      truth = meta$y[z$test_idx],
      score = ans$pred,
      stringsAsFactors = FALSE
    )
    detail <- data.frame(
      comparison = comparison,
      pair_id = z$pair_id,
      held_out = z$held_out,
      n_train = length(z$train_idx),
      n_test = length(z$test_idx),
      n_feature = ans$best$n_feature,
      kernel = ans$best$kernel,
      cost = ans$best$cost,
      gamma_mult = ans$best$gamma_mult,
      inner_auc = ans$best$inner_auc,
      test_auc = auc_value(meta$y[z$test_idx], ans$pred),
      stringsAsFactors = FALSE
    )
    list(pred = pred, detail = detail)
  }

  batch_ans <- parallel::mclapply(
    batch_task, run_batch_task, mc.cores = n_cores,
    mc.preschedule = FALSE
  )
  batch_pred <- rbindlist(lapply(batch_ans, `[[`, "pred"), fill = TRUE)
  batch_detail <- rbindlist(lapply(batch_ans, `[[`, "detail"), fill = TRUE)
  s2 <- summarize_patient_prediction(
    batch_pred, comparison, batch_validation,
    900000 + match(comparison, c("MDA5_vs_HC", "ARS_vs_HC", "MDA5_vs_ARS", "Pooled_severity"))
  )

  stability <- selected[, .(
    selected_folds = uniqueN(paste(repeat_id, fold, sep = "_")),
    mean_rank = mean(as.numeric(rank)),
    median_rank = as.numeric(median(rank))
  ), by = .(comparison, feature_id, gene_symbol, biotype)]
  stability[, selection_frequency := selected_folds / (n_repeat * n_outer)]
  setorder(stability, -selection_frequency, mean_rank, gene_symbol)

  list(
    pred = pred,
    detail = detail,
    selected = selected,
    repeat_auc = repeat_auc,
    performance = rbindlist(list(s1$performance, s2$performance), fill = TRUE),
    patient = rbindlist(list(s1$patient, s2$patient), fill = TRUE),
    stability = stability,
    batch_pred = batch_pred,
    batch_detail = batch_detail
  )
}

comparison_set <- strsplit(
  "Pooled_severity",
  ",", fixed = TRUE
)[[1]]
all_result <- lapply(comparison_set, function(comparison) {

  ans <- run_comparison(comparison)

  ans
})
names(all_result) <- comparison_set

write_one <- function(name, x) {
  fwrite(x$pred, file.path(result_dir, paste0(name, "_outer_predictions.tsv")), sep = "\t")
  fwrite(x$detail, file.path(result_dir, paste0(name, "_outer_fold_details.tsv")), sep = "\t")
  fwrite(x$selected, file.path(result_dir, paste0(name, "_selected_features_by_fold.tsv")), sep = "\t")
  fwrite(x$repeat_auc, file.path(result_dir, paste0(name, "_repeat_auc.tsv")), sep = "\t")
  fwrite(x$patient, file.path(result_dir, paste0(name, "_participant_mean_predictions.tsv")), sep = "\t")
  fwrite(x$stability, file.path(result_dir, paste0(name, "_feature_stability.tsv")), sep = "\t")
  fwrite(x$batch_pred, file.path(result_dir, paste0(name, "_batch_pair_predictions.tsv")), sep = "\t")
  fwrite(x$batch_detail, file.path(result_dir, paste0(name, "_batch_pair_details.tsv")), sep = "\t")
}

for (name in names(all_result)) write_one(name, all_result[[name]])
performance <- rbindlist(lapply(all_result, `[[`, "performance"), fill = TRUE)
fwrite(performance, file.path(result_dir, "Fig2_Fig3_model_performance.tsv"), sep = "\t")

input_summary <- data.frame(
  item = c(
    "input_mode", "count_matrix", "annotation", "usable_library_n", "MDA5_n", "ARS_n", "HC_n",
    "eligible_RNA_features", "outer_folds", "inner_folds", "repeats",
    "feature_number_grid", "kernel_grid", "cost_grid", "gamma_multiplier_grid",
    "bootstrap_iterations", "random_seed"
  ),
  value = c(
    input_mode, input_path, file.path(dm, "annotation.txt"),
    nrow(anno), sum(anno$group == "MDA5"), sum(anno$group == "ARS"),
    sum(anno$group == "HC"), nrow(count_all), n_outer, n_inner, n_repeat,
    paste(n_grid, collapse = ","), "linear,radial", paste(cost_grid, collapse = ","),
    paste(gamma_mult_grid, collapse = ","), n_boot, 20260826
  ),
  stringsAsFactors = FALSE
)
fwrite(input_summary, file.path(result_dir, "input_and_parameter_summary.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
md5 <- tools::md5sum(c(input_path, file.path(dm, "annotation.txt")))
fwrite(data.frame(path = names(md5), md5 = unname(md5)),
       file.path(log_dir, "input_md5.tsv"), sep = "\t")

print(performance)

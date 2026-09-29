# RStudio工作目录设为本文件夹
# 先运行biomarker.R
library(data.table)
library(e1071)


options(stringsAsFactors = FALSE)

dm <- "inputs"
root <- "generated/biomarker"
manifest_path <- file.path(root, "meta", "global_discovery_validation_manifest.tsv")
primary_result_dir <- file.path(root, "results")
result_dir <- file.path(primary_result_dir, "sensitivity_delta_auc_0.01")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

n_repeat <- 10L
n_fold <- 5L
n_boot <- 2000L
delta_auc <- 0.01

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

bootstrap_auc <- function(y, score, B, seed) {
  set.seed(seed)
  pos <- which(y == 1)
  neg <- which(y == 0)
  val <- replicate(B, {
    idx <- c(
      sample(pos, length(pos), replace = TRUE),
      sample(neg, length(neg), replace = TRUE)
    )
    auc_value(y[idx], score[idx])
  })
  unname(quantile(val, c(0.025, 0.975), na.rm = TRUE))
}

log_cpm <- function(count, library_size) {
  stopifnot(all(library_size > 0))
  log2(sweep(count, 2, library_size / 1e6, "/") + 1)
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

fit_svm_probability <- function(x_train, y_train, x_test, cost, seed) {
  set.seed(seed)
  fit <- svm(
    x = x_train,
    y = factor(y_train, levels = c(0, 1)),
    type = "C-classification",
    kernel = "linear",
    cost = cost,
    class.weights = class_weight(y_train),
    probability = TRUE,
    scale = FALSE
  )
  pred <- predict(fit, x_test, probability = TRUE)
  as.numeric(attr(pred, "probabilities")[, "1"])
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

best_threshold <- function(y, score) {
  candidate <- sort(unique(c(0.5, score)))
  stat <- rbindlist(lapply(candidate, function(th) {
    cls <- as.integer(score >= th)
    sensitivity <- mean(cls[y == 1] == 1)
    specificity <- mean(cls[y == 0] == 0)
    data.frame(
      threshold = th,
      sensitivity = sensitivity,
      specificity = specificity,
      youden = sensitivity + specificity - 1
    )
  }))
  stat[order(-youden, abs(threshold - 0.5))][1]
}

performance_row <- function(comparison, set_name, y, score, threshold, seed) {
  cls <- as.integer(score >= threshold)
  ci <- bootstrap_auc(y, score, n_boot, seed)
  data.frame(
    comparison = comparison,
    set = set_name,
    n = length(y),
    AUC = auc_value(y, score),
    AUC_low = ci[1],
    AUC_high = ci[2],
    threshold = threshold,
    accuracy = mean(cls == y),
    sensitivity = mean(cls[y == 1] == 1),
    specificity = mean(cls[y == 0] == 0),
    Brier = mean((score - y)^2)
  )
}

meta_all <- fread(manifest_path, data.table = FALSE)

count_path <- file.path(dm, "1204", "gencode.txt")
raw <- fread(count_path, data.table = FALSE, check.names = FALSE)
feature_id <- raw[[1]]
raw <- raw[, -1, drop = FALSE]
raw <- raw[, -c(40:44, 50:54), drop = FALSE]
colnames(raw)[35:44] <- sub("-GP", "", colnames(raw)[35:44], fixed = TRUE)
rownames(raw) <- feature_id
raw <- raw[, meta_all$library_id, drop = FALSE]

part <- tstrsplit(rownames(raw), "|", fixed = TRUE, fill = "")
biotype <- part[[4]]
gene_symbol <- part[[3]]
rna_type <- c(
  "protein_coding", "lncRNA", "Mt_tRNA", "scaRNA", "sRNA",
  "IG_C_gene", "IG_D_gene", "IG_V_gene", "IG_J_gene",
  "TR_C_gene", "TR_D_gene", "TR_V_gene", "TR_J_gene",
  "miRNA", "snRNA", "snoRNA"
)
keep <- biotype %in% rna_type & gene_symbol != "Y_RNA"
count_all <- as.matrix(raw[keep, , drop = FALSE])
storage.mode(count_all) <- "double"
rm(raw)

comparison_spec <- data.frame(
  comparison = c("MDA5_vs_HC", "ARS_vs_HC", "MDA5_vs_ARS"),
  positive = c("MDA5", "ARS", "MDA5"),
  negative = c("HC", "HC", "ARS"),
  seed_offset = c(1000L, 2000L, 3000L)
)

run_sensitivity <- function(comparison, positive_group, negative_group, seed_offset) {
  curve <- fread(file.path(
    primary_result_dir, comparison, "panel_size_one_se_selection.tsv"
  ))
  best_auc <- max(curve$mean_AUC)
  selected <- curve[mean_AUC >= best_auc - delta_auc][order(k)][1]
  chosen_k <- selected$k
  chosen_cost <- selected$cost

  primary_panel <- fread(file.path(
    primary_result_dir, comparison, "fixed_biomarker_panel.tsv"
  ))
  if (nrow(primary_panel) < chosen_k) {
    stop(comparison, ": the primary panel is shorter than the sensitivity panel")
  }
  fixed_panel <- primary_panel[order(panel_rank)][seq_len(chosen_k)]
  fixed_panel[, sensitivity_panel_rank := seq_len(.N)]
  fixed_feature <- fixed_panel$feature_id

  meta <- meta_all[meta_all$group %in% c(positive_group, negative_group), ]
  meta$y <- as.integer(meta$group == positive_group)
  count <- count_all[, meta$library_id, drop = FALSE]
  discovery_idx <- which(meta$analysis_set == "discovery")
  validation_idx <- which(meta$analysis_set == "internal_validation")
  discovery_meta <- meta[discovery_idx, ]
  validation_meta <- meta[validation_idx, ]
  discovery_count <- count[, discovery_idx, drop = FALSE]
  validation_count <- count[, validation_idx, drop = FALSE]

  oof_record <- list()
  for (repeat_id in seq_len(n_repeat)) {
    foldid <- make_foldid(
      paste(discovery_meta$group, discovery_meta$batch, sep = "_"),
      n_fold,
      500000 + seed_offset + repeat_id
    )
    for (fold in seq_len(n_fold)) {
      train <- which(foldid != fold)
      test <- which(foldid == fold)
      xtr <- t(log_cpm(
        discovery_count[fixed_feature, train, drop = FALSE],
        colSums(discovery_count[, train, drop = FALSE])
      ))
      xte <- t(log_cpm(
        discovery_count[fixed_feature, test, drop = FALSE],
        colSums(discovery_count[, test, drop = FALSE])
      ))
      ss <- scale_train_test(xtr, xte)
      score <- fit_svm_probability(
        ss$train, discovery_meta$y[train], ss$test, chosen_cost,
        600000 + seed_offset + repeat_id * 100 + fold
      )
      oof_record[[length(oof_record) + 1L]] <- data.frame(
        comparison = comparison,
        repeat_id = repeat_id,
        fold = fold,
        sample_id = discovery_meta$sample_id[test],
        truth = discovery_meta$y[test],
        score = score
      )
    }
  }
  oof <- rbindlist(oof_record)
  oof_mean <- oof[, .(truth = unique(truth), score = mean(score)),
                  by = .(comparison, sample_id)]
  threshold_stat <- best_threshold(oof_mean$truth, oof_mean$score)
  threshold <- threshold_stat$threshold

  x_discovery <- t(log_cpm(
    discovery_count[fixed_feature, , drop = FALSE], colSums(discovery_count)
  ))
  x_validation <- t(log_cpm(
    validation_count[fixed_feature, , drop = FALSE], colSums(validation_count)
  ))
  ss <- scale_train_test(x_discovery, x_validation)
  validation_score <- fit_svm_probability(
    ss$train, discovery_meta$y, ss$test, chosen_cost, 700000 + seed_offset
  )
  validation_pred <- data.frame(
    comparison = comparison,
    sample_id = validation_meta$sample_id,
    library_id = validation_meta$library_id,
    group = validation_meta$group,
    batch = validation_meta$batch,
    truth = validation_meta$y,
    score = validation_score,
    predicted = as.integer(validation_score >= threshold)
  )

  performance <- rbind(
    performance_row(
      comparison, "development_fixed_panel_resampling",
      oof_mean$truth, oof_mean$score, threshold, 800000 + seed_offset
    ),
    performance_row(
      comparison, "internal_validation_fixed_panel",
      validation_pred$truth, validation_pred$score, threshold,
      900000 + seed_offset
    )
  )
  confusion <- as.data.frame.matrix(table(
    true = factor(validation_pred$truth, levels = c(1, 0)),
    predicted = factor(validation_pred$predicted, levels = c(1, 0))
  ))
  confusion$true <- c(positive_group, negative_group)
  names(confusion)[1:2] <- c(
    paste0("predicted_", positive_group),
    paste0("predicted_", negative_group)
  )

  comparison_dir <- file.path(result_dir, comparison)
  dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(data.frame(
    comparison = comparison,
    best_discovery_mean_AUC = best_auc,
    delta_AUC = delta_auc,
    selected_panel_size = chosen_k,
    selected_cost = chosen_cost
  ), file.path(comparison_dir, "selection_rule.tsv"), sep = "\t")
  fwrite(fixed_panel, file.path(comparison_dir, "fixed_biomarker_panel.tsv"), sep = "\t")
  fwrite(oof_mean, file.path(comparison_dir, "development_fixed_panel_predictions.tsv"), sep = "\t")
  fwrite(threshold_stat, file.path(comparison_dir, "development_threshold.tsv"), sep = "\t")
  fwrite(validation_pred, file.path(comparison_dir, "internal_validation_predictions.tsv"), sep = "\t")
  fwrite(performance, file.path(comparison_dir, "fixed_panel_performance.tsv"), sep = "\t")
  fwrite(confusion, file.path(comparison_dir, "internal_validation_confusion_matrix.tsv"), sep = "\t")

  performance
}

performance <- rbindlist(lapply(seq_len(nrow(comparison_spec)), function(i) {
  run_sensitivity(
    comparison_spec$comparison[i],
    comparison_spec$positive[i],
    comparison_spec$negative[i],
    comparison_spec$seed_offset[i]
  )
}))
fwrite(performance, file.path(result_dir, "fixed_panel_performance_all.tsv"), sep = "\t")
print(performance)

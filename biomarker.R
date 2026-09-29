# RStudio工作目录设为本文件夹
library(data.table)
library(edgeR)
library(e1071)
library(openxlsx)


options(stringsAsFactors = FALSE)
set.seed(20260830)

dm <- "inputs"
root <- "generated/biomarker"

script_dir <- file.path(root, "scripts")
meta_dir <- file.path(root, "meta")
result_dir <- file.path(root, "results")
log_dir <- file.path(root, "logs")
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

n_repeat <- 10L
n_fold <- 5L
n_boot <- 2000L
max_panel <- 8L
fdr_cutoff <- 0.05
logfc_cutoff <- 1
cost_grid <- c(0.01, 0.1, 1, 10, 100)

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

log_cpm <- function(count, library_size = NULL) {
  lib <- if (is.null(library_size)) colSums(count) else library_size
  stopifnot(all(lib > 0))
  log2(sweep(count, 2, lib / 1e6, "/") + 1)
}

scale_train_test <- function(x_train, x_test) {
  center <- colMeans(x_train)
  scale <- apply(x_train, 2, sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  list(
    train = sweep(sweep(x_train, 2, center, "-"), 2, scale, "/"),
    test = sweep(sweep(x_test, 2, center, "-"), 2, scale, "/"),
    center = center,
    scale = scale
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
  prob <- attr(pred, "probabilities")
  if (is.null(prob) || !"1" %in% colnames(prob)) {
    stop("SVM probability output is incomplete")
  }
  as.numeric(prob[, "1"])
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

rank_de <- function(count, y) {
  group <- factor(y, levels = c(0, 1))
  dge <- DGEList(counts = count, group = group)
  keep <- filterByExpr(dge, group = group, min.count = 2, min.prop = 0.2)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  design <- model.matrix(~group)
  dge <- estimateDisp(dge, design, robust = TRUE)
  fit <- glmQLFit(dge, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = 2)
  tab <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  tab$feature_id <- rownames(tab)
  tab$pass <- tab$FDR < fdr_cutoff & abs(tab$logFC) >= logfc_cutoff
  tab <- tab[order(!tab$pass, tab$FDR, -abs(tab$logFC)), , drop = FALSE]
  tab$rank_all <- seq_len(nrow(tab))
  tab$rank_pass <- NA_integer_
  tab$rank_pass[tab$pass] <- seq_len(sum(tab$pass))
  rownames(tab) <- NULL
  tab
}

best_threshold <- function(y, score) {
  candidate <- sort(unique(c(0.5, score)))
  stat <- lapply(candidate, function(th) {
    cls <- as.integer(score >= th)
    sensitivity <- mean(cls[y == 1] == 1)
    specificity <- mean(cls[y == 0] == 0)
    data.frame(
      threshold = th,
      sensitivity = sensitivity,
      specificity = specificity,
      youden = sensitivity + specificity - 1
    )
  })
  stat <- rbindlist(stat)
  stat <- stat[order(-youden, abs(threshold - 0.5))]
  stat[1]
}

performance_row <- function(comparison, set_name, y, score, threshold, seed) {
  cls <- as.integer(score >= threshold)
  ci <- bootstrap_auc(y, score, n_boot, seed)
  eps <- 1e-6
  p <- pmin(pmax(score, eps), 1 - eps)
  cal <- try(glm(y ~ qlogis(p), family = binomial()), silent = TRUE)
  cal_intercept <- NA_real_
  cal_slope <- NA_real_
  if (!inherits(cal, "try-error")) {
    cc <- coef(cal)
    cal_intercept <- unname(cc[1])
    cal_slope <- unname(cc[2])
  }
  data.frame(
    comparison = comparison,
    set = set_name,
    n = length(y),
    positive_n = sum(y == 1),
    negative_n = sum(y == 0),
    AUC = auc_value(y, score),
    AUC_low = ci[1],
    AUC_high = ci[2],
    threshold = threshold,
    accuracy = mean(cls == y),
    sensitivity = mean(cls[y == 1] == 1),
    specificity = mean(cls[y == 0] == 0),
    Brier = mean((score - y)^2),
    calibration_intercept = cal_intercept,
    calibration_slope = cal_slope,
    stringsAsFactors = FALSE
  )
}

make_split_once <- function(meta, validation_fraction = 0.30, seed = 20260830) {
  set.seed(seed)
  meta$split_stratum <- paste(meta$group, meta$batch, sep = "_")
  validation <- rep(FALSE, nrow(meta))
  for (s in unique(meta$split_stratum)) {
    idx <- which(meta$split_stratum == s)
    n_val <- max(1L, round(length(idx) * validation_fraction))
    sex <- ifelse(is.na(meta$sex[idx]) | meta$sex[idx] == "", "Unknown", meta$sex[idx])
    sex_count <- table(sex)
    quota_raw <- n_val * as.numeric(sex_count) / sum(sex_count)
    quota <- floor(quota_raw)
    remainder <- n_val - sum(quota)
    if (remainder > 0) {
      add <- order(quota_raw - quota, decreasing = TRUE)[seq_len(remainder)]
      quota[add] <- quota[add] + 1L
    }
    names(quota) <- names(sex_count)
    chosen <- integer()
    for (sx in names(quota)) {
      sx_idx <- idx[sex == sx]
      take <- min(quota[sx], length(sx_idx))
      if (take > 0) chosen <- c(chosen, sample(sx_idx, take))
    }
    if (length(chosen) < n_val) {
      chosen <- c(chosen, sample(setdiff(idx, chosen), n_val - length(chosen)))
    }
    validation[chosen] <- TRUE
  }
  meta$analysis_set <- ifelse(validation, "internal_validation", "discovery")
  meta
}

split_balance_score <- function(meta) {
  score <- 0
  for (g in unique(meta$group)) {
    z <- meta[meta$group == g, ]
    age_discovery <- z$age[z$analysis_set == "discovery"]
    age_validation <- z$age[z$analysis_set == "internal_validation"]
    pooled_sd <- sd(z$age, na.rm = TRUE)
    if (is.finite(pooled_sd) && pooled_sd > 0) {
      score <- score + abs(
        mean(age_discovery, na.rm = TRUE) - mean(age_validation, na.rm = TRUE)
      ) / pooled_sd
    }
    sex <- as.character(z$sex)
    sex[is.na(sex) | sex == ""] <- "Unknown"
    level <- sort(unique(sex))
    p_discovery <- prop.table(table(factor(
      sex[z$analysis_set == "discovery"], levels = level
    )))
    p_validation <- prop.table(table(factor(
      sex[z$analysis_set == "internal_validation"], levels = level
    )))
    score <- score + sum(abs(p_discovery - p_validation))
  }
  as.numeric(score)
}

make_split <- function(meta, validation_fraction = 0.30,
                       base_seed = 20260830, n_candidate = 2000L) {
  candidate <- vector("list", n_candidate)
  audit <- vector("list", n_candidate)
  for (i in seq_len(n_candidate)) {
    seed <- base_seed + i - 1L
    z <- make_split_once(meta, validation_fraction, seed)
    score <- split_balance_score(z)
    candidate[[i]] <- z
    audit[[i]] <- data.frame(candidate_id = i, seed = seed, balance_score = score)
  }
  audit <- rbindlist(audit)
  setorder(audit, balance_score, seed)
  best <- audit$candidate_id[1]
  list(meta = candidate[[best]], audit = audit)
}

anno <- read.delim(file.path(dm, "annotation.txt"), check.names = FALSE)
colnames(anno)[2] <- "library_id"
anno <- anno[anno$mistake == 0, c("sample_id", "library_id", "group")]
anno$batch <- substr(anno$library_id, 1, 6)
stopifnot(nrow(anno) == 155L)
stopifnot(!anyDuplicated(anno$sample_id), !anyDuplicated(anno$library_id))

# Identical demographics supplied without source patient identifiers.
clinical_path <- file.path(dm, "biomarker_manifest.tsv")
supp_path <- clinical_path
meta_all <- read.delim(clinical_path, check.names = FALSE)
stopifnot(setequal(meta_all$sample_id, anno$sample_id))
# Use the published participant allocation; do not draw a new holdout split.
split_object <- list(meta = meta_all, audit = read.delim(file.path(dm, "split_covariate_balance_candidates.tsv")))
selected_split_seed <- split_object$audit$seed[1]

count_path <- file.path(dm, "1204", "gencode.txt")
raw <- fread(count_path, data.table = FALSE, check.names = FALSE)
feature_id <- raw[[1]]
raw <- raw[, -1, drop = FALSE]
raw <- raw[, -c(40:44, 50:54), drop = FALSE]
colnames(raw)[35:44] <- sub("-GP", "", colnames(raw)[35:44], fixed = TRUE)
rownames(raw) <- feature_id
stopifnot(setequal(colnames(raw), meta_all$library_id))
raw <- raw[, meta_all$library_id, drop = FALSE]

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

fwrite(meta_all, file.path(meta_dir, "global_discovery_validation_manifest.tsv"), sep = "\t")
fwrite(split_object$audit, file.path(meta_dir, "split_covariate_balance_candidates.tsv"), sep = "\t")

split_balance <- meta_all[, c(
  "sample_id", "group", "batch", "sex", "sex_code", "age", "analysis_set"
)]
split_summary <- as.data.table(split_balance)[, .(
  n = .N,
  sex_level_1_n = sum(sex_code == 1, na.rm = TRUE),
  sex_level_2_n = sum(sex_code == 2, na.rm = TRUE),
  sex_missing_n = sum(is.na(sex) | sex == ""),
  age_mean = mean(age, na.rm = TRUE),
  age_sd = sd(age, na.rm = TRUE),
  age_missing_n = sum(is.na(age))
), by = .(group, batch, analysis_set)]
fwrite(split_summary, file.path(meta_dir, "split_balance_by_group_batch.tsv"), sep = "\t")

run_comparison <- function(comparison, positive_group, negative_group, seed_offset) {

  meta <- meta_all[meta_all$group %in% c(positive_group, negative_group), ]
  meta$y <- as.integer(meta$group == positive_group)
  count <- count_all[, meta$library_id, drop = FALSE]
  discovery_idx <- which(meta$analysis_set == "discovery")
  validation_idx <- which(meta$analysis_set == "internal_validation")
  discovery_meta <- meta[discovery_idx, ]
  discovery_count <- count[, discovery_idx, drop = FALSE]
  validation_meta <- meta[validation_idx, ]
  validation_count <- count[, validation_idx, drop = FALSE]

  fold_record <- list()
  pred_record <- list()
  stability_record <- list()
  z <- 0L
  for (repeat_id in seq_len(n_repeat)) {
    foldid <- make_foldid(
      paste(discovery_meta$group, discovery_meta$batch, sep = "_"),
      n_fold,
      300000 + seed_offset + repeat_id
    )
    for (fold in seq_len(n_fold)) {
      train <- which(foldid != fold)
      test <- which(foldid == fold)
      de <- rank_de(discovery_count[, train, drop = FALSE], discovery_meta$y[train])
      selected <- head(de$feature_id, max_panel)
      stability_record[[length(stability_record) + 1L]] <- data.frame(
        comparison = comparison,
        repeat_id = repeat_id,
        fold = fold,
        feature_id = selected,
        rank = seq_along(selected),
        logFC = de$logFC[match(selected, de$feature_id)],
        FDR = de$FDR[match(selected, de$feature_id)],
        pass_threshold = de$pass[match(selected, de$feature_id)],
        stringsAsFactors = FALSE
      )

      x_train_all <- log_cpm(
        discovery_count[selected, train, drop = FALSE],
        colSums(discovery_count[, train, drop = FALSE])
      )
      x_test_all <- log_cpm(
        discovery_count[selected, test, drop = FALSE],
        colSums(discovery_count[, test, drop = FALSE])
      )
      for (k in seq_len(max_panel)) {
        n_use <- min(k, length(selected))
        use <- selected[seq_len(n_use)]
        xtr <- t(x_train_all[use, , drop = FALSE])
        xte <- t(x_test_all[use, , drop = FALSE])
        ss <- scale_train_test(xtr, xte)
        for (cost in cost_grid) {
          z <- z + 1L
          score <- fit_svm_probability(
            ss$train, discovery_meta$y[train], ss$test,
            cost,
            400000 + seed_offset + repeat_id * 10000 + fold * 100 + z
          )
          pred_record[[length(pred_record) + 1L]] <- data.frame(
            comparison = comparison,
            repeat_id = repeat_id,
            fold = fold,
            k = k,
            cost = cost,
            sample_id = discovery_meta$sample_id[test],
            truth = discovery_meta$y[test],
            score = score,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  cv_pred <- rbindlist(pred_record)
  dir.create(file.path(result_dir, comparison), recursive = TRUE, showWarnings = FALSE)
  fwrite(unique(cv_pred[, .(comparison, repeat_id, fold, sample_id)]),
         file.path(result_dir, comparison, "screening_fold_assignments.tsv"), sep = "\t")

  cv_repeat <- cv_pred[, .(AUC = auc_value(truth, score)),
                       by = .(comparison, k, cost, repeat_id)]
  cv_summary <- cv_repeat[, .(
    mean_AUC = mean(AUC, na.rm = TRUE),
    sd_AUC = sd(AUC, na.rm = TRUE),
    se_AUC = sd(AUC, na.rm = TRUE) / sqrt(.N),
    min_AUC = min(AUC, na.rm = TRUE),
    max_AUC = max(AUC, na.rm = TRUE)
  ), by = .(comparison, k, cost)]
  best_by_k <- cv_summary[order(k, -mean_AUC, se_AUC, cost), .SD[1], by = k]
  best_row <- best_by_k[which.max(mean_AUC)]
  one_se_limit <- best_row$mean_AUC - best_row$se_AUC
  chosen <- best_by_k[mean_AUC >= one_se_limit][order(k, -mean_AUC)][1]
  chosen_k <- chosen$k
  chosen_cost <- chosen$cost
  best_by_k$one_se_limit <- one_se_limit
  best_by_k$selected_panel_size <- best_by_k$k == chosen_k

  stability_long <- rbindlist(stability_record)
  stability <- stability_long[, .(
    selected_folds = uniqueN(paste(repeat_id, fold, sep = "_")),
    selection_frequency = uniqueN(paste(repeat_id, fold, sep = "_")) / (n_repeat * n_fold),
    median_rank = as.numeric(median(rank)),
    mean_rank = as.numeric(mean(rank)),
    threshold_pass_frequency = mean(pass_threshold),
    direction_positive_fraction = mean(logFC > 0),
    median_logFC = as.numeric(median(logFC)),
    median_FDR = as.numeric(median(FDR))
  ), by = .(comparison, feature_id)]
  stability$gene_symbol <- feature[stability$feature_id, "gene_symbol"]
  stability$biotype <- feature[stability$feature_id, "biotype"]

  full_de <- rank_de(discovery_count, discovery_meta$y)
  full_de$gene_symbol <- feature[full_de$feature_id, "gene_symbol"]
  full_de$biotype <- feature[full_de$feature_id, "biotype"]
  full_de$discovery_detection_positive <- rowMeans(
    discovery_count[full_de$feature_id, discovery_meta$y == 1, drop = FALSE] > 0
  )
  full_de$discovery_detection_negative <- rowMeans(
    discovery_count[full_de$feature_id, discovery_meta$y == 0, drop = FALSE] > 0
  )

  panel_rank <- merge(
    full_de[full_de$pass, ],
    stability,
    by = c("feature_id", "gene_symbol", "biotype"),
    all.x = TRUE,
    suffixes = c("_discovery", "_resampling")
  )
  panel_rank$selection_frequency[is.na(panel_rank$selection_frequency)] <- 0
  panel_rank$threshold_pass_frequency[is.na(panel_rank$threshold_pass_frequency)] <- 0
  panel_rank$median_rank[is.na(panel_rank$median_rank)] <- Inf
  panel_rank <- panel_rank[order(
    -panel_rank$selection_frequency,
    -panel_rank$threshold_pass_frequency,
    panel_rank$median_rank,
    panel_rank$FDR,
    -abs(panel_rank$logFC)
  ), ]
  fixed_panel <- head(panel_rank, chosen_k)
  fixed_panel$panel_rank <- seq_len(nrow(fixed_panel))
  fixed_feature <- fixed_panel$feature_id
  if (length(fixed_feature) != chosen_k) {
    stop(comparison, ": fixed panel is shorter than selected panel size")
  }

  fixed_oof_record <- list()
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
        ss$train, discovery_meta$y[train], ss$test,
        chosen_cost,
        600000 + seed_offset + repeat_id * 100 + fold
      )
      fixed_oof_record[[length(fixed_oof_record) + 1L]] <- data.frame(
        comparison = comparison,
        repeat_id = repeat_id,
        fold = fold,
        sample_id = discovery_meta$sample_id[test],
        truth = discovery_meta$y[test],
        score = score,
        stringsAsFactors = FALSE
      )
    }
  }
  fixed_oof <- rbindlist(fixed_oof_record)
  fixed_oof_mean <- fixed_oof[, .(
    truth = unique(truth),
    score = mean(score)
  ), by = .(comparison, sample_id)]
  threshold_stat <- best_threshold(fixed_oof_mean$truth, fixed_oof_mean$score)
  threshold <- threshold_stat$threshold

  x_discovery <- t(log_cpm(
    discovery_count[fixed_feature, , drop = FALSE],
    colSums(discovery_count)
  ))
  x_validation <- t(log_cpm(
    validation_count[fixed_feature, , drop = FALSE],
    colSums(validation_count)
  ))
  ss <- scale_train_test(x_discovery, x_validation)
  validation_score <- fit_svm_probability(
    ss$train, discovery_meta$y, ss$test,
    chosen_cost,
    700000 + seed_offset
  )
  validation_pred <- data.frame(
    comparison = comparison,
    sample_id = validation_meta$sample_id,
    library_id = validation_meta$library_id,
    group = validation_meta$group,
    batch = validation_meta$batch,
    truth = validation_meta$y,
    score = validation_score,
    predicted = as.integer(validation_score >= threshold),
    stringsAsFactors = FALSE
  )

  perf_discovery <- performance_row(
    comparison, "discovery_fixed_panel_repeated_oof",
    fixed_oof_mean$truth, fixed_oof_mean$score, threshold,
    800000 + seed_offset
  )
  perf_validation <- performance_row(
    comparison, "internal_validation_fixed_panel",
    validation_pred$truth, validation_pred$score, threshold,
    900000 + seed_offset
  )
  performance <- rbind(perf_discovery, perf_validation)

  confusion <- as.data.frame.matrix(table(
    true = factor(validation_pred$truth, levels = c(1, 0)),
    predicted = factor(validation_pred$predicted, levels = c(1, 0))
  ))
  confusion$true <- c(positive_group, negative_group)
  names(confusion)[1:2] <- c(paste0("predicted_", positive_group),
                             paste0("predicted_", negative_group))

  comparison_dir <- file.path(result_dir, comparison)
  dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(full_de, file.path(comparison_dir, "discovery_differential_expression.tsv"), sep = "\t")
  fwrite(stability, file.path(comparison_dir, "feature_stability.tsv"), sep = "\t")
  fwrite(stability_long, file.path(comparison_dir, "feature_selection_by_resample.tsv"), sep = "\t")
  fwrite(best_by_k, file.path(comparison_dir, "panel_size_one_se_selection.tsv"), sep = "\t")
  fwrite(cv_repeat, file.path(comparison_dir, "panel_size_repeat_auc.tsv"), sep = "\t")
  fwrite(fixed_panel, file.path(comparison_dir, "fixed_biomarker_panel.tsv"), sep = "\t")
  fwrite(fixed_oof_mean, file.path(comparison_dir, "discovery_fixed_panel_oof_predictions.tsv"), sep = "\t")
  fwrite(threshold_stat, file.path(comparison_dir, "discovery_threshold.tsv"), sep = "\t")
  fwrite(validation_pred, file.path(comparison_dir, "internal_validation_predictions.tsv"), sep = "\t")
  fwrite(performance, file.path(comparison_dir, "fixed_panel_performance.tsv"), sep = "\t")
  fwrite(confusion, file.path(comparison_dir, "internal_validation_confusion_matrix.tsv"), sep = "\t")


  list(
    panel = fixed_panel,
    performance = performance,
    validation = validation_pred,
    panel_size = best_by_k,
    threshold = threshold_stat
  )
}

comparison_spec <- data.frame(
  comparison = c("MDA5_vs_HC", "ARS_vs_HC", "MDA5_vs_ARS"),
  positive = c("MDA5", "ARS", "MDA5"),
  negative = c("HC", "HC", "ARS"),
  seed_offset = c(1000, 2000, 3000),
  stringsAsFactors = FALSE
)
requested_comparison <- strsplit(
  "MDA5_vs_HC,ARS_vs_HC,MDA5_vs_ARS",
  ",", fixed = TRUE
)[[1]]
comparison_spec <- comparison_spec[comparison_spec$comparison %in% requested_comparison, ]

all_result <- lapply(seq_len(nrow(comparison_spec)), function(i) {
  run_comparison(
    comparison_spec$comparison[i],
    comparison_spec$positive[i],
    comparison_spec$negative[i],
    comparison_spec$seed_offset[i]
  )
})
names(all_result) <- comparison_spec$comparison

available_comparison <- c("MDA5_vs_HC", "ARS_vs_HC", "MDA5_vs_ARS")
available_comparison <- available_comparison[file.exists(file.path(
  result_dir, available_comparison, "fixed_biomarker_panel.tsv"
))]
panel_summary <- rbindlist(lapply(available_comparison, function(x) {
  p <- fread(file.path(result_dir, x, "fixed_biomarker_panel.tsv"), data.table = FALSE)
  data.frame(
    comparison = x,
    panel_rank = p$panel_rank,
    feature_id = p$feature_id,
    gene_symbol = p$gene_symbol,
    biotype = p$biotype,
    logFC = p$logFC,
    FDR = p$FDR,
    selection_frequency = p$selection_frequency,
    stringsAsFactors = FALSE
  )
}), fill = TRUE)
performance_summary <- rbindlist(lapply(available_comparison, function(x) {
  fread(file.path(result_dir, x, "fixed_panel_performance.tsv"))
}), fill = TRUE)
fwrite(panel_summary, file.path(result_dir, "fixed_biomarker_panels_all.tsv"), sep = "\t")
fwrite(performance_summary, file.path(result_dir, "fixed_panel_performance_all.tsv"), sep = "\t")

parameter <- data.frame(
  item = c(
    "raw_count_matrix", "annotation", "clinical_meta", "HC_demo_source",
    "usable_libraries", "discovery_validation_split", "split_seed",
    "DE_FDR_cutoff", "DE_abs_log2FC_cutoff", "resampling",
    "maximum_panel_size", "panel_size_rule", "classifier",
    "SVM_cost_grid", "CI_method", "bootstrap_iterations"
  ),
  value = c(
    count_path, file.path(dm, "annotation.txt"), clinical_path, supp_path,
    nrow(meta_all), "70% discovery / 30% internal validation within group x batch",
    paste0("base=20260830; selected=", selected_split_seed,
           "; covariate balance only"), fdr_cutoff, logfc_cutoff,
    paste0(n_repeat, " repeats x ", n_fold, " folds in discovery"),
    max_panel, "smallest panel within one SE of best discovery CV AUC",
    "class-weighted linear SVM", paste(cost_grid, collapse = ","),
    "participant-level stratified bootstrap", n_boot
  ),
  stringsAsFactors = FALSE
)
fwrite(parameter, file.path(result_dir, "analysis_parameters.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
md5 <- tools::md5sum(c(count_path, file.path(dm, "annotation.txt"), clinical_path, supp_path))
fwrite(data.frame(path = names(md5), md5 = unname(md5)),
       file.path(log_dir, "input_md5.tsv"), sep = "\t")

print(panel_summary)
print(performance_summary)

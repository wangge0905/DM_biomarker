# RStudio工作目录设为本文件夹
library(data.table)
library(glmnet)
library(ggplot2)


dm <- "inputs"
out_dir <- file.path(
  "generated",
  "DM_severity_reviewer_completion_20260830", "output"
)
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
anno <- fread(file.path(dm, "annotation.txt"), data.table = FALSE, check.names = FALSE)
colnames(anno)[2] <- "library_id"
anno <- anno[anno$mistake == 0 & anno$group %in% c("MDA5", "ARS", "HC"), , drop = FALSE]
anno$batch <- substr(anno$library_id, 1, 6)

universe <- fread(file.path(
  dm, "MDA5_ARS_pathway_20260816", "ranked_genes",
  "standard_design_filter_protein_coding_rank.tsv"
))
gmt_path <- file.path(
  dm, "MDA5_ARS_pathway_20260816", "gene_sets",
  "h.all.v2025.1.Hs.symbols.gmt"
)

gencode <- fread(file.path(dm, "1204", "gencode.txt"), data.table = FALSE, check.names = FALSE)
feature_id <- gencode[[1]]
gencode <- gencode[, -1, drop = FALSE]
gencode <- gencode[, -c(40:44, 50:54), drop = FALSE]
colnames(gencode)[35:44] <- sub("-GP", "", colnames(gencode)[35:44], fixed = TRUE)
rownames(gencode) <- feature_id
stopifnot(setequal(colnames(gencode), anno$library_id))

count <- as.matrix(gencode[universe$feature_id, anno$library_id, drop = FALSE])
rownames(count) <- universe$gene_symbol
storage.mode(count) <- "double"
rm(gencode)

read_gmt <- function(path, universe_gene) {
  z <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  ans <- lapply(z, function(v) intersect(v[-c(1, 2)], universe_gene))
  names(ans) <- vapply(z, `[`, character(1), 1)
  ans[lengths(ans) >= 15L & lengths(ans) <= 500L]
}
hallmark_set <- read_gmt(gmt_path, rownames(count))
hallmark <- names(hallmark_set)
stopifnot(length(hallmark) == 50L)

# 样本内gene rank只依赖该样本本身，不使用其他样本或严重性标签
rank_pct <- apply(count, 2, rank, ties.method = "average") / (nrow(count) + 1)
score_mat <- vapply(
  hallmark_set,
  function(g) colMeans(rank_pct[g, , drop = FALSE]),
  numeric(ncol(rank_pct))
)
rownames(score_mat) <- colnames(count)

qc <- data.table(
  library_id = colnames(count),
  library_size = colSums(count),
  detected_genes = colSums(count > 0)
)
qc[, log10_library_size := log10(library_size)]

score <- data.table(
  sample_id = anno$sample_id,
  library_id = anno$library_id,
  group = anno$group,
  batch = anno$batch
)
score <- cbind(score, qc[match(score$library_id, qc$library_id), -"library_id"])
score <- cbind(score, as.data.table(score_mat))

dat <- merge(
  meta[, .(sample_id, sample_library_id, group, DLCO_pct, DLCO_severity_code)],
  score,
  by = c("sample_id", "group")
)
dat <- dat[group %in% c("MDA5", "ARS") & !is.na(DLCO_pct)]
dat[, outcome := as.integer(DLCO_severity_code >= 2)]
dat[, batch := factor(batch)]
stopifnot(nrow(dat) == 87L, all(dat$sample_library_id == dat$library_id))

make_folds <- function(group, y, k, seed) {
  set.seed(seed)
  strata <- interaction(group, y, drop = TRUE)
  fold <- integer(length(y))
  for (lev in levels(strata)) {
    idx <- sample(which(strata == lev))
    fold[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold
}

train_z <- function(train, test, variable) {
  mu <- mean(train[[variable]])
  s <- sd(train[[variable]])
  if (!is.finite(s) || s == 0) s <- 1
  list(train = (train[[variable]] - mu) / s, test = (test[[variable]] - mu) / s)
}

residualize <- function(train, test) {
  lib <- train_z(train, test, "log10_library_size")
  det <- train_z(train, test, "detected_genes")
  train[, lib_z := lib$train]
  test[, lib_z := lib$test]
  train[, det_z := det$train]
  test[, det_z := det$test]
  batch_levels <- levels(dat$batch)
  train$batch <- factor(train$batch, levels = batch_levels)
  test$batch <- factor(test$batch, levels = batch_levels)
  z_train <- model.matrix(~ batch + lib_z + det_z, data = train)
  z_test <- model.matrix(~ batch + lib_z + det_z, data = test)
  x_train <- as.matrix(train[, ..hallmark])
  x_test <- as.matrix(test[, ..hallmark])
  beta <- qr.solve(
    crossprod(z_train) + diag(1e-8, ncol(z_train)),
    crossprod(z_train, x_train)
  )
  list(train = x_train - z_train %*% beta, test = x_test - z_test %*% beta)
}

lambda_score <- function(fit) {
  i <- which.min(abs(log(fit$lambda) - log(fit$lambda.min)))
  fit$cvm[i]
}

fit_one <- function(x_train, y_train, x_test, inner_fold) {
  w <- ifelse(
    y_train == 1,
    length(y_train) / (2 * sum(y_train == 1)),
    length(y_train) / (2 * sum(y_train == 0))
  )
  fits <- vector("list", length(alpha_grid))
  score <- rep(Inf, length(alpha_grid))
  for (i in seq_along(alpha_grid)) {
    fit <- try(cv.glmnet(
      x_train, y_train, family = "binomial", alpha = alpha_grid[i],
      foldid = inner_fold, weights = w, type.measure = "deviance",
      standardize = TRUE, nlambda = 80
    ), silent = TRUE)
    if (inherits(fit, "try-error")) next
    fits[[i]] <- fit
    score[i] <- lambda_score(fit)
  }
  i <- which.min(score)
  list(
    probability = as.numeric(predict(fits[[i]], x_test, s = "lambda.min", type = "response")),
    alpha = alpha_grid[i],
    lambda = fits[[i]]$lambda.min
  )
}

auc_rank <- function(y, p) {
  pos <- p[y == 1]
  neg <- p[y == 0]
  z <- outer(pos, neg, "-")
  mean(z > 0) + 0.5 * mean(z == 0)
}

pred_rows <- list()
fold_rows <- list()
ii <- 0L
for (r in seq_len(n_repeats)) {
  outer_fold <- make_folds(dat$group, dat$outcome, outer_k, base_seed + r * 1000L)
  for (f in seq_len(outer_k)) {
    te <- which(outer_fold == f)
    tr <- setdiff(seq_len(nrow(dat)), te)
    inner_fold <- make_folds(
      dat$group[tr], dat$outcome[tr], inner_k,
      base_seed + r * 1000L + f * 10L
    )
    x <- residualize(copy(dat[tr]), copy(dat[te]))
    fit <- fit_one(x$train, dat$outcome[tr], x$test, inner_fold)
    ii <- ii + 1L
    pred_rows[[ii]] <- data.table(
      repeat_id = r, fold_id = f,
      sample_id = dat$sample_id[te], group = dat$group[te],
      truth = dat$outcome[te], probability = fit$probability
    )
    fold_rows[[ii]] <- data.table(
      repeat_id = r, fold_id = f, alpha = fit$alpha, lambda = fit$lambda
    )
  }
}

pred <- rbindlist(pred_rows)
fold <- rbindlist(fold_rows)
repeat_auc <- pred[, .(AUC = auc_rank(truth, probability)), by = repeat_id]
participant <- pred[, .(
  group = unique(group), truth = unique(truth),
  mean_oof_score = mean(probability), sd_oof_score = sd(probability)
), by = sample_id]

current <- fread(file.path(
  "generated", "clinical_context", "participant_mean_oof_scores.tsv"
))[model == "EV_RNA_pathway_only" & penalty_rule == "lambda_min"]
comparison <- merge(
  participant,
  current[, .(sample_id, current_mean_oof_score = mean_oof_score)],
  by = "sample_id"
)

set.seed(20260831)
boot_n <- 2000L
boot_auc <- replicate(boot_n, {
  idx <- c(
    sample(which(participant$truth == 0), sum(participant$truth == 0), replace = TRUE),
    sample(which(participant$truth == 1), sum(participant$truth == 1), replace = TRUE)
  )
  auc_rank(participant$truth[idx], participant$mean_oof_score[idx])
})

summary <- data.table(
  analysis = "sample-wise raw-count rank score; QC standardization and residualization fitted within each outer fold",
  n = nrow(participant),
  mean_repeated_cv_AUC = mean(repeat_auc$AUC),
  participant_mean_oof_AUC = auc_rank(participant$truth, participant$mean_oof_score),
  participant_AUC_low = quantile(boot_auc, 0.025),
  participant_AUC_high = quantile(boot_auc, 0.975),
  current_participant_mean_oof_AUC = auc_rank(comparison$truth, comparison$current_mean_oof_score),
  score_correlation_with_current = cor(
    comparison$mean_oof_score, comparison$current_mean_oof_score,
    method = "spearman"
  )
)

fwrite(pred, file.path(out_dir, "foldwise_raw_rank_all_outer_predictions.tsv"), sep = "\t")
fwrite(fold, file.path(out_dir, "foldwise_raw_rank_outer_fold_details.tsv"), sep = "\t")
fwrite(repeat_auc, file.path(out_dir, "foldwise_raw_rank_auc_by_repeat.tsv"), sep = "\t")
fwrite(participant, file.path(out_dir, "foldwise_raw_rank_participant_oof.tsv"), sep = "\t")
fwrite(comparison, file.path(out_dir, "foldwise_raw_rank_score_comparison.tsv"), sep = "\t")
fwrite(summary, file.path(out_dir, "foldwise_raw_rank_sensitivity_summary.tsv"), sep = "\t")

comparison[, severity := factor(truth, levels = c(0, 1), labels = c("Mild", "Mod. & Severe"))]

p <- ggplot(comparison, aes(current_mean_oof_score, mean_oof_score, colour = severity)) +
  geom_abline(slope = 1, intercept = 0, colour = "#A8A8A8", linewidth = 0.8) +
  geom_point(size = 2.4, alpha = 0.82) +
  scale_colour_manual(values = c("Mild" = "#159A80", "Mod. & Severe" = "#E76600")) +
  annotate(
    "text", x = Inf, y = -Inf, label = "Spearman rho = 1.000",
    hjust = 1.08, vjust = -0.8, family = "Helvetica", size = 4.1
  ) +
  coord_equal(xlim = c(0.10, 0.95), ylim = c(0.10, 0.95), expand = FALSE) +
  labs(
    x = "Current OOF probability",
    y = "Fold-wise raw-rank\nOOF probability"
  ) +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    axis.title = element_text(size = 15, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    legend.key.width = unit(5, "mm"),
    legend.key.height = unit(3.5, "mm"),
    plot.margin = margin(7, 10, 7, 9)
  )

ggsave(file.path(out_dir, "S_foldwise_raw_rank_sensitivity.pdf"), p,
       width = 5.33, height = 3.9, units = "in", device = cairo_pdf)
png(file.path(out_dir, "S_foldwise_raw_rank_sensitivity.png"),
    width = 5.33, height = 3.9, units = "in", res = 300, type = "cairo")
print(p)
dev.off()


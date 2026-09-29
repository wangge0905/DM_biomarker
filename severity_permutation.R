# RStudio工作目录设为本文件夹
library(data.table)
library(glmnet)
library(ggplot2)
library(parallel)


dm <- "inputs"
out_dir <- file.path(
  "generated",
  "DM_severity_reviewer_completion_20260830", "output"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_perm <- 200L
n_core <- 4L
base_seed <- 20260830L
outer_k <- 5L
inner_k <- 4L
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

dat <- merge(
  meta[, .(
    sample_id, sample_library_id, group, DLCO_pct, DLCO_severity_code
  )],
  score[, c(
    "sample_id", "library_id", "batch", "log10_library_size_z", "detected_genes_z",
    hallmark
  ), with = FALSE],
  by = "sample_id"
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

permute_within_group <- function(y, group, seed) {
  set.seed(seed)
  ans <- y
  for (g in unique(group)) {
    idx <- which(group == g)
    ans[idx] <- sample(y[idx])
  }
  ans
}

residualize <- function(train, test) {
  batch_levels <- levels(dat$batch)
  train$batch <- factor(train$batch, levels = batch_levels)
  test$batch <- factor(test$batch, levels = batch_levels)
  z_train <- model.matrix(~ batch + log10_library_size_z + detected_genes_z, data = train)
  z_test <- model.matrix(~ batch + log10_library_size_z + detected_genes_z, data = test)
  x_train <- as.matrix(train[, ..hallmark])
  x_test <- as.matrix(test[, ..hallmark])
  beta <- qr.solve(
    crossprod(z_train) + diag(1e-8, ncol(z_train)),
    crossprod(z_train, x_train)
  )
  list(train = x_train - z_train %*% beta, test = x_test - z_test %*% beta)
}

lambda_score <- function(fit) {
  idx <- which.min(abs(log(fit$lambda) - log(fit$lambda.min)))
  fit$cvm[idx]
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
      x_train, y_train,
      family = "binomial", alpha = alpha_grid[i], foldid = inner_fold,
      weights = w, type.measure = "deviance", standardize = TRUE, nlambda = 80
    ), silent = TRUE)
    if (inherits(fit, "try-error")) next
    fits[[i]] <- fit
    score[i] <- lambda_score(fit)
  }
  i <- which.min(score)
  as.numeric(predict(fits[[i]], x_test, s = "lambda.min", type = "response"))
}

auc_rank <- function(y, p) {
  pos <- p[y == 1]
  neg <- p[y == 0]
  z <- outer(pos, neg, "-")
  mean(z > 0) + 0.5 * mean(z == 0)
}

run_nested <- function(y, seed) {
  outer_fold <- make_folds(dat$group, y, outer_k, seed)
  p <- rep(NA_real_, nrow(dat))
  for (f in seq_len(outer_k)) {
    te <- which(outer_fold == f)
    tr <- setdiff(seq_len(nrow(dat)), te)
    inner_fold <- make_folds(dat$group[tr], y[tr], inner_k, seed + f * 100L)
    x <- residualize(dat[tr], dat[te])
    p[te] <- fit_one(x$train, y[tr], x$test, inner_fold)
  }
  auc_rank(y, p)
}

observed_auc <- run_nested(dat$outcome, base_seed)
cat("Observed 5-fold nested-CV AUC:", observed_auc, "\n")

perm_auc <- unlist(mclapply(seq_len(n_perm), function(i) {
  y_perm <- permute_within_group(dat$outcome, dat$group, base_seed + i * 10000L)
  a <- run_nested(y_perm, base_seed + i * 10000L + 500L)
  a
}, mc.cores = n_core, mc.preschedule = FALSE))

empirical_p <- (1 + sum(perm_auc >= observed_auc)) / (n_perm + 1)
result <- data.table(
  permutation = 0:n_perm,
  AUC = c(observed_auc, perm_auc),
  observed = c(TRUE, rep(FALSE, n_perm))
)
summary <- data.table(
  analysis = "EV-RNA Hallmark severity model",
  permutation_scheme = "severity labels permuted within antibody subtype; all outer folds, inner tuning, residualization and model fitting repeated",
  observed_nested_cv_AUC = observed_auc,
  permutations = n_perm,
  null_AUC_mean = mean(perm_auc),
  null_AUC_sd = sd(perm_auc),
  null_AUC_95_low = quantile(perm_auc, 0.025),
  null_AUC_95_high = quantile(perm_auc, 0.975),
  empirical_P = empirical_p
)
fwrite(result, file.path(out_dir, "severity_full_pipeline_permutation_auc.tsv"), sep = "\t")
fwrite(summary, file.path(out_dir, "severity_full_pipeline_permutation_summary.tsv"), sep = "\t")

p <- ggplot(data.table(AUC = perm_auc), aes(AUC)) +
  geom_histogram(binwidth = 0.025, boundary = 0.5, fill = "#B7B7B7", colour = "white") +
  geom_vline(xintercept = observed_auc, colour = "#C64B40", linewidth = 1.15) +
  annotate(
    "text", x = observed_auc, y = Inf,
    label = sprintf("Observed AUC = %.3f\nEmpirical P = %.3f", observed_auc, empirical_p),
    hjust = 1.05, vjust = 1.25, family = "Helvetica", size = 5.2
  ) +
  labs(x = "Nested-CV AUC under permuted labels", y = "Number of permutations") +
  theme_classic(base_family = "Helvetica", base_size = 17) +
  theme(
    axis.title = element_text(size = 19, colour = "black"),
    axis.text = element_text(size = 16, colour = "black"),
    plot.margin = margin(8, 10, 8, 8)
  )

ggsave(file.path(out_dir, "S_severity_full_pipeline_permutation.pdf"), p,
       width = 5.5, height = 4.4, units = "in", device = cairo_pdf)
png(file.path(out_dir, "S_severity_full_pipeline_permutation.png"),
    width = 5.5, height = 4.4, units = "in", res = 300, type = "cairo")
print(p)
dev.off()


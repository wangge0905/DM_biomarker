# RStudio工作目录设为本文件夹
library(data.table)
library(edgeR)
library(e1071)
library(parallel)


options(stringsAsFactors = FALSE)

dm <- "inputs"
root <- "generated/biomarker"
out_dir <- file.path(root, "results", "permutation_MDA5_vs_ARS")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_perm_fixed <- 1000L
n_perm_nested <- 200L
n_core <- 4L
n_fold <- 5L
max_panel <- 8L
cost_grid <- c(0.01, 0.1, 1, 10, 100)
fdr_cutoff <- 0.05
logfc_cutoff <- 1

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

log_cpm <- function(count, library_size = colSums(count)) {
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

permute_within <- function(y, strata, seed) {
  set.seed(seed)
  out <- y
  for (s in unique(strata)) {
    idx <- which(strata == s)
    out[idx] <- sample(y[idx], length(idx), replace = FALSE)
  }
  out
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
  rownames(tab) <- NULL
  tab
}

manifest <- fread(file.path(root, "meta", "global_discovery_validation_manifest.tsv"))
meta <- manifest[group %in% c("MDA5", "ARS")]
meta[, y := as.integer(group == "MDA5")]

raw <- fread(file.path(dm, "1204", "gencode.txt"), data.table = FALSE,
             check.names = FALSE)
feature_id <- raw[[1]]
raw <- raw[, -1, drop = FALSE]
raw <- raw[, -c(40:44, 50:54), drop = FALSE]
colnames(raw)[35:44] <- sub("-GP", "", colnames(raw)[35:44], fixed = TRUE)
rownames(raw) <- feature_id
raw <- raw[, meta$library_id, drop = FALSE]

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
count <- as.matrix(raw[keep, , drop = FALSE])
storage.mode(count) <- "double"
rm(raw)

discovery_idx <- which(meta$analysis_set == "discovery")
validation_idx <- which(meta$analysis_set == "internal_validation")
discovery_meta <- meta[discovery_idx]
validation_meta <- meta[validation_idx]
discovery_count <- count[, discovery_idx, drop = FALSE]
validation_count <- count[, validation_idx, drop = FALSE]

panel <- fread(file.path(
  root, "results", "sensitivity_delta_auc_0.01", "MDA5_vs_ARS",
  "fixed_biomarker_panel.tsv"
))
fixed_feature <- panel[order(sensitivity_panel_rank)]$feature_id
selection <- fread(file.path(
  root, "results", "sensitivity_delta_auc_0.01", "MDA5_vs_ARS",
  "selection_rule.tsv"
))
fixed_cost <- selection$selected_cost[1]

x_discovery_fixed <- t(log_cpm(
  discovery_count[fixed_feature, , drop = FALSE], colSums(discovery_count)
))
x_validation_fixed <- t(log_cpm(
  validation_count[fixed_feature, , drop = FALSE], colSums(validation_count)
))
scaled_fixed <- scale_train_test(x_discovery_fixed, x_validation_fixed)

fixed_stat <- function(i, observed = FALSE) {
  if (observed) {
    y_discovery <- discovery_meta$y
    y_validation <- validation_meta$y
  } else {
    y_discovery <- permute_within(
      discovery_meta$y, discovery_meta$batch, 1000000L + i
    )
    y_validation <- permute_within(
      validation_meta$y, validation_meta$batch, 2000000L + i
    )
  }
  score <- fit_svm_probability(
    scaled_fixed$train, y_discovery, scaled_fixed$test,
    fixed_cost, 3000000L + i
  )
  data.frame(
    permutation = if (observed) 0L else i,
    AUC = auc_value(y_validation, score),
    stringsAsFactors = FALSE
  )
}


observed_fixed <- fixed_stat(0L, observed = TRUE)
fixed_null <- rbindlist(mclapply(
  seq_len(n_perm_fixed), fixed_stat, mc.cores = n_core,
  mc.preschedule = TRUE
))
fixed_all <- rbind(observed_fixed, fixed_null)
fixed_p <- (1 + sum(fixed_null$AUC >= observed_fixed$AUC, na.rm = TRUE)) /
  (n_perm_fixed + 1)
fwrite(fixed_all, file.path(out_dir, "fixed_panel_permutation_auc.tsv"), sep = "\t")

nested_stat <- function(i, observed = FALSE) {
  if (observed) {
    y <- discovery_meta$y
  } else {
    y <- permute_within(discovery_meta$y, discovery_meta$batch, 4000000L + i)
  }
  foldid <- make_foldid(
    paste(discovery_meta$batch, y, sep = "_"), n_fold, 5000000L + i
  )
  pred <- vector("list", n_fold * max_panel * length(cost_grid))
  z <- 0L
  for (fold in seq_len(n_fold)) {
    train <- which(foldid != fold)
    test <- which(foldid == fold)
    de <- rank_de(discovery_count[, train, drop = FALSE], y[train])
    selected <- head(de$feature_id, max_panel)
    xtr_all <- log_cpm(
      discovery_count[selected, train, drop = FALSE],
      colSums(discovery_count[, train, drop = FALSE])
    )
    xte_all <- log_cpm(
      discovery_count[selected, test, drop = FALSE],
      colSums(discovery_count[, test, drop = FALSE])
    )
    for (k in seq_len(max_panel)) {
      use <- selected[seq_len(k)]
      ss <- scale_train_test(
        t(xtr_all[use, , drop = FALSE]),
        t(xte_all[use, , drop = FALSE])
      )
      for (cost in cost_grid) {
        z <- z + 1L
        score <- fit_svm_probability(
          ss$train, y[train], ss$test, cost,
          6000000L + i * 1000L + z
        )
        pred[[z]] <- data.frame(
          k = k, cost = cost, truth = y[test], score = score
        )
      }
    }
  }
  pred <- rbindlist(pred)
  stat <- pred[, .(AUC = auc_value(truth, score)), by = .(k, cost)]
  best <- stat[order(-AUC, k, cost)][1]
  data.frame(
    permutation = if (observed) 0L else i,
    AUC = best$AUC,
    selected_k = best$k,
    selected_cost = best$cost,
    stringsAsFactors = FALSE
  )
}


observed_nested <- nested_stat(0L, observed = TRUE)
nested_null <- rbindlist(mclapply(
  seq_len(n_perm_nested), nested_stat, mc.cores = n_core,
  mc.preschedule = TRUE
))
nested_all <- rbind(observed_nested, nested_null)
nested_p <- (1 + sum(nested_null$AUC >= observed_nested$AUC, na.rm = TRUE)) /
  (n_perm_nested + 1)
fwrite(nested_all, file.path(out_dir, "full_pipeline_permutation_auc.tsv"), sep = "\t")

summary <- data.frame(
  test = c(
    "locked_6_feature_panel_internal_validation",
    "full_discovery_pipeline_5_fold_sensitivity"
  ),
  permutation_scheme = c(
    "labels permuted within batch separately in discovery and internal validation; fixed features and cost; model refitted",
    "discovery labels permuted within batch; DE feature selection and k/cost tuning repeated inside 5-fold CV"
  ),
  observed_AUC = c(observed_fixed$AUC, observed_nested$AUC),
  permutations = c(n_perm_fixed, n_perm_nested),
  null_AUC_mean = c(mean(fixed_null$AUC), mean(nested_null$AUC)),
  null_AUC_sd = c(sd(fixed_null$AUC), sd(nested_null$AUC)),
  null_AUC_95_low = c(
    quantile(fixed_null$AUC, 0.025), quantile(nested_null$AUC, 0.025)
  ),
  null_AUC_95_high = c(
    quantile(fixed_null$AUC, 0.975), quantile(nested_null$AUC, 0.975)
  ),
  empirical_P = c(fixed_p, nested_p),
  stringsAsFactors = FALSE
)
fwrite(summary, file.path(out_dir, "permutation_test_summary.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)

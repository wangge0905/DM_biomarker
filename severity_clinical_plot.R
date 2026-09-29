# RStudio工作目录设为本文件夹
# 先运行severity_pathway.R
library(data.table)
library(ggplot2)


root <- "inputs"
input_dir <- file.path("generated", "clinical_context")
out_dir <- file.path(
  "generated",
  "DM_severity_reviewer_completion_20260830", "output"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260830)
boot_n <- 2000L

score <- fread(file.path(input_dir, "participant_mean_oof_scores.tsv"))
score <- score[penalty_rule == "lambda_min"]

model_order <- c("autoantibody_only", "clinical_only", "EV_RNA_pathway_only", "combined")
model_label <- c(
  autoantibody_only = "Antibody",
  clinical_only = "Clinical",
  EV_RNA_pathway_only = "EV-RNA",
  combined = "Combined"
)
model_color <- c(
  "Antibody" = "#4B4B4B",
  "Clinical" = "#7F7F7F",
  "EV-RNA" = "#3B83BD",
  "Combined" = "#2E665C"
)

auc_rank <- function(y, p) {
  pos <- p[y == 1]
  neg <- p[y == 0]
  z <- outer(pos, neg, "-")
  mean(z > 0) + 0.5 * mean(z == 0)
}

pr_points <- function(y, p) {
  o <- order(p, decreasing = TRUE)
  y <- y[o]
  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)
  data.table(
    recall = c(0, tp / sum(y == 1)),
    precision = c(1, tp / (tp + fp))
  )
}

average_precision <- function(y, p) {
  o <- order(p, decreasing = TRUE)
  y <- y[o]
  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)
  recall <- tp / sum(y == 1)
  precision <- tp / (tp + fp)
  sum(c(recall[1], diff(recall)) * precision)
}

roc_points <- function(y, p) {
  cutoff <- c(Inf, sort(unique(p), decreasing = TRUE), -Inf)
  rbindlist(lapply(cutoff, function(z) {
    pred <- p >= z
    data.table(
      fpr = sum(pred & y == 0) / sum(y == 0),
      tpr = sum(pred & y == 1) / sum(y == 1)
    )
  }))
}

stratified_index <- function(y) {
  c(
    sample(which(y == 0), sum(y == 0), replace = TRUE),
    sample(which(y == 1), sum(y == 1), replace = TRUE)
  )
}

metric_rows <- list()
pr_rows <- list()
roc_rows <- list()
cal_rows <- list()

for (m in model_order) {
  z <- score[model == m]
  stopifnot(nrow(z) == 87L, !anyDuplicated(z$sample_id))
  y <- z$truth
  p <- z$mean_oof_score
  ap_boot <- replicate(boot_n, {
    idx <- stratified_index(y)
    average_precision(y[idx], p[idx])
  })
  auc_boot <- replicate(boot_n, {
    idx <- stratified_index(y)
    auc_rank(y[idx], p[idx])
  })
  metric_rows[[m]] <- data.table(
    model = m,
    label = model_label[m],
    n = length(y),
    positive_n = sum(y == 1),
    prevalence = mean(y),
    AUC = auc_rank(y, p),
    AUC_low = quantile(auc_boot, 0.025),
    AUC_high = quantile(auc_boot, 0.975),
    PR_AUC = average_precision(y, p),
    PR_AUC_low = quantile(ap_boot, 0.025),
    PR_AUC_high = quantile(ap_boot, 0.975),
    Brier = mean((p - y)^2)
  )
  pr_rows[[m]] <- cbind(pr_points(y, p), model = model_label[m])
  roc_rows[[m]] <- cbind(roc_points(y, p), model = model_label[m])

  # 每组至少有8例，避免过细分箱造成不稳定
  q <- unique(quantile(p, probs = seq(0, 1, length.out = 7), na.rm = TRUE))
  bin <- cut(p, breaks = q, include.lowest = TRUE, labels = FALSE)
  cal_rows[[m]] <- data.table(y = y, p = p, bin = bin)[, .(
    mean_probability = mean(p),
    observed_proportion = mean(y),
    bin_n = .N
  ), by = bin][, model := model_label[m]]
}

metric <- rbindlist(metric_rows)
pr <- rbindlist(pr_rows)
roc <- rbindlist(roc_rows)
cal <- rbindlist(cal_rows)

paired_comparison <- rbindlist(lapply(list(
  c("EV_RNA_pathway_only", "autoantibody_only"),
  c("EV_RNA_pathway_only", "clinical_only"),
  c("combined", "autoantibody_only"),
  c("combined", "clinical_only"),
  c("combined", "EV_RNA_pathway_only")
), function(pair) {
  a <- score[model == pair[1], .(sample_id, truth, p_a = mean_oof_score)]
  b <- score[model == pair[2], .(sample_id, p_b = mean_oof_score)]
  z <- merge(a, b, by = "sample_id")
  boot <- replicate(boot_n, {
    idx <- stratified_index(z$truth)
    average_precision(z$truth[idx], z$p_a[idx]) -
      average_precision(z$truth[idx], z$p_b[idx])
  })
  data.table(
    model_a = pair[1], model_b = pair[2],
    PR_AUC_difference_a_minus_b =
      average_precision(z$truth, z$p_a) - average_precision(z$truth, z$p_b),
    CI_low = quantile(boot, 0.025),
    CI_high = quantile(boot, 0.975),
    bootstrap_iterations = boot_n
  )
}))

thresholds <- seq(0.05, 0.80, by = 0.01)
dca <- rbindlist(lapply(model_order, function(m) {
  z <- score[model == m]
  rbindlist(lapply(thresholds, function(pt) {
    pred <- z$mean_oof_score >= pt
    n <- nrow(z)
    tp <- sum(pred & z$truth == 1)
    fp <- sum(pred & z$truth == 0)
    data.table(
      threshold = pt,
      net_benefit = tp / n - fp / n * pt / (1 - pt),
      model = model_label[m]
    )
  }))
}))
prevalence <- unique(metric$prevalence)
dca_reference <- rbind(
  data.table(threshold = thresholds, net_benefit = 0, model = "Treat none"),
  data.table(
    threshold = thresholds,
    net_benefit = prevalence - (1 - prevalence) * thresholds / (1 - thresholds),
    model = "Treat all"
  )
)

fwrite(metric, file.path(out_dir, "severity_pr_auc_auc_brier.tsv"), sep = "\t")
fwrite(paired_comparison, file.path(out_dir, "severity_paired_pr_auc_comparisons.tsv"), sep = "\t")
fwrite(pr, file.path(out_dir, "severity_pr_curve_points.tsv"), sep = "\t")
fwrite(roc, file.path(out_dir, "severity_roc_clinical_context_points.tsv"), sep = "\t")
fwrite(cal, file.path(out_dir, "severity_calibration_points.tsv"), sep = "\t")
fwrite(rbind(dca, dca_reference), file.path(out_dir, "severity_decision_curve_points.tsv"), sep = "\t")

theme_original <- function(base_size = 14) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.title = element_text(size = base_size + 1, colour = "black"),
      axis.text = element_text(size = base_size, colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1),
      plot.margin = margin(7, 10, 7, 9)
    )
}

metric_label <- setNames(
  sprintf("%s: %.3f", metric$label, metric$PR_AUC), metric$label
)
p_pr <- ggplot(pr, aes(recall, precision, colour = model)) +
  geom_step(linewidth = 1.15, direction = "vh") +
  geom_hline(yintercept = prevalence, linetype = 2, colour = "#A8A8A8", linewidth = 0.7) +
  scale_colour_manual(values = model_color, labels = metric_label) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(x = "Recall", y = "Precision") + theme_original(14) +
  theme(
    legend.position = c(0.02, 0.03),
    legend.justification = c(0, 0),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key.width = unit(5, "mm"),
    legend.key.height = unit(3.5, "mm")
  )

auc_label <- setNames(
  sprintf("%s: %.3f", metric$label, metric$AUC), metric$label
)
p_roc <- ggplot(roc, aes(fpr, tpr, colour = model)) +
  geom_abline(slope = 1, intercept = 0, colour = "black", linewidth = 0.6) +
  geom_step(linewidth = 1.15, direction = "vh") +
  scale_colour_manual(values = model_color, labels = auc_label) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(x = "False positive rate", y = "True positive rate") + theme_original(14) +
  theme(
    legend.position = c(0.98, 0.03),
    legend.justification = c(1, 0),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key.width = unit(5, "mm"),
    legend.key.height = unit(3.5, "mm")
  )

p_cal <- ggplot(cal, aes(mean_probability, observed_proportion, colour = model)) +
  geom_abline(slope = 1, intercept = 0, colour = "#A8A8A8", linewidth = 0.7) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.7) +
  scale_colour_manual(values = model_color) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(x = "Mean predicted probability", y = "Observed proportion") + theme_original(14) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key.width = unit(5, "mm"),
    legend.key.height = unit(3.5, "mm")
  )

dca_color <- c(model_color, "Treat all" = "#B8B8B8", "Treat none" = "#333333")
dca_all <- rbind(dca, dca_reference)
p_dca <- ggplot(dca_all, aes(threshold, net_benefit, colour = model, linetype = model)) +
  geom_line(linewidth = 1.0) +
  scale_colour_manual(values = dca_color) +
  scale_linetype_manual(values = c(
    "Antibody" = 1, "Clinical" = 1, "EV-RNA" = 1, "Combined" = 1,
    "Treat all" = 2, "Treat none" = 3
  )) +
  coord_cartesian(xlim = c(0.05, 0.80), ylim = c(-0.10, 0.35), expand = FALSE) +
  labs(x = "Threshold probability", y = "Net benefit") + theme_original(14) +
  guides(
    colour = guide_legend(ncol = 1, override.aes = list(linewidth = 1)),
    linetype = "none"
  ) +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key.width = unit(5, "mm"),
    legend.key.height = unit(3.5, "mm")
  )

save_plot <- function(p, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p,
         width = width, height = height, units = "in", device = cairo_pdf)
  png(file.path(out_dir, paste0(stem, ".png")),
      width = width, height = height, units = "in", res = 300, type = "cairo")
  print(p)
  dev.off()
}

save_plot(p_pr, "S_severity_PR_curve", 5.33, 3.9)
save_plot(p_roc, "S_severity_ROC_clinical_context", 5.33, 3.9)
save_plot(p_cal, "S_severity_calibration", 5.33, 3.9)
save_plot(p_dca, "S_severity_decision_curve", 5.33, 3.9)


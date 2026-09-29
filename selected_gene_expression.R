# RStudio工作目录设为本文件夹
library(data.table)
library(ggplot2)
library(patchwork)


options(stringsAsFactors = FALSE)

dm <- "inputs"
root <- "generated/biomarker"
result_root <- file.path(root, "results", "sensitivity_delta_auc_0.01")
out_dir <- "generated/biomarker/figures/fixed_panel_DEF_20260830"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

col_mda5 <- "#B36A1B"
col_ars <- "#E0C377"
col_hc <- "#078E78"
col_development <- "#E31A1C"
col_validation <- "#2C7FB8"
font_family <- "Helvetica"

theme_manuscript <- function(base_size = 13) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(size = base_size),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      plot.margin = margin(4, 5, 4, 5)
    )
}

save_plot <- function(p, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p,
         width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(stem, ".png")), p,
         width = width, height = height, units = "in", dpi = 300,
         bg = "white", device = grDevices::png, type = "cairo")
}

auc_value <- function(y, score) {
  r <- rank(score, ties.method = "average")
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

roc_frame <- function(y, score) {
  threshold <- c(Inf, sort(unique(score), decreasing = TRUE), -Inf)
  out <- rbindlist(lapply(threshold, function(th) {
    pred <- as.integer(score >= th)
    data.frame(
      FPR = mean(pred[y == 0] == 1),
      TPR = mean(pred[y == 1] == 1)
    )
  }))
  unique(out[order(FPR, TPR)])
}

raw <- fread(file.path(dm, "1204", "gencode.txt"), data.table = FALSE,
             check.names = FALSE)
feature_id <- raw[[1]]
raw <- raw[, -1, drop = FALSE]
raw <- raw[, -c(40:44, 50:54), drop = FALSE]
colnames(raw)[35:44] <- sub("-GP", "", colnames(raw)[35:44], fixed = TRUE)
rownames(raw) <- feature_id
meta <- fread(file.path(root, "meta", "global_discovery_validation_manifest.tsv"))
raw <- raw[, meta$library_id, drop = FALSE]
log_expr <- log2(sweep(as.matrix(raw), 2, colSums(raw) / 1e6, "/") + 1)
rm(raw)

comparison_info <- list(
  MDA5_vs_HC = list(
    positive = "MDA5", negative = "HC",
    positive_label = "Anti-MDA5+", negative_label = "HCs",
    colors = c("Anti-MDA5+" = col_mda5, "HCs" = col_hc),
    box_features = c("MIR146B", "MIR3142HG", "MIR223", "DNM3OS"),
    box_width = 4.2, box_height = 4.7, box_ncol = 2
  ),
  ARS_vs_HC = list(
    positive = "ARS", negative = "HC",
    positive_label = "Anti-ARS+", negative_label = "HCs",
    colors = c("Anti-ARS+" = col_ars, "HCs" = col_hc),
    box_features = c("MEX3D", "C19orf81"),
    box_width = 3.7, box_height = 3.0, box_ncol = 2
  ),
  MDA5_vs_ARS = list(
    positive = "MDA5", negative = "ARS",
    positive_label = "Anti-MDA5+", negative_label = "Anti-ARS+",
    colors = c("Anti-MDA5+" = col_mda5, "Anti-ARS+" = col_ars),
    # The final PDF template repositions cropped instances to put CCN1 above CXCL8.
    box_features = c("EGR1", "H2BC17", "CXCL8", "CCN1", "C2CD4B", "MAMLD1"),
    box_width = 6.2, box_height = 4.7, box_ncol = 3
  )
)

display_group <- function(x) {
  ifelse(
    x == "Anti-MDA5+", "Anti-\nMDA5+",
    ifelse(x == "Anti-ARS+", "Anti-\nARS+", x)
  )
}

make_boxplot <- function(comparison, info) {
  panel <- fread(file.path(result_root, comparison, "fixed_biomarker_panel.tsv"))
  panel <- panel[gene_symbol %in% info$box_features]
  panel <- panel[match(info$box_features, gene_symbol)]
  zmeta <- meta[group %in% c(info$positive, info$negative)]
  plot_list <- lapply(seq_len(nrow(panel)), function(i) {
    z <- data.frame(
      group = ifelse(
        zmeta$group == info$positive,
        info$positive_label,
        info$negative_label
      ),
      value = as.numeric(log_expr[panel$feature_id[i], zmeta$library_id])
    )
    z$group <- factor(z$group, levels = c(info$positive_label, info$negative_label))
    fdr <- panel$FDR[i]
    star <- if (fdr < 0.001) "***" else if (fdr < 0.01) "**" else if (fdr < 0.05) "*" else "NS"
    ymax <- max(z$value, na.rm = TRUE)
    ymin <- min(z$value, na.rm = TRUE)
    span <- max(1, ymax - ymin)
    ggplot(z, aes(group, value, fill = group)) +
      geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.55) +
      geom_point(position = position_jitter(width = 0.14, height = 0, seed = 20260830), size = 1.35, colour = "black") +
      annotate("text", x = 1.5, y = ymax + 0.12 * span,
               label = star, size = 4.5, family = font_family) +
      scale_fill_manual(values = info$colors, guide = "none") +
      scale_x_discrete(labels = display_group) +
      scale_y_continuous(
        name = "log(CPM+1)",
        expand = expansion(mult = c(0.03, 0.17))
      ) +
      labs(x = NULL, title = ifelse(grepl("^MIR[0-9]+[A-Z]*$", panel$gene_symbol[i]) & panel$gene_symbol[i] != "MIR3142HG", sub("^mir", "hsa-mir-", tolower(panel$gene_symbol[i])), panel$gene_symbol[i])) +
      theme_manuscript(12) +
      theme(
        axis.text.x = element_text(size = 10.5, lineheight = 0.9),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(size = 13, face = "bold")
      )
  })
  wrap_plots(plot_list, ncol = info$box_ncol) + plot_layout(guides = "collect")
}

confusion_plot <- function(pred, threshold, info, title, palette) {
  pred <- copy(pred)
  if (!"predicted" %in% names(pred)) pred[, predicted := as.integer(score >= threshold)]
  cm <- as.data.table(table(
    truth = factor(pred$truth, levels = c(1, 0)),
    predicted = factor(pred$predicted, levels = c(1, 0))
  ))
  cm[, truth_label := factor(
    display_group(ifelse(truth == 1, info$positive_label, info$negative_label)),
    levels = display_group(c(info$negative_label, info$positive_label))
  )]
  cm[, pred_label := factor(
    display_group(ifelse(predicted == 1, info$positive_label, info$negative_label)),
    levels = display_group(c(info$positive_label, info$negative_label))
  )]
  sensitivity <- mean(pred$predicted[pred$truth == 1] == 1)
  specificity <- mean(pred$predicted[pred$truth == 0] == 0)
  accuracy <- mean(pred$predicted == pred$truth)
  metric <- sprintf(
    "Sensitivity = %.2f%%\nSpecificity = %.2f%%\nAccuracy = %.2f%%",
    100 * sensitivity, 100 * specificity, 100 * accuracy
  )
  ggplot(cm, aes(pred_label, truth_label, fill = N)) +
    geom_tile(colour = "white", linewidth = 1.3) +
    geom_text(aes(label = N), size = 5.1, family = font_family) +
    scale_fill_gradient(low = palette[1], high = palette[2], guide = "none") +
    coord_fixed(clip = "off") +
    labs(x = "Predicted class", y = "True class", title = title,
         caption = metric) +
    theme_void(base_family = font_family, base_size = 13) +
    theme(
      plot.title = element_text(size = 13.5, face = "bold", hjust = 0.5,
                                margin = margin(b = 5)),
      axis.text.x = element_text(size = 12, colour = "black", lineheight = 0.88),
      axis.text.y = element_text(size = 12, colour = "black", lineheight = 0.88),
      axis.title.x = element_text(size = 12.5, face = "bold", margin = margin(t = 6)),
      axis.title.y = element_text(size = 12.5, face = "bold", angle = 90,
                                  margin = margin(r = 6)),
      plot.caption = element_text(size = 10.5, hjust = 0.5, lineheight = 0.9,
                                  margin = margin(t = 8)),
      plot.margin = margin(5, 7, 5, 7)
    )
}

make_confusion_pair <- function(comparison, info) {
  dev <- fread(file.path(result_root, comparison, "development_fixed_panel_predictions.tsv"))
  val <- fread(file.path(result_root, comparison, "internal_validation_predictions.tsv"))
  threshold <- fread(file.path(result_root, comparison, "development_threshold.tsv"))$threshold[1]
  dev[, predicted := as.integer(score >= threshold)]
  p1 <- confusion_plot(dev, threshold, info, "Training set", c("#FDE0DD", "#B30000"))
  p2 <- confusion_plot(val, threshold, info, "Internal test set", c("#EFF3FF", "#08519C"))
  p1 + p2 + plot_layout(ncol = 2)
}

make_roc <- function(comparison, info) {
  dev <- fread(file.path(result_root, comparison, "development_fixed_panel_predictions.tsv"))
  val <- fread(file.path(result_root, comparison, "internal_validation_predictions.tsv"))
  perf <- fread(file.path(result_root, comparison, "fixed_panel_performance.tsv"))
  roc_dev <- roc_frame(dev$truth, dev$score)
  roc_dev$set <- "Training"
  roc_val <- roc_frame(val$truth, val$score)
  roc_val$set <- "Internal test"
  roc <- rbind(roc_dev, roc_val)
  roc$set <- factor(roc$set, levels = c("Training", "Internal test"))
  d <- perf[set == "development_fixed_panel_resampling"]
  v <- perf[set == "internal_validation_fixed_panel"]
  label_d <- sprintf(
    "Training: %.3f (%.3f-%.3f)",
    d$AUC, d$AUC_low, d$AUC_high
  )
  label_v <- sprintf(
    "Internal test: %.3f (%.3f-%.3f)",
    v$AUC, v$AUC_low, v$AUC_high
  )
  ggplot(roc, aes(FPR, TPR, colour = set)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.6, colour = "black") +
    geom_step(linewidth = 1.25, direction = "vh") +
    scale_colour_manual(
      values = c("Training" = col_development, "Internal test" = col_validation),
      breaks = c("Training", "Internal test"),
      labels = c(label_d, label_v),
      name = "AUC (95% CI)"
    ) +
    guides(colour = guide_legend(ncol = 1, byrow = TRUE,
                                 override.aes = list(linewidth = 1.25))) +
    coord_equal(xlim = c(0, 1.02), ylim = c(0, 1.02), clip = "off") +
    scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    labs(x = "False positive rate", y = "True positive rate") +
    theme_manuscript(13) +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12.5),
      legend.position = "top",
      legend.justification = "left",
      legend.direction = "vertical",
      legend.title = element_text(size = 10.5, face = "bold"),
      legend.text = element_text(size = 10.5),
      legend.key.width = grid::unit(18, "pt"),
      legend.key.height = grid::unit(8, "pt"),
      legend.margin = margin(0, 0, 3, 0),
      plot.margin = margin(3, 9, 5, 5)
    )
}

for (comparison in names(comparison_info)) {
  info <- comparison_info[[comparison]]
  p_d <- make_boxplot(comparison, info)
  p_e <- make_confusion_pair(comparison, info)
  p_f <- make_roc(comparison, info)
  save_plot(p_d, paste0(comparison, "_D_boxplots"), info$box_width, info$box_height)
  save_plot(p_e, paste0(comparison, "_E_confusion_development_internal"), 6.0, 3.2)
  save_plot(p_f, paste0(comparison, "_F_ROC_development_internal"), 4.35, 3.5)
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

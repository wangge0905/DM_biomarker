# RStudio工作目录设为本文件夹
library(data.table)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
library(gtable)
library(grid)


options(stringsAsFactors = FALSE)
set.seed(20260828)

dm <- "inputs"
root <- "generated/figure6"
main_dir <- file.path(root, "figures", "main_redesigned")
supp_dir <- file.path(root, "figures", "supplementary_redesigned")
data_dir <- file.path(root, "data_used", "fig5_redesigned")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Helvetica"
brbg <- brewer.pal(5, "BrBG")
group_color <- c(MDA5 = brbg[1], ARS = brbg[2])
severity_color <- c(Mild = brewer.pal(5, "Dark2")[1],
                    `Mod. & Severe` = brewer.pal(5, "Dark2")[2])
model_color <- c(`Transcript-level` = "#777777", `Pathway-level` = "#377EB8")

save_plot <- function(p, stem, width, height, folder = main_dir) {
  pdf_file <- file.path(folder, paste0(stem, ".pdf"))
  png_file <- file.path(folder, paste0(stem, ".png"))
  ggsave(pdf_file, p, width = width, height = height, units = "in",
         device = cairo_pdf, bg = "white")
  png(png_file, width = width, height = height, units = "in", res = 300,
      type = "cairo", bg = "white")
  print(p)
  dev.off()
  data.table(stem, pdf = pdf_file, png = png_file,
             width_in = width, height_in = height)
}

save_grid <- function(draw_fun, stem, width, height, folder = main_dir) {
  pdf_file <- file.path(folder, paste0(stem, ".pdf"))
  png_file <- file.path(folder, paste0(stem, ".png"))
  cairo_pdf(pdf_file, width = width, height = height, family = font_family)
  draw_fun()
  dev.off()
  png(png_file, width = width, height = height, units = "in", res = 300,
      type = "cairo", bg = "white")
  draw_fun()
  dev.off()
  data.table(stem, pdf = pdf_file, png = png_file,
             width_in = width, height_in = height)
}

theme_original <- function(base_size = 15) {
  theme_bw(base_family = font_family, base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.text = element_text(color = "black", size = base_size),
      axis.title = element_text(color = "black", size = base_size + 1),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_blank(),
      plot.margin = margin(6, 8, 6, 7)
    )
}

auc_value <- function(y, score) {
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

roc_coordinates <- function(y, score) {
  ord <- order(score, decreasing = TRUE)
  y <- y[ord]
  data.table(
    FPR = c(0, cumsum(y == 0) / sum(y == 0), 1),
    TPR = c(0, cumsum(y == 1) / sum(y == 1), 1)
  )
}

manifest <- list()

## 读取87例共同样本和pathway score -------------------------------------------

manifest87 <- fread(file.path(
  "generated", "clinical_context",
  "analysis_manifest_87.tsv"
))
score <- fread(file.path(
  "generated", "pathway_scores",
  "Hallmark_single_sample_rank_scores_wide.tsv.gz"
))
score87 <- merge(
  manifest87,
  score,
  by = c("sample_id", "library_id", "group", "batch"),
  all.x = TRUE,
  sort = FALSE
)
stopifnot(nrow(score87) == 87L)

pathway_show <- c(
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_BILE_ACID_METABOLISM",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_ESTROGEN_RESPONSE_LATE"
)
pathway_label <- c(
  "Cholesterol homeostasis",
  "IL-2/STAT5",
  "Mitotic spindle",
  "Bile acid metabolism",
  "Complement",
  "IL-6/JAK/STAT3",
  "Angiogenesis",
  "Late estrogen response"
)

setorder(score87, -DLCO_pct, group, sample_id)
mat <- as.matrix(score87[, ..pathway_show])
mat <- t(scale(mat))
mat[mat > 2] <- 2
mat[mat < -2] <- -2
rownames(mat) <- pathway_label
colnames(mat) <- score87$sample_id

ann_col <- data.frame(
  Severity = factor(
    ifelse(score87$outcome == 1, "Mod. & Severe", "Mild"),
    levels = c("Mild", "Mod. & Severe")
  ),
  Subtype = factor(score87$group, levels = c("MDA5", "ARS")),
  DLCO = score87$DLCO_pct,
  row.names = score87$sample_id
)
ann_colors <- list(
  Severity = severity_color,
  Subtype = group_color,
  DLCO = colorRampPalette(c("#D95F02", "white", "#1B9E77"))(100)
)

draw_pathway_heatmap <- function() {
  ph <- pheatmap(
    mat,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    breaks = seq(-2, 2, length.out = 51),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize = 13,
    fontsize_row = 13,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    annotation_names_col = TRUE,
    border_color = NA,
    legend = TRUE,
    annotation_legend = TRUE,
    fontfamily = font_family,
    silent = TRUE
  )
  gt <- ph$gtable
  legend_pos <- min(gt$layout$l[gt$layout$name %in% c("legend", "annotation_legend")])
  gt <- gtable_add_cols(gt, unit(2, "mm"), pos = legend_pos - 1)
  grid.newpage()
  grid.draw(gt)
}
manifest[[length(manifest) + 1L]] <- save_grid(
  draw_pathway_heatmap, "Figure6A_pathway_score_heatmap", 6.6, 3.8
)

fwrite(data.table(
  pathway = rep(pathway_label, ncol(mat)),
  sample_id = rep(colnames(mat), each = nrow(mat)),
  z_score = as.vector(mat)
), file.path(data_dir, "Figure6A_pathway_heatmap_data.tsv"), sep = "\t")
fwrite(cbind(data.table(sample_id = rownames(ann_col)), as.data.table(ann_col)),
       file.path(data_dir, "Figure6A_pathway_heatmap_annotations.tsv"), sep = "\t")

## 转录本与pathway模型的participant-level OOF结果 -----------------------------

gene_patient <- fread(file.path(
  "generated", "transcript", "results",
  "Pooled_severity_participant_mean_predictions.tsv"
))[validation == "repeated_nested_5fold"]

path_patient <- fread(file.path(
  "generated", "clinical_context",
  "participant_mean_oof_scores.tsv"
))[model == "EV_RNA_pathway_only" & penalty_rule == "lambda_min"]

patient <- merge(
  gene_patient[, .(sample_id, library_id, group, truth, transcript_score = score)],
  path_patient[, .(sample_id, pathway_score = mean_oof_score)],
  by = "sample_id",
  all = FALSE
)
stopifnot(nrow(patient) == 87L)
patient[, Severity := factor(
  ifelse(truth == 1, "Mod. & Severe", "Mild"),
  levels = c("Mild", "Mod. & Severe")
)]

score_long <- melt(
  patient,
  id.vars = c("sample_id", "library_id", "group", "truth", "Severity"),
  measure.vars = c("transcript_score", "pathway_score"),
  variable.name = "Model",
  value.name = "OOF_probability"
)
score_long[, Model := factor(
  Model,
  levels = c("transcript_score", "pathway_score"),
  labels = c("Transcript-level", "Pathway-level")
)]

score_test <- score_long[, .(
  p_value = wilcox.test(OOF_probability ~ Severity, exact = FALSE)$p.value,
  mild_median = median(OOF_probability[Severity == "Mild"]),
  mod_severe_median = median(OOF_probability[Severity == "Mod. & Severe"])
), by = Model]
score_test[, significance := fifelse(p_value < 0.001, "***",
                              fifelse(p_value < 0.01, "**",
                              fifelse(p_value < 0.05, "*", "NS")))]

p_score <- ggplot(score_long, aes(Severity, OOF_probability, fill = Severity)) +
  geom_boxplot(outlier.shape = NA, width = 0.62) +
  geom_jitter(shape = 16, width = 0.20, size = 1.65) +
  geom_text(
    data = score_test,
    aes(x = 1.5, y = 1.02, label = significance),
    inherit.aes = FALSE,
    family = font_family,
    size = 5
  ) +
  facet_wrap(~ Model, nrow = 1) +
  scale_fill_manual(values = severity_color) +
  scale_x_discrete(labels = c("Mild", "Mod. &\nSevere")) +
  scale_y_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "", y = "Out-of-fold probability") +
  theme_original(14) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 15, face = "bold"),
    panel.spacing.x = unit(8, "mm")
  )
manifest[[length(manifest) + 1L]] <- save_plot(
  p_score, "Figure6B_OOF_probability_model_comparison", 5.2, 3.8
)
fwrite(score_long, file.path(data_dir, "Figure6B_OOF_probability_data.tsv"), sep = "\t")
fwrite(score_test, file.path(data_dir, "Figure6B_OOF_probability_tests.tsv"), sep = "\t")

## C. 两类模型的阈值指标 ------------------------------------------------------

gene_perf <- fread(file.path(
  "generated", "transcript", "results", "Fig2_Fig3_model_performance.tsv"
))[comparison == "Pooled_severity" & validation == "repeated_nested_5fold"]
path_ci <- fread(file.path("generated",
  "DM_severity_reviewer_completion_20260830", "output", "severity_pr_auc_auc_brier.tsv"))[
  model == "EV_RNA_pathway_only"]
path_perf <- data.table(secondary_ci_low=path_ci$AUC_low,secondary_ci_high=path_ci$AUC_high,
  accuracy=mean((patient$pathway_score>=0.5)==patient$truth),
  sensitivity=mean(patient$pathway_score[patient$truth==1]>=0.5),
  specificity=mean(patient$pathway_score[patient$truth==0]<0.5))

metric <- rbind(
  data.table(
    Model = "Transcript-level",
    Metric = c("Accuracy", "Sensitivity", "Specificity"),
    Value = c(gene_perf$accuracy_0.5, gene_perf$sensitivity_0.5, gene_perf$specificity_0.5)
  ),
  data.table(
    Model = "Pathway-level",
    Metric = c("Accuracy", "Sensitivity", "Specificity"),
    Value = c(path_perf$accuracy, path_perf$sensitivity, path_perf$specificity)
  )
)
metric[, Model := factor(Model, levels = c("Transcript-level", "Pathway-level"))]
metric[, Metric := factor(Metric, levels = c("Accuracy", "Sensitivity", "Specificity"))]

p_metric <- ggplot(metric, aes(Metric, Value, fill = Model)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65) +
  geom_text(
    aes(label = sprintf("%.2f", Value)),
    position = position_dodge(width = 0.72), vjust = -0.35,
    family = font_family, size = 4.2
  ) +
  scale_fill_manual(values = model_color) +
  scale_y_continuous(limits = c(0, 0.86), expand = c(0, 0)) +
  labs(x = "", y = "Value") +
  theme_original(14) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.key.width = unit(6, "mm"),
    axis.text.x = element_text(size = 13)
  )
manifest[[length(manifest) + 1L]] <- save_plot(
  p_metric, "Figure6D_threshold_metric_comparison", 4.1, 3.8
)
fwrite(metric, file.path(data_dir, "Figure6D_threshold_metric_data.tsv"), sep = "\t")

## E. pathway模型的患者级OOF混淆矩阵 ---------------------------------------

patient[, pathway_predicted := as.integer(pathway_score >= 0.5)]
severity_labels <- c("Mod. &\nSev.", "Mild")
cm_pathway <- table(
  factor(patient$truth, levels = c(1, 0), labels = severity_labels),
  factor(patient$pathway_predicted, levels = c(1, 0), labels = severity_labels)
)
dimnames(cm_pathway) <- list(
  "True class" = rownames(cm_pathway),
  "Predicted class" = colnames(cm_pathway)
)

draw_pathway_confusion <- function() {
  pheatmap(
    cm_pathway,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    angle_col = 0,
    number_format = "%.0f",
    number_color = "black",
    fontsize_number = 28,
    fontsize = 24,
    cellwidth = 88,
    cellheight = 88,
    border_color = NA,
    legend = FALSE,
    color = colorRampPalette(brewer.pal(6, "Blues"))(11),
    fontfamily = font_family,
    silent = FALSE
  )
  grid.text(
    "Out-of-fold predictions", x = 0.48, y = 0.95,
    gp = gpar(fontfamily = font_family, fontsize = 20, fontface = "bold")
  )
  grid.text(
    "Predicted class", x = 0.48, y = 0.035,
    gp = gpar(fontfamily = font_family, fontsize = 20, fontface = "bold")
  )
  grid.text(
    "True class", x = 0.035, y = 0.53, rot = 90,
    gp = gpar(fontfamily = font_family, fontsize = 20, fontface = "bold")
  )
}
manifest[[length(manifest) + 1L]] <- save_grid(
  draw_pathway_confusion, "Figure6C_pathway_pathway_OOF_confusion_matrix", 5.4, 4.4
)
fwrite(
  as.data.table(as.table(cm_pathway)),
  file.path(data_dir, "Figure6C_pathway_pathway_OOF_confusion_matrix.tsv"),
  sep = "\t"
)

fwrite(rbindlist(manifest), file.path(root,"figure_manifest.tsv"),sep="\t")
writeLines(capture.output(sessionInfo()),file.path(root,"sessionInfo.txt"))

# RStudio工作目录设为本文件夹
library(data.table)
library(pheatmap)
library(RColorBrewer)
library(grid)


dm <- "inputs"
output_root <- "generated/figure6"
output_dir <- file.path(output_root, "figures", "main_redesigned")
data_dir <- file.path(output_root, "data_used", "fig5_redesigned")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Helvetica"
prediction_file <- file.path(
  "generated", "transcript", "results",
  "Pooled_severity_participant_mean_predictions.tsv"
)

patient <- fread(prediction_file)[validation == "repeated_nested_5fold"]
stopifnot(nrow(patient) == 87L)
patient[, predicted := as.integer(score >= 0.5)]

severity_labels <- c("Mod. &\nSev.", "Mild")
cm_transcript <- table(
  factor(patient$truth, levels = c(1, 0), labels = severity_labels),
  factor(patient$predicted, levels = c(1, 0), labels = severity_labels)
)
dimnames(cm_transcript) <- list(
  "True class" = rownames(cm_transcript),
  "Predicted class" = colnames(cm_transcript)
)

tp <- unname(cm_transcript[1, 1])
fn <- unname(cm_transcript[1, 2])
fp <- unname(cm_transcript[2, 1])
tn <- unname(cm_transcript[2, 2])

metrics <- data.table(
  model = "Transcript-level",
  threshold = 0.5,
  TP = tp,
  FN = fn,
  FP = fp,
  TN = tn,
  sensitivity = tp / (tp + fn),
  specificity = tn / (tn + fp),
  accuracy = (tp + tn) / sum(cm_transcript)
)

draw_transcript_confusion <- function() {
  pheatmap(
    cm_transcript,
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
    color = colorRampPalette(c("#F2F2F2", "#7F7F7F"))(11),
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

pdf_file <- file.path(output_dir, "Figure6C_transcript_transcript_OOF_confusion_matrix.pdf")
png_file <- file.path(output_dir, "Figure6C_transcript_transcript_OOF_confusion_matrix.png")

cairo_pdf(pdf_file, width = 5.4, height = 4.4, family = font_family)
draw_transcript_confusion()
dev.off()

png(
  png_file, width = 5.4, height = 4.4, units = "in",
  res = 300, type = "cairo", bg = "white"
)
draw_transcript_confusion()
dev.off()

fwrite(
  as.data.table(as.table(cm_transcript)),
  file.path(data_dir, "Figure6C_transcript_transcript_OOF_confusion_matrix.tsv"),
  sep = "\t"
)
fwrite(
  metrics,
  file.path(data_dir, "Figure6C_transcript_transcript_OOF_metrics.tsv"),
  sep = "\t"
)



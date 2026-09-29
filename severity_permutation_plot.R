# RStudio工作目录设为本文件夹
library(data.table)
library(ggplot2)


out_dir <- file.path(
  "generated",
  "DM_severity_reviewer_completion_20260830", "output"
)
z <- fread(file.path(out_dir, "severity_full_pipeline_permutation_auc.tsv"))
s <- fread(file.path(out_dir, "severity_full_pipeline_permutation_summary.tsv"))
perm_auc <- z[observed == FALSE, AUC]
observed_auc <- s$observed_nested_cv_AUC
empirical_p <- s$empirical_P

p <- ggplot(data.table(AUC = perm_auc), aes(AUC)) +
  geom_histogram(binwidth = 0.025, boundary = 0.5, fill = "#B7B7B7", colour = "white") +
  geom_vline(xintercept = observed_auc, colour = "#C64B40", linewidth = 1.15) +
  annotate(
    "text", x = observed_auc, y = Inf,
    label = sprintf("Observed AUC = %.3f\nEmpirical P = %.3f", observed_auc, empirical_p),
    hjust = 1.05, vjust = 1.25, family = "Helvetica", size = 4.1
  ) +
  labs(x = "Nested-CV AUC under permuted labels", y = "Number of permutations") +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    axis.title = element_text(size = 15, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.margin = margin(7, 10, 7, 9)
  )

ggsave(file.path(out_dir, "S_severity_full_pipeline_permutation.pdf"), p,
       width = 5.33, height = 3.9, units = "in", device = cairo_pdf)
png(file.path(out_dir, "S_severity_full_pipeline_permutation.png"),
    width = 5.33, height = 3.9, units = "in", res = 300, type = "cairo")
print(p)
dev.off()

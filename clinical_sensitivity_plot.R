# RStudio工作目录设为本文件夹
library(data.table)
library(ggplot2)


source_dir <- file.path("generated", "clinical_sensitivity", "results")
out_dir <- file.path(
  "generated",
  "DM_severity_reviewer_completion_20260830", "output"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_order <- c(
  "primary_105", "age_sex_105",
  "duration_subset_base_93", "duration_adjusted_93",
  "steroid_subset_base_73", "steroid_adjusted_73",
  "full_subset_base_65", "full_adjusted_65"
)

model_labels <- c(
  primary_105 = "Primary",
  age_sex_105 = "Age/sex-adjusted",
  duration_subset_base_93 = "Duration subset",
  duration_adjusted_93 = "Duration-adjusted",
  steroid_subset_base_73 = "Steroid subset",
  steroid_adjusted_73 = "Steroid-adjusted",
  full_subset_base_65 = "Complete-case subset",
  full_adjusted_65 = "Fully adjusted"
)

panel_genes <- c("EGR1", "H2BC17", "CXCL8", "CCN1", "C2CD4B", "MAMLD1")

read_model <- function(model) {
  x <- fread(file.path(source_dir, paste0(model, "_full.tsv.gz")))
  x[gene_symbol %in% panel_genes, .(
    model, feature_id, gene_symbol, biotype, logFC, FDR
  )]
}

panel <- rbindlist(lapply(model_order, read_model), use.names = TRUE)
if (panel[, uniqueN(gene_symbol)] != length(panel_genes)) {
  stop("Not all current MDA5-vs-ARS panel genes were found.")
}
if (panel[, .N] != length(model_order) * length(panel_genes)) {
  stop("Unexpected number of current-panel sensitivity rows.")
}

panel[, significant := FDR < 0.05]
panel[, model_x := c(
  primary_105 = 1, age_sex_105 = 2,
  duration_subset_base_93 = 3.25, duration_adjusted_93 = 4.25,
  steroid_subset_base_73 = 5.5, steroid_adjusted_73 = 6.5,
  full_subset_base_65 = 7.75, full_adjusted_65 = 8.75
)[model]]
panel[, gene_symbol := factor(gene_symbol, levels = rev(panel_genes))]
panel[, support := factor(
  fifelse(significant, "BH FDR < 0.05", "BH FDR >= 0.05"),
  levels = c("BH FDR < 0.05", "BH FDR >= 0.05")
)]

base_theme <- theme_minimal(base_family = "Helvetica", base_size = 11.5) +
  theme(
    plot.title = element_text(size = 15.5, face = "bold", hjust = 0, margin = margin(b = 3)),
    plot.subtitle = element_text(size = 10.4, colour = "#333333", hjust = 0,
                                 lineheight = 1.02, margin = margin(b = 8)),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 9.0, colour = "#1A1A1A", angle = 90,
                               hjust = 1, vjust = 0.5, margin = margin(t = 5)),
    axis.text.y = element_text(size = 10.6, face = "italic", colour = "#1A1A1A",
                               margin = margin(r = 5)),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10.2, face = "bold", lineheight = 0.95),
    legend.text = element_text(size = 9.2),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(7, 8, 7, 7)
  )

p <- ggplot(panel, aes(model_x, gene_symbol)) +
  geom_tile(
    aes(fill = logFC, alpha = support), width = 0.92, height = 0.86,
    colour = "white", linewidth = 0.65
  ) +
  scale_fill_gradient2(
    name = "log2 fold\nchange",
    low = "#35689A", mid = "#F5F5F2", high = "#B54A42", midpoint = 0,
    limits = c(-3.3, 3.3), breaks = c(-3, 0, 3), labels = c("-3", "0", "+3"),
    oob = scales::squish,
    guide = guide_colorbar(
      title.position = "top", title.hjust = 0,
      barheight = grid::unit(3.7, "cm"), barwidth = grid::unit(0.38, "cm"),
      ticks.colour = "#4D4D4D", frame.colour = "#666666"
    )
  ) +
  scale_alpha_manual(
    values = c("BH FDR < 0.05" = 1, "BH FDR >= 0.05" = 0.22), guide = "none"
  ) +
  scale_x_continuous(
    breaks = c(1, 2, 3.25, 4.25, 5.5, 6.5, 7.75, 8.75),
    labels = unname(model_labels[model_order]),
    expand = expansion(add = c(0.1, 0.1))
  ) +
  scale_y_discrete(drop = FALSE, expand = expansion(add = c(0.05, 0.05))) +
  labs(title = NULL, subtitle = NULL) +
  base_theme

fwrite(
  panel[, .(
    model, model_label = unname(model_labels[model]), gene_symbol = as.character(gene_symbol),
    feature_id, logFC, FDR, significant
  )],
  file.path(out_dir, "current_MDA5_ARS_panel_clinical_sensitivity.tsv"), sep = "\t"
)

ggsave(
  file.path(out_dir, "S_current_MDA5_ARS_panel_clinical_sensitivity.pdf"), p,
  width = 4.95, height = 4.60, units = "in", device = cairo_pdf, bg = "white"
)
png(
  file.path(out_dir, "S_current_MDA5_ARS_panel_clinical_sensitivity.png"),
  width = 4.95, height = 4.60, units = "in", res = 400, type = "cairo"
)
print(p)
dev.off()

cat("Current fixed-panel clinical sensitivity outputs written to", out_dir, "\n")

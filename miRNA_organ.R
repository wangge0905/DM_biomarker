# RStudio工作目录设为本文件夹
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ragg)

root <- "generated/miRNA"
tissue_file <- "inputs/miRNA_tissue_selected.csv"

result_dir <- file.path(root, "results")
figure_dir <- file.path(root, "figures")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

tissue_mirnas <- c(
  "hsa-miR-146b-3p",
  "hsa-miR-146b-5p",
  "hsa-miR-223-3p",
  "hsa-miR-223-5p"
)


mirna_labels <- setNames(tolower(tissue_mirnas), tissue_mirnas)

font_family <- "Helvetica"

theme_manuscript <- theme_bw(base_size = 12, base_family = font_family) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.margin = margin(6, 8, 6, 6)
  )

save_plot <- function(plot, name, width, height) {
  ggsave(
    file.path(figure_dir, paste0(name, ".pdf")),
    plot = plot, width = width, height = height,
    units = "in", device = cairo_pdf
  )
  png_tmp <- file.path(tempdir(), paste0(name, ".png"))
  ggsave(
    png_tmp,
    plot = plot, width = width, height = height,
    units = "in", dpi = 300, bg = "white",
    device = ragg::agg_png
  )
  file.copy(
    png_tmp,
    file.path(figure_dir, paste0(name, ".png")),
    overwrite = TRUE
  )
  unlink(png_tmp)
}


# Figure S6C: TissueAtlas expression
tissue_raw <- fread(tissue_file)

tissue_plot_data <- tissue_raw %>%
  filter(type == "mirna", acc %in% tissue_mirnas) %>%
  group_by(organ, acc) %>%
  summarise(mean_expression = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  group_by(acc) %>%
  mutate(relative_expression = mean_expression / sum(mean_expression, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    organ_label = str_to_title(str_replace_all(organ, "_", " ")),
    organ_label = factor(organ_label, levels = sort(unique(organ_label))),
    miRNA = factor(
      acc,
      levels = rev(tissue_mirnas),
      labels = rev(unname(mirna_labels[tissue_mirnas]))
    )
  )

fwrite(tissue_plot_data, file.path(result_dir, "FigureS6C_tissue_expression.csv"))

p_tissue <- ggplot(
  tissue_plot_data,
  aes(x = organ_label, y = miRNA, size = relative_expression)
) +
  geom_point(shape = 21, fill = "#A50F15", colour = "#A50F15", stroke = 0.25) +
  scale_size_continuous(
    name = "Relative\nexpression",
    range = c(1.3, 6.0),
    breaks = c(0.01, 0.05, 0.10, 0.20),
    labels = c("0.01", "0.05", "0.10", "0.20")
  ) +
  labs(x = "Human tissue", y = NULL) +
  theme_manuscript +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 12, margin = margin(t = 5)),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.height = unit(0.36, "cm"),
    legend.position = "right"
  )

save_plot(p_tissue, "FigureS6C_updated_tissue_expression", 8.6, 3.55)


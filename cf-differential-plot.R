# RStudio工作目录设为本文件夹
library(data.table);library(ggplot2)
de_file <- "generated/paired_DE/results/paired_DE_all_primary.tsv"
figure_dir <- result_dir <- "generated/figure2E"
dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)
direction_col <- c(Plasma_higher = "#4C78A8", EV_higher = "#E07A5F",
                   Not_significant = "#C8C8C8")
biotype_col <- c(
  mRNA = "#4E79A7", lncRNA = "#A0CBE8", pseudogene = "#F28E2B",
  miRNA = "#E15759", snRNA = "#76B7B2", snoRNA = "#59A14F",
  `misc RNA` = "#EDC948", `mitochondrial RNA` = "#B07AA1",
  Other = "#BAB0AC"
)

theme_manuscript <- function(base_size = 15) {
  theme_classic(base_family = "Helvetica", base_size = base_size) +
    theme(
      axis.text = element_text(color = "black", size = base_size - 1.5),
      axis.title = element_text(color = "black", size = base_size),
      plot.title = element_text(face = "bold", size = base_size + 1,
                                hjust = 0),
      plot.subtitle = element_text(size = base_size - 2, color = "#444444"),
      legend.title = element_text(face = "bold", size = base_size - 1.5),
      legend.text = element_text(size = base_size - 2.5),
      legend.key.height = grid::unit(0.42, "cm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size - 2),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
}

png_cairo <- function(filename, width, height, units = "in", res = 400, ...) {
  grDevices::png(filename = filename, width = width, height = height,
                 units = units, res = res, type = "cairo",
                 family = "Helvetica", ...)
}

save_plot <- function(p, name, width, height) {
  ggsave(file.path(figure_dir, paste0(name, ".pdf")), p,
         width = width, height = height, device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(name, ".png")), p,
         width = width, height = height, dpi = 400, bg = "white",
         device = png_cairo)
}

collapse_biotype <- function(x) {
  out <- rep("Other", length(x))
  out[x == "protein_coding"] <- "mRNA"
  out[x == "lncRNA"] <- "lncRNA"
  out[grepl("pseudogene", x, ignore.case = TRUE)] <- "pseudogene"
  out[x == "miRNA"] <- "miRNA"
  out[x == "snRNA"] <- "snRNA"
  out[x == "snoRNA"] <- "snoRNA"
  out[x %in% c("misc_RNA", "scaRNA", "vault_RNA", "ribozyme",
               "sRNA", "scRNA", "Y_RNA")] <- "misc RNA"
  out[x %in% c("Mt_rRNA", "Mt_tRNA")] <- "mitochondrial RNA"
  factor(out, levels = names(biotype_col))
}

de <- read.delim(de_file, check.names = FALSE)
de$biotype_group <- collapse_biotype(de$biotype)
de$plot_status <- "Not_significant"
de$plot_status[de$FDR < 0.05 & de$logFC > 0] <- "EV_higher"
de$plot_status[de$FDR < 0.05 & de$logFC < 0] <- "Plasma_higher"
de$plot_status <- factor(de$plot_status,
                         levels = c("Plasma_higher", "Not_significant", "EV_higher"))

biotype_de <- as.data.frame(table(
  de$biotype_group[de$FDR < 0.05],
  ifelse(de$logFC[de$FDR < 0.05] > 0, "EV_higher", "Plasma_higher")
))
names(biotype_de) <- c("biotype", "direction", "n_features")
biotype_de$signed_n <- ifelse(biotype_de$direction == "EV_higher",
                              biotype_de$n_features, -biotype_de$n_features)
biotype_de$biotype <- factor(biotype_de$biotype,
                             levels = rev(names(biotype_col)))
fwrite(biotype_de, file.path(result_dir,
                             "FDR05_feature_counts_by_biotype.tsv"), sep = "\t")

pD <- ggplot(biotype_de, aes(signed_n, biotype, fill = direction)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, color = "#555555", linewidth = 0.45) +
  scale_fill_manual(values = direction_col[c("Plasma_higher", "EV_higher")],
                    labels = c(Plasma_higher = "Plasma higher",
                               EV_higher = "EV higher")) +
  scale_x_continuous(labels = abs, expand = expansion(mult = c(0.12, 0.12))) +
  labs(title = "Differential features by RNA class", subtitle = "BH FDR < 0.05",
       x = "Number of differential features", y = NULL, fill = NULL) +
  theme_manuscript(15) +
  theme(legend.position = "top")
save_plot(pD, "Figure2E_differential_feature_biotypes", 6.1, 4.25)

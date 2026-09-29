# RStudio工作目录设为本文件夹
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

root <- "inputs"
out <- file.path("generated", "figureS4")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

group_cols <- c(MDA5 = "#A6611A", ARS = "#DFC27D", HC = "#018571")
font_family <- "Helvetica"

qc_raw <- read.delim(file.path(root, "1204/QC.txt"), check.names = FALSE)
names(qc_raw)[36:55] <- paste0(substr(names(qc_raw)[36:55], 1, 12),
                               substr(names(qc_raw)[36:55], 16, 20))
qc_raw <- qc_raw[, -c(41:45, 51:55)]
anno <- read.delim(file.path(root, "annotation.txt"), check.names = FALSE)
names(anno)[names(anno) == "sample_library_id"] <- "library_id"

qc_mat <- as.data.frame(t(qc_raw[, -1, drop = FALSE]), check.names = FALSE)
names(qc_mat) <- make.unique(qc_raw[[1]], sep = ".")
qc_mat$library_id <- rownames(qc_mat)
qc <- anno %>%
  select(sample_id, library_id, group) %>%
  inner_join(qc_mat, by = "library_id") %>%
  mutate(group = factor(group, levels = c("MDA5", "ARS", "HC")))

num <- function(x) suppressWarnings(as.numeric(x))
metric <- function(x) if (x %in% names(qc)) num(qc[[x]]) else rep(0, nrow(qc))

qc2 <- qc %>%
  mutate(
    clean_reads = metric("clean"),
    usable_reads = metric("star_hg38_v38_dedup"),
    rRNA_reads = metric("star_rRNA"),
    exonic_reads = metric("MT_mRNA") + metric("MT_tRNA") + metric("chrM.all") +
      metric("mRNA") + metric("lncRNA") + metric("snoRNA") + metric("snRNA") +
      metric("srpRNA") + metric("tRNA") + metric("tucpRNA") + metric("Y_RNA") +
      metric("misc_RNA") + metric("pseudogene") + metric("exon"),
    intronic_reads = metric("intron"),
    intergenic_reads = metric("intergenic"),
    mRNA_reads = metric("mRNA"),
    lncRNA_reads = metric("lncRNA"),
    miRNA_reads = metric("miRNA"),
    sn_sno_reads = metric("snRNA") + metric("snoRNA"),
    other_small_reads = metric("tRNA") + metric("Y_RNA") + metric("srpRNA") +
      metric("misc_RNA") + metric("pseudogene") + metric("tucpRNA")
  ) %>%
  filter(!is.na(clean_reads), clean_reads > 0)
qc_analyzable <- qc2 %>% filter(clean_reads >= 1e6, usable_reads >= 3e5)

theme_dm <- function(base_size = 15) {
  theme_classic(base_family = font_family, base_size = base_size) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 13, colour = "#444444"),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13, colour = "black"),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12),
      plot.margin = margin(8, 10, 8, 8)
    )
}

save_pair <- function(p, stem, width, height) {
  ggsave(file.path(out, paste0(stem, ".pdf")), p, width = width, height = height,
         device = cairo_pdf, family = font_family)
  ggsave(file.path(out, paste0(stem, ".png")), p, width = width, height = height,
         dpi = 400, bg = "white", device = function(filename, ...) png(filename, type = "cairo", units="in", res=400, ...))
}

comp_cols <- c(Exonic="#D07A5F", Intronic="#7BA6B7", Intergenic="#D9B55F")
## S06: sample-level genomic-region composition
region_sample <- qc_analyzable %>%
  transmute(sample_id, group, Exonic = exonic_reads, Intronic = intronic_reads,
            Intergenic = intergenic_reads) %>%
  pivot_longer(c(Exonic, Intronic, Intergenic), names_to = "component", values_to = "reads") %>%
  group_by(sample_id) %>%
  mutate(fraction = reads / sum(reads, na.rm = TRUE)) %>%
  ungroup()
sample_order <- region_sample %>% filter(component == "Exonic") %>%
  arrange(group, desc(fraction)) %>% pull(sample_id)
region_sample$sample_id <- factor(region_sample$sample_id, levels = sample_order)
p06 <- ggplot(region_sample, aes(sample_id, fraction, fill = component)) +
  geom_col(width = 0.92) +
  facet_grid(~group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = comp_cols[c("Exonic", "Intronic", "Intergenic")]) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(title = "Genomic-region composition across libraries", x = NULL, y = "Read fraction", fill = NULL) +
  theme_dm() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "top", strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"), panel.spacing.x = unit(0.12, "cm"))
save_pair(p06, "exon_ratio", 7.2, 4.0)


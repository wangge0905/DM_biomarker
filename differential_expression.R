# RStudio工作目录设为本文件夹
library(data.table)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(grid)


options(stringsAsFactors = FALSE)
set.seed(20260831)

## 路径 -----------------------------------------------------------------------

dm <- "inputs"
output_root <- "generated/figure345_ABC"

figure_dir <- file.path(output_root, "figures")
data_dir <- file.path(output_root, "data_used")
log_dir <- file.path(output_root, "logs")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

combat_count_path <- file.path("generated", "descriptive_counts", "gencode_rmbatch.all.txt")
annotation_path <- file.path(dm, "annotation.txt")
metainfo_path <- file.path(dm, "metainfo.txt")

## 全局作图参数 ----------------------------------------------------------------

font_family <- "Helvetica"
heatmap_each_direction_n <- 20L
fdr_cutoff <- 0.05
logfc_cutoff <- 1

figure_size <- list(
  heatmap = c(7.0, 7.0),
  biotype = c(6.1, 4.25),
  volcano = c(4.6, 3.6)
)

brbg <- brewer.pal(5, "BrBG")
group_color <- c(
  `Anti-MDA5+` = brbg[1],
  `Anti-ARS+` = brbg[2],
  HCs = brbg[5]
)
display_rna <- function(x) {
  selected <- grepl("^MIR[0-9A-Z-]+$", x) & !grepl("HG$",x)
  x[selected] <- sub("^mir", "hsa-mir-", tolower(x[selected]))
  x
}
set1 <- brewer.pal(3, "Set1")
volcano_color <- c(UP = set1[1], DOWN = set1[2], Not = "grey")

comparison_info <- list(
  MDA5_vs_HC = list(
    figure = "Figure3",
    positive = "MDA5", negative = "HC",
    positive_label = "Anti-MDA5+", negative_label = "HCs",
    heatmap_size = c(7.0, 7.0)
  ),
  ARS_vs_HC = list(
    figure = "Figure4",
    positive = "ARS", negative = "HC",
    positive_label = "Anti-ARS+", negative_label = "HCs",
    heatmap_size = c(7.0, 5.0)
  ),
  MDA5_vs_ARS = list(
    figure = "Figure5",
    positive = "MDA5", negative = "ARS",
    positive_label = "Anti-MDA5+", negative_label = "Anti-ARS+",
    heatmap_size = c(7.0, 6.0)
  )
)

rna_class_levels <- c(
  "mRNA", "lncRNA", "pseudogene", "miRNA", "snRNA",
  "snoRNA", "misc RNA", "mitochondrial RNA", "Other"
)

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
  factor(out, levels = rna_class_levels)
}

theme_volcano <- function() {
  theme_bw(base_family = font_family, base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 1, fill = "transparent"),
      legend.title = element_blank(),
      legend.text = element_text(color = "black", size = 10),
      legend.key.size = unit(0.30, "cm"),
      axis.text = element_text(color = "black", size = 15),
      axis.title = element_text(color = "black", size = 15),
      plot.margin = margin(7, 9, 7, 7)
    )
}

theme_biotype <- function() {
  theme_classic(base_family = font_family, base_size = 15) +
    theme(
      axis.text = element_text(color = "black", size = 13.5),
      axis.title = element_text(color = "black", size = 15),
      legend.title = element_blank(),
      legend.text = element_text(size = 12.5),
      legend.key.height = unit(0.42, "cm"),
      legend.position = "top",
      plot.margin = margin(6, 12, 6, 9)
    )
}

save_ggplot <- function(p, stem, size) {
  pdf_file <- file.path(figure_dir, paste0(stem, ".pdf"))
  png_file <- file.path(figure_dir, paste0(stem, ".png"))
  ggsave(
    pdf_file, p, width = size[1], height = size[2], units = "in",
    device = cairo_pdf, bg = "white"
  )
  png(
    png_file, width = size[1], height = size[2], units = "in",
    res = 300, type = "cairo", bg = "white"
  )
  print(p)
  dev.off()
  data.table(
    stem = stem,
    pdf = file.path("figures", basename(pdf_file)),
    png = file.path("figures", basename(png_file)),
    width_in = size[1], height_in = size[2]
  )
}

save_heatmap <- function(draw_fun, stem, size) {
  pdf_file <- file.path(figure_dir, paste0(stem, ".pdf"))
  png_file <- file.path(figure_dir, paste0(stem, ".png"))

  cairo_pdf(pdf_file, width = size[1], height = size[2], family = font_family)
  draw_fun()
  dev.off()

  png(
    png_file, width = size[1], height = size[2], units = "in",
    res = 300, type = "cairo", bg = "white"
  )
  draw_fun()
  dev.off()

  data.table(
    stem = stem,
    pdf = file.path("figures", basename(pdf_file)),
    png = file.path("figures", basename(png_file)),
    width_in = size[1], height_in = size[2]
  )
}

## 读取ComBat-seq矩阵和样本信息 ------------------------------------------------


combat_count <- as.matrix(read.table(
  combat_count_path,
  sep = "\t", header = TRUE, row.names = 1,
  check.names = FALSE, quote = "", comment.char = ""
))
storage.mode(combat_count) <- "double"

anno <- fread(annotation_path)
setnames(anno, "sample_library_id", "library_id")
anno <- anno[mistake == 0]
anno[, original_order := .I]

missing_library <- setdiff(anno$library_id, colnames(combat_count))
if (length(missing_library) > 0) {
  stop("Libraries absent from ComBat-seq matrix: ", paste(missing_library, collapse = ", "))
}
combat_count <- combat_count[, anno$library_id, drop = FALSE]

metainfo <- fread(metainfo_path)
setnames(metainfo, "sample_library_id", "library_id")
anno <- merge(
  anno,
  metainfo[, .(sample_id, library_id, gender, age)],
  by = c("sample_id", "library_id"), all.x = TRUE, sort = FALSE
)

feature_split <- tstrsplit(rownames(combat_count), "|", fixed = TRUE)
feature_meta <- data.table(
  feature_id = rownames(combat_count),
  gene_symbol = feature_split[[3]],
  biotype = feature_split[[4]]
)

## 沿用原稿的RNA feature universe
original_nc_biotypes <- c(
  "lncRNA", "Mt_tRNA", "scaRNA", "sRNA",
  "IG_C_gene", "IG_D_gene", "IG_V_gene", "IG_J_gene",
  "TR_C_gene", "TR_D_gene", "TR_V_gene", "TR_J_gene",
  "miRNA", "snRNA", "snoRNA"
)
feature_keep <- feature_meta$biotype %in% c(original_nc_biotypes, "protein_coding") &
  feature_meta$gene_symbol != "Y_RNA"
combat_count <- combat_count[feature_keep, , drop = FALSE]
feature_meta <- feature_meta[feature_keep]

stopifnot(identical(rownames(combat_count), feature_meta$feature_id))

## edgeR：沿用原稿的基本流程，只把显著性统一为FDR ------------------------------

run_edger <- function(count_matrix, sample_table, positive, negative) {
  des <- copy(sample_table[group %in% c(positive, negative)])
  des[, analysis_group := fifelse(group == positive, "positive", "negative")]
  des[, analysis_group := factor(analysis_group, levels = c("negative", "positive"))]
  setorder(des, analysis_group, original_order)

  y <- DGEList(
    counts = count_matrix[, des$library_id, drop = FALSE],
    group = des$analysis_group
  )
  keep <- filterByExpr(y, group = des$analysis_group, min.count = 2, min.prop = 0.2)
  y <- y[keep, , keep.lib.sizes = TRUE]
  y <- calcNormFactors(y, method = "TMM")
  design <- model.matrix(~ analysis_group, data = des)
  rownames(design) <- des$library_id
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  test <- glmQLFTest(fit, coef = 2)
  result <- as.data.table(topTags(test, n = Inf, adjust.method = "BH")$table, keep.rownames = "feature_id")
  result <- merge(result, feature_meta, by = "feature_id", all.x = TRUE, sort = FALSE)
  result[, significant := FDR < fdr_cutoff & abs(logFC) >= logfc_cutoff]

  list(result = result, samples = des, y = y)
}

manifest <- list()
summary_list <- list()

for (comparison in names(comparison_info)) {
  info <- comparison_info[[comparison]]

  fit <- run_edger(combat_count, anno, info$positive, info$negative)
  de <- fit$result
  des <- fit$samples

  fwrite(
    de,
    file.path(data_dir, paste0(comparison, "_ComBatSeq_edgeR_full.tsv.gz")),
    sep = "\t"
  )

  ## A. 热图：每个方向最多20个RNA，按FDR排序 -------------------------------
  hm_up <- de[FDR < fdr_cutoff & logFC >= logfc_cutoff][order(FDR, -logFC)]
  hm_down <- de[FDR < fdr_cutoff & logFC <= -logfc_cutoff][order(FDR, logFC)]
  hm_feature <- rbind(
    head(hm_up, heatmap_each_direction_n),
    head(hm_down, heatmap_each_direction_n)
  )
  if (nrow(hm_feature) == 0) stop("No heatmap features for ", comparison)

  display_group <- setNames(
    c(info$positive_label, info$negative_label),
    c(info$positive, info$negative)
  )
  des[, display_group := factor(
    unname(display_group[as.character(group)]),
    levels = c(info$positive_label, info$negative_label)
  )]
  setorder(des, display_group, original_order)

  y_display <- DGEList(counts = combat_count[, des$library_id, drop = FALSE])
  y_display <- calcNormFactors(y_display, method = "TMM")
  log_cpm <- cpm(y_display, log = TRUE, prior.count = 1)
  hm_matrix <- log_cpm[hm_feature$feature_id, , drop = FALSE]
  hm_matrix <- t(scale(t(hm_matrix), center = TRUE, scale = TRUE))
  hm_matrix[!is.finite(hm_matrix)] <- 0
  hm_matrix[hm_matrix > 2] <- 2
  hm_matrix[hm_matrix < -2] <- -2
  rownames(hm_matrix) <- make.unique(ifelse(
    is.na(hm_feature$gene_symbol) | hm_feature$gene_symbol == "",
    hm_feature$feature_id,
    display_rna(hm_feature$gene_symbol)
  ))

  age_levels <- paste(seq(10, 80, by = 10), seq(19, 89, by = 10), sep = "-")
  ann_col <- data.frame(
    class = as.character(des$display_group),
    gender = as.character(des$gender),
    age = cut(
      as.numeric(des$age), breaks = seq(10, 90, by = 10),
      labels = age_levels, include.lowest = TRUE
    )
  )
  rownames(ann_col) <- des$library_id
  ann_row <- data.frame(gene_type = hm_feature$biotype)
  rownames(ann_row) <- rownames(hm_matrix)

  gene_type_levels <- unique(hm_feature$biotype)
  gene_type_colors <- brewer.pal(max(3, length(gene_type_levels)), "Set3")
  ann_colors <- list(
    class = group_color[c(info$positive_label, info$negative_label)],
    gender = c(M = "black", F = "white"),
    age = setNames(
      colorRampPalette(c("white", "darkblue"))(length(age_levels)),
      age_levels
    ),
    gene_type = setNames(
      gene_type_colors[seq_along(gene_type_levels)], gene_type_levels
    )
  )

  draw_heatmap <- function() {
    ph <- pheatmap(
      hm_matrix,
      color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
      breaks = seq(-2, 2, length.out = 51),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      show_colnames = FALSE,
      fontsize = 12,
      fontsize_row = 12,
      annotation_col = ann_col,
      annotation_row = ann_row,
      annotation_colors = ann_colors,
      annotation_names_row = FALSE,
      border_color = NA,
      legend = FALSE,
      annotation_legend = FALSE,
      fontfamily = font_family,
      silent = TRUE
    )

    grid.newpage()
    pushViewport(viewport(
      x = 0, y = 0.50, width = 0.77, height = 0.98,
      just = c("left", "center"), clip = "on"
    ))
    grid.draw(ph$gtable)
    popViewport()

    pushViewport(viewport(
      x = 0.78, y = 0.50, width = 0.22, height = 0.96,
      just = c("left", "center"), clip = "off"
    ))

    legend_text <- gpar(fontfamily = font_family, fontsize = 12, col = "black")
    legend_title <- gpar(
      fontfamily = font_family, fontsize = 12,
      fontface = "bold", col = "black"
    )
    swatch_x <- 0.08
    label_x <- 0.20
    swatch_w <- 0.10
    swatch_h <- 0.030

    z_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
    grid.raster(
      as.raster(matrix(rev(z_colors), ncol = 1)),
      x = swatch_x, y = 0.91, width = swatch_w, height = 0.16,
      interpolate = TRUE
    )
    grid.text(c("2", "0", "-2"), x = label_x,
              y = c(0.99, 0.91, 0.83), just = "left", gp = legend_text)

    draw_items <- function(title, labels, colors, title_y, first_y, step = 0.040) {
      grid.text(title, x = 0.02, y = title_y, just = "left", gp = legend_title)
      yy <- first_y - (seq_along(labels) - 1) * step
      grid.rect(
        x = swatch_x, y = yy, width = swatch_w, height = swatch_h,
        gp = gpar(fill = colors, col = NA)
      )
      grid.text(labels, x = label_x, y = yy, just = "left", gp = legend_text)
    }

    draw_items("age", names(ann_colors$age), unname(ann_colors$age),
               title_y = 0.79, first_y = 0.755)
    if (comparison == "ARS_vs_HC") {
      draw_items("gender", names(ann_colors$gender), unname(ann_colors$gender),
                 title_y = 0.44, first_y = 0.400)
      draw_items("class", names(ann_colors$class), unname(ann_colors$class),
                 title_y = 0.32, first_y = 0.280)
      draw_items("gene_type", names(ann_colors$gene_type), unname(ann_colors$gene_type),
                 title_y = 0.205, first_y = 0.165, step = 0.034)
    } else {
      draw_items("gender", names(ann_colors$gender), unname(ann_colors$gender),
                 title_y = 0.43, first_y = 0.395)
      draw_items("class", names(ann_colors$class), unname(ann_colors$class),
                 title_y = 0.30, first_y = 0.265)
      draw_items("gene_type", names(ann_colors$gene_type), unname(ann_colors$gene_type),
                 title_y = 0.18, first_y = 0.145, step = 0.029)
    }
    popViewport()
  }

  stem_a <- paste0(info$figure, "A_ComBatSeq_FDR_heatmap")
  manifest[[length(manifest) + 1L]] <- save_heatmap(
    draw_heatmap, stem_a, info$heatmap_size
  )
  fwrite(
    hm_feature,
    file.path(data_dir, paste0(comparison, "_heatmap_features.tsv")),
    sep = "\t"
  )

  ## B. RNA类别及方向的绝对数量 --------------------------------------------
  count_data <- de[FDR < fdr_cutoff & abs(logFC) >= logfc_cutoff]
  count_data[, rna_class := collapse_biotype(biotype)]
  count_data[, direction := ifelse(
    logFC > 0, info$positive_label, info$negative_label
  )]
  count_data <- count_data[, .(n_features = .N), by = .(rna_class, direction)]

  complete_grid <- CJ(
    rna_class = factor(rna_class_levels, levels = rna_class_levels),
    direction = c(info$negative_label, info$positive_label),
    unique = TRUE
  )
  count_data <- merge(
    complete_grid, count_data,
    by = c("rna_class", "direction"), all.x = TRUE, sort = FALSE
  )
  count_data[is.na(n_features), n_features := 0L]
  count_data[, signed_n := ifelse(
    direction == info$positive_label, n_features, -n_features
  )]
  count_data[, rna_class := factor(rna_class, levels = rev(rna_class_levels))]
  count_data[, comparison := comparison]

  label_data <- count_data[n_features > 0]
  label_data[, hjust_value := ifelse(signed_n < 0, 1.18, -0.18)]
  fill_values <- setNames(
    group_color[c(info$negative_label, info$positive_label)],
    c(info$negative_label, info$positive_label)
  )
  legend_labels <- setNames(
    c(
      paste(info$negative_label, "higher"),
      paste(info$positive_label, "higher")
    ),
    c(info$negative_label, info$positive_label)
  )

  p_biotype <- ggplot(count_data, aes(signed_n, rna_class, fill = direction)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, color = "#555555", linewidth = 0.45) +
    geom_text(
      data = label_data,
      aes(label = comma(n_features), hjust = hjust_value),
      family = font_family, size = 4.1, color = "black",
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = fill_values,
      breaks = c(info$negative_label, info$positive_label),
      labels = legend_labels
    ) +
    scale_x_continuous(
      labels = function(x) comma(abs(x)),
      expand = expansion(mult = c(0.16, 0.16))
    ) +
    labs(x = "Number of differential RNAs", y = NULL, fill = NULL) +
    coord_cartesian(clip = "off") +
    theme_biotype()

  stem_b <- paste0(info$figure, "B_ComBatSeq_FDR_RNA_class_counts")
  manifest[[length(manifest) + 1L]] <- save_ggplot(
    p_biotype, stem_b, figure_size$biotype
  )
  fwrite(
    count_data,
    file.path(data_dir, paste0(comparison, "_RNA_class_counts.tsv")),
    sep = "\t"
  )

  ## C. 火山图：纵轴为BH-adjusted P ----------------------------------------
  de[, threshold := factor(
    ifelse(
      FDR < fdr_cutoff & abs(logFC) >= logfc_cutoff,
      ifelse(logFC > 0, "UP", "DOWN"), "Not"
    ),
    levels = c("UP", "DOWN", "Not")
  )]
  de[, plot_y := -log10(pmax(FDR, .Machine$double.xmin))]
  top_gene <- rbind(
    head(de[threshold == "UP"][order(FDR, -logFC)], 5),
    head(de[threshold == "DOWN"][order(FDR, logFC)], 5)
  )

  p_volcano <- ggplot(de, aes(logFC, plot_y, color = threshold)) +
    geom_point(size = 2, alpha = 0.9, shape = 16) +
    scale_color_manual(values = volcano_color) +
    geom_vline(
      xintercept = c(-logfc_cutoff, logfc_cutoff),
      linetype = 4, color = "grey", linewidth = 0.6
    ) +
    geom_hline(
      yintercept = -log10(fdr_cutoff),
      linetype = 4, color = "grey", linewidth = 0.6
    ) +
    geom_text_repel(
      data = top_gene,
      aes(label = display_rna(gene_symbol)),
      color = "black", size = 2.5,
      box.padding = 0.5, point.padding = 0.5,
      max.overlaps = Inf, show.legend = FALSE
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.07, 0.08))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
    labs(
      x = "log2 fold change",
      y = expression(-log[10] * "(BH-adjusted " * italic(P) * ")"),
      title = ""
    ) +
    coord_cartesian(clip = "off") +
    theme_volcano()

  stem_c <- paste0(info$figure, "C_ComBatSeq_FDR_volcano")
  manifest[[length(manifest) + 1L]] <- save_ggplot(
    p_volcano, stem_c, figure_size$volcano
  )

  summary_list[[comparison]] <- data.table(
    comparison = comparison,
    positive_group = info$positive_label,
    negative_group = info$negative_label,
    tested_RNAs = nrow(de),
    positive_higher_FDR05_abs_logFC1 = sum(
      de$FDR < fdr_cutoff & de$logFC >= logfc_cutoff
    ),
    negative_higher_FDR05_abs_logFC1 = sum(
      de$FDR < fdr_cutoff & de$logFC <= -logfc_cutoff
    ),
    total_FDR05_abs_logFC1 = sum(de$FDR < fdr_cutoff & abs(de$logFC) >= logfc_cutoff)
  )
}

manifest_table <- rbindlist(manifest, fill = TRUE)
manifest_table[, threshold := "FDR < 0.05 and |log2FC| >= 1"]
manifest_table[, volcano_y_axis := "-log10(BH-adjusted P)"]
fwrite(manifest_table, file.path(output_root, "figure_manifest.tsv"), sep = "\t")
fwrite(rbindlist(summary_list), file.path(output_root, "DE_summary.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))



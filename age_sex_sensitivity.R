# RStudio工作目录设为本文件夹
library(data.table)
library(edgeR)
library(ggplot2)
library(openxlsx)


dm <- "inputs"
panel_root <- file.path("generated", "biomarker")
out_dir <- file.path(
  "generated",
  "DM_severity_reviewer_completion_20260830", "output"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta <- fread(file.path(panel_root, "meta", "global_discovery_validation_manifest.tsv"))
meta <- meta[!is.na(age) & !is.na(sex) & sex != ""]
meta[, batch := factor(batch)]
meta[, sex := factor(sex)]
meta[, age_z := as.numeric(scale(age))]

raw <- fread(file.path(dm, "1204", "gencode.txt"), data.table = FALSE, check.names = FALSE)
feature_id <- raw[[1]]
raw <- raw[, -1, drop = FALSE]
raw <- raw[, -c(40:44, 50:54), drop = FALSE]
colnames(raw)[35:44] <- sub("-GP", "", colnames(raw)[35:44], fixed = TRUE)
rownames(raw) <- feature_id
stopifnot(all(meta$library_id %in% colnames(raw)))
raw <- raw[, meta$library_id, drop = FALSE]

part <- tstrsplit(rownames(raw), "|", fixed = TRUE, fill = "")
annotation <- data.table(
  feature_id = rownames(raw), gene_symbol = part[[3]], biotype = part[[4]]
)
rna_type <- c(
  "protein_coding", "lncRNA", "Mt_tRNA", "scaRNA", "sRNA",
  "IG_C_gene", "IG_D_gene", "IG_V_gene", "IG_J_gene",
  "TR_C_gene", "TR_D_gene", "TR_V_gene", "TR_J_gene",
  "miRNA", "snRNA", "snoRNA"
)
keep <- annotation$biotype %in% rna_type & annotation$gene_symbol != "Y_RNA"
count <- as.matrix(raw[keep, , drop = FALSE])
storage.mode(count) <- "double"
annotation <- annotation[keep]
rm(raw)

comparison_def <- list(
  MDA5_vs_HC = c("MDA5", "HC"),
  ARS_vs_HC = c("ARS", "HC"),
  MDA5_vs_ARS = c("MDA5", "ARS")
)

panels <- rbindlist(lapply(names(comparison_def), function(comparison) {
  x <- fread(file.path(
    panel_root, "results", "sensitivity_delta_auc_0.01",
    comparison, "fixed_biomarker_panel.tsv"
  ))
  x[, comparison := comparison]
  if (!"panel_rank" %in% names(x)) x[, panel_rank := seq_len(.N)]
  x
}), use.names = TRUE, fill = TRUE)

fit_de <- function(comparison, positive, negative, adjusted) {

  z <- copy(meta[group %in% c(positive, negative)])
  z[, batch := droplevels(batch)]
  z[, sex := droplevels(sex)]
  z[, group2 := factor(group, levels = c(negative, positive))]
  y <- DGEList(counts = count[, z$library_id, drop = FALSE], group = z$group2)
  keep_gene <- filterByExpr(y, group = z$group2, min.count = 2, min.prop = 0.2)
  y <- y[keep_gene, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  batch_estimable <- comparison == "MDA5_vs_ARS"
  design <- if (adjusted && batch_estimable) {
    model.matrix(~batch + age_z + sex + group2, data = z)
  } else if (adjusted) {
    model.matrix(~age_z + sex + group2, data = z)
  } else if (batch_estimable) {
    model.matrix(~batch + group2, data = z)
  } else {
    model.matrix(~group2, data = z)
  }
  stopifnot(qr(design)$rank == ncol(design))
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  coef_group <- grep("^group2", colnames(design))
  stopifnot(length(coef_group) == 1L)
  qlf <- glmQLFTest(fit, coef = coef_group)
  tab <- as.data.table(topTags(qlf, n = Inf, sort.by = "none")$table, keep.rownames = "feature_id")
  tab <- merge(tab, annotation, by = "feature_id", all.x = TRUE)
  tab[, `:=`(
    comparison = comparison,
    model = if (adjusted) "Age/sex-adjusted" else
      if (batch_estimable) "Batch-adjusted" else "Unadjusted",
    batch_estimable = batch_estimable,
    n = nrow(z), positive_n = sum(z$group2 == positive), negative_n = sum(z$group2 == negative)
  )]

  tab
}

cache_file <- file.path(out_dir, "all_current_panels_age_sex_full_de_cache.tsv.gz")
if (file.exists(cache_file)) {
  all_de <- fread(cache_file)
} else {
  all_de <- list()
  for (comparison in names(comparison_def)) {
    g <- comparison_def[[comparison]]
    all_de[[paste0(comparison, "_base")]] <- fit_de(comparison, g[1], g[2], FALSE)
    all_de[[paste0(comparison, "_adjusted")]] <- fit_de(comparison, g[1], g[2], TRUE)
  }
  all_de <- rbindlist(all_de, use.names = TRUE, fill = TRUE)
  fwrite(all_de, cache_file, sep = "\t")
}

panel_effect <- merge(
  all_de,
  panels[, .(comparison, feature_id, panel_rank, panel_gene = gene_symbol)],
  by = c("comparison", "feature_id")
)
panel_effect[, significant := FDR < 0.05]
setorder(panel_effect, comparison, panel_rank, model)

fwrite(
  panel_effect[, .(
    comparison, model, panel_rank, feature_id, gene_symbol, biotype,
    logFC, logCPM, F, PValue, FDR, significant, batch_estimable,
    n, positive_n, negative_n
  )],
  file.path(out_dir, "all_current_panels_age_sex_sensitivity.tsv"), sep = "\t"
)

theme_panel <- theme_minimal(base_family = "Helvetica", base_size = 11.5) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 10.5, colour = "#1A1A1A", angle = 90,
                               hjust = 1, vjust = 0.5, margin = margin(t = 5)),
    axis.text.y = element_text(size = 11.5, face = "italic", colour = "#1A1A1A",
                               margin = margin(r = 5)),
    axis.ticks = element_blank(), panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10.2, face = "bold", lineheight = 0.95),
    legend.text = element_text(size = 9.2),
    plot.margin = margin(7, 8, 7, 7)
  )

draw_panel <- function(comparison, width, height) {
  comparison_name <- comparison
  z <- copy(panel_effect[comparison == comparison_name])
  gene_order <- panels[comparison == comparison_name][order(panel_rank), gene_symbol]
  z[, gene_symbol := factor(gene_symbol, levels = rev(gene_order))]
  first_model <- if (comparison == "MDA5_vs_ARS") "Batch-adjusted" else "Unadjusted"
  z[, model := factor(model, levels = c(first_model, "Age/sex-adjusted"))]
  z[, support := factor(
    fifelse(significant, "BH FDR < 0.05", "BH FDR >= 0.05"),
    levels = c("BH FDR < 0.05", "BH FDR >= 0.05")
  )]

  p <- ggplot(z, aes(model, gene_symbol)) +
    geom_tile(
      aes(fill = logFC, alpha = support), width = 0.92, height = 0.86,
      colour = "white", linewidth = 0.65
    ) +
    scale_fill_gradient2(
      name = "log2 fold\nchange", low = "#35689A", mid = "#F5F5F2", high = "#B54A42",
      midpoint = 0, limits = c(-5, 5), breaks = c(-5, 0, 5), labels = c("-5", "0", "+5"),
      oob = scales::squish,
      guide = guide_colorbar(
        title.position = "top", title.hjust = 0,
        barheight = grid::unit(3.0, "cm"), barwidth = grid::unit(0.38, "cm"),
        ticks.colour = "#4D4D4D", frame.colour = "#666666"
      )
    ) +
    scale_alpha_manual(
      values = c("BH FDR < 0.05" = 1, "BH FDR >= 0.05" = 0.22), guide = "none"
    ) +
    scale_x_discrete(drop = FALSE, expand = expansion(add = c(0.05, 0.05))) +
    scale_y_discrete(drop = FALSE, expand = expansion(add = c(0.05, 0.05)),
      labels = function(x) {x[x=="MIR146B"] <- "hsa-mir-146b"; x[x=="MIR223"] <- "hsa-mir-223"; x}) +
    theme_panel

  stem <- paste0("S_", comparison, "_current_panel_age_sex_sensitivity")
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p,
         width = width, height = height, units = "in", device = cairo_pdf, bg = "white")
  png(file.path(out_dir, paste0(stem, ".png")), width = width, height = height,
      units = "in", res = 400, type = "cairo")
  print(p)
  dev.off()
}

draw_panel("MDA5_vs_HC", 3.65, 3.80)
draw_panel("ARS_vs_HC", 3.65, 3.20)
draw_panel("MDA5_vs_ARS", 3.65, 4.25)

cat("Age/sex sensitivity for all current fixed panels written to", out_dir, "\n")

# RStudio工作目录设为本文件夹
options(stringsAsFactors = FALSE)


library(edgeR)


dm <- "inputs"
out <- "generated/paired_DE"
dir.create(file.path(out, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "logs"), recursive = TRUE, showWarnings = FALSE)

meta_file <- file.path(dm, "cfRNA", "metadata.txt")
count_file <- file.path(dm, "cfRNA", "paired-gencode.txt")
old_filter_file <- file.path(dm, "cfRNA", "paired-count_matrix_data_filter.txt")

meta0 <- read.delim(meta_file, check.names = FALSE)
counts0 <- read.delim(count_file, check.names = FALSE, row.names = 1)

meta <- meta0[, c("library_id", "sample_id", "patient", "subtype", "class")]
names(meta)[3] <- "donor_id"
meta$class <- factor(meta$class, levels = c("CF", "EV"))
meta$subtype <- factor(meta$subtype, levels = c("MDA5", "ARS", "HC"))

stopifnot(nrow(meta) == 18)
stopifnot(length(unique(meta$donor_id)) == 9)
stopifnot(all(table(meta$donor_id) == 2))
stopifnot(all(table(meta$donor_id, meta$class) == 1))
stopifnot(setequal(meta$sample_id, colnames(counts0)))

counts0 <- counts0[, meta$sample_id]
stopifnot(!anyNA(counts0), all(counts0 >= 0), all(counts0 == round(counts0)))

parts <- strsplit(rownames(counts0), "|", fixed = TRUE)
anno <- data.frame(
  feature_id = rownames(counts0),
  gene_id = vapply(parts, function(x) x[1], character(1)),
  gene_length = suppressWarnings(as.numeric(vapply(parts, function(x) x[2], character(1)))),
  gene_symbol = vapply(parts, function(x) if (length(x) >= 3) x[3] else NA_character_, character(1)),
  biotype = vapply(parts, function(x) if (length(x) >= 4) x[4] else NA_character_, character(1)),
  stringsAsFactors = FALSE
)

design <- model.matrix(~ donor_id + class, data = meta)
rownames(design) <- meta$sample_id
stopifnot("classEV" %in% colnames(design), qr(design)$rank == ncol(design))

run_ql <- function(counts, keep, robust = TRUE) {
  y <- DGEList(counts = counts)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, design, robust = robust)
  fit <- glmQLFit(y, design, robust = robust)
  qlf <- glmQLFTest(fit, coef = "classEV")
  tab <- topTags(qlf, n = Inf, sort.by = "none")$table
  tab$feature_id <- rownames(tab)
  list(y = y, tab = tab)
}

# 主分析：复现原分析的subtype表达过滤，配对donor固定效应
y0 <- DGEList(counts = counts0)
keep_subtype <- filterByExpr(y0, group = meta$subtype, min.count = 5, min.prop = 0.2)
fit_primary <- run_ql(counts0, keep_subtype, robust = FALSE)
fit_primary_robust <- run_ql(counts0, keep_subtype, robust = TRUE)

# 直接按CF/EV过滤作敏感性分析
keep_class_default <- filterByExpr(y0, group = meta$class)
keep_class_loose <- filterByExpr(y0, group = meta$class, min.count = 5, min.total.count = 10)
fit_class_default <- run_ql(counts0, keep_class_default, robust = TRUE)
fit_class_loose <- run_ql(counts0, keep_class_loose, robust = TRUE)

# 存档的14439-feature输入，用来核对旧224/155结果
counts_old <- read.delim(old_filter_file, check.names = FALSE, row.names = 1)
counts_old <- counts_old[, meta$sample_id]
fit_old <- run_ql(counts_old, rep(TRUE, nrow(counts_old)), robust = FALSE)

tab <- merge(anno, fit_primary$tab, by = "feature_id", all.y = TRUE, sort = FALSE)
tab <- tab[match(fit_primary$tab$feature_id, tab$feature_id), ]

logcpm <- cpm(fit_primary$y, log = TRUE, prior.count = 0.5)
pair_diff <- sapply(levels(factor(meta$donor_id)), function(id) {
  ev <- meta$sample_id[meta$donor_id == id & meta$class == "EV"]
  cf <- meta$sample_id[meta$donor_id == id & meta$class == "CF"]
  logcpm[, ev] - logcpm[, cf]
})
colnames(pair_diff) <- levels(factor(meta$donor_id))

tab$n_donors_EV_higher <- rowSums(pair_diff[tab$feature_id, , drop = FALSE] > 0)
tab$n_donors_CF_higher <- 9 - tab$n_donors_EV_higher
tab$n_donors_same_as_model <- ifelse(tab$logFC >= 0, tab$n_donors_EV_higher, tab$n_donors_CF_higher)
tab$median_paired_logCPM_difference <- apply(pair_diff[tab$feature_id, , drop = FALSE], 1, median)
tab$direction <- ifelse(tab$logFC >= 0, "EV_higher", "plasma_CF_higher")

tab <- tab[order(tab$PValue), ]
write.table(tab, file.path(out, "results", "paired_DE_all_primary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tab[tab$FDR < 0.05, ], file.path(out, "results", "paired_DE_FDR_lt_0.05.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tab[tab$FDR < 0.10, ], file.path(out, "results", "paired_DE_exploratory_FDR_lt_0.10.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tab[tab$FDR < 0.10 & abs(tab$logFC) >= 1, ],
            file.path(out, "results", "paired_DE_historical_cutoff_FDR_lt_0.10_abs_log2FC_ge_1.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

old <- fit_old$tab
old_ev <- sum(old$FDR <= 0.10 & old$logFC >= 1)
old_cf <- sum(old$FDR <= 0.10 & old$logFC <= -1)

summary_tab <- data.frame(
  analysis = c("primary_subtype_filter_nonrobustQL", "same_filter_robustQL",
               "class_default_filter_robustQL", "class_loose_filter_robustQL",
               "archived_14439_nonrobustQL"),
  n_features_tested = c(nrow(fit_primary$tab), nrow(fit_primary_robust$tab),
                        nrow(fit_class_default$tab), nrow(fit_class_loose$tab), nrow(fit_old$tab)),
  FDR_lt_0.05 = c(sum(fit_primary$tab$FDR < 0.05), sum(fit_primary_robust$tab$FDR < 0.05),
                  sum(fit_class_default$tab$FDR < 0.05), sum(fit_class_loose$tab$FDR < 0.05),
                  sum(fit_old$tab$FDR < 0.05)),
  FDR_lt_0.10 = c(sum(fit_primary$tab$FDR < 0.10), sum(fit_primary_robust$tab$FDR < 0.10),
                  sum(fit_class_default$tab$FDR < 0.10), sum(fit_class_loose$tab$FDR < 0.10),
                  sum(fit_old$tab$FDR < 0.10)),
  FDR_lt_0.10_abs_log2FC_ge_1 = c(
    sum(fit_primary$tab$FDR < 0.10 & abs(fit_primary$tab$logFC) >= 1),
    sum(fit_primary_robust$tab$FDR < 0.10 & abs(fit_primary_robust$tab$logFC) >= 1),
    sum(fit_class_default$tab$FDR < 0.10 & abs(fit_class_default$tab$logFC) >= 1),
    sum(fit_class_loose$tab$FDR < 0.10 & abs(fit_class_loose$tab$logFC) >= 1),
    sum(fit_old$tab$FDR <= 0.10 & abs(fit_old$tab$logFC) >= 1)
  ),
  EV_higher_historical_cutoff = c(
    sum(fit_primary$tab$FDR < 0.10 & fit_primary$tab$logFC >= 1),
    sum(fit_primary_robust$tab$FDR < 0.10 & fit_primary_robust$tab$logFC >= 1),
    sum(fit_class_default$tab$FDR < 0.10 & fit_class_default$tab$logFC >= 1),
    sum(fit_class_loose$tab$FDR < 0.10 & fit_class_loose$tab$logFC >= 1),
    old_ev
  ),
  plasma_CF_higher_historical_cutoff = c(
    sum(fit_primary$tab$FDR < 0.10 & fit_primary$tab$logFC <= -1),
    sum(fit_primary_robust$tab$FDR < 0.10 & fit_primary_robust$tab$logFC <= -1),
    sum(fit_class_default$tab$FDR < 0.10 & fit_class_default$tab$logFC <= -1),
    sum(fit_class_loose$tab$FDR < 0.10 & fit_class_loose$tab$logFC <= -1),
    old_cf
  )
)

write.table(meta, file.path(out, "paired_9donor_18library_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_tab, file.path(out, "results", "paired_DE_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample_id = colnames(fit_primary$y),
                       library_size = fit_primary$y$samples$lib.size,
                       norm_factor = fit_primary$y$samples$norm.factors),
            file.path(out, "results", "sample_library_size_and_TMM_factor.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(check = c("n_donors", "n_libraries", "two_libraries_per_donor",
                                  "one_CF_and_one_EV_per_donor", "matrix_columns_match_manifest",
                                  "subtype_filter_matches_archived_14439",
                                  "archived_result_224_CF_155_EV_reproduced"),
                       value = c(9, 18, all(table(meta$donor_id) == 2),
                                 all(table(meta$donor_id, meta$class) == 1),
                                 identical(colnames(counts0), meta$sample_id),
                                 setequal(rownames(counts0)[keep_subtype], rownames(counts_old)),
                                 old_cf == 224 & old_ev == 155)),
            file.path(out, "results", "validation_checks.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out, "logs", "sessionInfo.txt"))
writeLines(capture.output(tools::md5sum(c(meta_file, count_file, old_filter_file))),
           file.path(out, "logs", "input_md5.txt"))

print(summary_tab)

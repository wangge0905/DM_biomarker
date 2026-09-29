# RStudio工作目录设为本文件夹
library(data.table)
library(edgeR)


dm <- "inputs"
out_dir <- "generated/pathway_scores"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

set.seed(20260817)
bootstrap_n <- 2000L

raw_path <- file.path(dm, "1204", "gencode.txt")
anno_path <- file.path(dm, "annotation.txt")
rank_path <- file.path(
  dm, "MDA5_ARS_pathway_20260816", "ranked_genes",
  "standard_design_filter_protein_coding_rank.tsv"
)
gmt_path <- file.path(
  dm, "MDA5_ARS_pathway_20260816", "gene_sets",
  "h.all.v2025.1.Hs.symbols.gmt"
)

# 155个可用library
anno <- fread(anno_path, data.table = FALSE, check.names = FALSE)
colnames(anno)[2] <- "library_id"
anno <- anno[anno$mistake == 0 & anno$group %in% c("MDA5", "ARS", "HC"), , drop = FALSE]
anno$batch <- substr(anno$library_id, 1, 6)
anno$group <- factor(anno$group, levels = c("HC", "ARS", "MDA5"))
anno$iim_group <- ifelse(anno$group == "HC", "HC", "IIM")
stopifnot(nrow(anno) == 155L, !anyDuplicated(anno$library_id))

# count矩阵的整理方式与前面DE一致
gencode <- fread(raw_path, data.table = FALSE, check.names = FALSE)
feature_id <- gencode[[1]]
gencode <- gencode[, -1, drop = FALSE]
gencode <- gencode[, -c(40:44, 50:54), drop = FALSE]
colnames(gencode)[35:44] <- sub("-GP", "", colnames(gencode)[35:44], fixed = TRUE)
rownames(gencode) <- feature_id
stopifnot(setequal(colnames(gencode), anno$library_id))
gencode <- gencode[, anno$library_id, drop = FALSE]

# 固定使用MDA5-vs-ARS通路分析中的15,252个protein-coding gene universe
universe <- fread(rank_path)
stopifnot(nrow(universe) == 15252L, !anyDuplicated(universe$gene_symbol))
stopifnot(!anyDuplicated(universe$feature_id), all(universe$feature_id %in% rownames(gencode)))
count <- as.matrix(gencode[universe$feature_id, anno$library_id, drop = FALSE])
rownames(count) <- universe$gene_symbol
storage.mode(count) <- "double"
rm(gencode)
gc()

sample_qc <- data.table(
  library_id = colnames(count),
  library_size_fixed_universe = colSums(count),
  detected_genes_fixed_universe = colSums(count > 0)
)
sample_qc[, log10_library_size := log10(library_size_fixed_universe)]
sample_qc[, detected_gene_fraction := detected_genes_fixed_universe / nrow(count)]
sample_qc[, log10_library_size_z := as.numeric(scale(log10_library_size))]
sample_qc[, detected_genes_z := as.numeric(scale(detected_genes_fixed_universe))]

# TMM-logCPM用于排名；后续score只使用每个样本内的gene rank
y <- DGEList(counts = count)
y <- calcNormFactors(y, method = "TMM")
log_cpm <- cpm(y, log = TRUE, prior.count = 0.5)
rank_pct <- apply(log_cpm, 2, rank, ties.method = "average") / (nrow(log_cpm) + 1)
rownames(rank_pct) <- rownames(log_cpm)

read_gmt <- function(path, universe_gene) {
  z <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  ans <- lapply(z, function(v) intersect(v[-c(1, 2)], universe_gene))
  names(ans) <- vapply(z, `[`, character(1), 1)
  ans[lengths(ans) >= 15L & lengths(ans) <= 500L]
}

hallmark <- read_gmt(gmt_path, rownames(rank_pct))
stopifnot(length(hallmark) == 50L)

# 单样本通路分数：通路gene在该样本全部gene中的平均百分位秩
score_mat <- vapply(hallmark, function(g) colMeans(rank_pct[g, , drop = FALSE]), numeric(ncol(rank_pct)))
rownames(score_mat) <- colnames(rank_pct)

score_wide <- data.table(
  sample_id = anno$sample_id,
  library_id = anno$library_id,
  group = as.character(anno$group),
  iim_group = anno$iim_group,
  batch = anno$batch
)
sample_qc_ordered <- sample_qc[match(score_wide$library_id, sample_qc$library_id)]
stopifnot(identical(score_wide$library_id, sample_qc_ordered$library_id))
score_wide <- cbind(score_wide, sample_qc_ordered[, -"library_id"])
score_wide <- cbind(score_wide, as.data.table(score_mat))
score_long <- melt(
  score_wide,
  id.vars = c(
    "sample_id", "library_id", "group", "iim_group", "batch",
    "library_size_fixed_universe", "detected_genes_fixed_universe",
    "log10_library_size", "detected_gene_fraction",
    "log10_library_size_z", "detected_genes_z"
  ),
  variable.name = "pathway", value.name = "rank_score"
)
score_long[, pathway := as.character(pathway)]


fwrite(score_wide, file.path(out_dir, "Hallmark_single_sample_rank_scores_wide.tsv.gz"), sep="\t")
fwrite(score_long, file.path(out_dir, "Hallmark_single_sample_rank_scores_long.tsv.gz"), sep="\t")
fwrite(data.table(pathway=names(hallmark), n_genes=lengths(hallmark)), file.path(out_dir,"gene_set_membership_counts.tsv"),sep="\t")

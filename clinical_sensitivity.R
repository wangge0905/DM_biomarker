# RStudio工作目录设为本文件夹
library(data.table)
library(edgeR)


dm <- "inputs"
out_dir <- file.path("generated", "clinical_sensitivity")
for (p in c("results", "meta", "logs")) dir.create(file.path(out_dir,p),recursive=TRUE,showWarnings=FALSE)
result_dir <- file.path(out_dir, "results")
meta_dir <- file.path(out_dir, "meta")
log_dir <- file.path(out_dir, "logs")

raw_path <- file.path(dm, "1204", "gencode.txt")
anno_path <- file.path(dm, "annotation.txt")
clinical_path <- file.path(
  dm, "refine-logs", "clinical_rna_meta_20260815",
  "analysis_ready_clinical_rna_meta.tsv"
)
base_result_path <- "inputs/clinical_sensitivity_features.txt"

fdr_cutoff <- 0.05

# 105例患者及临床变量
anno <- read.delim(anno_path, check.names = FALSE, stringsAsFactors = FALSE)
colnames(anno)[2] <- "library_id"
anno <- anno[anno$mistake == 0 & anno$group %in% c("MDA5", "ARS"), , drop = FALSE]
anno$batch <- substr(anno$library_id, 1, 6)

clinical <- fread(clinical_path, data.table = FALSE, check.names = FALSE)
clinical <- clinical[clinical$group %in% c("MDA5", "ARS"), , drop = FALSE]
clinical <- clinical[, c(
  "sample_id", "sample_library_id", "age", "sex",
  "disease_duration_years", "steroid_daily_mg"
)]
colnames(clinical)[2] <- "library_id"

stopifnot(!anyDuplicated(anno$library_id), !anyDuplicated(clinical$library_id))
stopifnot(setequal(anno$library_id, clinical$library_id))
meta <- merge(anno, clinical, by = c("sample_id", "library_id"), all.x = TRUE, sort = FALSE)
stopifnot(nrow(meta) == 105, !anyNA(meta$age), !anyNA(meta$sex), all(meta$sex != ""))
meta <- meta[order(meta$batch, meta$group, meta$library_id), , drop = FALSE]

covariates <- c("age", "sex", "disease_duration_years", "steroid_daily_mg")
completeness <- rbindlist(lapply(covariates, function(v) {
  rbindlist(lapply(c("All", "MDA5", "ARS"), function(g) {
    z <- if (g == "All") meta else meta[meta$group == g, , drop = FALSE]
    present <- !is.na(z[[v]]) & z[[v]] != ""
    data.table(
      variable = v, group = g, n_total = nrow(z),
      n_available = sum(present), n_missing = sum(!present),
      percent_available = round(100 * mean(present), 1)
    )
  }))
}))
fwrite(completeness, file.path(meta_dir, "clinical_covariate_completeness.tsv"), sep = "\t")

missing_pattern <- copy(as.data.table(meta[, c("sample_id", "library_id", "group", "batch")]))
for (v in covariates) missing_pattern[, (paste0(v, "_available")) := !is.na(meta[[v]]) & meta[[v]] != ""]
fwrite(missing_pattern, file.path(meta_dir, "clinical_covariate_availability_by_sample.tsv"), sep = "\t")

# count矩阵，处理方式与原差异分析一致
raw <- fread(raw_path, data.table = FALSE, check.names = FALSE)
feature_id <- raw[[1]]
raw <- raw[, -1, drop = FALSE]
raw <- raw[, -c(40:44, 50:54), drop = FALSE]
colnames(raw)[35:44] <- sub("-GP", "", colnames(raw)[35:44], fixed = TRUE)
rownames(raw) <- feature_id

stopifnot(setequal(colnames(raw), read.delim(anno_path, check.names = FALSE)[read.delim(anno_path, check.names = FALSE)$mistake == 0, 2]))
raw <- raw[, meta$library_id, drop = FALSE]

parts <- tstrsplit(rownames(raw), "|", fixed = TRUE, fill = "")
feature_info <- data.frame(
  feature_id = rownames(raw), gene_symbol = parts[[3]], biotype = parts[[4]],
  stringsAsFactors = FALSE
)
rownames(feature_info) <- feature_info$feature_id

base_result <- fread(base_result_path, data.table = FALSE)
fixed_feature_id <- base_result$feature_id
stopifnot(!anyDuplicated(fixed_feature_id), all(fixed_feature_id %in% rownames(raw)))
count_all <- as.matrix(raw[fixed_feature_id, , drop = FALSE])
storage.mode(count_all) <- "double"
rm(raw)
gc()

model_def <- list(
  primary_105 = list(vars = character(), formula = "~ batch + group", note = "105例主模型"),
  age_sex_105 = list(vars = c("age", "sex"), formula = "~ batch + age + sex + group", note = "105例年龄性别调整"),
  duration_subset_base_93 = list(vars = "disease_duration_years", formula = "~ batch + group", note = "病程完整子集，未加协变量"),
  duration_adjusted_93 = list(vars = c("age", "sex", "disease_duration_years"), formula = "~ batch + age + sex + disease duration + group", note = "病程完整子集，病程敏感性模型"),
  steroid_subset_base_73 = list(vars = "steroid_daily_mg", formula = "~ batch + group", note = "激素剂量完整子集，未加协变量"),
  steroid_adjusted_73 = list(vars = c("age", "sex", "steroid_daily_mg"), formula = "~ batch + age + sex + steroid dose + group", note = "激素剂量完整子集，治疗敏感性模型"),
  full_subset_base_65 = list(vars = c("disease_duration_years", "steroid_daily_mg"), formula = "~ batch + group", note = "病程和激素均完整子集，未加协变量"),
  full_adjusted_65 = list(vars = c("age", "sex", "disease_duration_years", "steroid_daily_mg"), formula = "~ batch + age + sex + disease duration + steroid dose + group", note = "病程和激素均完整的探索性完整模型")
)

make_design <- function(z, name) {
  z$group_factor <- factor(z$group, levels = c("ARS", "MDA5"))
  z$batch <- factor(z$batch)
  z$sex <- factor(z$sex)
  z$age_z <- as.numeric(scale(z$age))
  z$duration_z <- as.numeric(scale(z$disease_duration_years))
  z$steroid_z <- as.numeric(scale(z$steroid_daily_mg))

  if (name %in% c("primary_105", "duration_subset_base_93", "steroid_subset_base_73", "full_subset_base_65")) {
    design <- model.matrix(~ batch + group_factor, data = z)
  } else if (name == "age_sex_105") {
    design <- model.matrix(~ batch + age_z + sex + group_factor, data = z)
  } else if (name == "duration_adjusted_93") {
    design <- model.matrix(~ batch + age_z + sex + duration_z + group_factor, data = z)
  } else if (name == "steroid_adjusted_73") {
    design <- model.matrix(~ batch + age_z + sex + steroid_z + group_factor, data = z)
  } else if (name == "full_adjusted_65") {
    design <- model.matrix(~ batch + age_z + sex + duration_z + steroid_z + group_factor, data = z)
  } else {
    stop("unknown model: ", name)
  }
  rownames(design) <- z$library_id
  stopifnot(qr(design)$rank == ncol(design))
  list(meta = z, design = design)
}

fit_one <- function(name, def) {
  cc_vars <- unique(def$vars)
  use <- if (!length(cc_vars)) rep(TRUE, nrow(meta)) else complete.cases(meta[, cc_vars, drop = FALSE])
  z <- meta[use, , drop = FALSE]
  z <- z[order(z$batch, z$group, z$library_id), , drop = FALSE]
  built <- make_design(z, name)
  z <- built$meta
  design <- built$design

  y <- DGEList(counts = count_all[, z$library_id, drop = FALSE], group = z$group_factor)
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  coef_name <- "group_factorMDA5"
  qlf <- glmQLFTest(fit, coef = match(coef_name, colnames(design)))
  de <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  de$feature_id <- rownames(de)
  de$gene_symbol <- feature_info[de$feature_id, "gene_symbol"]
  de$biotype <- feature_info[de$feature_id, "biotype"]
  de$model <- name
  de <- de[, c("model", "feature_id", "gene_symbol", "biotype", "logFC", "logCPM", "F", "PValue", "FDR")]

  fwrite(de, file.path(result_dir, paste0(name, "_full.tsv.gz")), sep = "\t")
  fwrite(de[de$FDR < fdr_cutoff, , drop = FALSE], file.path(result_dir, paste0(name, "_FDR_lt_0.05.tsv")), sep = "\t")
  fwrite(data.frame(library_id = rownames(design), design, check.names = FALSE), file.path(meta_dir, paste0(name, "_design_matrix.tsv")), sep = "\t")
  fwrite(z[, c("sample_id", "library_id", "group", "batch", "age", "sex", "disease_duration_years", "steroid_daily_mg")], file.path(meta_dir, paste0(name, "_sample_manifest.tsv")), sep = "\t")

  data.frame(
    model = name, n = nrow(z), n_MDA5 = sum(z$group == "MDA5"), n_ARS = sum(z$group == "ARS"),
    design_columns = ncol(design), design_rank = qr(design)$rank,
    feature_universe = nrow(de), FDR_lt_0.05 = sum(de$FDR < fdr_cutoff),
    FDR_lt_0.05_abs_log2FC_ge_1 = sum(de$FDR < fdr_cutoff & abs(de$logFC) >= 1),
    formula = def$formula, note = def$note, stringsAsFactors = FALSE
  )
}

model_summary <- rbindlist(lapply(names(model_def), function(nm) fit_one(nm, model_def[[nm]])), fill = TRUE)
fwrite(model_summary, file.path(out_dir, "clinical_sensitivity_model_summary.tsv"), sep = "\t")

model_dictionary <- copy(model_summary)
model_dictionary[, `:=`(
  comparison = "anti-MDA5-positive versus anti-ARS-positive",
  model_family = "edgeR quasi-likelihood negative-binomial GLM",
  normalization = "TMM",
  feature_filter = "fixed features retained by the original current filter (count >=2 in >=20% of a group)",
  effect = "log2FC; positive values indicate higher abundance in anti-MDA5-positive",
  multiple_testing = "BH FDR across the fixed feature universe",
  FDR_cutoff = 0.05
)]
fwrite(model_dictionary, file.path(out_dir, "analysis_data_dictionary.tsv"), sep = "\t")

md5 <- tools::md5sum(c(raw_path, anno_path, clinical_path, base_result_path))
fwrite(data.frame(path = names(md5), md5 = unname(md5)), file.path(log_dir, "input_md5.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(log_dir, "clinical_sensitivity_sessionInfo.txt"))

print(completeness)
print(model_summary)

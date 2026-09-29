# RStudio工作目录设为本文件夹
library(dplyr)
library(stringr)
library(edgeR)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
dir.create("generated/TE",recursive=TRUE,showWarnings=FALSE)
{
  count <- read.table("inputs/TE/TEtranscript_count_matrix_heatmap.txt",sep="\t",header=T,check.names=F,row.names=1)
  filtered_vector <- colnames(count)[!grepl("GW_", colnames(count))]
  cleaned_vector <- gsub("GP_", "", filtered_vector)
  cleaned_vector <- gsub("-BC", "_BC", cleaned_vector)
  count <- count[,filtered_vector]
  colnames(count) <- cleaned_vector
  metadata <- read.table("inputs/annotation.txt", sep = "\t", header = T,check.names = F)
  metadata <- metadata[which(metadata$mistake == 0),]
  
  overlap <- intersect(metadata$sample_library_id,colnames(count))
  diff <- setdiff(metadata$sample_library_id,colnames(count))
  metadata <- metadata[metadata$sample_library_id %in% overlap,]
  count <- count[,metadata$sample_library_id]
  colnames(count) <- metadata$sample_id

  type <- metadata$group
  #class <- metadata$class
  
}

tecount=count[substr(rownames(count),1,4)!="ENSG",]
#heatmap
library(pheatmap)
library(RColorBrewer)
class = brewer.pal(5,'BrBG')[c(1,2,5)]
anno <- read.table("inputs/annotation.txt", sep = "\t", header = T)
colnames(anno)[2] <- "library_id"
anno <- anno[-c(117,155:158),]
#排序，可以注释掉
group_order <- c("MDA5", "ARS", "HC")
anno <- anno %>%
  arrange(factor(group, levels = group_order))
group <- anno$group
anno$group <- as.factor(anno$group)
anno$date <- paste0(substr(anno$library_id,1,6))
anno$date <- as.factor(anno$date)
str(anno)
levels(anno$group)
levels(anno$date)
sup_anno <- read.table("inputs/metainfo.txt", sep = "\t", header = T)
anno$gender <- sup_anno$gender
anno$age <- sup_anno$age
anno$gender <- as.factor(anno$gender)
anno$age <- as.numeric(anno$age)


matrix <- tecount
matrix <- matrix[,anno$sample_id]
subtype <-  anno$group
subtype <- factor(subtype, levels = c("MDA5", "ARS", "HC"))

group2 <- group
group2[group2 == "MDA5" | group2 == "ARS"] <- "DM"
group2 <- factor(group2, levels = c("DM", "HC"))

annotation <- anno
class_colors <- class[c(1,2,3)]
class_col <- c("MDA5","ARS","HC")
#matrix <- mt[,c(56:155)]
#subtype <- group[group != "MDA5"]
#annotation <- anno[c(56:155),]
#class_colors <- class[c(2,3)]
#class_col <- c("ARS","HC")
deg_heatmap <- function(matrix, subtype, annotation, class_colors, class_col) {
  y <- DGEList(counts = matrix, group = group2)
  keep <- filterByExpr(y, group = group2, min.count = 2, min.prop = 0.2)
  y <- y[keep, , keep.lib.size = TRUE]
  y <- calcNormFactors(y, method = "TMM")
  design <- model.matrix(~group2)
  
  rownames(design) <- colnames(y)
  y <- estimateDisp(y, design)
  y$common.dispersion
  
  cpm <- edgeR::cpm(y)
  

  
  fit.ql <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit.ql, coef=2)
  de <- topTags(qlf, n=Inf)$table
  de.pvalue <- filter(de, PValue < 0.01)
  
  de.dw <- filter(de.pvalue, logFC > 1)
  de.up <- filter(de.pvalue, logFC < -1)
  
  
  # 对 de.up 按 logFC 从大到小排序并取前 20 行
  if (nrow(de.up) > 20) {
    de.up <- de.up[order(de.up$logFC, decreasing = FALSE), ][1:20, ]
  }
  
  # 对 de.dw 按 logFC 从小到大排序并取前 20 行
  if (nrow(de.dw) > 20) {
    de.dw <- de.dw[order(de.dw$logFC, decreasing = TRUE), ][1:20, ]
  }
  
  logRPM <- edgeR::cpm(y, log = TRUE)
  logRPM <- logRPM[c(rownames(de.up), rownames(de.dw)), ]
  logRPM.scale <- scale(t(logRPM), center = TRUE, scale = TRUE)
  logRPM.scale <- t(logRPM.scale)
  
  # 去除重复的行，保留第一个
  tename <- unlist(lapply(strsplit(rownames(logRPM.scale), ":", fixed = TRUE), function(x) x[1]))
  tefamily <- unlist(lapply(strsplit(rownames(logRPM.scale), ":", fixed = TRUE), function(x) x[2]))
  teclass <- unlist(lapply(strsplit(rownames(logRPM.scale), ":", fixed = TRUE), function(x) x[3]))
  
  #logRPM.scale <- logRPM.scale[!duplicated(temp_rownames),]
  #genename <- rownames(logRPM.scale)
  #rownames(logRPM.scale) <- unlist(lapply(strsplit(rownames(logRPM.scale), "|", fixed = TRUE), function(x) x[3]))
  rownames(logRPM.scale) <- tename
  #gene_types <- sapply(strsplit(genename, "\\|"), function(x) x[4])
  
  ann_col <- data.frame(
    class = as.character(subtype), 
    gender = as.character(anno$gender),
    age = cut(anno$age, breaks = seq(10, 90, by = 10), labels = paste(seq(10, 80, by = 10), seq(19, 89, by = 10), sep = "-"))
  )
  rownames(ann_col) <- colnames(logRPM.scale)
  
  ann_colors <- list()
  ann_colors$class <- setNames(class_colors, class_col)
  ann_colors$gender <- c("M" = "black", "F" = "white")
  age_groups <- paste(seq(10, 80, by = 10), seq(19, 89, by = 10), sep = "-")
  age_colors <- colorRampPalette(c("white", "darkblue"))(length(age_groups))
  ann_colors$age <- setNames(age_colors, age_groups)
  
  ann_row <- data.frame(gene_type = teclass)
  rownames(ann_row) <- rownames(logRPM.scale)
  gene_type_colors <- brewer.pal(length(unique(teclass)), "Set3")
  ann_colors$gene_type <- setNames(brewer.pal(7,"Set3"), c("LINE","SINE","LTR","Retroposon","DNA","Satellite","Unknown"))
  
  col <- brewer.pal(9, "Set1")
  
  heatmap <- pheatmap(
    logRPM.scale, 
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    breaks = seq(-2, 2, length.out = 51),
    cutree_col = 0,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 12,
    angle_col = 90,
    annotation_col = ann_col,
    annotation_row = ann_row,
    annotation_colors = ann_colors,
    annotation_names_row = FALSE,
    border = FALSE
  )
  
  return(heatmap)
}

heatmap=deg_heatmap(matrix,subtype,annotation,class_colors,class_col)
ggsave("generated/TE/TE_DEG4.pdf",heatmap,width=9,height=6)

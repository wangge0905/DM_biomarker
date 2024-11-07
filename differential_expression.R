setwd("/Users/wangge/Documents/DM/")
library(tidyverse)
library(edgeR)
anno <- read.table("annotation.txt", sep = "\t", header = T)
colnames(anno)[2] <- "library_id"
anno <- anno[-c(117,155:158),]
#排序，可以注释掉
group_order <- c("MDA5", "ARS", "HC")
anno <- anno %>%
  arrange(factor(group, levels = group_order))

anno$group <- as.factor(anno$group)
anno$date <- paste0(substr(anno$library_id,1,6))
anno$date <- as.factor(anno$date)
str(anno)
levels(anno$group)
levels(anno$date)

ann <- anno[,-c(1,4)]
rownames(ann) = ann$library_id
ann <- ann[,-1]

sup_anno <- read.table("metainfo.txt", sep = "\t", header = T)
anno$gender <- sup_anno$gender
anno$age <- sup_anno$age
anno$gender <- as.factor(anno$gender)
anno$age <- as.numeric(anno$age)
{
  count_matrix.adjusted.all <- read.table("gencode_rmbatch.all.txt",sep="\t", header = TRUE, row.names = 1,check.names = FALSE)
  table(str_split(rownames(count_matrix.adjusted.all),"\\|",simplify=T)[,4])
  gn.type <- c("protein_coding","lncRNA",
               #"IG_C_pseudogene","IG_J_pseudogene","IG_pseudogene","IG_V_pseudogene",
               #"polymorphic_pseudogene","processed_pseudogene","pseudogene",
               #"TR_J_pseudogene","TR_V_pseudogene",
               #"transcribed_processed_pseudogene","transcribed_unitary_pseudogene","transcribed_unprocessed_pseudogene","translated_processed_pseudogene","translated_unprocessed_pseudogene","unitary_pseudogene","unprocessed_pseudogene",
               "Mt_tRNA",
               "scaRNA",
               "sRNA",
               "IG_C_gene","IG_D_gene","IG_V_gene","IG_J_gene",
               "TR_C_gene","TR_D_gene","TR_V_gene","TR_J_gene",
               "miRNA",
               "snRNA","snoRNA")
  gn.nc = count_matrix.adjusted.all[str_split(rownames(count_matrix.adjusted.all),"\\|",simplify=T)[,4] %in% gn.type,]       
  
  #gn.srp = count_matrix.adjusted.all[substr(str_split(rownames(count_matrix.adjusted.all),"\\|",simplify=T)[,3],1,5)=="RN7SL",]
  #gn.y = count_matrix.adjusted.all[str_split(rownames(count_matrix.adjusted.all),"\\|",simplify=T)[,3]=="Y_RNA",]
  
  gn.nc <- rbind(gn.nc, gn.srp, gn.y)
  mt <- gn.nc
  dim(mt)
  mt[1:4,1:2]
  # 137047
  colSums(mt)
  # 57642
}



library(pheatmap)
library(RColorBrewer)
class = brewer.pal(5,'BrBG')[c(1,2,5)]

matrix <- mt
subtype <- group2
annotation <- anno
class_colors <- class[c(1,2,3)]
class_col <- c("MDA5","ARS","HC")
matrix <- mt[,c(56:155)]
subtype <- group[group != "MDA5"]
annotation <- anno[c(56:155),]
class_colors <- class[c(2,3)]
class_col <- c("ARS","HC")
p5 = deg_heatmap(mt5 ,group5, anno5, class_colors, class_col)
p5
deg_heatmap <- function(matrix, subtype, annotation, class_colors, class_col) {
  y <- DGEList(counts = matrix, group = subtype)
  keep <- filterByExpr(y, group = subtype, min.count = 2, min.prop = 0.2)
  y <- y[keep, , keep.lib.size = TRUE]
  y <- calcNormFactors(y, method = "TMM")
  design <- model.matrix(~subtype)
  
  rownames(design) <- colnames(y)
  y <- estimateDisp(y, design)
  y$common.dispersion
  
  cpm <- edgeR::cpm(y)
  
  if(FALSE){
    fit.lr <- glmFit(y, design)
    lrt <- glmLRT(fit.lr, coef=2)
    #lrt <- glmLRT(fit.lr, contrast = c(-1,1))
    de <- topTags(lrt, n=Inf)$table
  }
  
  fit.ql <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit.ql, coef=2)
  de <- topTags(qlf, n=Inf)$table
  de.pvalue <- filter(de, PValue < 0.01)
  
  de.dw <- filter(de.pvalue, logFC > 1)
  de.up <- filter(de.pvalue, logFC < -1)

  
  if (nrow(de.up) > 20) {
    de.up <- de.up[1:20,]
  }
  
  if (nrow(de.dw) > 20) {
    de.dw <- de.dw[1:20,]
  }
  
  logRPM <- edgeR::cpm(y, log = TRUE)
  logRPM <- logRPM[c(rownames(de.up), rownames(de.dw)), ]
  logRPM.scale <- scale(t(logRPM), center = TRUE, scale = TRUE)
  logRPM.scale <- t(logRPM.scale)
  
  # 去除重复的行，保留第一个
  temp_rownames <- unlist(lapply(strsplit(rownames(logRPM.scale), "|", fixed = TRUE), function(x) x[3]))
  logRPM.scale <- logRPM.scale[!duplicated(temp_rownames),]
  genename <- rownames(logRPM.scale)
  rownames(logRPM.scale) <- unlist(lapply(strsplit(rownames(logRPM.scale), "|", fixed = TRUE), function(x) x[3]))
  gene_types <- sapply(strsplit(genename, "\\|"), function(x) x[4])
  
  if (nrow(annotation) > 150) {
    subtype <- annotation$group
  }
  subtype <- as.character(subtype)
  ann_col <- data.frame(
    class = as.character(subtype), 
    gender = as.character(annotation$gender),
    age = cut(annotation$age, breaks = seq(10, 90, by = 10), labels = paste(seq(10, 80, by = 10), seq(19, 89, by = 10), sep = "-"))
  )
  rownames(ann_col) <- colnames(logRPM.scale)
  
  ann_colors <- list()
  ann_colors$class <- setNames(class_colors, class_col)
  ann_colors$gender <- c("M" = "black", "F" = "white")
  age_groups <- paste(seq(10, 80, by = 10), seq(19, 89, by = 10), sep = "-")
  age_colors <- colorRampPalette(c("white", "darkblue"))(length(age_groups))
  ann_colors$age <- setNames(age_colors, age_groups)
  
  ann_row <- data.frame(gene_type = gene_types)
  rownames(ann_row) <- rownames(logRPM.scale)
  gene_type_colors <- brewer.pal(length(unique(gene_types)), "Set3")
  ann_colors$gene_type <- setNames(gene_type_colors, unique(gene_types))
  
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

group <- anno$group
group <- as.character(group)
class = brewer.pal(5,'BrBG')[c(1,2,5)]

group2 <- group
group2[group2 == "MDA5" | group2 == "ARS"] <- "DM"
class_colors <- class[c(1,2,3)]
class_col <- c("MDA5","ARS","HC")
p2 = deg_heatmap(mt ,group2, anno, class_colors, class_col)

mt3 <- mt[,c(1:105)]
group3 <- group[group != "HC"]
anno3 <- anno[c(1:105),]
class_colors <- class[c(1,2)]
class_col <- c("MDA5","ARS")
p3 = deg_heatmap(mt3 ,group3, anno3, class_colors, class_col)
#ggsave("plot/MDA5vsARS_total-heatmap.rmbatch-all.new.pdf",p3,width=7,height=5)#dev.off()

mt4 <- mt[,c(1:55,106:155)]
group4 <- group[group != "ARS"]
anno4 <- anno[c(1:55,106:155),]
class_colors <- class[c(1,3)]
class_col <- c("MDA5","HC")
p4 = deg_heatmap(mt4 ,group4, anno4, class_colors, class_col)
#ggsave("plot/MDA5vsHC_total-heatmap.rmbatch-all.new.pdf",p4,width=7,height=5.4)


mt5 <- mt[,c(56:155)]
group5 <- group[group != "MDA5"]
anno5 <- anno[c(56:155),]
class_colors <- class[c(2,3)]
class_col <- c("ARS","HC")
p5 = deg_heatmap(mt5 ,group5, anno5, class_colors, class_col)
p5
ggsave("plot/ARSvsHC_total-heatmap.rmbatch-all.new2.pdf",p5,width=7,height=6)


## ComplexHeatmap 暂时没用，作为备用
{
  library(ComplexHeatmap)
  library(circlize)
  {
    library(ggsci)
    col1 = brewer.pal(8,"Blues")
    col2 = brewer.pal(8,"YlOrRd")
    
    ann_colors = list(cell = c(brewer.pal(5,"Blues"), #B
                               brewer.pal(9,"Reds"),  #CD4
                               brewer.pal(4,"Greens"), #CD8
                               pal_igv()(8)[4],  #NK
                               brewer.pal(8,"PuBu")[c(5,8)],    #DC
                               brewer.pal(8,"Greys")[c(3,5,6,8)], #Monocyte
                               brewer.pal(8,"Purples")[c(5,8)], #Neutrophil
                               pal_igv()(10)[8], #platelet
                               pal_igv()(10)[9], #EV
                               pal_igv()(10)[10] # plasma
                               #brewer.pal(8,"YlOrRd")[3],
    ),
    Type = brewer.pal(8,"Set2")[1:3],
    group = brewer.pal(12,"Paired")[1:5])
    names(ann_colors$cell) = celltype
    names(ann_colors$Type) = unique(res.row$Type)
    names(ann_colors$group) = unique(res.row$Fun_V34)
    
    res.row <- data.frame(id=unique(select)) %>% left_join(., unique(res[,2:6]), by="id")
    res.row[1:4,]
    table(res.row$Type)
    table(res.row$Fun_V34)
    dim(res.row)
    
  }
  
  p <- Heatmap(logRPM.scale,
               col <- colorRampPalette(c("navy", "white", "firebrick3"))(50),
               #col <- colorRamp2(c(-4.3,0,4.3),c(col1[8],"white",col2[8])),
               heatmap_legend_param = list(grid_height = unit(10,'mm')),  #图例高度设置
               name = "-log(pvalue)",
               show_row_names = T,
               show_column_names = T,
               cluster_rows = T,
               cluster_columns = F,
               #column_order = celltype,
               #row_order = unique(select),
               #top_annotation = HeatmapAnnotation(simple_anno_size = unit(5, 'mm'),
               #                        Gender = adar.col$Gender
               top_annotation = HeatmapAnnotation(#simple_anno_size = unit(5, 'mm'),
                 Group=anno_text(ann_col$class,location=1, just="right"),
                 Type=res.row$Type,
                 #Repeat=res.row$Repeat,
                 #Gene=res.row$Gene_V34)
                 Gene=anno_text(res.row$Gene_V34,location=1, just="right"),
                 group=res.row$Fun_V34,
                 col = ann_colors
               ))
  p
}
####以前的脚本
{{
  ##DMvsHC
  {
    group <- anno$group
    group <- as.character(group)
    group2 <- group
    group2[group2 == "MDA5" | group2 == "ARS"] <- "DM"
    
    #mt1 <- mt[!grepl("^piR-hsa|protein_coding|pseudogene|lncRNA", rownames(mt)), ]
    #mt_pc <- mt[grep("protein_coding", rownames(mt)), ]
    #lncRNA_rows <- grepl("\\|lncRNA", rownames(mt))
    #mt_lncRNA <- mt[lncRNA_rows, ]
    #miRNA_rows <- grepl("\\|miRNA", rownames(mt))
    #mt_miRNA <- mt[miRNA_rows, ]
    #group4 <- group[c(1:44,56:105)]
    
    y <- DGEList(counts=mt, group=group2)
    y
    keep <- filterByExpr(y, group=group2, min.count=2, min.prop=0.2) # min.count=5, min.prop=0.5
    table(keep) #32670
    y <- y[keep, ,keep.lib.size=TRUE]
    y <- calcNormFactors(y, method="TMM")
    #design <- model.matrix(~group)
    design <- model.matrix(~group2)
  }
  ##MDA5vsARS
  {
    group <- anno$group
    group <- as.character(group)
    group3 <- group[group != "HC"]
    anno3 <- anno[c(1:105),]
    #mt1 <- mt[!grepl("^piR-hsa|protein_coding|pseudogene|lncRNA", rownames(mt)), ]
    #mt_pc <- mt[grep("protein_coding", rownames(mt)), ]
    #lncRNA_rows <- grepl("\\|lncRNA", rownames(mt))
    #mt_lncRNA <- mt[lncRNA_rows, ]
    #miRNA_rows <- grepl("\\|miRNA", rownames(mt))
    #mt_miRNA <- mt[miRNA_rows, ]
    #group4 <- group[c(1:44,56:105)]
    
    y <- DGEList(counts=mt[,1:105], group=group3)
    y
    keep <- filterByExpr(y, group=group3, min.count=2, min.prop=0.2) # min.count=5, min.prop=0.5
    table(keep) #30827
    y <- y[keep, ,keep.lib.size=TRUE]
    y <- calcNormFactors(y, method="TMM")
    #design <- model.matrix(~group)
    design <- model.matrix(~group3)
  }
  
  rownames(design) <- colnames(y)
  y <- estimateDisp(y, design)
  y$common.dispersion
  
  cpm <- edgeR::cpm(y)
  
  if(FALSE){
    fit.lr <- glmFit(y, design)
    lrt <- glmLRT(fit.lr, coef=2)
    #lrt <- glmLRT(fit.lr, contrast = c(-1,1))
    de <- topTags(lrt, n=Inf)$table
  }
  
  fit.ql <- glmQLFit(y, design)
  #####
  qlf <- glmQLFTest(fit.ql, coef=2)
  #anov <- glmQLFTest(fit.ql, contrast=con)
  #qlf <- glmQLFTest(fit.ql, coef=2)
  de <- topTags(qlf, n=Inf)$table
  de.pvalue <- filter(de, PValue< 0.01)#0.01
  #de.FDR <- filter(de.pvalue,FDR< 0.05)#0.05
  
  de.dw <- filter(de.pvalue, logFC> 1)#0 #group2 4.5 #group3 2        #19
  de.up <- filter(de.pvalue, logFC< -1)#0 #group2 -3.5 #group3 -1.5   #731
  #de.up <- de.up[c(1:30),]#画DMvsHC用到
  #注意这里logFC为负的是上调，因为默认后面的组是实验组，根据实际数据要调整这个数值
}
  # heatmap
  {
    logRPM <- edgeR::cpm(y, log=T)
    logRPM <- logRPM[c(rownames(de.up),rownames(de.dw)), ]
    
    dim(logRPM) # 48,42
    logRPM.scale <- scale(t(logRPM), center = T, scale = T)
    logRPM.scale <- t(logRPM.scale)
    #logRPM.scale <- logRPM.scale[-21,]#画DMvsHC用到，去掉重复的行：U3
    genename <- rownames(logRPM.scale)
    #rownames(logRPM.scale)[c(26:32,34,35,39:41,43) ] <- paste0("||", rownames(logRPM.scale)[c(26:32,34,35,39:41,43)], "|tRNA")
    rownames(logRPM.scale) <- unlist(lapply(strsplit(rownames(logRPM.scale),"|",fixed=T),function(x) x[3])) #paste(x[3], x[4], sep = "|")
    gene_types <- sapply(strsplit(genename, "\\|"), function(x) x[4])
    #group. <- group3[order(group3,decreasing=F)]
    #logRPM.scale. <- logRPM.scale[,order(group2,decreasing=F)]
    #group1 <- factor(group,levels=c("MDA5","ARS","HC"))
    #group. <- group[order(group1)]
    #logRPM.scale. <- logRPM.scale[,order(group3)]
    
    # pheatmap
    {
      library(pheatmap)
      library(RColorBrewer)
      #ann_col = data.frame(class = as.character(group3))######
      #rownames(ann_col) = colnames(logRPM.scale)
      #ann_colors = list(class = brewer.pal(5,'BrBG')[c(1,2)])
      #names(ann_colors$class) = c("MDA5","ARS")
      
      #######DMvsHC
      {
        ann_col = data.frame(
          class = as.character(group), 
          gender = as.character(anno$gender),
          age = cut(anno$age, breaks = seq(10, 90, by = 10), labels = paste(seq(10, 80, by=10), seq(19, 89, by=10), sep = "-"))
        )
        rownames(ann_col) = colnames(logRPM.scale)
        ann_colors = list(class = brewer.pal(5,'BrBG')[c(1,2,5)])
        names(ann_colors$class) = c("MDA5","ARS","HC")
        ann_colors$gender = c("M" = "black", "F" = "white")
        age_groups <- paste(seq(10, 80, by=10), seq(19, 89, by=10), sep = "-")
        age_colors <- colorRampPalette(c("white", "darkblue"))(length(age_groups))
        ann_colors$age = setNames(age_colors, age_groups)
      }
      #######MDA5vsARS
      {
        ann_col = data.frame(
          class = as.character(group3), 
          gender = as.character(anno3$gender),
          age = cut(anno3$age, breaks = seq(10, 90, by = 10), labels = paste(seq(10, 80, by=10), seq(19, 89, by=10), sep = "-"))
        )
        rownames(ann_col) = colnames(logRPM.scale)
        ann_colors = list(class = brewer.pal(5,'BrBG')[c(1,2)])
        names(ann_colors$class) = c("MDA5","ARS")
        ann_colors$gender = c("M" = "black", "F" = "white")
        age_groups <- paste(seq(10, 80, by=10), seq(19, 89, by=10), sep = "-")
        age_colors <- colorRampPalette(c("white", "darkblue"))(length(age_groups))
        ann_colors$age = setNames(age_colors, age_groups)
      }
      ann_row <- data.frame(
        gene_type = gene_types
      )
      rownames(ann_row) <- rownames(logRPM.scale)
      gene_type_colors <- brewer.pal(length(unique(gene_types)), "Set3")
      ann_colors$gene_type <- setNames(gene_type_colors, unique(gene_types))
      #GeneClass = c(Path1 = "#807DBA", Path2 = "#9E9AC8", Path3 = "#BCBDDC"))
      #"#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF" "#999999"
      col = brewer.pal(9,"Set1")
      #pdf("plot/pheatmap-pvalue0.001.pdf", width = 20, height = 16)
      par(mar=c(2,2,2,2))
      heatmap = pheatmap(logRPM.scale, 
                         #color = colorRampPalette(c(col[2],"white",col[1]))(1000),
                         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                         breaks = seq(-2, 2, length.out = 51),
                         cutree_col = 0, 
                         #cutree_row = 3, #break up the heatmap by clusters you define
                         cluster_rows = F,
                         cluster_cols = F, #by default, pheatmap clusters by both row and col
                         show_rownames = T,
                         show_colnames = F,
                         fontsize_row = 12,
                         angle_col = 90,
                         annotation_col = ann_col,
                         annotation_row = ann_row,
                         annotation_colors = ann_colors,
                         annotation_names_row = F,
                         #legend_breaks=-4:3,
                         border=F)
      heatmap
      
      #ggsave("plot/DMvsHC_total-heatmap.rmbatch-all.pdf",heatmap,width=10,height=12)#dev.off()
      #ggsave("plot/MDA5vsARS_total-heatmap.rmbatch-all.new.pdf",heatmap,width=7,height=5)#dev.off()
    }
  }
}
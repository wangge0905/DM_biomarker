library(edgeR)
library(ggplot2)
library(dplyr)
library(ggsignif)
library(RColorBrewer)

setwd("/Users/wangge/Documents/DM/")
mt <- read.table("totalRNA.matrix_rmbatch.all.txt", sep="\t", header=T, check.names=F)

# 计算 CPM
y <- DGEList(counts = mt)
subtype <- unlist(lapply(strsplit(colnames(mt), "-", fixed=T), function(x) x[1]))
subtype <- factor(subtype, levels = c("MDA5", "ARS", "HC"))
cpm <- edgeR::cpm(y)
cpm <- as.data.frame(cpm)

# 固定颜色映射
color_map <- brewer.pal(5, 'BrBG')[c(1, 2, 5)]
names(color_map) <- c("MDA5", "ARS", "HC")

# 绘图函数
group_gene_expression <- function(gene_ids, group = subtype, compare_groups, gene_order = NULL) {
  all_gene_expression <- data.frame()
  
  for (gene_id in gene_ids) {
    gene_expression <- as.data.frame(t(cpm[gene_id, ]))
    #gene_expression$Sample <- rownames(gene_expression) 
    colnames(gene_expression) <- "CPM"
    gene_expression$Gene <- unlist(lapply(strsplit(gene_id, "|", fixed = TRUE), function(x) x[3]))
    gene_expression$group <- subtype
    rownames(gene_expression) <- NULL
    
    all_gene_expression <- rbind(all_gene_expression, gene_expression)
  }
  
  if (!is.null(gene_order)) {
    all_gene_expression$Gene <- factor(all_gene_expression$Gene, levels = gene_order)
  }
  
  all_gene_expression <- all_gene_expression[all_gene_expression$group %in% compare_groups, ]
  
  
  p <- ggplot(all_gene_expression, aes(x = group, y = log(CPM + 1), fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    labs(y = "log(CPM + 1)") +
    xlab("") +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      panel.border = element_rect(colour = NA),
      axis.title.y = element_text(angle = 90, vjust = 2, size = 12),
      axis.text = element_text(size = 12),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      legend.position = "none",
      legend.title = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_blank()
    ) +
    scale_fill_manual(values = color_map) +
    geom_jitter(shape = 16, position = position_jitter(0.25)) +
    geom_signif(comparisons = list(compare_groups), step_increase = 0.1, vjust = 0.5, map_signif_level = TRUE) +
    facet_wrap(~ Gene, scales = "free_y", nrow = 1)
  p
  return(p)
}

# 根据对比组绘制组合图
# HC vs ARS
genes_hc_vs_ars <- c("ENSG00000206047.3|498|DEFA1|protein_coding", 
                     "ENSG00000207741.1|97|MIR590|miRNA", 
                     "ENSG00000237346.2|803|RP3-448I9.2|lncRNA", 
                     "ENSG00000207759.1|110|MIR181A1|miRNA")
gene_order_hc_vs_ars <- c("DEFA1", "MIR590", "MIR181A1", "RP3-448I9.2")
plot_hc_vs_ars <- group_gene_expression(genes_hc_vs_ars, compare_groups = c("HC", "ARS"), gene_order = gene_order_hc_vs_ars)
ggsave("plot/HC_vs_ARS_boxplot.pdf", plot_hc_vs_ars, width = 7, height = 2.5)

# HC vs MDA5
genes_hc_vs_mda5 <- c("ENSG00000206047.3|498|DEFA1|protein_coding", 
                      "ENSG00000207741.1|97|MIR590|miRNA", 
                      "ENSG00000199107.3|79|MIR409|miRNA", 
                      "ENSG00000207651.1|86|MIR28|miRNA", 
                      "ENSG00000202569.4|73|MIR146B|miRNA")
gene_order_hc_vs_mda5 <- c("DEFA1", "MIR590", "MIR28", "MIR409", "MIR146B")
plot_hc_vs_mda5 <- group_gene_expression(genes_hc_vs_mda5, compare_groups = c("HC", "MDA5"), gene_order = gene_order_hc_vs_mda5)
ggsave("plot/HC_vs_MDA5_boxplot.pdf", plot_hc_vs_mda5, width = 8, height = 2.5)

# MDA5 vs ARS
genes_mda5_vs_ars <- c("ENSG00000274641.2|467|H2BC17|protein_coding", 
                       "ENSG00000169429.11|2239|CXCL8|protein_coding", 
                       "ENSG00000187608.10|867|ISG15|protein_coding", 
                       "ENSG00000227121.2|1615|LINC02672|lncRNA", 
                       "ENSG00000231128.6|1582|RP5-1073O3.2|lncRNA",
                       "ENSG00000254251.1|833|RP11-662G23.1|lncRNA", 
                       "ENSG00000201271.1|164|RNU1-112P|snRNA")
gene_order_mda5_vs_ars <- c("H2BC17", "ISG15", "CXCL8", "LINC02672", "RP5-1073O3.2", "RP11-662G23.1", "RNU1-112P")
plot_mda5_vs_ars <- group_gene_expression(genes_mda5_vs_ars, compare_groups = c("MDA5", "ARS"), gene_order = gene_order_mda5_vs_ars)
ggsave("plot/MDA5_vs_ARS_boxplot.pdf", plot_mda5_vs_ars, width = 10.7, height = 2.8)


##### 没使用的
library(edgeR)
library(ggplot2)
library(dplyr)
library(ggsignif)
library(RColorBrewer)
setwd("/Users/wangge/Documents/DM/")
mt <- read.table("totalRNA.matrix_rmbatch.all.txt",sep="\t",header=T,check.names=F)

y <- DGEList(counts = mt)
libsize = y$samples
subtype <- unlist(lapply(strsplit(colnames(mt),"-",fixed=T),function(x) x[1]))
subtype <- factor(subtype, levels = c("MDA5","ARS","HC"))
cpm <- edgeR::cpm(y)
cpm <- as.data.frame(cpm)
###DEFA1 expression
{
  DEFA1_cpm <- as.data.frame(t(cpm["ENSG00000206047.3|498|DEFA1|protein_coding",]))
  #rownames(DEFA1_cpm) <- c("CPM")
  DEFA1_cpm$group <- subtype
  names(DEFA1_cpm)[names(DEFA1_cpm) == "ENSG00000206047.3|498|DEFA1|protein_coding"] <- "CPM"
  
  {
    a=ggplot(DEFA1_cpm,aes(x=group,y=log(CPM+1), fill=group)) +
      geom_boxplot(outlier.shape = NA) +
      labs(title="DEFA1",y="log(CPM+1)") +
      xlab("")+
      theme_bw()+
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
            #panel.background = element_rect(colour = NA),
            #plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title.y = element_text(angle = 90, vjust = 2, size = 16),
            #axis.title.x = element_text(vjust = -0.2, size = base_size),
            axis.text = element_text(size = 12),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line())+
      scale_fill_manual(values = brewer.pal(5,'BrBG')[c(1,2,5)])+
      scale_y_continuous(
        limits = c(-0.2,9),
        breaks = c(0,2,4,6,8),
        labels = c(0,2,4,6,8))+
      #geom_hline(aes(yintercept=.4, color="red"), linetype="dashed",show.legend = FALSE)+
      geom_jitter(shape=16, position = position_jitter(0.25),)+
      geom_signif(comparisons = list( c("ARS", "HC"),c("MDA5", "ARS"), c("MDA5", "HC")),step_increase = 0.1,
                  map_signif_level = TRUE)
    a
  }
  #ggsave("plot/DEFA1_expression.pdf",a,width = 4, height = 4)
  
}

plot_gene_expression <- function(gene_id, 
                                 expression_matrix = cpm, 
                                 group = subtype, 
                                 #save_dir = "plot", 
                                 #width = 4, 
                                 #height = 4, 
                                 compare_groups = c("MDA5", "ARS")) {  # 新增 compare_groups 参数
  
  gene_expression <- as.data.frame(t(cpm[gene_id, ]))
  colnames(gene_expression) <- "CPM"
  gene <- unlist(lapply(strsplit(gene_id, "|", fixed = TRUE), function(x) x[3]))
  gene_expression$group <- group
  
  # 筛选出需要比较的两组
  gene_expression <- gene_expression[gene_expression$group %in% compare_groups, ]
  
  # 设定固定的颜色映射
  color_map <- brewer.pal(5, 'BrBG')[c(1, 2, 5)]
  names(color_map) <- c("MDA5", "ARS", "HC")
  
  p <- ggplot(gene_expression, aes(x = group, y = log(CPM + 1), fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = gene, y = "log(CPM+1)") +
    xlab("") +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      panel.border = element_rect(colour = NA),
      axis.title.y = element_text(angle = 90, vjust = 2, size = 16),
      axis.text = element_text(size = 12),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line()
    ) +
    scale_fill_manual(values = color_map) +  # 使用固定的颜色映射
    scale_y_continuous(
      limits = c(-0.2, 5.5),
      breaks = c(0, 1, 2, 3, 4, 5),
      labels = c(0, 1, 2, 3, 4, 5)
    ) +
    geom_jitter(shape = 16, position = position_jitter(0.25)) +
    geom_signif(comparisons = list(compare_groups), 
                step_increase = 0.1, map_signif_level = TRUE)
  
  return(p)
}

# Plot for all pairwise comparisons for DEFA1

# Plot for HC vs ARS
plot1 <- plot_gene_expression("ENSG00000207741.1|97|MIR590|miRNA", compare_groups = c("HC", "ARS"))
plot2 <- plot_gene_expression("ENSG00000237346.2|803|RP3-448I9.2|lncRNA", compare_groups = c("HC", "ARS"))
plot13 <- plot_gene_expression("ENSG00000207759.1|110|MIR181A1|miRNA", compare_groups = c("HC", "ARS"))
plot18 <- plot_gene_expression("ENSG00000206047.3|498|DEFA1|protein_coding", compare_groups = c("HC", "ARS"))

# Plot for HC vs MDA5
plot3 <- plot_gene_expression("ENSG00000207651.1|86|MIR28|miRNA", compare_groups = c("HC", "MDA5"))
plot4 <- plot_gene_expression("ENSG00000199107.3|79|MIR409|miRNA", compare_groups = c("HC", "MDA5"))
plot5 <- plot_gene_expression("ENSG00000202569.4|73|MIR146B|miRNA", compare_groups = c("HC", "MDA5"))
plot14 <- plot_gene_expression("ENSG00000207741.1|97|MIR590|miRNA", compare_groups = c("HC", "MDA5"))
plot19 <- plot_gene_expression("ENSG00000206047.3|498|DEFA1|protein_coding", compare_groups = c("HC", "MDA5"))

# Plot for MDA5 vs ARS
plot9 <- plot_gene_expression("ENSG00000227121.2|1615|LINC02672|lncRNA", compare_groups = c("MDA5", "ARS"))
plot10 <- plot_gene_expression("ENSG00000201271.1|164|RNU1-112P|snRNA", compare_groups = c("MDA5", "ARS"))
plot11 <- plot_gene_expression("ENSG00000231128.6|1582|RP5-1073O3.2|lncRNA", compare_groups = c("MDA5", "ARS"))
plot12 <- plot_gene_expression("ENSG00000254251.1|833|RP11-662G23.1|lncRNA", compare_groups = c("MDA5", "ARS"))
plot6 <- plot_gene_expression("ENSG00000274641.2|467|H2BC17|protein_coding", compare_groups = c("MDA5", "ARS"))
plot7 <- plot_gene_expression("ENSG00000169429.11|2239|CXCL8|protein_coding", compare_groups = c("MDA5", "ARS"))
plot8 <- plot_gene_expression("ENSG00000187608.10|867|ISG15|protein_coding", compare_groups = c("MDA5", "ARS"))


###DEG barplot function wilcox.test
plot_gene_expression <- function(gene_id, expression_matrix = cpm, group = subtype, save_dir = "plot", width = 4, height = 4) {
  gene_expression <- as.data.frame(t(cpm[gene_id, ]))
  colnames(gene_expression) <- "CPM"
  gene <- unlist(lapply(strsplit(gene_id, "|", fixed = TRUE), function(x) x[3]))
  gene_expression$group <- group
  
  p <- ggplot(gene_expression, aes(x = group, y = log(CPM + 1), fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = gene, y = "log(CPM+1)") +
    xlab("") +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      panel.border = element_rect(colour = NA),
      axis.title.y = element_text(angle = 90, vjust = 2, size = 16),
      axis.text = element_text(size = 12),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line()
    ) +
    scale_fill_manual(values = brewer.pal(5, 'BrBG')[c(1, 2, 5)]) +
    scale_y_continuous(
      limits = c(-0.2, 4.5),
      breaks = c(0, 1, 2, 3, 4),
      labels = c(0, 1, 2, 3, 4)
    ) +
    geom_jitter(shape = 16, position = position_jitter(0.25)) +
    geom_signif(comparisons = list(c("ARS", "HC"), c("MDA5", "ARS"), c("MDA5", "HC")),
                step_increase = 0.1, map_signif_level = TRUE)
  
  #if (!is.null(save_dir)) {
  #  save_path <- sprintf("%s/%s_expression.pdf", save_dir, gene)
  #  ggsave(save_path, plot = p, width = width, height = height)
  #}
  
  return(p)
}
###函数使用
{
  #plot <- plot_gene_expression(gene_id, expression_matrix, subtype)
  #####
  #HCvsARS: DEFA1, MIR590, RP3-44819.2, MIR181A1
  #HCvsMDA5: DEFA1, MIR590, MIR28, MIR409, MIR146B
  #MDA5vsARS: H2BC17, CXCL8, ISG15, LINC02672, RNU1-112P, RP5-1073O3.2, RP11-662G23.1
  #ncRNA parameter:
  #scale_y_continuous(
  #  limits = c(-0.2, 4.5),
  #  breaks = c(0, 1, 2, 3, 4),
  #  labels = c(0, 1, 2, 3, 4)
  #  )
  plot1 <- plot_gene_expression("ENSG00000207741.1|97|MIR590|miRNA")
  plot2 <- plot_gene_expression("ENSG00000237346.2|803|RP3-448I9.2|lncRNA")
  plot3 <- plot_gene_expression("ENSG00000207651.1|86|MIR28|miRNA")
  plot4 <- plot_gene_expression("ENSG00000199107.3|79|MIR409|miRNA")
  plot5 <- plot_gene_expression("ENSG00000202569.4|73|MIR146B|miRNA")
  plot9 <- plot_gene_expression("ENSG00000227121.2|1615|LINC02672|lncRNA")
  plot10 <- plot_gene_expression("ENSG00000201271.1|164|RNU1-112P|snRNA")
  plot11 <- plot_gene_expression("ENSG00000231128.6|1582|RP5-1073O3.2|lncRNA")
  plot12 <- plot_gene_expression("ENSG00000254251.1|833|RP11-662G23.1|lncRNA")
  plot13 <- plot_gene_expression("ENSG00000207759.1|110|MIR181A1|miRNA")
  #protein_coding RNA parameter:
  #scale_y_continuous(
  #  limits = c(-0.2, 6.7),
  #  breaks = c(0, 2, 4, 6),
  #  labels = c(0, 2, 4, 6)
  #  )
  plot6 <- plot_gene_expression("ENSG00000274641.2|467|H2BC17|protein_coding")
  plot7 <- plot_gene_expression("ENSG00000169429.11|2239|CXCL8|protein_coding")
  plot8 <- plot_gene_expression("ENSG00000187608.10|867|ISG15|protein_coding")
}

plot13 <- plot_gene_expression("ENSG00000150593.18|4547|PDCD4|protein_coding")
###组合图
{
  p1 <- ggpubr::ggarrange(a,plot1,plot2,plot13, nrow = 1, ncol = 4, #labels = c('A'), 
                          font.label = list(color = 'black'), 
                          common.legend = T, 
                          legend = "right")
  p1
  ggsave("plot/HCvsARS_ML_gene.pdf", plot = p1, width = 9, height = 3)
  
  p2 <- ggpubr::ggarrange(a,plot1,plot3,plot4,plot5, nrow = 1, ncol = 5, #labels = c('A'), 
                          font.label = list(color = 'black'), 
                          common.legend = T, 
                          legend = "right")
  p2
  #ggsave("plot/HCvsMDA5_ML_gene.pdf", plot = p2, width = 11.5, height = 3)
  
  p3 <- ggpubr::ggarrange(plot6,plot7,plot8,plot9,plot10,plot11,plot12, nrow = 2, ncol = 5, #labels = c('A'), 
                          font.label = list(color = 'black'), 
                          common.legend = T, 
                          legend = "right")
  p3
  ggsave("plot/MDA5vsARS_ML_gene3.pdf", plot = p3, width = 11.5, height = 5.5)
  
}
##### 没使用的小提琴图
{
  summary_stats <- DEFA1_cpm %>%
    group_by(subtype) %>%
    summarize(mean_y = mean(log(CPM+1)),
              q25 = quantile(log(CPM+1), 0.25),
              q75 = quantile(log(CPM+1), 0.75))
  p3 <- ggplot(DEFA1_cpm,aes(x=subtype,y=log(CPM+1), fill=subtype)) +
    geom_violin(show.legend = TRUE, width = 0.8, scale = "width") +
    geom_segment(data = summary_stats, aes(x = subtype, xend=subtype, y = q25, yend = q75), color = "black") +
    # 在竖线上画点表示均值
    geom_point(data = summary_stats, aes(x = subtype, y = mean_y), color = "black", size = 3) +
    labs(title="miRNA species",y="Species number") +
    xlab("")+
    theme_bw()+
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
          #panel.background = element_rect(colour = NA),
          #plot.background = element_rect(colour = NA),
          panel.border = element_rect(colour = NA),
          axis.title.y = element_text(angle = 90, vjust = 2, size = 16),
          #axis.title.x = element_text(vjust = -0.2, size = base_size),
          axis.text = element_text(size = 12),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line())+
    scale_fill_manual(values = brewer.pal(5,'BrBG')[c(1,2,5)])+
    scale_y_continuous(
      limits = c(0,6),
      breaks = c(1,3,5),
      labels = c(1,3,5))
  #geom_hline(aes(yintercept=6, color="red"), linetype="dashed",show.legend = FALSE)
  p3
  
}



#### GO & KEGG
KEGG_GO_Expression <- function(Differential_result,
                               output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                               pvalue_cutoff,log2Foldchange_cutoff){
  #all_gene <- row.names(Differential_result)
  #strsplit is highly depend on the colnames of alteration matrix
  #all_gene <- as.character(lapply(strsplit(all_gene,"\\."),function(x) x[1]))
  #background <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
  #                    filters = "ensembl_gene_id",
  #                    values=all_gene, mart= mart,useCache = FALSE)
  #background <- background[-which(background$entrezgene_id=="")]
  #background <- background$entrezgene_id
  
  filtered_down <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange < -log2Foldchange_cutoff),])
  filtered_down <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
  forenrich_down <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                          filters = "ensembl_gene_id",
                          values=filtered_down, mart= mart,useCache = FALSE)
  if(length(which(forenrich_down$entrezgene_id==""))==0){
    forenrich_down <- forenrich_down
  } else {
    forenrich_down <- forenrich_down[-which(forenrich_down$entrezgene_id==""),]
  }
  forenrich_down <- forenrich_down$entrezgene_id
  
  filtered_up <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange > log2Foldchange_cutoff),])
  filtered_up <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
  forenrich_up <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                        filters = "ensembl_gene_id",
                        values=filtered_up, mart= mart,useCache = FALSE)
  if(length(which(forenrich_up$entrezgene_id==""))==0){
    forenrich_up <- forenrich_up
  } else {
    forenrich_up <- forenrich_up[-which(forenrich_up$entrezgene_id==""),]
  }
  forenrich_up <- forenrich_up$entrezgene_id
  
  #KEGG
  {
    KEGG_res_down <- enrichKEGG(
      forenrich_down,
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1,
      use_internal_data = FALSE)
    KEGG_output_down <- KEGG_res_down@result
    KEGG_output_down$GeneEnrichedIn <- "Down regulated"
    
    KEGG_res_up <- enrichKEGG(
      forenrich_up,
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1,
      use_internal_data = FALSE)
    KEGG_output_up <- KEGG_res_up@result
    KEGG_output_up$GeneEnrichedIn <- "Up regulated"
    
    KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
    write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
  }
  #GO_BP
  {
    library(clusterProfiler)
    GO_res_down <- enrichGO(
      forenrich_down,
      'org.Hs.eg.db',
      ont = "BP",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1)
    GO_output_down <- GO_res_down@result
    GO_output_down$GeneEnrichedIn <- "Down regulated"
    
    GO_res_up <- enrichGO(
      forenrich_up,
      'org.Hs.eg.db',
      ont = "BP",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1)
    GO_output_up <- GO_res_up@result
    GO_output_up$GeneEnrichedIn <- "Up regulated"
    
    GO_output <- rbind(GO_output_up,GO_output_down)
    
    write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
  }
  #GO_CC
  {
    library(clusterProfiler)
    GO_res_down <- enrichGO(
      forenrich_down,
      'org.Hs.eg.db',
      ont = "CC",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1)
    GO_output_down <- GO_res_down@result
    GO_output_down$GeneEnrichedIn <- "Down regulated"
    
    GO_res_up <- enrichGO(
      forenrich_up,
      'org.Hs.eg.db',
      ont = "CC",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1)
    GO_output_up <- GO_res_up@result
    GO_output_up$GeneEnrichedIn <- "Up regulated"
    
    GO_output <- rbind(GO_output_up,GO_output_down)
    
    write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
  }
  #GO_MF
  {
    library(clusterProfiler)
    GO_res_down <- enrichGO(
      forenrich_down,
      'org.Hs.eg.db',
      ont = "MF",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1)
    GO_output_down <- GO_res_down@result
    GO_output_down$GeneEnrichedIn <- "Down regulated"
    
    GO_res_up <- enrichGO(
      forenrich_up,
      'org.Hs.eg.db',
      ont = "MF",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1)
    GO_output_up <- GO_res_up@result
    GO_output_up$GeneEnrichedIn <- "Up regulated"
    
    GO_output <- rbind(GO_output_up,GO_output_down)
    
    write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
  }
}

#function for enrichment analysis for miRNA: clusterprofiler
KEGG_GO_miRNA_target <- function(Differential_result,
                                 output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                                 pvalue_cutoff,log2Foldchange_cutoff){
  #all_tx <- row.names(Differential_result)
  #strsplit is highly depend on the colnames of alteration matrix
  #all_tx <- as.character(lapply(strsplit(all_tx,"\\."),function(x) x[1]))
  #all_mir <- getBM(attributes=c("ensembl_transcript_id", "mirbase_id"),
  #                         filters = "ensembl_transcript_id",
  #                         values=all_tx, mart= mart,useCache = FALSE)
  #all_mir <- all_mir[-which(all_mir$mirbase_id=="")]
  #background <- get_multimir(mirna = all_mir$mirbase_id, summary = TRUE)
  #background <- unique(background@data$target_entrez)
  
  filtered_down <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange < -log2Foldchange_cutoff),])
  filtered_down <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
  down_tx <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
  down_mir <- getBM(attributes=c("ensembl_transcript_id", "mirbase_id"),
                    filters = "ensembl_transcript_id",
                    values=down_tx, mart= mart,useCache = FALSE)
  if(length(which(down_mir$entrezgene_id==""))==0){
    down_mir <- down_mir
  } else {
    down_mir <- down_mir[-which(down_mir$mirbase_id==""),]
  }
  forenrich_down <- get_multimir(mirna = down_mir$mirbase_id, summary = TRUE)
  forenrich_down <- unique(forenrich_down@data$target_entrez)
  
  filtered_up <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange > log2Foldchange_cutoff),])
  filtered_up <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
  up_tx <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
  up_mir <- getBM(attributes=c("ensembl_transcript_id", "mirbase_id"),
                  filters = "ensembl_transcript_id",
                  values=up_tx, mart= mart,useCache = FALSE)
  if(length(which(up_mir$entrezgene_id==""))==0){
    up_mir <- up_mir
  } else {
    up_mir <- up_mir[-which(up_mir$mirbase_id==""),]
  }
  forenrich_up <- get_multimir(mirna = up_mir$mirbase_id, summary = TRUE)
  forenrich_up <- unique(forenrich_up@data$target_entrez)
  
  #KEGG
  {
    KEGG_res_down <- enrichKEGG(
      forenrich_down,
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1,
      use_internal_data = FALSE)
    KEGG_output_down <- KEGG_res_down@result
    KEGG_output_down$GeneEnrichedIn <- "Down regulated"
    
    KEGG_res_up <- enrichKEGG(
      forenrich_up,
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1,
      use_internal_data = FALSE)
    KEGG_output_up <- KEGG_res_up@result
    KEGG_output_up$GeneEnrichedIn <- "Up regulated"
    
    KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
    write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
  }
  #GO_BP
  {
    library(clusterProfiler)
    GO_res_down <- enrichGO(
      forenrich_down,
      'org.Hs.eg.db',
      ont = "BP",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1)
    GO_output_down <- GO_res_down@result
    GO_output_down$GeneEnrichedIn <- "Down regulated"
    
    GO_res_up <- enrichGO(
      forenrich_up,
      'org.Hs.eg.db',
      ont = "BP",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1)
    GO_output_up <- GO_res_up@result
    GO_output_up$GeneEnrichedIn <- "Up regulated"
    
    GO_output <- rbind(GO_output_up,GO_output_down)
    
    write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
  }
  #GO_CC
  {
    library(clusterProfiler)
    GO_res_down <- enrichGO(
      forenrich_down,
      'org.Hs.eg.db',
      ont = "CC",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1)
    GO_output_down <- GO_res_down@result
    GO_output_down$GeneEnrichedIn <- "Down regulated"
    
    GO_res_up <- enrichGO(
      forenrich_up,
      'org.Hs.eg.db',
      ont = "CC",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1)
    GO_output_up <- GO_res_up@result
    GO_output_up$GeneEnrichedIn <- "Up regulated"
    
    GO_output <- rbind(GO_output_up,GO_output_down)
    
    write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
  }
  #GO_MF
  {
    library(clusterProfiler)
    GO_res_down <- enrichGO(
      forenrich_down,
      'org.Hs.eg.db',
      ont = "MF",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1)
    GO_output_down <- GO_res_down@result
    GO_output_down$GeneEnrichedIn <- "Down regulated"
    
    GO_res_up <- enrichGO(
      forenrich_up,
      'org.Hs.eg.db',
      ont = "MF",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1)
    GO_output_up <- GO_res_up@result
    GO_output_up$GeneEnrichedIn <- "Up regulated"
    
    GO_output <- rbind(GO_output_up,GO_output_down)
    
    write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
  }
}

KEGG_GO_miRNA <- function(Differential_result,
                          output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                          pvalue_cutoff,log2Foldchange_cutoff){
  #all_tx <- row.names(Differential_result)
  #strsplit is highly depend on the colnames of alteration matrix
  #all_tx <- as.character(lapply(strsplit(all_tx,"\\."),function(x) x[1]))
  #all_mir <- getBM(attributes=c("ensembl_transcript_id", "mirbase_id"),
  #                         filters = "ensembl_transcript_id",
  #                         values=all_tx, mart= mart,useCache = FALSE)
  #all_mir <- all_mir[-which(all_mir$mirbase_id=="")]
  #background <- get_multimir(mirna = all_mir$mirbase_id, summary = TRUE)
  #background <- unique(background@data$target_entrez)
  
  filtered_down <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange < -log2Foldchange_cutoff),])
  filtered_down <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
  down_tx <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
  down_mir <- getBM(attributes=c("ensembl_transcript_id", "entrezgene_id"),
                    filters = "ensembl_transcript_id",
                    values=down_tx, mart= mart,useCache = FALSE)
  forenrich_down <- down_mir$entrezgene_id
  
  filtered_up <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange > log2Foldchange_cutoff),])
  filtered_up <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
  up_tx <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
  up_mir <- getBM(attributes=c("ensembl_transcript_id", "entrezgene_id"),
                  filters = "ensembl_transcript_id",
                  values=up_tx, mart= mart,useCache = FALSE)
  forenrich_up <- up_mir$entrezgene_id
  
  #KEGG
  {
    KEGG_res_down <- enrichKEGG(
      forenrich_down,
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1,
      use_internal_data = FALSE)
    KEGG_output_down <- KEGG_res_down@result
    KEGG_output_down$GeneEnrichedIn <- "Down regulated"
    
    KEGG_res_up <- enrichKEGG(
      forenrich_up,
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1,
      use_internal_data = FALSE)
    KEGG_output_up <- KEGG_res_up@result
    KEGG_output_up$GeneEnrichedIn <- "Up regulated"
    
    KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
    write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
  }
  #GO_BP
  {
    library(clusterProfiler)
    GO_res_down <- enrichGO(
      forenrich_down,
      'org.Hs.eg.db',
      ont = "BP",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1)
    GO_output_down <- GO_res_down@result
    GO_output_down$GeneEnrichedIn <- "Down regulated"
    
    GO_res_up <- enrichGO(
      forenrich_up,
      'org.Hs.eg.db',
      ont = "BP",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1)
    GO_output_up <- GO_res_up@result
    GO_output_up$GeneEnrichedIn <- "Up regulated"
    
    GO_output <- rbind(GO_output_up,GO_output_down)
    
    write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
  }
  #GO_CC
  {
    library(clusterProfiler)
    GO_res_down <- enrichGO(
      forenrich_down,
      'org.Hs.eg.db',
      ont = "CC",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1)
    GO_output_down <- GO_res_down@result
    GO_output_down$GeneEnrichedIn <- "Down regulated"
    
    GO_res_up <- enrichGO(
      forenrich_up,
      'org.Hs.eg.db',
      ont = "CC",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1)
    GO_output_up <- GO_res_up@result
    GO_output_up$GeneEnrichedIn <- "Up regulated"
    
    GO_output <- rbind(GO_output_up,GO_output_down)
    
    write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
  }
  #GO_MF
  {
    library(clusterProfiler)
    GO_res_down <- enrichGO(
      forenrich_down,
      'org.Hs.eg.db',
      ont = "MF",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1)
    GO_output_down <- GO_res_down@result
    GO_output_down$GeneEnrichedIn <- "Down regulated"
    
    GO_res_up <- enrichGO(
      forenrich_up,
      'org.Hs.eg.db',
      ont = "MF",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      #universe=as.character(background),
      minGSSize = 0,
      maxGSSize = 500,
      qvalueCutoff = 1)
    GO_output_up <- GO_res_up@result
    GO_output_up$GeneEnrichedIn <- "Up regulated"
    
    GO_output <- rbind(GO_output_up,GO_output_down)
    
    write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
  }
}
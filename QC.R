# RStudio工作目录设为本文件夹
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(RColorBrewer)
library(scales)
library(patchwork)
library(edgeR)
set.seed(20260929)
dir.create("generated/QC",recursive=TRUE,showWarnings=FALSE)
{
  qc <- read.table("inputs/1204/QC.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F)
  #colnames(qc)[2:21] <- paste0("230821-",colnames(qc)[2:21])
  colnames(qc)[36:55] <- paste0(substr(colnames(qc)[36:55],1,12),substr(colnames(qc)[36:55],16,20))
  qc <- qc[, -c(41:45, 51:55)]
  
  anno <- read.table("inputs/annotation.txt", sep = "\t", header = T)
  colnames(anno)[2] <- "library_id"
  anno <- filter(anno, library_id %in% colnames(qc))
  reshapee <- function(qc){
    library_id = colnames(qc)[-1]
    qc <- qc[-c(32:35),]#去掉32到35行
    qc <- data.frame(t(qc), stringsAsFactors = F)
    colnames(qc) <- qc[1,]
    qc <- qc[-1,]
    qc <- as.data.frame(lapply(qc,as.numeric))
    rownames(qc) <- anno[match(library_id,anno$library_id),"sample_id"]
    return(qc)
  }
  #注意这里如果library_id和anno$library_id不匹配的话，无法执行
  qc <- reshapee(qc)
  group_order <- c("MDA5", "ARS", "HC")
  # 使用 arrange 函数按照 group 列的顺序*排序*数据框
  anno <- anno %>%
    arrange(factor(group, levels = group_order))
  qc <- qc[anno$sample_id,]
  #相当于sort了一遍，按anno顺序排序
  qc$group <- anno$group
  
  
  # 将分组变量转换为有序因子，指定新的顺序
  #a <- factor(qc$group, levels = c("MDA5", "ARS", "HC"))
  #qc$group <- qc$group[order(a)]
  #qc$sample_id <- rownames(qc)
  #qc$sample_id <- factor(qc$sample_id, levels=c(paste0("MDA5-",c(1:55)),paste0("ARS-",c(1:50)),paste0("HC-",c(1:11,13:55))))
  #qc <- qc[order(factor(qc$sample_id,levels = c(paste("HC", 1:11, sep = "-"), paste("HC", 12:55, sep = "-"), paste("MDA5", 1:55, sep = "-"), paste("ARS", 1:50, sep = "-")))),]
  
}

{
  qc_sample <- qc %>%
    transmute(sample_id = rownames(qc),
              clean = clean,
              `genome` = star_hg38_v38,  
              `genome_nondup` = star_hg38_v38_dedup,
              rRNA = star_rRNA,
              group=group
    )
  qc_sample$log_clean <- log10(qc_sample$clean)
  p1 <- ggplot(qc_sample,aes(x=group,y=log_clean, fill=group)) +
    geom_boxplot(outlier.shape = NA,show.legend = FALSE) +
    labs(title="Clean reads > 1M",y="Reads number") +
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
      limits = c(5,8.2),
      breaks = c(5,6,7,8),
      labels = c(expression(10^5),expression(10^6), expression(10^7), expression(10^8)))+
    geom_hline(aes(yintercept=6, color="red"), linetype="dashed",show.legend = FALSE)+
    geom_jitter(shape=16, position = position_jitter(width=0.25,height=0,seed=20260929),show.legend = FALSE) 
  p1
  
#=========================================================================  
  qc_sample$log_usable <- log10(qc_sample$`genome_nondup`)
  p2 <- ggplot(qc_sample,aes(x=group,y=log_usable, fill=group)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title="Usable reads > 0.3M",y="Reads number") +
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
      limits = c(5,8.2),
      breaks = c(5,6,7,8),
      labels = c(expression(10^5),expression(10^6), expression(10^7), expression(10^8)))+
    geom_hline(aes(yintercept=log10(300000), color="red"), linetype="dashed",show.legend = FALSE)+
    geom_jitter(shape=16, position = position_jitter(width=0.25,height=0,seed=20260929),) 
  p2
  
  # 将两张图拼在一起，共用相同的y轴
  combined_plot <- p1 + p2 + plot_layout(nrow = 1, widths = c(0.5, 0.5))
  print(combined_plot)
  ggsave("generated/QC/clean&usable.pdf", combined_plot,width=7,height=4)
}

qc_filtered=qc[!rownames(qc) %in% c("HC-50","HC-51","HC-52","HC-53"),]
qc_filtered$group=factor(qc_filtered$group,levels=c("MDA5","ARS","HC"))
qc_filtered$rRNA_ratio <- qc_filtered$star_rRNA / qc_filtered$clean
ggplot(qc_filtered,aes(x=group,y=rRNA_ratio, fill=group)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title="rRNA ratio",y="Reads ratio(%)") +
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
    limits = c(0,0.6),
    breaks = c(0.1,0.2,0.3,0.4,0.5),
    labels = c('10%','20%','30%','40%','50%'))+
  geom_hline(aes(yintercept=0.4, color="red"), linetype="dashed",show.legend = FALSE)+
  geom_jitter(shape=16, position = position_jitter(width=0.25,height=0,seed=20260929),) 
ggsave("generated/QC/rRNA_ratio_barplot.pdf",width =4 ,height = 4)

mt=read.table("inputs/1204/totalRNA.matrix.txt",sep="\t",header=TRUE,check.names=FALSE)
mt=mt[rowSums(mt)!=0,]
rna.type=read.table("inputs/1204/totalRNA.matrix.type.txt",sep="\t",header=TRUE)
rownames(rna.type)=rna.type[[1]]
cpm=as.data.frame(edgeR::cpm(mt))
  all <- data.frame(genename=rownames(cpm))
  all$type <- rna.type[all$genename,]$type2
  cpm_miRNA <- cpm[all$type == "miRNA",]#也可以是其他RNA
  # 创建一个新的数据框来存储结果
  count_miRNA <- data.frame(rownames = character(0), number = numeric(0))
  # 遍历数据框的每一列，计算大于一的数值个数，并将结果存储在新数据框中
  for (col in colnames(cpm_miRNA)) {
    count <- sum(cpm_miRNA[, col] > 1)
    count_miRNA <- rbind(count_miRNA, data.frame(rownames = col, number = count))
  }
  count_miRNA$group <- sub("^([^_]+)-.*", "\\1", count_miRNA$rownames)
  count_miRNA$group <- factor(count_miRNA$group,levels = c("MDA5","ARS","HC"))
  count_miRNA <- count_miRNA[order(count_miRNA$group),]
  summary_stats <- count_miRNA %>%
    group_by(group) %>%
    summarize(mean_y = mean(number),
              q25 = quantile(number, 0.25),
              q75 = quantile(number, 0.75))
  p3 <- ggplot(count_miRNA,aes(x=group,y=number, fill=group)) +
    geom_violin(show.legend = TRUE) +
    geom_segment(data = summary_stats, aes(x = group, xend=group, y = q25, yend = q75), color = "black") +
    # 在竖线上画点表示均值
    geom_point(data = summary_stats, aes(x = group, y = mean_y), color = "black", size = 3) +
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
      limits = c(0,400),
      breaks = c(100,200,300),
      labels = c(100,200,300))
    #geom_hline(aes(yintercept=6, color="red"), linetype="dashed",show.legend = FALSE)
  p3
  ggsave("generated/QC/miRNA-species.pdf",p3,width = 4.3, height = 3.7)
ggsave("generated/QC/clean_reads.pdf",p1,width=3.5,height=4)
ggsave("generated/QC/usable_reads.pdf",p2,width=3.5,height=4)

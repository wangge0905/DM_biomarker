### wangge 
### 20240110

library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(RColorBrewer)
library(scales)
library(grid)
library(cowplot)
library(tidyverse)

library(patchwork)

options(scipen = 2)

setwd("/Users/wangge/Documents/DM/")


######################################
# 1.1 total QC
{
  qc <- read.table("1204/QC.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F)
  #colnames(qc)[2:21] <- paste0("230821-",colnames(qc)[2:21])
  colnames(qc)[36:55] <- paste0(substr(colnames(qc)[36:55],1,12),substr(colnames(qc)[36:55],16,20))
  qc <- qc[, -c(41:45, 51:55)]
  
  anno <- read.table("annotation.txt", sep = "\t", header = T)
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


##################################### 
# 1.2 clean/genome/rna/stacked plot
{
  #clean reads and ratio
  clean_reads <- qc %>%
    transmute(sample_id = rownames(qc),
              clean = clean,
              genome = star_hg38_v38 + star_circRNA,
              rRNA = star_rRNA,
              microbiome = unmapped - unclassified,
              unmapped = clean-genome-rRNA-microbiome,
              group=group)
  clean_ratio <- clean_reads %>%
    transmute(sample_id = sample_id,
              genome = genome / clean,
              rRNA = rRNA / clean,
              microbiome = microbiome / clean,
              unmapped = unmapped / clean,
              group=group)
  
  clean_ratio$genome <- as.numeric(clean_ratio$genome)
  clean_ratio <- clean_ratio %>%
    arrange(desc(genome))
  clean_ratio <- melt(clean_ratio,id=c("sample_id","group"))
  clean_ratio$sample_id <- factor(clean_ratio$sample_id, levels=anno$sample_id)
  #clean_ratio$sample_id <- factor(clean_ratio$sample_id,levels = unique(clean_ratio$sample_id))
  {
    col1 = pal_jco()(8)
    col2 = pal_npg()(8)
    col4=pal_igv()(8)
    a=ggplot(clean_ratio, aes(x=sample_id, y=value, fill=variable))+
      geom_bar(position = "fill", stat = "identity", width=0.8)+
      #scale_fill_manual(values=c(col1[5],col1[6],col2[7],col2[5],col1[4],col1[3]))+
      scale_fill_manual(values=c(col1[5],col1[4],col2[5],col1[3]))+
      scale_y_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1))+
      theme_bw()+
      guides(fill=guide_legend(title=NULL, ncol=1))+
      theme(
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        plot.title = element_text(color="black", size=12, hjust=0.5),
        #legend.position="top",
        #legend.direction = "horizontal",
        #legend.margin = margin(c(50,0,50,0)),
        legend.text= element_text(color="black", size=12),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length.y = unit(0.2,"cm"),
        axis.ticks.length.x = unit(0.1,"cm"),
        axis.text.x = element_text(color="black", size=5,angle=90,hjust=1,vjust=0.5),
        axis.text.y = element_text(color="black", size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
      labs(x="",y="",title="", face="bold")#+
    a
  }
  
  genome_reads <- qc %>%
    transmute(sample_id = rownames(qc),
              `genome_non-dup` = star_hg38_v38_dedup,  #+ circRNA_dedup,
              exonic = (MT_mRNA+MT_tRNA+chrM.all+mRNA+lncRNA+snoRNA+snRNA+srpRNA+tRNA+tucpRNA+Y_RNA+misc_RNA+pseudogene+exon), #+circRNA_dedup
              intronic = intron,
              intergenic = intergenic)
  genome_ratio <- genome_reads %>%
    transmute(sample_id = sample_id,
              exonic = exonic / `genome_non-dup`,
              intronic = intronic / `genome_non-dup`,
              intergenic = intergenic /`genome_non-dup`)
  genome_ratio <- melt(genome_ratio,id=c("sample_id"))
  genome_ratio$sample_id <- factor(genome_ratio$sample_id, levels=anno$sample_id)
  
  {
    col1 = pal_jco()(8)
    b=ggplot(genome_ratio, aes(x=sample_id, y=value, fill=variable))+
      geom_bar(position = "fill", stat = "identity",width=0.8)+
      theme_bw()+
      scale_fill_manual(values=c(col1[1],"#EEB4B4","#8B6969"))+
      scale_y_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1))+
      guides(fill=guide_legend(title=NULL, ncol=1))+      
      theme(
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        plot.title = element_text(color="black", size=12, hjust=0.5),
        legend.margin = margin(c(30,12,30,10)),
        #legend.position="top",
        #legend.box = "horizontal",
        legend.text= element_text(color="black", size=12),
        panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.ticks.length.y = unit(0.2,"cm"),
        axis.ticks.length.x = unit(0.1,"cm"),
        axis.text.x = element_text(color="black", size=5,angle=90,hjust=1,vjust=0.5),
        axis.text.y = element_text(color="black", size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
      labs(x="",y="",title="", face="bold")
    
    b
  }
  
  p=plot_grid(a,b,ncol=1,nrow=2,axis="l",rel_heights = c(1,1),align="v")
  p
  ggsave("plot/QC/qc_clean_genome_ratio.pdf",p,width=12,height=7)
  
}
{
  #rna reads and ratio
  rna_reads <- qc %>%
    transmute(sample_id=rownames(qc),
              mRNA=mRNA,
              lncRNA=lncRNA,
              #circRNA=circRNA_dedup,
              pseudogene=pseudogene,
              misc=misc_RNA,
              sncRNA=(snoRNA + snRNA + srpRNA + tRNA + Y_RNA),
              tucpRNA=tucpRNA,
              exon=exon,  
              total=mRNA+lncRNA+pseudogene+misc+sncRNA+tucpRNA+exon)
  #total=total)
  rna_ratio <- rna_reads %>%
    transmute(sample_id=sample_id,
              mRNA=mRNA / total,
              lncRNA=lncRNA / total,
              #circRNA=circRNA / total,
              pseudogene=pseudogene / total,
              misc=misc / total,
              sncRNA=sncRNA / total,
              tucpRNA=tucpRNA /total,
              other_exon=exon / total)
  
  rna_ratio <- melt(rna_ratio,id=c("sample_id"))
  rna_ratio$sample_id <- factor(rna_ratio$sample_id, levels=metadata$sample_id)
  
  {
    c=ggplot(rna_ratio, aes(x=sample_id, y=value, fill=variable))+
      geom_bar(position = "fill", stat = "identity",width=0.8)+
      theme_bw()+
      scale_fill_igv()+
      scale_y_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1))+
      guides(fill=guide_legend(title=NULL, ncol=1))+      
      theme(
        plot.title = element_text(color="black", size=12, hjust=0.5),
        strip.background = element_rect(size=1, fill="transparent"),
        strip.text = element_text(color="black", size=12),
        #legend.position="top",
        legend.text= element_text(color="black", size=12),
        legend.key.size = unit(0.2,"inches"),
        panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.ticks.length.y = unit(0.2,"cm"),
        axis.ticks.length.x = unit(0,"cm"),
        axis.text.x = element_text(color="black", size=7,angle=90,hjust=1,vjust=0.5),
        axis.text.y = element_text(color="black", size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
      labs(x="",y="",title="", face="bold")
    c
    ggsave("plot/cfev_qc_exon_ratio.pdf",c,width=12,height=4)
  }
}

####################################
# 1.3 boxplot
{
  #qc summary
  qc_reads <- qc %>%
    transmute(sample_id = rownames(qc),
              raw = bcumi,
              clean = clean,
              `genome` = star_hg38_v38,  
              `genome_nondup` = star_hg38_v38_dedup,
              exonic = (mRNA+lncRNA+snoRNA+snRNA+srpRNA+tRNA+tucpRNA+Y_RNA+misc_RNA+pseudogene+exon), #+circRNA_dedup
              microbiome = unmapped - unclassified,
              rRNA = star_rRNA
              #group=group
    )
  qc_reads <- melt(qc_reads,id.vars = c("sample_id"))#,"group"))
  
  qc_ratio <- qc %>%
    transmute(sample_id = rownames(qc),
              clean = clean / (bcumi*2),
              `genome` =  star_hg38_v38 / qc$clean, 
              `genome_nondup` = star_hg38_v38_dedup / star_hg38_v38,
              exonic = (mRNA+lncRNA+snoRNA+snRNA+srpRNA+tRNA+tucpRNA+Y_RNA+misc_RNA+pseudogene+exon) / star_hg38_v38_dedup, #+circRNA_dedup
              microbiome = (unmapped - unclassified) / qc$clean,
              rRNA = star_rRNA / qc$clean
              #group=group
    )
  qc_ratio <- melt(qc_ratio,id.vars = c("sample_id")) #,"group"))
  
  
  
  # plot
  {
    summary_stats <- qc_reads %>%
      group_by(variable) %>%
      summarize(mean_y = mean(value),
                q25 = quantile(value, 0.25),
                q75 = quantile(value, 0.75))
    p <- ggplot(qc_reads, aes(variable, value)) +
      geom_violin(trim = FALSE, alpha = 0.9, scale = "width", position = position_dodge(0.9),aes(fill = variable)) +
      geom_segment(data = summary_stats, aes(x = variable, xend=variable, y = q25, yend = q75), color = "black") +
      # 在竖线上画点表示均值
      geom_point(data = summary_stats, aes(x = variable, y = mean_y), color = "black", size = 3) +
      #geom_point(aes(color = variable), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9), alpha = 0.7, size = 2) +
      scale_y_log10(limits = c(1 * 10^4, 5 * 10^7), labels = trans_format("log10", math_format(10^.x))) +
      theme_bw() +
      scale_fill_nejm() +
      theme(
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5),
        strip.text = element_text(size = 12, color = "black"),
        strip.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "none",
        legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12)
      ) +
      ylab("number of reads") +
      xlab("")
    p
    ggsave("plot/QC/qc-reads.pdf",p,width=5.5,height=3.5)
  }
  
  
  #####暂时不用    
  {  
    p=ggplot(qc_reads, aes(variable, value, fill=variable)) +
      #facet_grid( ~ variable)+
      geom_boxplot(notch = FALSE, alpha = 0.9, size=0.4, outlier.shape = NA, position=position_dodge(0.9)) +
      geom_point()+
      #geom_jitter(position = position_jitterdodge(jitter.width = 0, dodge.width = 0), alpha = 0.7, size=0.8)+
      scale_y_log10(limits = c(1*10^4, 5*10^7), labels = trans_format("log10",math_format(10^.x)))+
      #geom_line(aes(group=patient), color="grey60", size=0.3)+
      theme_bw()+
      scale_fill_nejm()+
      theme_classic()+
      theme(panel.background=element_rect(fill="white", colour="black", size=0.5),
            strip.text = element_text(size=12, color="black"),
            strip.background = element_blank(),
            axis.line=element_blank(),
            axis.title=element_text(size=15, color="black"),
            axis.text.x = element_text(size=10, color="black",angle=45,hjust=1,vjust=1),
            axis.text.y = element_text(size=12, color="black"),
            legend.position="none",
            legend.title= element_text(color="black", size=12),
            legend.text= element_text(color="black", size=12))+
      ylab("number of reads")+
      xlab("")
    p
  }
  ggsave("plot/QC/qc-reads_barplot.pdf",width=5.5,height=3.5,p)
  
  
  #####qc-ratio
  {
    q=ggplot(qc_ratio,  aes(variable, value, fill=variable)) +
      #facet_grid(~ variable)+
      geom_boxplot(notch = FALSE, alpha = 0.9, size=0.4, outlier.shape = NA, position=position_dodge(0.9)) +
      geom_jitter(position = position_jitterdodge(jitter.width = 0, dodge.width = 0), alpha = 0.7, size=0.8)+
      #geom_line(aes(group=patient), color="grey60", size=0.3)+
      scale_y_continuous(limits = c(0,1))+
      scale_fill_nejm()+
      theme_classic()+
      theme(panel.background=element_rect(fill="white", colour="black", size=0.5),
            strip.text = element_text(size=12, color="black"),
            strip.background = element_blank(),
            axis.line=element_blank(),
            axis.title=element_text(size=15, color="black"),
            axis.text.x = element_text(size=10, color="black",angle=45,hjust=1,vjust=1),
            axis.text.y = element_text(size=12, color="black"),
            legend.position="none",
            legend.title= element_text(color="black", size=12),
            legend.text= element_text(color="black", size=12))+
      ylab("ratio")+
      xlab("")
    q
    
  }
  ggsave("plot/QC/qc-ratio.pdf",width=5.5,height=3.5,q)
  
  
}


#===================================
# 2--sample filter
# 2.1 clean reads > 1M & hg38_dedup reads > 0.3M
#===================================
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
    geom_jitter(shape=16, position = position_jitter(0.25),show.legend = FALSE) 
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
    geom_jitter(shape=16, position = position_jitter(0.25),) 
  p2
  
  # 将两张图拼在一起，共用相同的y轴
  combined_plot <- p1 + p2 + plot_layout(nrow = 1, widths = c(0.5, 0.5))
  print(combined_plot)
  #ggsave("/Users/wangge/Documents/DM/plot/QC/clean&usable.pdf", combined_plot,width =7 ,height = 4)
}


#===================================
# 2.2 filter low quality sample
#####
#qc <- qc[-c(154:157),]
rows_to_keep <- !(rownames(qc) %in% c("HC-50","HC-51","HC-52","HC-53"))
qc <- qc[rows_to_keep,]
rows_to_keep <- !(anno$sample_id %in% c("HC-50","HC-51","HC-52","HC-53"))
anno <- anno[rows_to_keep,]


######
ratio <- function(x){
  as.numeric(x)
  (x/sum(x))*100
}

#=======================
#RNA_category piechat
#####
#还没画

########### 2.clean/genome/rna/stacked plot
#
{
  clean_reads <- qc %>%
    transmute(sample_id = rownames(qc),
              clean = clean,
              genome = star_hg38_v38 + star_circRNA,
              rRNA = star_rRNA,
              microbiome = unmapped - unclassified,
              unmapped = clean-genome-rRNA-microbiome,
              group=group)
  clean_ratio <- clean_reads %>%
    transmute(sample_id = sample_id,
              genome = genome / clean,
              rRNA = rRNA / clean,
              microbiome = microbiome / clean,
              unmapped = unmapped / clean,
              group=group)
  
  clean_ratio$genome <- as.numeric(clean_ratio$genome)
  clean_ratio <- clean_ratio %>%
    arrange(desc(genome))
  clean_ratio <- melt(clean_ratio,id=c("sample_id","group"))
  clean_ratio$sample_id <- factor(clean_ratio$sample_id, levels=anno$sample_id)
  #clean_ratio$sample_id <- factor(clean_ratio$sample_id,levels = unique(clean_ratio$sample_id))
  {
    col1 = pal_jco()(8)
    col2 = pal_npg()(8)
    col4=pal_igv()(8)
    a=ggplot(clean_ratio, aes(x=sample_id, y=value, fill=variable))+
      geom_bar(position = "fill", stat = "identity", width=0.8)+
      #scale_fill_manual(values=c(col1[5],col1[6],col2[7],col2[5],col1[4],col1[3]))+
      scale_fill_manual(values=c(col1[5],col1[4],col2[5],col1[3]))+
      scale_y_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1))+
      theme_bw()+
      guides(fill=guide_legend(title=NULL, ncol=1))+
      theme(
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        plot.title = element_text(color="black", size=12, hjust=0.5),
        #legend.position="top",
        #legend.direction = "horizontal",
        #legend.margin = margin(c(50,0,50,0)),
        legend.text= element_text(color="black", size=12),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length.y = unit(0.2,"cm"),
        axis.ticks.length.x = unit(0.1,"cm"),
        axis.text.x = element_text(color="black", size=7,angle=90,hjust=1,vjust=0.5),
        axis.text.y = element_text(color="black", size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
      labs(x="",y="",title="", face="bold")#+
    a
  }
  
  genome_reads <- qc %>%
    transmute(sample_id = rownames(qc),
              `genome_non-dup` = star_hg38_v38_dedup,  #+ circRNA_dedup,
              exonic = (MT_mRNA+MT_tRNA+chrM.all+mRNA+lncRNA+snoRNA+snRNA+srpRNA+tRNA+tucpRNA+Y_RNA+misc_RNA+pseudogene+exon), #+circRNA_dedup
              intronic = intron,
              intergenic = intergenic)
  genome_ratio <- genome_reads %>%
    transmute(sample_id = sample_id,
              exonic = exonic / `genome_non-dup`,
              intronic = intronic / `genome_non-dup`,
              intergenic = intergenic /`genome_non-dup`)
  genome_ratio <- melt(genome_ratio,id=c("sample_id"))
  genome_ratio$sample_id <- factor(genome_ratio$sample_id, levels=anno$sample_id)
  
  {
    col1 = pal_jco()(8)
    b=ggplot(genome_ratio, aes(x=sample_id, y=value, fill=variable))+
      geom_bar(position = "fill", stat = "identity",width=0.8)+
      theme_bw()+
      scale_fill_manual(values=c(col1[1],"#EEB4B4","#8B6969"))+
      scale_y_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1))+
      guides(fill=guide_legend(title=NULL, ncol=1))+      
      theme(
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        plot.title = element_text(color="black", size=12, hjust=0.5),
        legend.margin = margin(c(30,12,30,10)),
        #legend.position="top",
        #legend.box = "horizontal",
        legend.text= element_text(color="black", size=12),
        panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.ticks.length.y = unit(0.2,"cm"),
        axis.ticks.length.x = unit(0.1,"cm"),
        axis.text.x = element_text(color="black", size=7,angle=90,hjust=1,vjust=0.5),
        axis.text.y = element_text(color="black", size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
      labs(x="",y="",title="", face="bold")
    
    b
  }
  
  p=plot_grid(a,b,ncol=1,nrow=2,axis="l",rel_heights = c(1,1),align="v")
  p
  ggsave("plot/qc_clean_genome_ratio.pdf",p,width=12,height=7)
  
}
{
  rna_reads <- qc %>%
    transmute(sample_id=rownames(qc),
              mRNA=mRNA,
              lncRNA=lncRNA,
              #circRNA=circRNA_dedup,
              pseudogene=pseudogene,
              misc=misc_RNA,
              sncRNA=(snoRNA + snRNA + srpRNA + tRNA + Y_RNA),
              tucpRNA=tucpRNA,
              exon=exon,  
              total=mRNA+lncRNA+pseudogene+misc+sncRNA+tucpRNA+exon)
  #total=total)
  rna_ratio <- rna_reads %>%
    transmute(sample_id=sample_id,
              mRNA=mRNA / total,
              lncRNA=lncRNA / total,
              #circRNA=circRNA / total,
              pseudogene=pseudogene / total,
              misc=misc / total,
              sncRNA=sncRNA / total,
              tucpRNA=tucpRNA /total,
              other_exon=exon / total)
  
  rna_ratio <- melt(rna_ratio,id=c("sample_id"))
  rna_ratio$sample_id <- factor(rna_ratio$sample_id, levels=anno$sample_id)
  
  {
    c=ggplot(rna_ratio, aes(x=sample_id, y=value, fill=variable))+
      geom_bar(position = "fill", stat = "identity",width=0.8)+
      theme_bw()+
      scale_fill_igv()+
      scale_y_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1))+
      guides(fill=guide_legend(title=NULL, ncol=1))+      
      theme(
        plot.title = element_text(color="black", size=12, hjust=0.5),
        strip.background = element_rect(size=1, fill="transparent"),
        strip.text = element_text(color="black", size=12),
        #legend.position="top",
        legend.text= element_text(color="black", size=12),
        legend.key.size = unit(0.2,"inches"),
        panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.ticks.length.y = unit(0.2,"cm"),
        axis.ticks.length.x = unit(0,"cm"),
        axis.text.x = element_text(color="black", size=7,angle=90,hjust=1,vjust=0.5),
        axis.text.y = element_text(color="black", size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
      labs(x="",y="",title="", face="bold")
    c
    ggsave("plot/qc_exon_ratio.pdf",c,width=12,height=4)
  }
}
#=======================
#rRNA_ratio
#####
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
  geom_jitter(shape=16, position = position_jitter(0.25),) 
ggsave("/Users/wangge/Documents/DM/QC/rRNA_ratio_barplot.pdf",width =4 ,height = 4)


ggplot(qc_filtered,aes(x=1,y=rRNA_ratio)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title="rRNA ratio",y="Reads ratio(%)") +
  xlab("")+
  theme_bw()+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
        panel.grid = element_blank(),
        #panel.background = element_rect(colour = NA),
        #plot.background = element_rect(colour = NA),
        panel.border = element_rect(colour = NA),
        axis.title.y = element_text(angle = 90, vjust = 2, size = 16),
        #axis.title.x = element_text(vjust = -0.2, size = base_size),
        axis.text.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line())+
  scale_fill_manual(values = brewer.pal(5,'BrBG')[c(1,2,5)])+
  scale_y_continuous(
    limits = c(0,0.6),
    breaks = c(0.1,0.2,0.3,0.4,0.5),
    labels = c('10%','20%','30%','40%','50%'))+
  geom_jitter(shape=16, position = position_jitter(0.25)) +
  geom_hline(aes(yintercept=0.4, color="red"), linetype="dashed",show.legend = FALSE)
  
#=======================
#错误的exon-ratio barplot
#####
qc_exon <- qc_filtered %>% 
  transmute(#group = group,
            #raw_reads = raw,
            #clean_reads = clean,
            #genome_align_reads = hg38 + circRNA,   # circRNA
            #`De-duplicated RNA reads`= hg38_dedup,
            #libsize = libsize[rownames(qc),"lib.size"], ############# 换成matrix中的reads总数/libsize
            # longRNA = mRNA + lncRNA,
            #`intron-spanning` = hg38_intron_spanning,
            exonic = (mRNA+lncRNA+snoRNA+snRNA+srpRNA+tRNA+tucpRNA+Y_RNA+misc_RNA+pseudogene+exon),
            intronic = intron,
            intergenic = intergenic,
            #spikeIn = spikein_long,
            #rRNA = rRNA,
            #mtRNA = chrM,
            #microbiome = unmapped - unclassified,
            #unmapped_reads = clean_reads-genome_align_reads-mtRNA-rRNA-spikeIn-microbiome,
            # clean_ratio = trimGC / raw,
            # genome_ratio = (hg38+circRNA) / trimGC,   # circRNA    
            # `non-dup_ratio` = hg38_dedup / trimGC,  ### changed     
            # `longRNA_ratio` = (mRNA + lncRNA) / hg38_dedup,
            # `intron-spanning_ratio` = hg38_intron_spanning / hg38_dedup,
            # exonic_ratio = (mRNA+lncRNA+snoRNA+snRNA+srpRNA+tRNA+tucpRNA+Y_RNA+misc_RNA+pseudogene+exon) / hg38_dedup,
            # intronic_ratio = intron / hg38_dedup,
            # intergenic_ratio = intergenic / hg38_dedup,
            # rRNA_ratio = rRNA / trimGC,
            # mtRNA_ratio = chrM / trimGC,
            # cancer = cancer,
            # class = class
  )%>%
  t() %>% 
  as.data.frame() %>% 
  mutate(across(contains("-"), ratio)) %>%
  t() %>% 
  as.data.frame()
qc_exon$group <- sub("-.*","",rownames(qc_exon))
#qc_exon <- qc_exon %>%
#  mutate(group = factor(group, levels = c("MDA5", "ARS", "HC")))%>%
#  arrange(group)
exon_level <- qc_exon %>%
  mutate(sample_id=rownames(.)) %>%
  group_by(group) %>%
  arrange(group, desc(exonic)) %>%
  ungroup() 
rownames(exon_level) <- exon_level$sample_id
exon_level <- exon_level %>%
  t()%>% 
  as.data.frame() 
exon_level <- exon_level[-c(4:5),]%>%
  mutate(parameters=rownames(.)) %>% 
  pivot_longer(cols=colnames(.)[colnames(.) != "parameters"],names_to="Sample",values_to = "reads_ratio") %>% 
  as.data.frame()
######
#合适的exon-ratio barplot
#####
qc_reads_dist <- qc_filtered %>%
  transmute(#sample_id = sample_id,
            exonic = (mRNA+lncRNA+snoRNA+snRNA+srpRNA+tRNA+tucpRNA+Y_RNA+misc_RNA+pseudogene+exon),
            intronic = intron,
            intergenic = intergenic
            
  )%>% 
  #column_to_rownames(.,"sample_id") %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(parameters=rownames(.)) %>% 
  mutate(across(contains("-"),ratio)) %>%
  pivot_longer(cols=colnames(.)[colnames(.) != "parameters"],names_to="Sample",values_to = "reads_ratio") %>% 
  as.data.frame()
s_level <- qc_reads_dist %>% filter(parameters=="exonic") %>% 
  arrange(-reads_ratio) 
qc_reads_dist$Sample <- factor(qc_reads_dist$Sample,levels = s_level$Sample)
qc_reads_dist$group <- sub("-.*","",qc_reads_dist$Sample)


ggplot(qc_reads_dist, aes(x = Sample, y = reads_ratio, fill = parameters)) + 
  geom_bar(stat = "identity", aes(fill=parameters)) + 
  labs(x="", y="reads ratio (%)") +
  theme_classic() + 
  # scale_fill_brewer(palette = "Paired")+
  scale_fill_discrete_sequential(palette = "BrwnYl", 
                                 nmax = 6, 
                                 rev = T, 
                                 order = 3:5) +
  theme(
    legend.position="top",
    legend.title = element_text(face="bold", color="black", size=12),
    legend.text= element_text(face="bold", color="black", size=12),
    plot.title = element_text(hjust = 0.5,size=15,face="bold"),
    axis.text.x = element_text(face="bold", color="black", size=4,angle = 90,hjust = 1),
    axis.text.y = element_text(face="bold",  color="black", size=12),
    axis.title.y = element_text(face="bold",color="black", size=12))
ggsave("/Users/wangge/Documents/DM/QC/exon_ratio.pdf",width = 8,height = 4)






# clean reads > 1M
rawreads <- QC_reads_number %>% 
  filter(parameters %in% "raw_reads") %>% 
  filter(values < 4000000) %>% 
  mutate(status = "raw_reads-fail") #20221124-SLE-1-BC7 fail
QC_reads_number %>% 
  filter(parameters %in% "raw_reads") %>% 
  ggplot(aes(x=parameters, y=values)) + 
  # geom_violin(fill="gray")+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=2,alpha=0.7, color="transparent")+
  # geom_jitter(position=position_jitter(0.1), cex=2)+
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="crossbar", width=0.5)+
  scale_y_log10(limits = c(10^6, 10^8), labels = trans_format("log10",math_format(10^.x)))+
  # labs(title="Sample QC",x="", y = "Reads number")+
  labs(title="",x="raw reads", y = "Reads number")+
  geom_hline(aes(yintercept=4000000, color="red"), linetype="dashed",)+ # raw reads 4M
  # geom_hline(aes(yintercept=3800000, color="red"), linetype="dashed")+ # clean reads 3.8M 
  # geom_hline(aes(yintercept=2000000, color="Blue"), linetype="dashed")+ # genome reads 2M
  theme_classic()+
  # facet_wrap(~parameters,scales = "free_y",ncol = 3)
  theme(legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=12),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=12),
        plot.title = element_text(hjust = 0.5,size=15,face="bold"),
        # axis.text.x = element_text(face="bold", color="black", size=12,angle = 45,hjust = 1),
        axis.text.x = element_text(face="bold", color="black", size=0),
        axis.text.y = element_text(face="bold",  color="black", size=12),
        axis.title.x = element_text(face="bold", color="black", size=12),
        axis.title.y = element_text(face="bold",color="black", size=12))
ggsave("./cfRNA/results/QC/raw_reads.pdf",width = 2,height=3)




qc.1 <- qc %>% transmute(sample_id=sample_id,
                         #spikein_small=spikein_small,
                         rRNA=rRNA,
                         mtRNA=mtRNA,
                         sncRNA=(miRNA+piRNA+snoRNA+snRNA+tRNA+Y_RNA),
                         lncRNA=lncRNA,
                         mRNA=mRNA,
                         srpRNA=srpRNA,
                         tucpRNA=tucpRNA,
                         pseudogene=pseudogene,
                         genome=genome,
                         #unmapped=trimGC-rRNA-mtRNA-lncRNA-mRNA-sncRNA-tucpRNA-pseudogene-genome,
                         unmapped_bt2=unmapped_bt2)

qc.1 <- qc %>% transmute(sample_id=sample_id,
                         miRNA=miRNA,
                         piRNA=piRNA,
                         snoRNA=snoRNA,
                         snRNA=snRNA,
                         tRNA=tRNA,
                         Y_RNA=Y_RNA)

qc.1 <- qc %>% transmute(sample_id=sample_id,
                         #rRNA=rRNA_dedup,  #### rRNA
                         lncRNA=lncRNA_dedup,
                         mRNA=mRNA_dedup,
                         sncRNA=miRNA_dedup+piRNA_dedup+snoRNA_dedup+snRNA_dedup+tRNA_dedup+Y_RNA_dedup,
                         srpRNA=srpRNA_dedup,
                         tucpRNA=mit_tucpRNA_dedup,
                         pseudogene=pseudogene_dedup)
#genome=genome_dedup)

qc.1 <- qc %>% transmute(sample_id=sample_id,
                         miRNA=miRNA_dedup,
                         piRNA=piRNA_dedup,
                         snoRNA=snoRNA_dedup,
                         snRNA=snRNA_dedup,
                         tRNA=tRNA_dedup,
                         Y_RNA=Y_RNA_dedup)

rt<-function(x){x/rowSums(qc.1[,-1])}
qc.ratio <- apply(qc.1[,-1],2,rt)




qc.1 <- melt(qc.1)
qc.1$sample_id <- factor(qc.1$sample_id, levels=anno$sample_id)








######################################
## non-coding genes matrix
{
  # mRNA, lncRNA, pseudogene, mtRNA-tRNA 
  # 不用bowtie2的结果，从gencode中提取
  {
    gencode <- read.table("1204/gencode.txt",sep="\t",header=T,check.names=F,row.names = 1)
    gencode <- gencode[,-c(40:44, 50:54)]
    colnames(gencode)[35:44] <- paste0(substr(colnames(gencode)[35:44],1,12),substr(colnames(gencode)[35:44],16,20))
    #colnames(gencode)
    colnames(gencode) <- anno[match(colnames(gencode),anno$library_id),]$sample_id
    #anno <- anno[-c(154:157),]#
    gencode <- gencode[,anno$sample_id]
    group = anno$group

    #### gencode
    table(str_split(rownames(gencode),"\\|",simplify=T)[,4])
    
    gn.type <- c("protein_coding","lncRNA",
                 "IG_C_pseudogene","IG_J_pseudogene","IG_pseudogene","IG_V_pseudogene",
                 "polymorphic_pseudogene","processed_pseudogene","pseudogene",
                 "TR_J_pseudogene","TR_V_pseudogene",
                 "transcribed_processed_pseudogene","transcribed_unitary_pseudogene","transcribed_unprocessed_pseudogene","translated_processed_pseudogene","translated_unprocessed_pseudogene","unitary_pseudogene","unprocessed_pseudogene",
                 "Mt_tRNA",
                 "scaRNA",
                 "sRNA",
                 "IG_C_gene","IG_D_gene","IG_V_gene","IG_J_gene",
                 "TR_C_gene","TR_D_gene","TR_V_gene","TR_J_gene",
                 "miRNA",
                 "snRNA","snoRNA")
    gn.nc = gencode[str_split(rownames(gencode),"\\|",simplify=T)[,4] %in% gn.type,]       
    
    gn.srp = gencode[substr(str_split(rownames(gencode),"\\|",simplify=T)[,3],1,5)=="RN7SL",]
    gn.y = gencode[str_split(rownames(gencode),"\\|",simplify=T)[,3]=="Y_RNA",]
    
    gn.nc <- rbind(gn.nc, gn.srp, gn.y)
    # 38306
    # 58261
    gn.code <- gencode[str_split(rownames(gencode),"\\|",simplify=T)[,4]=="protein_coding",]
    write.table(gn.code,"1204/protein-coding.matrix.txt",sep="\t",quote=F)
  }
  
  mt.pi <- read.table("1204/count_bt2_matrix/piRNA.matrix",sep="\t",header=T,check.names=F,row.names = 1)
  mt.t <- read.table("1204/count_bt2_matrix/tRNA.matrix",sep="\t",header=T,check.names=F,row.names = 1)
  mt.tucp <- read.table("1204/count_bt2_matrix/mit_tucpRNA.matrix",sep="\t",header=T,check.names=F,row.names = 1)
  mt1 <- rbind(mt.pi, mt.t, mt.tucp)
  
  mt1 <- mt1[,-c(40:44, 50:54,165:168)]
  colnames(mt1)[1:20] <- paste0("230821-",colnames(mt1)[1:20])
  colnames(mt1)[35:44] <- paste0(substr(colnames(mt1)[35:44],1,12),substr(colnames(mt1)[35:44],16,20))
  #colnames(gencode)
  colnames(mt1) <- anno[match(colnames(mt1),anno$library_id),]$sample_id
  #anno <- anno[-c(154:157),]
  mt1 <- mt1[,anno$sample_id]
  #检查是否一致
  colnames(gencode)==colnames(mt1)

  mt <- rbind(gn.nc, mt1)
  dim(mt)
  mt[1:4,1:2]
  # 137047
  colSums(mt)
  write.table(mt,"1204/totalRNA.matrix.txt",sep="\t",quote=F)
  
  
  df1 <- data.frame(ncrnaid=c(rownames(gn.nc),rownames(mt.pi),rownames(mt.t),rownames(mt.tucp)))
  df1$type <- c(str_split(rownames(gn.nc),"\\|",simplify=T)[,4], rep("piRNA",nrow(mt.pi)), rep("tRNA",nrow(mt.t)), rep("tucpRNA",nrow(mt.tucp)))
  
  df1$type2 <- "pseudogene"
  
  df1$type2[df1$type %in% c("IG_C_gene","IG_J_gene","IG_D_gene","IG_V_gene","TR_C_gene","TR_J_gene","TR_D_gene","TR_V_gene")] <- "IGTR"
  df1$type2[df1$type == "protein_coding"] = "protein_coding"
  df1$type2[df1$type == "lncRNA"] = "lncRNA"
  df1$type2[df1$type == "miRNA"] = "miRNA"
  df1$type2[df1$type == "piRNA"] = "piRNA"
  df1$type2[df1$type == "tRNA"] = "tRNA"
  df1$type2[df1$type == "Mt_RNA"] = "Mt_RNA"
  df1$type2[df1$type == "snRNA"] = "snRNA"
  df1$type2[df1$type == "snoRNA"] = "snoRNA"
  df1$type2[df1$type == "tucpRNA"] = "tucpRNA"
  df1$type2[df1$type == "sRNA"] = "sRNA"
  df1$type2[df1$type == "scaRNA"] = "scaRNA"
  df1$type2[str_split(df1$ncrnaid,"\\|",simplify=T)[,3] == "Y_RNA"] = "Y_RNA"
  df1$type2[substr(str_split(df1$ncrnaid,"\\|",simplify=T)[,3],1,5) == "RN7SL"] = "srpRNA"
  
  table(df1$type)
  table(df1$type2)
  write.table(df1,"1204/totalRNA.matrix.type.txt",sep="\t",quote=F,row.names = F)
  
}
###figure 2B

libsize <- colSums(mt)
rtype <- data.frame(libsize=libsize)
df_list <- split(df1, df1$type2)
list2env(setNames(df_list, paste0("df_", names(df_list))), envir = .GlobalEnv)
df_vars <- ls(pattern = "^df_")
rm(list = df_vars)

identical(rownames(mt), df1$ncrnaid)
mt$type2 <- df1$type2
mt_sums <- aggregate(. ~ type2, data = mt, FUN = sum)
######################################
## Detected gene number 
{
  mt <- read.table("1204/totalRNA.matrix.txt",sep="\t",header=T,check.names = F)
  #anno <- anno[-c(154:157),]#
  group <- anno$group

  mt <- mt[rowSums(mt)!=0,] 
  mt[1:4,1:4]
  # 106227

  # libsize 这是在干嘛？？？？？
  table(rownames(qc)==colnames(mt))
  libsize = qc$star_hg38_v38_dedup
  colSums(mt)
  libsize
  
  # rnatypes
  rna.type <- read.table("1204/totalRNA.matrix.type.txt",sep="\t",header=T)
  colnames(rna.type)[1] = "id"
  rownames(rna.type) <- rna.type$id
  
  library(edgeR)
  y <- DGEList(counts=mt, group=anno$group)
  y
  ########## filter noncoding genes?
  #keep <- filterByExpr(y, group=cell, min.count=2, min.prop=0.2) # 34877
  #keep <- filterByExpr(y, group=cell, min.count=5, min.prop=0.2) # 5194
  #table(keep)
  #y <- y[keep, ,keep.lib.size=TRUE]
  
  cpm <- as.data.frame(edgeR::cpm(y))
  test.type <- rna.type[rownames(cpm),]$type2
#####
  # PCA for all data
  {
    #y <- DGEList(counts=data, group=group)
    keep <- filterByExpr(y, group=group, min.count=5, min.prop=0.2)
    table(keep)
    y <- y[keep, ,keep.lib.size=TRUE]
    
    
    library("preprocessCore")
    library("Rtsne")
    X <- t(edgeR::cpm(y))
    X_norm.4 <- scale(X, center = T, scale = T)
    
    # PCA
    {
      pca <- prcomp(X_norm.4, center = F, scale = F, rank. = 2)  
      pca.var <- pca$sdev^2  ## sdev是标准偏差，十个样本，就有十个标准偏差，平方是避免负数的干扰
      pca.var.per <- round(pca.var/sum(pca.var)*100, 1)  ##求每个样本的variation
      pca.var.per
      barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")  ##用柱状图可视化
      summary(pca)
      
      {
        pca.data <- data.frame(pca$x, group=group)
        percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
        percentage <- paste(colnames(pca.data)[1:9],"(", paste(as.character(percentage)), "%", ")", sep="")
        p=ggplot(pca.data, aes(x=PC1,y=PC2,color=group)) + 
          geom_point(size=2,alpha=1) +
          #stat_ellipse(aes(fill=group), type = "t", geom = "polygon", level = 0.95, size=0.6, alpha=0.5)+
          scale_color_brewer(palette = "Set1") +
          theme_bw()+
          guides(shape=guide_legend(NULL),
                 fill=guide_legend(NULL),
                 color=guide_legend(title=NULL, ncol=1))+
          theme(
            plot.margin = unit(c(1, 1, 1, 1),"cm"),
            panel.border=element_rect(size=1),
            panel.grid=element_blank(),
            legend.background = element_blank(),
            legend.key.size = unit(0,"cm"),
            legend.position = "top",
            legend.title = element_blank(),
            legend.spacing.y = unit(0,"cm"),
            legend.box = "horizontal",
            legend.box.margin = unit(c(0,0,0,0),"cm"),
            legend.text= element_text(color="black", size=10),
            axis.text.x = element_text(color="black", size=15),
            axis.text.y = element_text( color="black", size=15),
            axis.title.x = element_text(color="black", size=15),
            axis.title.y = element_text(color="black", size=15))+
          xlab(percentage[1]) + ylab(percentage[2])
        p
        ggsave("plot/gencode-PCA.pdf",p,width=7,height=7)
      }
      
      {
        #BiocManager::install("ggbiplot")
        library(ggbiplot)
        ggbiplot(pca, choices = 2,
                 obs.scale = 1, var.scale = 1,
                 ellipse = F, groups = pca.data[, 3], # 注意这里即可，刚刚是定义了数据集第一列的分组信息，这里保持一致即可，其他参数可以直接使用
                 ellipse.prob = 0.95, circle = F,
                 varname.size = 5, varname.adjust = 1.5,
                 var.axes = T) 
      }
    }
    # tSNE
    {
      set.seed(100)  ## 每次结果不一样 
      tsne <- Rtsne(X_norm.4, dims=2, check_duplicates = FALSE, perplexity = 3) 
      # dims参数设置降维之后的维度，默认值为2
      # perplexity参数的取值必须小于(nrow(data) - 1 )/ 3
      # theta参数取值越大，结果的准确度越低，默认值为0.5
      # max_iter参数设置最大迭代次数。
      #plot(tsne$Y, col = colors_plot, asp = 3, pch = 20, xlab = "tSNE_1", ylab = "tSNE_2", main = "")
      #legend("topright", inset=0.02, legend=unique(Y), col = a)
      
      tsne.data <- as.data.frame(tsne$Y)
      col3=pal_nejm()(8)
      show_col(col3)
      {
        p=ggplot(tsne.data, aes(x=V1,y=V2,color=group)) + 
          geom_point(size=2,alpha=1) +
          stat_ellipse(level = 0.95, size=0.7, alpha=0.7)+
          scale_color_manual(values=c(col3[3],col3[2],col3[1]))+
          #scale_x_continuous(limits=c(-20,20), breaks=c(-15,0,15))+
          #scale_y_continuous(limits=c(-20,20), breaks=c(-20,0,20))+
          theme_bw()+
          theme(
            plot.margin = unit(c(1, 1, 1, 1),"cm"),
            panel.border=element_rect(size=1),
            panel.grid=element_blank(),
            legend.background = element_rect(fill="transparent"),
            legend.key.size = unit(0.1,"cm"),
            legend.position = c(0.8,0.9),
            legend.title = element_blank(),
            legend.text= element_text(color="black", size=12),
            axis.text.x = element_text(color="black", size=15),
            axis.text.y = element_text( color="black", size=15),
            axis.title.x = element_text(color="black", size=15),
            axis.title.y = element_text(color="black", size=15))+
          labs(x="tSNE_1",y="tSNE_2")
        p
      }
      ggsave("plot/tsne2.pdf",p,width = 6.5, height = 6)
      
    }
  }  
###### 
  #coding+ncRNA type number and count
#####
  #查看检测到的基因数目 每个样本中的miRNA
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
    geom_violin(outlier.shape = NA,show.legend = TRUE) +
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
  ggsave("plot/miRNA-species.pdf",p3,width = 4.3, height = 3.7)
  ######
  {
    grouptest="ARS"
    dfall = data.frame()
    for (grouptest in unique(group)){
      
      test <- cpm[, group==grouptest]
      
      #### gene number
      n = ncol(test) * 0.2
      gn <- data.frame(detected = rownames(test)[rowSums(test > 1) > n]) 
      gn$type <- rna.type[gn$detected,]$type2
      gntype <- as.data.frame(table(gn$type))
      
      ### gene count
      test2 <- aggregate(test, by=list(test.type), sum)
      rownames(test2) <- test2$Group.1
      test2 <- test2[,-1]
      gncount <- as.data.frame(rowMeans(test2))
      
      
      ### merge
      df <- data.frame(type = gntype$Var1,
                       number = gntype$Freq,
                       count = gncount[gntype$Var1,1])
      df$celltype = grouptest
      
      dfall <- rbind(dfall, df)
    }
    write.table(dfall, "1204/genenumber.txt", sep="\t", quote=F)
    
  }
#######  
  dfall <- read.table("1204/genenumber.txt",sep="\t",header=T)
  # 
  dfall. <- dfall
  dfall.$type2 <- dfall.$type
  dfall. <- dfall.[dfall.$type!="IGTR",]
  dfall.$type2[dfall.$type2 %in% c("miRNA","piRNA","snoRNA","snRNA","scaRNA","sRNA","tRNA","Y_RNA")] = "sncRNA"
  dfall.$type2 <- factor(dfall.$type2, levels=c("protein_coding","lncRNA","sncRNA","pseudogene","srpRNA","tucpRNA"))
  dfall.$celltype <- factor(dfall.$celltype, levels = c("MDA5", "ARS", "HC"))
  write.table(dfall., "1204/genecount.txt", sep="\t", quote=F)
  
  # barplot
  selected_colors <- brewer.pal(11, "Spectral")[c(3,4,6,8,10,11)]
  {
    a=ggplot(dfall., aes(x=celltype, y=number, fill=type2))+
      #coord_polar(theta = 'y')+
      #coord_flip()+
      geom_bar(position = "fill", stat = "identity",width=0.8)+
      scale_fill_manual(values = selected_colors) +
      #scale_y_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1))+
      theme_bw()+
      guides(fill=guide_legend(title=NULL, ncol=1))+
      theme(
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position="right",
        #legend.margin = margin(c(50,0,50,0)),
        legend.text= element_text(color="black", size=10),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length.y = unit(0,"cm"),
        axis.ticks.length.x = unit(0.1,"cm"),
        axis.text.x = element_text(color="black", size=8),
        axis.text.y = element_text(color="black", size=8, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
      labs(x="",y="",title="", face="bold")
    a
  }
  dfall2 <- filter(dfall., dfall.$type2=="sncRNA")
  dfall2 <- dfall2[dfall2$type != "piRNA",]
  dfall2$type <- factor(dfall2$type, levels=c("miRNA","snoRNA","snRNA","scaRNA","sRNA","tRNA","Y_RNA"))#"piRNA",
  write.table(dfall2, "1204/genecount2.txt", sep="\t", quote=F)
  # barplot
  {
    b=ggplot(dfall2, aes(x=celltype, y=count, fill=type))+
      #coord_polar(theta = 'y')+
      #coord_flip()+
      geom_bar(position = "fill", stat = "identity",width=0.8)+
      scale_fill_nejm()+
      #scale_y_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1))+
      theme_bw()+
      guides(fill=guide_legend(title=NULL, ncol=1))+
      theme(
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position="right",
        #legend.margin = margin(c(50,0,50,0)),
        legend.text= element_text(color="black", size=10),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length.y = unit(0,"cm"),
        axis.ticks.length.x = unit(0.1,"cm"),
        axis.text.x = element_text(color="black", size=8),
        axis.text.y = element_text(color="black", size=8, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
      labs(x="",y="",title="", face="bold")
    b
  }
  p=ggarrange(a,b,ncol=2,align="v",widths=c(1,1))
  p
  ggsave("plot/barplot-rnatype-count.pdf",p,width=9,height=4.6)
#######
}

######################################
## edgeR diffexpressed
{  
  group2 <- group
  group2[group2 == "MDA5" | group2 == "ARS"] <- "DM"
  mt1 <- mt[!grepl("^piR-hsa|protein_coding|pseudogene|lncRNA", rownames(mt)), ]
  mt_pc <- mt[grep("protein_coding", rownames(mt)), ]
  lncRNA_rows <- grepl("\\|lncRNA", rownames(mt))
  mt_lncRNA <- mt[lncRNA_rows, ]
  #miRNA_rows <- grepl("\\|miRNA", rownames(mt))
  #mt_miRNA <- mt[miRNA_rows, ]
  group4 <- group[c(1:44,56:105)]
  
  y <- DGEList(counts=mt1, group=group2)
  
  group3 <- group[group != "HC"]  #######
  y <- DGEList(counts=mt_lncRNA[,c(1:44,56:105)], group=group4) ######
  
  y <- DGEList(counts=mt_lncRNA[,c(1:105)], group=group3) ######
  y
  ########## filter noncoding genes?
  #keep <- filterByExpr(y, group=group, min.count=2, min.prop=0.2) 
  keep <- filterByExpr(y, group=group3, min.count=2, min.prop=0.2) # min.count=5, min.prop=0.5
  table(keep)
  y <- y[keep, ,keep.lib.size=TRUE]
  
  y <- calcNormFactors(y, method="TMM")
  #design <- model.matrix(~group)
  design <- model.matrix(~group2)
  design <- model.matrix(~group3)
  design <- model.matrix(~group4) #####
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
  con <- makeContrasts(
    ARS_vs_MDA5 = groupARS - groupMDA5,
    HC_vs_MDA5 = groupHC - groupMDA5,
    HC_vs_ARS = groupHC - groupARS,
    levels = colnames(design)
  )
#####
  qlf <- glmQLFTest(fit.ql, coef=2)
  #anov <- glmQLFTest(fit.ql, contrast=con)
  #qlf <- glmQLFTest(fit.ql, coef=2)
  de <- topTags(qlf, n=Inf)$table
  de.pvalue <- filter(de, PValue< 0.01)#0.01
  #de.FDR <- filter(de.pvalue,FDR< 0.05)#0.05
  
  de.up <- filter(de.pvalue, logFC> 1.2)#0
  de.dw <- filter(de.pvalue, logFC< -1)#0
  
}
# heatmap
{
  logRPM <- edgeR::cpm(y, log=T)
  logRPM <- logRPM[c(rownames(de.up),rownames(de.dw)), ]

  dim(logRPM) # 48,42
  logRPM.scale <- scale(t(logRPM), center = T, scale = T)
  logRPM.scale <- t(logRPM.scale)
  rownames(logRPM.scale)[c(26:32,34,35,39:41,43) ] <- paste0("||", rownames(logRPM.scale)[c(26:32,34,35,39:41,43)], "|tRNA")
  rownames(logRPM.scale) <- unlist(lapply(strsplit(rownames(logRPM.scale),"|",fixed=T),function(x) x[3])) #paste(x[3], x[4], sep = "|")
  
  #group. <- group3[order(group3,decreasing=F)]
  #logRPM.scale. <- logRPM.scale[,order(group2,decreasing=F)]
  #group1 <- factor(group,levels=c("MDA5","ARS","HC"))
  #group. <- group[order(group1)]
  #logRPM.scale. <- logRPM.scale[,order(group3)]
  
  # pheatmap
  {
    library(pheatmap)
    ann_col = data.frame(class = as.character(group3))######
    rownames(ann_col) = colnames(logRPM.scale)
    ann_colors = list(class = brewer.pal(5,'BrBG')[c(1,2)])
    names(ann_colors$class) = c("MDA5","ARS")
    
    ann_col = data.frame(class = as.character(group))
    rownames(ann_col) = colnames(logRPM.scale)
    ann_colors = list(class = brewer.pal(5,'BrBG')[c(1,2,5)])
    names(ann_colors$class) = c("MDA5","ARS","HC")
    #GeneClass = c(Path1 = "#807DBA", Path2 = "#9E9AC8", Path3 = "#BCBDDC"))
    
    col = brewer.pal(9,"Set1")
    #pdf("plot/pheatmap-pvalue0.001.pdf", width = 20, height = 16)
    par(mar=c(2,2,2,2))
    heatmap = pheatmap(logRPM.scale, 
             #color = colorRampPalette(c(col[2],"white",col[1]))(1000),
             color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
             breaks = seq(-3, 3, length.out = 51),
             cutree_col = 0, 
             #cutree_row = 3, #break up the heatmap by clusters you define
             cluster_rows = T,
             cluster_cols = F, #by default, pheatmap clusters by both row and col
             show_rownames = T,
             show_colnames = F,
             fontsize_col = 12,
             angle_col = 90,
             annotation_col = ann_col,
             #annotation_row = ann_row,
             annotation_colors = ann_colors,
             annotation_names_row = F,
             #legend_breaks=-4:3,
             border=F)
    heatmap
    #ggsave("plot/DMvsHC_lncRNA-heatmap.pdf",heatmap,width=6,height=5)#dev.off()
    ggsave("plot/MDA5vsARS_lncRNA-heatmap.pdf",heatmap,width=7,height=5)#dev.off()
  }
  
  
  # ggplot
  {
    dd <- as.data.frame(logRPM.scale) %>% mutate(feature=rownames(logRPM.scale))
    dd <- melt(dd,id.vars = "feature")
    dd$value[dd$value>3]=3
    dd$value[dd$value<(-3)]=-3
    dd$variable <- factor(dd$variable)#, levels=c(colnames(data)[order(class,decreasing = F)]))
    b=ggplot(dd, aes(variable, feature)) +
      geom_tile(aes(fill = value),colour = "white") + 
      scale_fill_gradient2(name="Z-score",limits=c(-3,3), breaks=c(-3,0,3), low = "dodgerblue4", mid="white",high = "#A52A2A")+
      theme_grey(base_size = 10) + 
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(limits=rev(rownames(logRPM.scale)),expand = c(0, 0)) + 
      theme(legend.position = "right",
            panel.background = element_blank(),
            #panel.border = element_rect(),
            axis.ticks = element_blank(), 
            axis.title = element_blank(),
            #axis.text.x = element_text(size = 8, angle=90, hjust = 1),
            axis.text.x = element_text(size=5,angle=90,hjust=1,vjust=0.5),
            axis.text.y = element_text(size = 8),
            axis.text = element_text(color='black'),
            legend.key.height = unit(0.3,"cm"),
            legend.key.width = unit(0.15,"cm"))
    
    b
    ggsave("output/OSCC_output_202308/count_matrix/heatmap2-pvalue0.01.pdf",b,width=5.5,height=6)
  }
  
  # pheatmap
  {
    logRPM. <- logRPM[,order(class,decreasing=F)]
    library(pheatmap)
    ann_col = data.frame(class = as.character(class.))
    rownames(ann_col) = colnames(logRPM.scale.)
    ann_colors = list(class = brewer.pal(4,"Dark2")[1:2])
    names(ann_colors$class) = c("Primary","Metastatic")
    #GeneClass = c(Path1 = "#807DBA", Path2 = "#9E9AC8", Path3 = "#BCBDDC"))
    
    col = brewer.pal(9,"Set1")
    pdf("output/OSCC_output_202308/count_matrix/pheatmap-pvalue0.01.pdf", width = 6, height = 6)
    par(mar=c(2,2,2,2))
    pheatmap(logRPM., 
             color = colorRampPalette(c(col[2],"white",col[1]))(1000),
             #color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
             cutree_col = 1, 
             #cutree_row = 3, #break up the heatmap by clusters you define
             cluster_rows = F,
             cluster_cols = F, #by default, pheatmap clusters by both row and col
             show_rownames = T,
             show_colnames = F,
             fontsize_col = 12,
             angle_col = 0,
             annotation_col = ann_col,
             #annotation_row = ann_row,
             annotation_colors = ann_colors,
             annotation_names_row = F,
             #legend_breaks=-4:4,
             border=F)
    dev.off()
  }
}
# volcano plot
{
  de.plt <- de
  de.plt$threshold <- factor(ifelse(de.plt$PValue < 0.01 & abs(de.plt$logFC) >=1, ifelse(de.plt$logFC > 1 ,'UP-regulated','DOWN-regulated'),'Not'), levels = c('UP-regulated','DOWN-regulated',"Not"))

  col = brewer.pal(3,"Set1")
  {
    p=ggplot(de.plt) + 
      #geom_point(data=de.plt, aes(x=logFC, y =-log10(FDR), color=threshold), size=2, alpha=0.9, shape=16)+
      geom_point(data=de.plt, aes(x=logFC, y =-log10(PValue), color=threshold), size=2, alpha=0.9, shape=16)+
      scale_color_manual(values=c(col[1],col[2],"grey"))+
      #scale_y_continuous(limits=c(0,2))+
      theme_bw(base_size = 12) +
      geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6)+
      geom_hline(yintercept = -log10(0.01),lty=4,col="grey",lwd=0.6)+
      guides(fill="none")+
      theme(legend.position=c(0.78,0.9),
            panel.grid=element_blank(),
            panel.border=element_rect(size=1, fill="transparent"),
            legend.title = element_blank(),
            legend.text= element_text(color="black", size=10),
            legend.key.size = unit(0.3,"cm"),
            plot.title = element_text(hjust = 0.5),
            strip.background = element_blank(),
            strip.text = element_text(color="black", size=15),
            axis.text.x = element_text(color="black", size=15),
            axis.text.y = element_text(color="black", size=15),
            axis.title.x = element_text(color="black", size=15),
            axis.title.y = element_text(color="black", size=15))+
      #labs(x="log2FoldChange",y="-log10 (FDR)",title="")
      labs(x="log2FoldChange",y="-log10 (PValue)",title="")
    
    p
  }
  ggsave("1001/volcano.pdf",p,width = 4.3,height = 4.5)
}











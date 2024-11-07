########################################
## TE for DM EV
## liuke modified from zhanqing
## 202408

setwd("/Users/keliu/Desktop/IM EV/")
figdir="/Users/keliu/Desktop/IM EV/TE/"

library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(edgeR)
library(RColorBrewer)

### ref for TEtranscripts
{
  ref <- read.table("reference/hg38_rmsk_TE_20200804.gtf",sep="\t",header=F)
  colnames(ref) <- c("chr","source","type","start","end","n","strand",".","info")
  ref[1:4,]
  te_info <- data.frame(unlist(str_split(ref$info, "[ ;]+", simplify = TRUE)))
  te_info <- te_info[,-c(1,3,5,7,9)]
  colnames(te_info) <- c("gene_id","transcript_id","family_id","class_id")
  dim(te_info)
  # 139188
  te_info$gene_id <- str_split(te_info$gene_id," ",simplify=T)[,2]
  te_info$transcript_id <- str_split(te_info$transcript_id," ",simplify=T)[,3]
  te_info$family_id <- str_split(te_info$family_id," ",simplify=T)[,3]
  te_info$class_id <- str_split(te_info$class_id," ",simplify=T)[,3]
  #write.table(te_info,"reference/hg38_rmsk_TE_ref1.txt",sep="\t",quote=F,row.names=F)
  
  te_info <- te_info[!duplicated(te_info$gene_id),-2]
  dim(te_info)
  # 1030
  #write.table(te_info,"reference/hg38_rmsk_TE_ref2.txt",sep="\t",quote=F,row.names=F)
  table(te_info$class_id)
  table(te_info$family_id)
}

################
# TEtranscripts —— Detector2-seq data
# 155 MDA5\ARS\HC EV data
{
  count <- read.table("TEtranscript_count_matrix.txt",sep="\t",header=T,check.names=F,row.names=1)
  filtered_vector <- colnames(count)[!grepl("GW_", colnames(count))]
  cleaned_vector <- gsub("GP_", "", filtered_vector)
  cleaned_vector <- gsub("-BC", "_BC", cleaned_vector)
  count <- count[,filtered_vector]
  colnames(count) <- cleaned_vector
  metadata <- read.table("annotation.txt", sep = "\t", header = T,check.names = F)
  metadata <- metadata[which(metadata$mistake == 0),]
  
  overlap <- intersect(metadata$sample_library_id,colnames(count))
  diff <- setdiff(metadata$sample_library_id,colnames(count))
  metadata <- metadata[metadata$sample_library_id %in% overlap,]
  count <- count[,metadata$sample_library_id]
  colnames(count) <- metadata$sample_id

  type <- metadata$group
  #class <- metadata$class
  
}
# TE class
{
  gene <- rownames(count)
  
  te <- gene[substr(gene,1,4)!="ENSG"]
  te_name <- str_split(te,":",simplify=T)[,1]
  te_family <- str_split(te,":",simplify=T)[,2]
  te_class <- str_split(te,":",simplify=T)[,3]
  
  tecount <- count[te,]
  
  tecount.sum <- aggregate(tecount, by=list(family=te_class), sum)
  rownames(tecount.sum) <- tecount.sum$family
  tecount.sum <- tecount.sum[,-1]
  colSums(count)
  colSums(tecount)
  
  tecount.sum.ratio <- as.data.frame(apply(tecount.sum, 1, function(x) x/colSums(count)))
  tecount.sum.ratio <- tecount.sum.ratio[,c("DNA","LINE","SINE","LTR","Retroposon","Satellite")]
  tecount.sum.ratio$Gencode <- 1-rowSums(tecount.sum.ratio)
  tecount.sum.ratio$sample_id <- rownames(tecount.sum.ratio)
  tecount.sum.ratio$type <- type
  #tecount.sum.ratio$class <- class
  #tecount.sum.ratio$class <- factor(tecount.sum.ratio$class, levels=c("CF","EV"))
  tecount.sum.ratio$type <- factor(tecount.sum.ratio$type, levels=c("HC","MDA5","ARS"))
  
  tecount.sum.ratio <- melt(tecount.sum.ratio, id.var=c("sample_id","type")) #,"class"
  
  # type
  {
    bb=pal_nejm()(8)
    custom_colors <- brewer.pal(5,'BrBG')[c(1,2,5)]
    names(custom_colors) <- c("MDA5","ARS","HC")
    tecount.sum.ratio$type <- factor(tecount.sum.ratio$type, levels = c("MDA5","ARS","HC"))
    test="DNA"
    tecount.sum.ratio.filer <- tecount.sum.ratio[tecount.sum.ratio$variable==test,]
    if (test == 'DNA'){
      tecount.sum.ratio.filer$variable = 'DNA transposon'
    }
    p=ggplot(tecount.sum.ratio.filer, aes(type, value, fill=type))+
      scale_fill_manual(values = custom_colors) +
      facet_grid(~ variable)+
      geom_boxplot(notch = FALSE, alpha = 0.9, size=0.4, outlier.shape = NA, position=position_dodge(0.9)) +
      geom_jitter(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9), alpha = 0.7, size=0.8)+
      #scale_y_continuous(limits = c(0,600))+
      theme_bw()+
      #scale_fill_igv()+
      theme_classic()+
      theme(panel.background=element_rect(fill="white", colour="black", size=0.5),
            plot.title = element_text(size=12, hjust=0.5),
            strip.text = element_text(size=12, color="black"),
            strip.background = element_blank(),
            axis.line=element_blank(),
            axis.title=element_text(size=15, color="black"),
            axis.text.x = element_text(size=12, color="black"),
            axis.text.y = element_text(size=12, color="black"),
            legend.position="none",
            legend.title= element_text(color="black", size=12),
            legend.text= element_text(color="black", size=12))+
      labs(y="Reads ratio",x="")
    p
    
    my_comparisons <- list(c("MDA5","ARS"),
                           c("HC","ARS"),
                           c("HC","MDA5"))
    p <- p+geom_signif(comparisons = my_comparisons,
                       map_signif_level = T,
                       test = wilcox.test,
                       step_increase = 0.07)
    p
  }
  
  ggsave(paste0(figdir,"EV_ratio-boxplot-",test,".pdf"),p,width=3,height=4)
  
  
  
}


######
#EVvsCF 9:9
####
{
  count <- read.table("EVvsCF_TEtranscript_count_matrix.txt",sep="\t",header=T,check.names=F,row.names=1)
  filtered_vector <- colnames(count)[!grepl("GW_", colnames(count))]
  cleaned_vector <- gsub("GP_", "", filtered_vector)
  cleaned_vector <- gsub("L1_", "", cleaned_vector)
  cleaned_vector <- gsub("2023", "23", cleaned_vector)
  cleaned_vector <- gsub("-BC", "_BC", cleaned_vector)
  count <- count[,filtered_vector]
  colnames(count) <- cleaned_vector
  metadata <- read.table("TE/EVvsCF_metadata.txt", sep = "\t", header = T,check.names = F)

  overlap <- intersect(metadata$library_id,colnames(count))
  diff <- setdiff(metadata$library_id,colnames(count))
  metadata <- metadata[metadata$library_id %in% overlap,]
  count <- count[,metadata$library_id]
  colnames(count) <- metadata$sample_id
  
  type <- metadata$subtype
  class <- metadata$class
  
}
# TE class
{
  gene <- rownames(count)
  
  te <- gene[substr(gene,1,4)!="ENSG"]
  te_name <- str_split(te,":",simplify=T)[,1]
  te_family <- str_split(te,":",simplify=T)[,2]
  te_class <- str_split(te,":",simplify=T)[,3]
  
  tecount <- count[te,]
  
  tecount.sum <- aggregate(tecount, by=list(family=te_class), sum)
  rownames(tecount.sum) <- tecount.sum$family
  tecount.sum <- tecount.sum[,-1]
  colSums(count)
  colSums(tecount)
  
  tecount.sum.ratio <- as.data.frame(apply(tecount.sum, 1, function(x) x/colSums(count)))
  tecount.sum.ratio <- tecount.sum.ratio[,c("DNA","LINE","SINE","LTR","Retroposon","Satellite")]
  tecount.sum.ratio$Gencode <- 1-rowSums(tecount.sum.ratio)
  tecount.sum.ratio$sample_id <- rownames(tecount.sum.ratio)
  tecount.sum.ratio$type <- type
  tecount.sum.ratio$class <- class
  tecount.sum.ratio$class <- factor(tecount.sum.ratio$class, levels=c("CF","EV"))
  tecount.sum.ratio$type <- factor(tecount.sum.ratio$type, levels=c("HC","MDA5","ARS"))
  
  tecount.sum.ratio <- melt(tecount.sum.ratio, id.var=c("sample_id","type","class")) #

}



# class - paired boxplot
{
  bb=pal_nejm()(8)
  test="Satellite"
  tecount.sum.ratio.filer <- tecount.sum.ratio[tecount.sum.ratio$variable==test,] #[,-c(2,4)]
  tecount.sum.ratio.filer$sample_id <- paste0(str_split_fixed(tecount.sum.ratio.filer$sample_id,'-',3)[,1],'-',
                                              str_split_fixed(tecount.sum.ratio.filer$sample_id,'-',3)[,2])
  if (test == 'DNA'){
    tecount.sum.ratio.filer$variable = 'DNA transposon'
  }
  p=ggplot(tecount.sum.ratio.filer, aes(class, value, fill=class)) +
    facet_grid(~ variable)+
    stat_boxplot(geom = "errorbar",width=0.4) + 
    geom_boxplot(notch = FALSE, size=0.4, outlier.shape = 16, outlier.size = 0.8,position=position_dodge(0.9)) +
    theme_bw()+
    scale_fill_manual(values=c(bb[2],bb[1]))+
    scale_x_discrete(labels=c("Plasma","EV"))+
    theme_classic()+
    theme(panel.background=element_rect(fill="white", colour="black", size=1),
          strip.background = element_blank(),
          strip.placement = "inside",
          strip.text = element_text(size=12,color="black"),
          strip.switch.pad.grid = unit(0, "cm"),
          strip.switch.pad.wrap = unit(0, "cm"),
          axis.line=element_blank(),
          axis.title=element_text(color="black",size=12),
          axis.text.y = element_text(color="black",size=12),
          axis.text.x = element_text(color="black",size=11),
          legend.position="none",
          legend.title= element_text(color="black", size=12),
          legend.text= element_text(color="black", size=12))+
    labs(y="Reads ratio",x="")
  p
  
  p <- p+ stat_compare_means(method = "wilcox.test", paired = TRUE, 
                       aes(group = class, label = after_stat(p.signif)),
                       group.by = "sample_id",
                       comparisons = list(c("CF", "EV"))) #+ geom_line(aes(group = sample_id), color = "black", size = 0.25, alpha = 0.5)
  p
  
}
#ggsave(paste0(figdir,"figure6/EVvsCF_TE-boxplit-",test,".pdf"),p,width=1.5,height=2)
ggsave(paste0(figdir,"EVvsCF_TE-boxplit-",test,".pdf"),p,width=3,height=4.2)


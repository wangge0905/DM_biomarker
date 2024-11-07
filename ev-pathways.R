########################################
###### differential expression for i-pico PhaseII multiplex
###### zhanqing
###### 20220322

.libPaths()
library(ggplot2)
library(gg.gap)
library(dplyr)
library(ggpubr)
library(ggsci)
library(tidyverse)
library(edgeR)
library(RUVSeq)
library(RColorBrewer)
library(VennDiagram)
library(GGally)
library(gplots)
library(pheatmap)
library(scales)
library(reshape2)
display.brewer.all()

setwd("/Users/wangge/Documents/DM/")

# gencode+circRNA
data <- read.delim2("gencode_rmbatch.all.txt", header = T, row.names = 1, check.names = F)
metainfo <- read.delim2("metainfo.txt", header = T, row.names = 1, check.names = F)
metainfo$group <- sapply(rownames(metainfo), function(x) str_split(x, "-")[[1]][1])
colnames(data) <- rownames(metainfo)

y <- DGEList(counts=data)
libsize = y$samples
#write.table(libsize, "cfRNA/paired-libsize.txt", sep="\t", quote=F)
### paired group
{
  # put cell-free as reference level
  class <- unlist(lapply(strsplit(colnames(data),"-",fixed=T),function(x) x[1]))
  class <- factor(subtype, levels = c("MDA5","ARS","HC"))
  #patient <- unlist(lapply(strsplit(colnames(data),"-EV",fixed=T),function(x) x[1]))
  #patient <- unlist(lapply(strsplit(patient,"-CF",fixed=T),function(x) x[1]))
  #libsize <- read.delim2("cfRNA/paired-libsize.txt", sep="\t", header=T, check.names=F)
}

data_MDA5vsHC <- data[,c(1:55,106:155)]
class_MDA5vsHC <- class[c(1:55,106:155)]
class_MDA5vsHC <- factor(class_MDA5vsHC, levels = c("HC","MDA5"))
data_ARSvsHC <- data[,c(56:155)]
class_ARSvsHC <- class[c(56:155)]
class_ARSvsHC <- factor(class_ARSvsHC, levels = c("HC","ARS"))
data_MDA5vsARS <- data[,c(1:105)]
class_MDA5vsARS <- class[c(1:105)]
class_MDA5vsARS <- factor(class_MDA5vsARS, levels = c("ARS","MDA5"))

########################################################
###### design: paired test in all samples
data <- data_MDA5vsHC
class <- class_MDA5vsHC
patient

fc.cutoff=1  
fdr.cutoff=0.1
difexp_cfEV_all <- function(data, class, my_coef, fc.cutoff, fdr.cutoff){
  design <- model.matrix(~class)
  # paired sample 和 EV cell-free class影响，不考虑癌症cancer之间的差异。
  
  y <- DGEList(counts=data,  group=class)
  y <- calcNormFactors(y, method="TMM")
  
  rownames(design) <- colnames(y)
  y <- estimateDisp(y, design)
  y$common.dispersion
  # 0.6840856
  fit.ql <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit.ql, coef= my_coef)
  de <- topTags(qlf, n=Inf)$table
  
  #edgeR::cpm(y)[rownames(de)[1:10],]
  #summary(decideTests(qlf))
  #plotMD(qlf)
  #abline(h=c(-1, 1), col="blue")
  
  de.down <- filter(de, logFC <= -fc.cutoff, FDR <= fdr.cutoff)
  de.up <- filter(de, logFC >= fc.cutoff, FDR <= fdr.cutoff)
  print(paste0("Control group enrichend RNA: ",nrow(de.down)))
  print(paste0("Experimental group enriched RNA: ",nrow(de.up)))
  return(list(de=de, 
              de.down=de.down, 
              de.up=de.up))
}
data_MDA5vsHC <- data[,c(1:55,106:155)]
class_MDA5vsHC <- class[c(1:55,106:155)]
class_MDA5vsHC <- factor(class_MDA5vsHC, levels = c("HC","MDA5"))
data_ARSvsHC <- data[,c(56:155)]
class_ARSvsHC <- class[c(56:155)]
class_ARSvsHC <- factor(class_ARSvsHC, levels = c("HC","ARS"))
data_MDA5vsARS <- data[,c(1:105)]
class_MDA5vsARS <- class[c(1:105)]
class_MDA5vsARS <- factor(class_MDA5vsARS, levels = c("ARS","MDA5"))

de.MDA5vsHC <- difexp_cfEV_all(data_MDA5vsHC, class_MDA5vsHC,"classMDA5", 1, 0.1)
#[1] "Control group enrichend RNA: 137"
#[1] "Experimental group enriched RNA: 5263"
de.ARSvsHC <- difexp_cfEV_all(data_ARSvsHC, class_ARSvsHC,"classARS", 1, 0.1)
#[1] "Control group enrichend RNA: 154"
#[1] "Experimental group enriched RNA: 7183"
de.MDA5vsARS <- difexp_cfEV_all(data_MDA5vsARS, class_MDA5vsARS,"classMDA5", 1, 0.1)
#[1] "Control group enrichend RNA: 188"
#[1] "Experimental group enriched RNA: 99"

### summary differential genes
{
  write.table(de.MDA5vsHC$de, "pathway_enrichment/de.MDA5vsHC_de.txt", sep="\t", quote=F)
  write.table(de.ARSvsHC$de, "pathway_enrichment/de.ARSvsHC_de.txt", sep="\t", quote=F)
  write.table(de.MDA5vsARS$de, "pathway_enrichment/de.MDA5vsARS_de.txt", sep="\t", quote=F)

  fc.cutoff=0   ### change logFC cutoff to 0!!!
  fdr.cutoff=0.1
  de.cf <- filter(de.all$de, logFC <= -fc.cutoff, FDR <= fdr.cutoff)
  de.ev <- filter(de.all$de, logFC >= fc.cutoff, FDR <= fdr.cutoff)
  de.neutral <- filter(de.all$de, FDR>fdr.cutoff)
  
  write.table(unlist(lapply(strsplit(rownames(de.cf),".",fixed=T),function(x) x[1])), 
              "result/paired/diffgene/diffgene-CFEV-FCcutoff0/de.all_de.cf.txt", sep="\t", quote=F, row.names = F, col.names = F)
  write.table(unlist(lapply(strsplit(rownames(de.ev),".",fixed=T),function(x) x[1])), 
              "result/paired/diffgene/diffgene-CFEV-FCcutoff0/de.all_de.ev.txt", sep="\t", quote=F, row.names = F, col.names = F)
  write.table(unlist(lapply(strsplit(rownames(de.neutral),".",fixed=T),function(x) x[1])), 
              "result/paired/diffgene/diffgene-CFEV-FCcutoff0/de.all_de.neutral.txt", sep="\t", quote=F, row.names = F, col.names = F)
  
}
### volcano plot
{
  de.plt <- de.all$de
  de.plt$threshold <- factor(ifelse(de.plt$FDR < 0.1 & abs(de.plt$logFC) >=1, ifelse(de.plt$logFC > 1 ,'Up','Down'),'Not'), levels = c("Down","Up","Not"))
  
  de.nc.plt <- de.nc$de
  de.nc.plt$threshold <- factor(ifelse(de.nc.plt$FDR < 0.1 & abs(de.nc.plt$logFC) >=1, ifelse(de.nc.plt$logFC > 1 ,'Up','Down'),'Not'), levels = c("Down","Up","Not"))
  de.crc.plt <- de.crc$de
  de.crc.plt$threshold <- factor(ifelse(de.crc.plt$FDR < 0.1 & abs(de.crc.plt$logFC) >=1, ifelse(de.crc.plt$logFC > 1 ,'Up','Down'),'Not'), levels = c("Down","Up","Not"))
  de.luad.plt <- de.luad$de
  de.luad.plt$threshold <- factor(ifelse(de.luad.plt$FDR < 0.1 & abs(de.luad.plt$logFC) >=1, ifelse(de.luad.plt$logFC > 1 ,'Up','Down'),'Not'), levels = c("Down","Up","Not"))

  col = brewer.pal(3,"Set1")
  {
    p=ggplot(de.plt, aes(x=logFC, y =-log10(FDR), color=threshold)) + 
      scale_color_manual(values=c(col[2],col[1],"grey"))+
      geom_point(size=2, alpha=0.8, shape=16, stroke = 0)+
      #scale_y_continuous(limits=c(0,2.8))+
      theme_bw(base_size = 12) +
      geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6)+
      geom_hline(yintercept = -log10(0.1),lty=4,col="grey",lwd=0.6)+
      theme(legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_rect(size=1, fill="transparent"),
            legend.title = element_blank(),
            legend.text= element_text(color="black", size=15),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(color="black", size=15),
            axis.text.y = element_text(color="black", size=15),
            axis.title.x = element_text(color="black", size=15),
            axis.title.y = element_text(color="black", size=15))+
      labs(x="log2FoldChange",y="-log10 (adjusted p-value)",title="")
    p
  }
  {
    p=ggplot(de.crc.plt, aes(x=logFC, y =-log10(FDR), color=threshold)) + 
      scale_color_manual(values=c(col[2],col[1],"grey"))+
      geom_point(size=2, alpha=0.8, shape=16, stroke = 0)+
      #scale_y_continuous(limits=c(0,2.8))+
      theme_bw(base_size = 12) +
      geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6)+
      geom_hline(yintercept = -log10(0.1),lty=4,col="grey",lwd=0.6)+
      theme(legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_rect(size=1, fill="transparent"),
            legend.title = element_blank(),
            legend.text= element_text(color="black", size=15),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(color="black", size=15),
            axis.text.y = element_text(color="black", size=15),
            axis.title.x = element_text(color="black", size=15),
            axis.title.y = element_text(color="black", size=15))+
      labs(x="log2FoldChange",y="-log10 (p-value)",title="")
    p
  }
  ggsave("result/paired/difexp_crc_volcano.pdf",p,width = 4.4,height = 3.8)
}
### heatmap
{
  {
    difgenes <- c(rownames(de.all$de.cf),rownames(de.all$de.ev))
    
    logRPM <- edgeR::cpm(y, log=T)
    logRPM <- logRPM[difgenes,]
    logRPM.scale <- scale(t(logRPM), center = T, scale = T)
    logRPM.scale <- t(logRPM.scale)
    
    class. <- class[order(class)]
    subtype. <- subtype[order(class,subtype)]
    logRPM.scale. <- logRPM.scale[,order(class,subtype)]
  }
  {
    difgenes <- c(rownames(de.nc$de.cf),rownames(de.nc$de.ev))
    difgenes <- c(rownames(de.crc$de.cf),rownames(de.crc$de.ev))
    difgenes <- c(rownames(de.luad$de.cf),rownames(de.luad$de.ev))

    logRPM <- edgeR::cpm(y, log=T)
    logRPM <- logRPM[difgenes,]
    logRPM.scale <- scale(t(logRPM), center = T, scale = T)
    logRPM.scale <- t(logRPM.scale)
    
    cancer. <- cancer[cancer=="NC"]
    class. <- class[cancer=="NC"]
    logRPM.scale. <- logRPM.scale[,cancer=="NC"]
    
    cancer. <- cancer[cancer=="CRC"]
    class. <- class[cancer=="CRC"]
    logRPM.scale. <- logRPM.scale[,cancer=="CRC"]
    
    cancer. <- cancer[cancer=="LUAD"]
    class. <- class[cancer=="LUAD"]
    logRPM.scale. <- logRPM.scale[,cancer=="LUAD"]
  
  }
  {
    nc.cf <- rownames(filter(de.all.nc$de, logFC<0, FDR<0.1))
    nc.ev <- rownames(filter(de.all.nc$de, logFC>0, FDR<0.1))
    crc.cf <- rownames(filter(de.all.crc$de, logFC<0, FDR<0.1))
    crc.ev <- rownames(filter(de.all.crc$de, logFC>0, FDR<0.1))   
    lc.cf <- rownames(filter(de.all.lc$de, logFC < -1, FDR<0.05))
    lc.ev <- rownames(filter(de.all.lc$de, logFC > 1, FDR<0.05))
    
    df = rbind(data.frame(features=Reduce(intersect, list(nc.cf,crc.cf)),group="Universal CF enrichment"),
               data.frame(features=intersect(setdiff(nc.cf,crc.cf),setdiff(nc.cf,lc.cf)),group="NC-specific CF enrichment"),
               data.frame(features=intersect(setdiff(nc.ev,crc.ev),setdiff(nc.ev,lc.ev)),group="NC-specific EV enrichment"),
               data.frame(features=intersect(setdiff(crc.cf,nc.cf),setdiff(crc.cf,lc.cf)),group="CRC-specific CF enrichment"),
               data.frame(features=intersect(setdiff(crc.ev,nc.ev),setdiff(crc.ev,lc.ev)),group="CRC-specific EV enrichment"),
               data.frame(features=intersect(setdiff(lc.cf,nc.cf),setdiff(lc.cf,crc.cf)),group="LC-specific CF enrichment"),
               data.frame(features=intersect(setdiff(lc.ev,nc.ev),setdiff(lc.ev,crc.ev)),group="LC-specific EV enrichment"))
    
    logRPM <- edgeR::cpm(y, log=T)
    logRPM <- logRPM[df$features,]
    logRPM.scale <- scale(t(logRPM), center = T, scale = T)
    logRPM.scale <- t(logRPM.scale)
    cancer.=cancer
    class.=class
    logRPM.scale.=logRPM.scale
  }
  ## pheatmap::pheatmap 
  {
    library(pheatmap)
    ann_col = data.frame(subtype. = as.character(subtype.),
                         class. = as.character(class.))
    rownames(ann_col) = colnames(logRPM.scale.)
    #ann_row = data.frame(GeneClass = factor(rep(c("Universal CF enrichment","NC-specific CF enrichment","NC-specific EV enrichment","CRC-specific CF enrichment","CRC-specific EV enrichment","LC-specific CF enrichment","LC-specific EV enrichment"), c(8,2,1,13,15,90,12))))
    #rownames(ann_row) = df$features
    #set colors of each group
    ann_colors = list(subtype. = brewer.pal(4,"Set2")[1:3], 
                      class. = brewer.pal(3,"Set1")[c(2,1)])
    names(ann_colors$subtype.) = c("MDA5","ARS","HC")
    names(ann_colors$class.) = c("CF","EV")
    #GeneClass = c(Path1 = "#807DBA", Path2 = "#9E9AC8", Path3 = "#BCBDDC"))
    
    col = brewer.pal(9,"Set1")
    #pdf("result/paired/heatmap/pheatmap_universal-specific.pdf", width = 7, height = 5)
    par(mar=c(2,2,2,2))
    a1=pheatmap(logRPM.scale., 
             #color = colorRampPalette(c(col[2],"white",col[1]))(1000),
             color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
             cutree_col = 1, 
             #cutree_row = 3, #break up the heatmap by clusters you define
             cluster_rows = F,
             cluster_cols = F, #by default, pheatmap clusters by both row and col
             show_rownames = F,
             show_colnames = F,
             fontsize_col = 12,
             angle_col = 0,
             annotation_col = ann_col,
             #annotation_row = ann_row,
             annotation_colors = ann_colors,
             annotation_names_row = F,
             legend_breaks=-4:4,
             #labels_row = labels_row,
             border=F)
    dev.off()
  }
  
}
### Venn plot
{
  rownames(de.crc$de.cf)
  library(VennDiagram)
  color = brewer.pal(6, "Set1")
  venn <- venn.diagram(
    x = list(#NC.cf=rownames(de.nc$de.cf),
      NC.EV=rownames(de.nc$de.ev),
      #CRC.cf=rownames(de.crc$de.cf),
      CRC.EV=rownames(de.crc$de.ev),
      #LUAD.cf=rownames(de.luad$de.cf)),
      LUAD.EV=rownames(de.luad$de.ev)),
    filename = NULL,
    fill = color[1:3],
    alpha = 0.4,
    lwd = 2,
    lty = "longdash",
    col = "black", #线条色
    cex = 1.5,
    fontfamily = "sans",
    cat.col = "black",
    cat.default.pos = "text", #位置, outer内 text外
    cat.cex = 1.5,
    cat.fontfamily = "sans",     
    cat.dist = c(0.08, 0.12, 0.08), #位置，用圆的距离
    cat.pos = c(0,0,90), #位置，用圆的度数
    main = "",
    main.col = "black",
    main.cex = 1.5
  );
  grid.draw(venn);
  ggsave("result/paired/diffgenes_venn_ev.pdf",venn,width = 4,height = 4)
}  
{
  venn=ggvenn(list(cf=cf.genes, EV=ev.genes),
              fill_color=color[1:2],
              stroke_linetype="longdash",
              set_name_size=8,
              text_size=6,
              show_percentage = F)
  ggsave("result/paired/gene_number_cf-EV_TPM>1_venn2.pdf",venn,width = 5,height = 5)
}
### diffgenes 
{
  ### genetype 
  {
    # FC.cutoff=1
    difgenes = rbind(data.frame(gene=rownames(de.all$de.cf))%>%mutate(group="ALL.CF"),
                     data.frame(gene=rownames(de.all$de.ev))%>%mutate(group="ALL.EV"),
                     data.frame(gene=rownames(de.all.nc$de.cf))%>%mutate(group="NC.CF"),
                     data.frame(gene=rownames(de.all.nc$de.ev))%>%mutate(group="NC.EV"),
                     data.frame(gene=rownames(de.all.crc$de.cf))%>%mutate(group="CRC.CF"),
                     data.frame(gene=rownames(de.all.crc$de.ev))%>%mutate(group="CRC.EV"),
                     data.frame(gene=rownames(de.all.lc$de.cf))%>%mutate(group="LC.CF"),
                     data.frame(gene=rownames(de.all.lc$de.ev))%>%mutate(group="LC.EV"))
    # FC.cutoff=0
    fc.cutoff=0
    fdr.cutoff=0.1
    difgenes = rbind(data.frame(gene=rownames(filter(de.all$de, logFC <= -fc.cutoff, FDR <= fdr.cutoff)))%>%mutate(group="ALL.CF"),
                     data.frame(gene=rownames(filter(de.all$de, logFC >= fc.cutoff, FDR <= fdr.cutoff)))%>%mutate(group="ALL.EV"),
                     data.frame(gene=rownames(filter(de.all.nc$de, logFC <= -fc.cutoff, FDR <= fdr.cutoff)))%>%mutate(group="NC.CF"),
                     data.frame(gene=rownames(filter(de.all.nc$de, logFC >= fc.cutoff, FDR <= fdr.cutoff)))%>%mutate(group="NC.EV"),
                     data.frame(gene=rownames(filter(de.all.crc$de, logFC <= -fc.cutoff, FDR <= fdr.cutoff)))%>%mutate(group="CRC.CF"),
                     data.frame(gene=rownames(filter(de.all.crc$de, logFC >= fc.cutoff, FDR <= fdr.cutoff)))%>%mutate(group="CRC.EV"),
                     data.frame(gene=rownames(filter(de.all.lc$de, logFC <= -fc.cutoff, FDR <= fdr.cutoff)))%>%mutate(group="LC.CF"),
                     data.frame(gene=rownames(filter(de.all.lc$de, logFC >= fc.cutoff, FDR <= fdr.cutoff)))%>%mutate(group="LC.EV"))
    difgenes$genetype = unlist(lapply(strsplit(difgenes$gene,"|",fixed=T),function(x) x[4]))
    difgenes$genetype[substr(difgenes$gene,1,3)=="hsa"] = "circRNA"
    table(difgenes$genetype)
    difgenes$genetype = factor(difgenes$genetype, levels=rev(c("protein_coding","lncRNA","circRNA","snRNA","snoRNA","other_ncRNA")))
    difgenes$genetype[is.na(difgenes$genetype)] = "other_ncRNA"
    difgenes$class =  unlist(lapply(strsplit(difgenes$group,".",fixed=T),function(x) x[1]))
    difgenes$group =  unlist(lapply(strsplit(difgenes$group,".",fixed=T),function(x) x[2]))
    difgenes$class = factor(difgenes$class, levels=c("ALL","NC","CRC","LC"))

    {
      aa=pal_igv()(8)
      show_col(aa)
      d=ggplot(difgenes[difgenes$class=="LC",], aes(x=group, fill=genetype))+
        #facet_grid(~class)+
        geom_bar(stat="count",position="stack",width=0.7, col="black")+
        scale_fill_manual(limits=rev(c("protein_coding","lncRNA","circRNA","snRNA","snoRNA","other_ncRNA")),
                          values=c("#5050FFFF","#CE3D32FF","#802268FF","#F0E685FF","#BA6338FF","#5DB1DDFF"))+
        scale_x_discrete(label=c("CF-enriched","EV-enriched"))+
        theme_bw()+
        guides(fill=guide_legend(title=NULL, ncol=1, reverse=T))+
        theme(
          plot.margin = unit(c(1, 1, 0, 1),"cm"),
          panel.border=element_blank(),
          panel.grid=element_blank(),
          legend.position="right",
          legend.background = element_blank(),
          legend.title = element_blank(),
          legend.text= element_text(color="black", size=11),
          legend.key.size = unit(0.4,"cm"),
          plot.title = element_text(color="black", size=15, hjust=0.5),
          axis.line = element_line(color="black", size=0.5),
          axis.text.x = element_text(color="black", size=12, angle=30,hjust=1,vjust=1),
          axis.text.y = element_text(color="black", size=15),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color="black", size=15))+
        labs(y="Gene numbers",title="LC")
      #geom_vline(aes(xintercept=6.5))+
     a
    }
    p=ggarrange(a,b,c,d,ncol=2,nrow=2,common.legend = T,align="h")
    p
    ggsave("result/paired/diffgene/diff_genetype_barplot_FCcutoff0.pdf",p,width = 5, height = 8)
    ggsave("result/paired/diffgene/diff_genetype_barplot_all_FCcutoff0.pdf.pdf",a,width = 4, height = 3.3)
    
  }
  ### gene
  {
    #sort by logFC
    aa <- de.all$de.cf[order(de.all$de.cf$logFC),]
    bb <- de.all$de.ev[order(-de.all$de.ev$logFC),]
    difgenes <- rbind(bb[1:10,]%>%mutate(group="EV"),aa[c(1:4,6:11),]%>%mutate(group="CF"))
    #NC
    bb <- filter(de.all.nc$de, logFC>0, FDR<0.1)
    aa <- filter(de.all.nc$de, logFC<0, FDR<0.1)
    difgenes <- rbind(aa[c(1:2,4:11),] %>% mutate(group="CF"),bb%>%mutate(group="EV"))
    #CRC
    bb <- filter(de.all.crc$de, logFC>0, FDR<0.1)
    aa <- filter(de.all.crc$de, logFC<0, FDR<0.1)
    difgenes <- rbind(aa[c(1:5,7:11),] %>% mutate(group="CF"),bb[c(1:10),]%>%mutate(group="EV"))
    #LC
    bb <- filter(de.all.lc$de, logFC>0, FDR<0.1)
    aa <- filter(de.all.lc$de, logFC<0, FDR<0.1)
    difgenes <- rbind(aa[c(1:2,4:11),] %>% mutate(group="CF"),bb[c(1:10),]%>%mutate(group="EV"))
    
    difgenes$gene <- unlist(lapply(strsplit(rownames(difgenes),"|",fixed=T),function(x) x[3]))
    difgenes$gene[is.na(difgenes$gene)] <- "hsa_circ_0048555"
    difgenes$gene[difgenes$gene=="Y_RNA"] <- unlist(lapply(strsplit(rownames(difgenes[difgenes$gene=="Y_RNA",]),".",fixed=T),function(x) x[1]))
    difgenes$gene <- factor(difgenes$gene, levels=difgenes$gene[order(difgenes$logFC)])
    difgenes$logp <- -log10(difgenes$FDR)
    difgenes <- difgenes[order(-difgenes$logFC),]
    
    {
      #sort by FDR
      aa <- de.all$de.cf
      bb <- de.all$de.ev
      difgenes <- rbind(bb[1:5,]%>%mutate(group="EV"),aa[c(1,3:6),]%>%mutate(group="CF"))
      difgenes$gene <- unlist(lapply(strsplit(rownames(difgenes),"|",fixed=T),function(x) x[3]))
      difgenes$gene[is.na(difgenes$gene)] <- "hsa_circ_0048555"
      difgenes$gene <- factor(difgenes$gene, levels=difgenes$gene[order(difgenes$logFC)])
      difgenes$logp <- -log10(difgenes$FDR)
    }
    
    #barplot
    {
      library(ggnewscale)
      col2 = pal_nejm()(8)
      a=ggplot()+
        coord_flip()+
        geom_bar(data=difgenes[difgenes$group=="EV",],aes(x=gene,y=logFC,fill=logp),
                 stat="identity", position="stack", width=0.7,color="black",size=0.2)+
        scale_fill_gradient2(name="EV -log10(FDR)",limits=c(0,13.2),breaks=c(2,7,12), high=col2[1], guide="colorbar")+
        new_scale("fill")+
        geom_bar(data=difgenes[difgenes$group=="CF",],aes(x=gene,y=logFC,fill=logp),
                 stat="identity", position="stack", width=0.7,color="black",size=0.2)+
        scale_fill_gradient2(name="CF -log10(FDR)",limits=c(0,13.2),breaks=c(2,7,12), high=col2[2], guide="colorbar")+
        scale_y_continuous(position="right",limits = c(-4,4), breaks=c(-4,0,4),
                           sec.axis = sec_axis(~.,breaks=c(-2,2),label=c("CF-enriched","EV-enriched"),guide=guide_axis(angle = 30)))+
        scale_x_discrete(limits=difgenes$gene[order(difgenes$logFC)])+
        theme_bw()+
        theme(
          plot.margin = unit(c(0.5, 1, 0.5, 0.5),"cm"),
          panel.grid = element_blank(),
          panel.border = element_rect(size=0.5, fill="transparent"),
          #panel.border = element_blank(),
          #axis.line = element_line(size=0.5, color="black"),
          plot.title = element_text(color="black", size=12, hjust=0.5),
          legend.position = "right",
          legend.title = element_text(color="black", size=10),
          legend.text = element_text(color="black", size=10),
          legend.key.height = unit(0.3,"cm"),
          legend.key.width = unit(0.2,"cm"),
          axis.text.x = element_text(color="black", size=12),
          axis.text.y = element_text(color="black", size=12),
          axis.title.x = element_text(color="black", size=12),
          axis.title.y = element_blank())
      a
    }
    {
      library(ggnewscale)
      col2 = pal_nejm()(8)
      a=ggplot()+
        coord_flip()+
        geom_bar(data=difgenes[difgenes$group=="EV",],aes(x=gene,y=logFC,fill=logp),
                 stat="identity", position="stack", width=0.7,color="black",size=0.2)+
        scale_fill_gradient2(name="-log10(FDR)",limits=c(1,10),breaks=c(1,8), high=col2[1], guide="colorbar")+
        scale_y_continuous(limits = c(0,5.1), breaks=c(0,5))+
        scale_x_discrete(limits=rev(difgenes$gene[difgenes$group=="EV"]))+
        theme_bw()+
        theme(
          plot.margin = unit(c(0.5, 1, 0.5, 0.5),"cm"),
          panel.grid = element_blank(),
          panel.border = element_rect(size=1, fill="transparent"),
          plot.title = element_text(color="black", size=15, hjust=0.5),
          legend.position = "right",
          legend.title = element_text(color="black", size=10),
          legend.text = element_text(color="black", size=10),
          legend.key.height = unit(0.3,"cm"),
          legend.key.width = unit(0.2,"cm"),
          axis.text.x = element_text(color="black", size=12),
          axis.text.y = element_text(color="black", size=12),
          axis.title.x = element_text(color="black", size=12),
          axis.title.y = element_blank())+
        labs(title="EV-enriched")
      a
    }
    {
      library(ggnewscale)
      col2 = pal_nejm()(8)
      b=ggplot()+
        coord_flip()+
        geom_bar(data=difgenes[difgenes$group=="CF",],aes(x=gene,y=-logFC,fill=logp),
                 stat="identity", position="stack", width=0.7,color="black",size=0.2)+
        scale_fill_gradient2(name="-log10(FDR)",limits=c(0,10),breaks=c(1,8), high=col2[2], guide="colorbar")+
        scale_y_continuous(limits = c(0,5.1), breaks=c(0,5))+
        scale_x_discrete(limits=difgenes$gene[difgenes$group=="CF"])+
        theme_bw()+
        theme(
          plot.margin = unit(c(0.5, 1, 0.5, 0.5),"cm"),
          panel.grid = element_blank(),
          panel.border = element_rect(size=1, fill="transparent"),
          plot.title = element_text(color="black", size=15, hjust=0.5),
          legend.position = "right",
          legend.title = element_text(color="black", size=10),
          legend.text = element_text(color="black", size=10),
          legend.key.height = unit(0.3,"cm"),
          legend.key.width = unit(0.2,"cm"),
          axis.text.x = element_text(color="black", size=12),
          axis.text.y = element_text(color="black", size=12),
          axis.title.x = element_text(color="black", size=12),
          axis.title.y = element_blank())+
        labs(title="CF-enriched")
      b
    }
    p=ggarrange(b,a,nrow=2,align="v",heights = c(2.5,1))
    p=ggarrange(b,a,nrow=2,align="v",heights = c(1,1))
    p
    #ggsave("result/paired/diffgene/diff_genebar_topFDR_NC.pdf",p,width = 4.2, height = 4.5)
    ggsave("result/paired/diffgene/diff_genebar_topFDR_LC.pdf",p,width = 4.2, height = 6)
    
  }
}

### enrichment analysis
{
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(DOSE)

  {
    ref <- read.table("/Users/wangge/Documents/DM/reference/gene.gtf",sep = "\t",header = F)
    colnames(ref) <- c("chr","source","region","start","end","score","strand","phase","meta")
    geneid <- unlist(lapply(strsplit(ref[,9],".",fixed=T), function(x) x[1]))
    ref$ensg <- unlist(lapply(strsplit(geneid," ",fixed=T), function(x) x[2]))
    genename <- unlist(lapply(strsplit(ref[,9],";",fixed=T), function(x) x[3]))
    genename <- unlist(lapply(strsplit(genename," ",fixed=T), function(x) x[3]))
    ref$genename <- genename
    library(gprofiler2)
    tmp <- gprofiler2::gconvert(query = ref$ensg, organism = "hsapiens", target = "ENTREZGENE_ACC" )
    tmp <- tmp[,c("input","target")]
    colnames(tmp) <- c("ensg","ENTREZID")
    ref$ENTREZID <- tmp$ENTREZID[match(ref$ensg,tmp$ensg)]
    ref <- ref[!duplicated(ref$ensg),]
    table(duplicated(ref$ENTREZID))
    # 35124
    #write.table(ref,"/Users/wangge/Documents/DM/reference/gencode.v38.txt",quote = F,sep = "\t",col.names = T,row.names = F)
  }
  
  ##### initialize 
  de.order <- de.MDA5vsHC$de
  de.order <- de.ARSvsHC$de
  de.order <- de.MDA5vsARS$de
  
  pValue=1
  qValue=1
  pAdjust="BH"
  
  cutoff.FDR = 0.1
  cutoff.logFC = 0
  ORA <- function(de.order, method, cutoff.FDR, cutoff.logFC, outdir){
    de.order$ENSG <- rownames(de.order)
    de.order$ENSG <- as.character(lapply(strsplit(de.order$ENSG,"|",fixed=TRUE), function(x) x[1]))
    de.order$ENSG <- as.character(lapply(strsplit(de.order$ENSG,".",fixed=TRUE), function(x) x[1]))
    de.order <- de.order[!duplicated(de.order$ENSG),]
    ##### get all geneList (for GSEA)
    {
      #gene.list.ncbiid <- bitr(de.order$name, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
      #gene.list.ncbiid <- read.table("/BioII/lulab_b/zhanqing/cfRNAseq/reference/geneset/gencode.v38.txt",sep = "\t",header = T)
      gene.list.ncbiid <- read.table("/Users/wangge/Documents/DM/reference/gencode.v38.txt",sep = "\t",header = T)
      gene.list.ncbiid <- gene.list.ncbiid[!is.na(gene.list.ncbiid$ENTREZID),]
      gene.list.ncbiid <- gene.list.ncbiid[,c("ensg","ENTREZID")]
      gene.list.ncbiid <- gene.list.ncbiid[!duplicated(gene.list.ncbiid$ensg) & !duplicated(gene.list.ncbiid$ENTREZID),]
      colnames(gene.list.ncbiid)[1] <- "ENSG"
      
      #de.order <- merge(x = de.order, y = gene.list.ncbiid, by.x="name", by.y="ENSG")
      de.order <- dplyr::left_join(x = de.order, y = gene.list.ncbiid)
      de.order <- de.order[order(de.order$logFC, -de.order$PValue, decreasing = T),] # order by logFC then PValue
      
      gene.list.ncbiid <- de.order$logFC[!duplicated(de.order$ENTREZID)]
      names(gene.list.ncbiid) <- de.order$ENTREZID[!duplicated(de.order$ENTREZID)]
      message(paste0("remaining ENTREZID geneList lenth: ",length(gene.list.ncbiid)))
      gene.list.ensgid <- de.order$logFC
      names(gene.list.ensgid) <- de.order$ENSG
      message(paste0("remaining ENSEMBL geneList lenth: ",length(gene.list.ensgid)))
    }
    ##### get top geneList (for ORA) over representation analysis
    {
      if (method=="cutoff.FDR"){
        message(paste0("FDR cutoff: ", cutoff.FDR, "  logFC cutoff: ", cutoff.logFC))
        gene.top.ncbiid <- as.character(na.omit(de.order$ENTREZID[abs(de.order$logFC) > cutoff.logFC & de.order$FDR <= cutoff.FDR]))   # FDR cutoff
        gene.top.ncbiid.up <- as.character(na.omit(de.order$ENTREZID[de.order$logFC > cutoff.logFC & de.order$FDR <= cutoff.FDR]))
        gene.top.ncbiid.down <- as.character(na.omit(de.order$ENTREZID[de.order$logFC < -cutoff.logFC & de.order$FDR <= cutoff.FDR]))
        message(paste0("remaining ENTREZID top gene lenth: ",length(gene.top.ncbiid)))
        message(paste0("remaining ENTREZID top up gene lenth: ",length(gene.top.ncbiid.up)))
        message(paste0("remaining ENTREZID top down gene lenth: ",length(gene.top.ncbiid.down)))
        gene.top.ensgid <- de.order$ENSG[abs(de.order$logFC) > cutoff.logFC & de.order$FDR <= cutoff.FDR]   # FDR cutoff
        gene.top.ensgid.up <- de.order$ENSG[de.order$logFC > cutoff.logFC & de.order$FDR <= cutoff.FDR] 
        gene.top.ensgid.down <- de.order$ENSG[de.order$logFC < -cutoff.logFC & de.order$FDR <= cutoff.FDR]
        message(paste0("remaining ENSEMBL top gene lenth: ",length(gene.top.ensgid)))
        message(paste0("remaining ENSEMBL top up gene lenth: ",length(gene.top.ensgid.up)))
        message(paste0("remaining ENSEMBL top down gene lenth: ",length(gene.top.ensgid.down)))
      }
      else if (method=="cutoff.pvalue"){
        cutoff.pvalue = 0.05
        cutoff.logFC = 1
        message(paste0("p-value cutoff: ", cutoff.pvalue, "  logFC cutoff: ", cutoff.logFC))
        gene.top.ncbiid <- as.character(na.omit(de.order$ENTREZID[abs(de.order$logFC) > cutoff.logFC & de.order$PValue <= cutoff.pvalue]))   # p-value cutoff
        gene.top.ncbiid.up <- as.character(na.omit(de.order$ENTREZID[de.order$logFC > cutoff.logFC & de.order$PValue <= cutoff.pvalue])) 
        gene.top.ncbiid.down <- as.character(na.omit(de.order$ENTREZID[de.order$logFC < -cutoff.logFC & de.order$PValue <= cutoff.pvalue]))
        message(paste0("remaining ENTREZID top gene lenth: ",length(gene.top.ncbiid)))
        message(paste0("remaining ENTREZID top up gene lenth: ",length(gene.top.ncbiid.up)))
        message(paste0("remaining ENTREZID top down gene lenth: ",length(gene.top.ncbiid.down)))
        gene.top.ensgid <- de.order$ENSG[abs(de.order$logFC) > cutoff.logFC & de.order$PValue <= cutoff.pvalue]   # FDR cutoff
        gene.top.ensgid.up <- de.order$ENSG[de.order$logFC > cutoff.logFC & de.order$PValue <= cutoff.pvalue] 
        gene.top.ensgid.down <- de.order$ENSG[de.order$logFC < -cutoff.logFC & de.order$PValue <= cutoff.pvalue]
        message(paste0("remaining ENSEMBL top gene lenth: ",length(gene.top.ensgid)))
        message(paste0("remaining ENSEMBL top up gene lenth: ",length(gene.top.ensgid.up)))
        message(paste0("remaining ENSEMBL top down gene lenth: ",length(gene.top.ensgid.down)))
      }
      else if (method=="FDR.top"){
        FDR.top = 100
        message(paste0("FDR top: ", FDR.top))
        de.order <- de.order[order(-de.order$FDR, decreasing = T),]  # order by FDR
        gene.top.ncbiid <- as.character(na.omit(de.order$ENTREZID[1:FDR.top]))
        gene.top.ncbiid.up <- as.character(na.omit(de.order$ENTREZID[de.order$logFC > 0][1:FDR.top]))
        gene.top.ncbiid.down <- as.character(na.omit(de.order$ENTREZID[de.order$logFC < 0][1:FDR.top]))
        message(paste0("remaining ENTREZID top gene lenth: ",length(gene.top.ncbiid)))
        message(paste0("remaining ENTREZID top up gene lenth: ",length(gene.top.ncbiid.up)))
        message(paste0("remaining ENTREZID top down gene lenth: ",length(gene.top.ncbiid.down)))
        gene.top.ensgid <- de.order$ENSG[1:FDR.top]   
        gene.top.ensgid.up <- de.order$ENSG[1:FDR.top]
        gene.top.ensgid.down <- de.order$ENSG[1:FDR.top]
        message(paste0("remaining ENSEMBL top gene lenth: ",length(gene.top.ensgid)))
        message(paste0("remaining ENSEMBL top up gene lenth: ",length(gene.top.ensgid.up)))
        message(paste0("remaining ENSEMBL top down gene lenth: ",length(gene.top.ensgid.down)))
      }else{
        message(paste0("No methods specified"))
      }
    }
    # 1. Over-Representation Analysis 超几何检验或Fisher精确检验
    message("start ORA enrich...")
    # 1.1 ORA all gene
    {
      if(FALSE){
        enrich.GO <- enrichGO(
          gene = gene.top.ensgid,  # gene.top.ncbiid
          #universe = names(gene.list.ncbiid),  #选择background gene?
          OrgDb = org.Hs.eg.db,
          ont = "ALL",
          keyType = "ENSEMBL",  # ENTREZID
          pvalueCutoff = pValue,
          pAdjustMethod = pAdjust,
          qvalueCutoff = qValue)
      }
      enrich.GO.all <- enrichGO(
        gene = gene.top.ncbiid,
        #universe = names(gene.list.ncbiid),  #选择background gene?
        OrgDb = org.Hs.eg.db,
        ont = "ALL",
        keyType = "ENTREZID",
        pvalueCutoff = pValue,
        pAdjustMethod = pAdjust,
        qvalueCutoff = qValue
      )
      enrich.DO.all <- enrichDO(
        gene = gene.top.ncbiid, 
        ont = "DO",
        pvalueCutoff = pValue,
        pAdjustMethod = pAdjust,
        qvalueCutoff = qValue,
        #universe = names(gene.list.ncbiid),
        minGSSize = 5,
        maxGSSize = 500,
        readable = FALSE
      )
      enrich.KEGG.all <- enrichKEGG(
        gene = gene.top.ncbiid,
        keyType = "ncbi-geneid",
        pvalueCutoff = pValue,
        pAdjustMethod = pAdjust,
        qvalueCutoff = qValue
      )
    }
    # 1.2 ORA up gene
    {
      if(FALSE){
        enrich.GO.up <- enrichGO(
          gene = gene.top.ncbiid.up,
          #universe = names(gene.list.ncbiid),  #选择background gene?
          OrgDb = org.Hs.eg.db,
          ont = "ALL",
          keyType = "ENTREZID",
          pvalueCutoff = pValue,
          pAdjustMethod = pAdjust,
          qvalueCutoff = qValue)
      }
      enrich.GO.up <- enrichGO(
        gene = gene.top.ensgid.up,
        #universe = names(gene.list.ncbiid),  #选择background gene?
        OrgDb = org.Hs.eg.db,
        ont = "ALL",
        keyType = "ENSEMBL",
        pvalueCutoff = pValue,
        pAdjustMethod = pAdjust,
        qvalueCutoff = qValue
      )
      enrich.DO.up <- enrichDO(
        gene = gene.top.ncbiid.up,
        ont = "DO",
        pvalueCutoff = pValue,
        pAdjustMethod = pAdjust,
        qvalueCutoff = qValue,
        #universe = names(gene.list.ncbiid),
        minGSSize = 5,
        maxGSSize = 500,
        readable = FALSE
      )
      enrich.KEGG.up <- enrichKEGG(
        gene = gene.top.ncbiid.up,
        keyType = "ncbi-geneid",
        pvalueCutoff = pValue,
        pAdjustMethod = pAdjust,
        qvalueCutoff = qValue
      )
    }
    # 1.3 ORA down gene
    {
      if(FALSE){
        enrich.GO.down <- enrichGO(
          gene = gene.top.ncbiid.down,
          #universe = names(gene.list.ncbiid),  #选择background gene?
          OrgDb = org.Hs.eg.db,
          ont = "ALL",
          keyType = "ENTREZID",
          pvalueCutoff = pValue,
          pAdjustMethod = pAdjust,
          qvalueCutoff = qValue)
      }
      enrich.GO.down <- enrichGO(
        gene = gene.top.ensgid.down,
        #universe = names(gene.list.ncbiid),  #选择background gene?
        OrgDb = org.Hs.eg.db,
        ont = "ALL",
        keyType = "ENSEMBL",
        pvalueCutoff = pValue,
        pAdjustMethod = pAdjust,
        qvalueCutoff = qValue
      )
      enrich.DO.down <- enrichDO(
        gene = gene.top.ncbiid.down,
        ont = "DO",
        pvalueCutoff = pValue,
        pAdjustMethod = pAdjust,
        qvalueCutoff = qValue,
        #universe = names(gene.list.ncbiid),
        minGSSize = 5,
        maxGSSize = 500,
        readable = FALSE
      )
      enrich.KEGG.down <- enrichKEGG(
        gene = gene.top.ncbiid.down,
        keyType = "ncbi-geneid",
        pvalueCutoff = pValue,
        pAdjustMethod = pAdjust,
        qvalueCutoff = qValue
      )
    }
    # output table
    {
      types <- c(enrich.GO.all,enrich.DO.all,enrich.KEGG.all,
                 enrich.GO.up,enrich.DO.up,enrich.KEGG.up,
                 enrich.GO.down,enrich.DO.down,enrich.KEGG.down)
      names(types) <-c("enrich.GO.all","enrich.DO.all","enrich.KEGG.all",
                       "enrich.GO.up","enrich.DO.up","enrich.KEGG.up",
                       "enrich.GO.down","enrich.DO.down","enrich.KEGG.down")
      dir.create(paste0(outdir), showWarnings = T)
      message("start writing tables...")
      for(i in 1:length(types)){
        write.table(types[[i]],paste0(outdir,"/", names(types)[i], ".txt"),row.names = F,col.names = T, quote = F,sep = "\t")}
    }
  }
  outdir <- '/Users/wangge/Documents/DM/pathway_enrichment'
  #ORA(de.all$de, "cutoff.FDR", "result/paired/enrich/all-sample-table") # remaining ENSEMBL top gene lenth: 166, remaining ENSEMBL top up gene lenth: 95, remaining ENSEMBL top down gene lenth: 71
  #ORA(de.nc$de, "cutoff.pvalue", "result/paired/enrich/nc-sample-table") # remaining ENSEMBL top up gene lenth: 111, remaining ENSEMBL top down gene lenth: 201
  #ORA(de.crc$de, "cutoff.pvalue", "result/paired/enrich/crc-sample-table") # remaining ENSEMBL top up gene lenth: 220, remaining ENSEMBL top down gene lenth: 54
  #ORA(de.luad$de, "cutoff.FDR", "result/paired/enrich/luad-sample-table") # remaining ENSEMBL top up gene lenth: 69, remaining ENSEMBL top down gene lenth: 131

  
  ORA(de.MDA5vsHC$de, "cutoff.FDR", 0.1, 0, "pathway_enrichment/enrich_MDA5vsHC") # remaining ENSEMBL top gene lenth: 545, remaining ENSEMBL top up gene lenth: 274, remaining ENSEMBL top down gene lenth: 271
  ORA(de.ARSvsHC$de, "cutoff.FDR", 0.1, 0, "pathway_enrichment/enrich_ARSvsHC") # remaining ENSEMBL top gene lenth: 545, remaining ENSEMBL top up gene lenth: 274, remaining ENSEMBL top down gene lenth: 271
  ORA(de.MDA5vsARS$de, "cutoff.FDR", 0.1, 0, "pathway_enrichment/enrich_MDA5vsARS") # remaining ENSEMBL top gene lenth: 545, remaining ENSEMBL top up gene lenth: 274, remaining ENSEMBL top down gene lenth: 271
  
  
  
  # ORA enrichment plot
  {
    options(scipen=1000000) 
    tabledir="pathway_enrichment/enrich_ARSvsHC/"
    outdir="pathway_enrichment/enrich_ARSvsHC/plot/"
    dir.create(outdir, showWarnings = T)
    col = brewer.pal(3,"Set1")
    
    ### enrichplot::dotplot
    {
      pdf(paste0(outdir,"dotplot-enrich.KEGG.up",".pdf"), width = 5, height = 4)
      enrichplot::dotplot(enrich.KEGG.up, font.size=12, title="enrich.KEGG.up", showCategory = 10) +
        scale_color_continuous(low=col[1],high=col[2]) +
        scale_x_continuous(limits = c(0,1), breaks = c(0,0.5,1))
      dev.off()
      pdf(paste0(outdir,"dotplot-enrich.KEGG.down",".pdf"), width = 5, height = 4)
      enrichplot::dotplot(enrich.KEGG.down, font.size=12, title="enrich.KEGG.down", showCategory = 10) +
        scale_color_continuous(low=col[1],high=col[2]) +
        scale_x_continuous(limits = c(0,0.5), breaks = c(0,0.25,0.5))
      dev.off()
      
      pdf(paste0(outdir,"dotplot-enrich.DO.up",".pdf"), width = 5.5, height = 4.5)
      enrichplot::dotplot(enrich.DO.up, font.size=12, title="enrich.DO.up", showCategory = 10) +
        scale_color_continuous(low=col[1],high=col[2]) +
        scale_x_continuous(limits = c(0,0.3), breaks = c(0.1,0.2))
      dev.off()
      
      pdf(paste0(outdir,"dotplot-enrich.DO.down",".pdf"), width = 5, height = 4)
      enrichplot::dotplot(enrich.DO.down, font.size=12, title="enrich.DO.down", showCategory = 10) +
        scale_color_continuous(low=col[1],high=col[2]) +
        scale_x_continuous(limits = c(0,0.5), breaks = c(0,0.25,0.5))
      dev.off()
    }
    ### GO
    {
      all.go.up <- read.table(paste0(tabledir,"enrich.GO.up.txt"),sep="\t",header=T,quote="")
      all.go.up <- all.go.up[order(all.go.up$p.adjust),]
      all.go.up <- all.go.up[1:10,]
      all.go.up$logp <- -log10(all.go.up$p.adjust)
      all.go.up$ONTOLOGY <- factor(all.go.up$ONTOLOGY, levels=c("CC","BP","MF"))
      all.go.up$Description <- factor(all.go.up$Description, levels=rev(all.go.up$Description[order(all.go.up$ONTOLOGY,all.go.up$p.adjust)]))
      
      all.go.down <- read.table(paste0(tabledir,"enrich.GO.down.txt"),sep="\t",header=T,quote="")
      all.go.down <- all.go.down[order(all.go.down$p.adjust),]
      all.go.down <- all.go.down[1:10,]
      all.go.down$logp <- -log10(all.go.down$p.adjust)
      all.go.down$ONTOLOGY <- factor(all.go.down$ONTOLOGY, levels=c("CC","BP","MF"))
      all.go.down$Description <- factor(all.go.down$Description, levels=rev(all.go.down$Description[order(all.go.down$ONTOLOGY,all.go.down$p.adjust)]))
    }
    ### KEGG
    {
      all.kegg.up <- read.table(paste0(tabledir,"enrich.KEGG.up.txt"),sep="\t",header=T,quote="")
      all.kegg.up <- all.kegg.up[order(all.kegg.up$p.adjust),]
      all.kegg.up <- all.kegg.up[1:10,]
      all.kegg.up$logp <- -log10(all.kegg.up$p.adjust)
      all.kegg.up$Description <- factor(all.kegg.up$Description, levels=rev(all.kegg.up$Description[order(all.kegg.up$p.adjust)]))
      
      all.kegg.down <- read.table(paste0(tabledir,"enrich.KEGG.down.txt"),sep="\t",header=T,quote="")
      all.kegg.down <- all.kegg.down[order(all.kegg.down$p.adjust),]
      all.kegg.down <- all.kegg.down[1:10,]
      all.kegg.down$logp <- -log10(all.kegg.down$p.adjust)
      all.kegg.down$Description <- factor(all.kegg.down$Description, levels=rev(all.kegg.down$Description[order(all.kegg.down$p.adjust)]))
    }
    ### ggplot
    {
      #GO
      p=ggplot(all.go.down, aes(logp, Description)) +
        geom_point(aes(color = ONTOLOGY, size = Count)) +
        scale_color_aaas()+
        scale_x_continuous(limits=c(1,4),breaks=c(1,2,3,4))+
        theme_bw()+
        guides(fill=guide_legend(title=NULL, ncol=1))+
        coord_fixed(ratio=1)+
        theme(
          plot.margin = unit(c(1, 1, 1, 1),"cm"),
          panel.border=element_rect(size=1, fill="transparent"),
          panel.grid=element_line(size=rel(0.5)),
          #legend.position=c(0.25,0.85),
          legend.position="right",
          legend.background = element_blank(),
          legend.title = element_text(color="black", size=10),
          legend.text= element_text(color="black", size=10),
          axis.text.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=10),
          axis.title.x = element_text(color="black", size=10),
          axis.title.y = element_blank())+
        labs(title="GO.down")+
        xlab("-log10(p.adjust)")
      p
      
      p=ggplot(all.go.up, aes(logp, Description)) +
        geom_point(aes(color = ONTOLOGY, size = Count)) +
        scale_color_aaas()+
        scale_x_continuous(limits=c(5,15),breaks=c(5,10,15))+
        theme_bw()+
        guides(fill=guide_legend(title=NULL, ncol=1))+
        coord_fixed(ratio=1)+
        theme(
          plot.margin = unit(c(1, 1, 1, 1),"cm"),
          panel.border=element_rect(size=1, fill="transparent"),
          panel.grid=element_line(size=rel(0.5)),
          #legend.position=c(0.25,0.85),
          legend.position="right",
          legend.background = element_blank(),
          legend.title = element_text(color="black", size=10),
          legend.text= element_text(color="black", size=10),
          axis.text.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=10),
          axis.title.x = element_text(color="black", size=10),
          axis.title.y = element_blank())+
        labs(title="GO.up")+
        xlab("-log10(p.adjust)")
      p
    }
    ggsave(paste0(outdir,"dotplot-enrich.go.up",".pdf"), width =15, height = 4)
    {
      #KEGG
      p=ggplot(all.kegg.down, aes(logp, Description)) +
        geom_point(aes(color=col[1],size = Count)) +
        scale_color_aaas()+
        #scale_x_continuous(limits=c(0.5, 1.8),breaks=c(0.8,1.6))+
        scale_x_continuous(limits=c(0,10),breaks=c(2,4,6,8,10))+
        theme_bw()+
        guides(color="none")+
        coord_fixed(ratio=0.6)+
        theme(
          plot.margin = unit(c(1, 1, 1, 1),"cm"),
          panel.border=element_rect(size=1, fill="transparent"),
          panel.grid=element_line(size=rel(0.5)),
          #legend.position=c(0.25,0.85),
          legend.position="right",
          legend.background = element_blank(),
          legend.title = element_text(color="black", size=10),
          legend.text= element_text(color="black", size=10),
          axis.text.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=10),
          axis.title.x = element_text(color="black", size=10),
          axis.title.y = element_blank())+
        labs(title="KEGG.up")+
        xlab("-log10(p.adjust)")
      p
    }
    ggsave(paste0(plotdir,"dotplot-enrich.kegg.up-2",".pdf"), width = 6, height = 5)
  }
  
  
  
  
  de.order <- de$de
  gsea <- function(de.order, outdir){
    de.order$ENSG <- rownames(de.order)
    de.order$ENSG <- as.character(lapply(strsplit(de.order$ENSG,"|",fixed=TRUE), function(x) x[1]))
    de.order$ENSG <- as.character(lapply(strsplit(de.order$ENSG,".",fixed=TRUE), function(x) x[1]))
    de.order <- de.order[!duplicated(de.order$ENSG),]
    ##### get all geneList (for GSEA)
    {
      #gene.list.ncbiid <- bitr(de.order$name, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
      #gene.list.ncbiid <- read.table("/BioII/lulab_b/zhanqing/cfRNAseq/reference/geneset/gencode.v38.txt",sep = "\t",header = T)
      gene.list.ncbiid <- read.table("/Users/wangge/Documents/DM/reference/gencode.v38.txt",sep = "\t",header = T)
      gene.list.ncbiid <- gene.list.ncbiid[!is.na(gene.list.ncbiid$ENTREZID),]
      gene.list.ncbiid <- gene.list.ncbiid[,c("ensg","ENTREZID")]
      gene.list.ncbiid <- gene.list.ncbiid[!duplicated(gene.list.ncbiid$ensg) & !duplicated(gene.list.ncbiid$ENTREZID),]
      colnames(gene.list.ncbiid)[1] <- "ENSG"
      
      #de.order <- merge(x = de.order, y = gene.list.ncbiid, by.x="name", by.y="ENSG")
      de.order <- dplyr::left_join(x = de.order, y = gene.list.ncbiid)
      de.order <- de.order[order(de.order$logFC, -de.order$PValue, decreasing = T),] # order by logFC then PValue
      
      gene.list.ncbiid <- de.order$logFC[!duplicated(de.order$ENTREZID)]
      names(gene.list.ncbiid) <- de.order$ENTREZID[!duplicated(de.order$ENTREZID)]
      message(paste0("remaining ENTREZID geneList lenth: ",length(gene.list.ncbiid)))
      gene.list.ensgid <- de.order$logFC
      names(gene.list.ensgid) <- de.order$ENSG
      message(paste0("remaining ENSEMBL geneList lenth: ",length(gene.list.ensgid)))
    }
    # 2.Gene Set Enrichment Analysis online
    message("start gsea...")
    {
      gse.GO <- gseGO(
        geneList = gene.list.ensgid, 
        OrgDb = org.Hs.eg.db,
        ont = "ALL", 
        keyType = "ENSEMBL",  
        pvalueCutoff = pValue,
        pAdjustMethod = pAdjust
        #by = "fgsea", #'DOSE'
      ) 
      if(FALSE){
        gse.GO <- gseGO(
          geneList = gene.list.ncbiid, 
          OrgDb = org.Hs.eg.db,
          ont = "ALL", 
          keyType = "ENTREZID",  
          pvalueCutoff = pValue,
          pAdjustMethod = pAdjust
          #by = "fgsea", #'DOSE'
        ) 
      }
      gse.KEGG <- gseKEGG(
        geneList = gene.list.ncbiid,
        keyType = 'kegg', 
        organism = 'hsa',
        nPerm  = 1000,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = pValue,
        pAdjustMethod = pAdjust
        #by = "fgsea",
        #seed = T  
      )
      gse.DO <- gseDO(
        geneList = gene.list.ncbiid,
        nPerm  = 1000,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = pValue,
        pAdjustMethod = pAdjust
        #by = "fgsea",
        #seed = T
      )
    }
    # output table
    {
      types <- c(gse.GO, gse.KEGG, gse.DO)
      names(types) <-c("gse.GO","gse.KEGG","gse.DO")
      dir.create(paste0(outdir), showWarnings = T)
      message("start writing tables...")
      for(i in 1:length(types)){
        write.table(types[[i]],paste0(outdir,"/", names(types)[i], ".txt"),row.names = F,col.names = T, quote = F,sep = "\t")}
    }
  }
  
  gsea(de.all$de, "cfRNA/enrich_all-sample-table")
  gsea(de.nc$de, "result/paired/enrich/nc-sample-table")
  gsea(de.crc$de, "result/paired/enrich/crc-sample-table")
  gsea(de.luad$de, "result/paired/enrich/luad-sample-table")

  
  # GSEA enrichment plot 
  {
    options(scipen=1000000) 
    enrich
    tabledir="result/paired/enrich/all-sample-table/"
    plotdir="result/paired/enrich/all-sample-plot/"
    dir.create(plotdir, showWarnings = T)
    col = brewer.pal(3,"Set1")
    
    ### enrichment::dotplot
    {
      outdir="result/enrich/all-sample-plot/"
      col = brewer.pal(8,"Set1")
      # 气泡图，展示geneset被激活还是抑制
      pdf(paste0(outdir,"dotplot-gse.GO",".pdf"), width = 8, height = 5.8)
      enrichplot::dotplot(gse.GO, split=".sign", font.size=10, title="gsea.GO", showCategory = 10) +
        facet_grid(~.sign) + 
        scale_color_continuous(low=col[1],high=col[2]) +
        scale_x_continuous(limits = c(0.18,0.85),breaks = c(0.2,0.4,0.6,0.8))
      dev.off()
      
      pdf(paste0(outdir,"dotplot-gse.KEGG",".pdf"), width = 6, height = 4)
      enrichplot::dotplot(gse.KEGG, split=".sign", font.size=10, title="gsea.KEGG", showCategory = 10) +
        facet_grid(~.sign) + 
        scale_color_continuous(low=col[1],high=col[2]) 
      scale_x_continuous(limits = c(0.25,0.55),breaks = c(0.3,0.4,0.5))
      dev.off()
    }
    ### GSEA.GO
    {
      all.gse.go <- read.table(paste0(tabledir,"gse.GO.txt"),sep="\t",header=T,quote="")
      all.gse.go <- all.gse.go %>% mutate(logp = -log10(p.adjust),
                                          Count = str_count(core_enrichment,"/")+1)
      all.gse.go <- all.gse.go[order(all.gse.go$p.adjust, -abs(all.gse.go$NES), -all.gse.go$Count),]
      
      # up
      up <- filter(all.gse.go, NES>1, p.adjust<0.1)

      up.BP <- filter(up, ONTOLOGY=="BP")
      up.BP <- up.BP[1:10,]
      up.BP$Description <- factor(up.BP$Description, levels=up.BP$Description[order(-up.BP$p.adjust)])
      up.CC <- filter(up, ONTOLOGY=="CC")
      up.CC <- up.CC[1:10,]
      up.CC$Description <- factor(up.CC$Description, levels=up.CC$Description[order(-up.CC$p.adjust)])
      up.MF <- filter(up, ONTOLOGY=="MF")
      up.MF <- up.MF[1:10,]
      up.MF$Description <- factor(up.MF$Description, levels=up.MF$Description[order(-up.MF$p.adjust)])
      
      # down
      down <- filter(all.gse.go, NES < -1, p.adjust<0.1)
      
      down.BP <- filter(down, ONTOLOGY=="BP")
      down.BP <- down.BP[1:10,]
      down.BP$Description <- factor(down.BP$Description, levels=down.BP$Description[order(-down.BP$p.adjust)])
      down.CC <- filter(down, ONTOLOGY=="CC")
      down.CC <- down.CC[1:10,]
      down.CC$Description <- factor(down.CC$Description, levels=down.CC$Description[order(-down.CC$p.adjust)])
      down.MF <- filter(down, ONTOLOGY=="MF")
      down.MF$Description <- factor(down.MF$Description, levels=down.MF$Description[order(-down.MF$p.adjust)])
      
      # union
      BP <- rbind(up.BP %>% mutate(group = "EV"), down.BP %>% mutate(group = "CF"))
      CC <- rbind(up.CC %>% mutate(group = "EV"), down.CC %>% mutate(group = "CF"))
      MF <- rbind(up.MF %>% mutate(group = "EV"), down.MF %>% mutate(group = "CF"))
    }
    ### GSEA.KEGG
    {
      all.gse.kegg <- read.table(paste0(tabledir,"gse.KEGG.txt"),sep="\t",header=T,quote="")
      all.gse.kegg <- all.gse.kegg %>% mutate(logp = -log10(p.adjust),
                                              Count = str_count(core_enrichment,"/")+1)
      all.gse.kegg <- all.gse.kegg[order(all.gse.kegg$p.adjust, -abs(all.gse.kegg$NES), -all.gse.kegg$Count),]
      # up
      up <- all.gse.kegg[all.gse.kegg$NES>1 & all.gse.kegg$p.adjust < 0.1,]
      up <- up[1:10,]
      up$Description <- factor(up$Description, levels=rev(up$Description[order(up$p.adjust)]))
      # down
      down <- all.gse.kegg[all.gse.kegg$NES < -1 & all.gse.kegg$p.adjust < 0.1,]
      #down <- down[1:10,]
      down$Description <- factor(down$Description, levels=rev(down$Description[order(down$p.adjust)]))
      # union
      kegg <- rbind(up %>% mutate(group = "EV"), down %>% mutate(group = "CF"))
    }
    ### ggplot
    {
      p=ggplot(kegg, aes(group, Description, color=logp)) +
        #geom_point(aes(color = logp),size=5) +
        geom_point(aes(color = logp, size = Count)) +
        scale_color_gradient2(name="log10(FDR)", limits=c(0.5,2), breaks=c(1,2), 
                             #mid="#87CEEB", high="#333333")+  # shyblue
                             high="#104E8B")+
        theme_bw()+
        guides(fill=guide_legend(title=NULL, ncol=1))+
        coord_fixed(ratio=0.7)+
        theme(
          plot.margin = unit(c(1, 1, 1, 1),"cm"),
          panel.border=element_rect(size=1, fill="transparent"),
          panel.grid=element_blank(),
          #legend.position=c(0.25,0.85),
          legend.position="right",
          legend.background = element_blank(),
          legend.title = element_text(color="black", size=8),
          legend.text= element_text(color="black", size=8),
          legend.key.size = unit(0.3,"cm"),
          axis.text.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=10),
          axis.title.x = element_text(color="black", size=10),
          axis.title.y = element_blank())+
        labs(title="")
      p
    }
    ggsave(paste0(plotdir,"dotplot-gse.goBP",".pdf"), p, width = 10, height = 5)
    ggsave(paste0(plotdir,"dotplot-gse.kegg",".pdf"), p, width = 10, height = 3.3)
    ### GSEA.GO.pickup
    {
      col2 = pal_nejm()(8)
      all.gse.go <- read.table(paste0(tabledir,"gse.GO.txt"),sep="\t",header=T,quote="")
      all.gse.go <- all.gse.go %>% mutate(logp = -log10(p.adjust),
                                          Count = str_count(core_enrichment,"/")+1)
      all.gse.go <- all.gse.go[order(all.gse.go$p.adjust, -abs(all.gse.go$NES), -all.gse.go$Count),]
      # up
      up <- filter(all.gse.go, NES>1, p.adjust<0.1)
      select <- c(#transcriprion activity
        "DNA-binding transcription factor activity", 
        "regulatory region nucleic acid binding",   
        #immune response
        "lymphocyte activation",
        "defense response to other organism",
        #"immune response-regulating signaling pathway", 
        "Fc receptor signaling pathway",
        "antigen receptor-mediated signaling pathway",
        #"activation of innate immune response",
        "T cell receptor signaling pathway",
        #"nuclear speck",
        #"ribonucleoprotein complex",
        #"spliceosomal complex",
        
        #cell migration related —— only LC???????
        #"anchoring junction",
        "focal adhesion",
        "cell-substrate junction")
      up <- up[up$Description %in% select,]
      up$Description <- factor(up$Description, levels=rev(select))
      {
        p=ggplot(up, aes(x=Description,y=NES,fill=logp))+
          geom_bar(stat="identity",width=0.7,size=0.2,color="black")+
          #geom_text(aes(label=signif),size=10,vjust=0.8,hjust=1)+
          #geom_text(aes(label=signif),size=8,vjust=0.8,hjust=0)+
          coord_flip()+
          #geom_point(aes(color=p.adjust,size = Count)) +
          scale_fill_gradient2(name="log10(FDR)",breaks=c(2,4,6,8), high=col2[1], guide="colorbar")+
          #scale_y_continuous(limits = c(0,220), breaks=c(0,150))+
          scale_y_continuous(limits = c(0,2.5), breaks=c(0,1,2))+
          theme_bw()+
          theme(
            plot.margin = unit(c(1, 1, 1, 1),"cm"),
            plot.title = element_text(color="black", size=15, hjust = 0.5),
            panel.border = element_rect(size=1, color="black"),
            panel.grid = element_blank(),
            strip.background = element_rect(size=1, color="black",fill="transparent"),
            strip.text = element_text(color="black", size=12),
            legend.position="right",
            legend.background = element_blank(),
            legend.title = element_text(color="black", size=8),
            legend.text= element_text(color="black", size=8),
            legend.key.size = unit(0.2,"cm"),
            axis.text.x = element_text(color="black", size=10),
            axis.text.y = element_text(color="black", size=10, margin=ggplot2::margin(1,0.5,0,0,"cm")),
            axis.title.x = element_text(color="black", size=10),
            axis.title.y = element_blank())+
          labs(title="EV-enriched",y="NES")
        p
      }
      # down
      down <- filter(all.gse.go, NES < -1, p.adjust<0.1)
      down$logp <- as.numeric(down$logp)
      down <- down[order(down$p.adjust, abs(down$NES), down$Count),]
      select <- c(#splicing and RNP 
        "mRNA 5'-splice site recognition",
        "U1 snRNP",
        "formation of quadruple SL/U4/U5/U6 snRNP",
        "spliceosomal snRNP complex", 
        "Sm-like protein family complex", 
        #immune response
        "antimicrobial humoral response",   # NC signif
        "innate immune response in mucosa",
        "organ or tissue specific immune response",
        "defense response to Gram-negative bacterium"
        #others
        )
      down <- down[down$Description %in% select,]
      down$Description <- factor(down$Description, levels=rev(select))
      {
        q=ggplot(down, aes(x=Description,y=abs(NES),fill=logp))+
          geom_bar(stat="identity",width=0.7,size=0.2,color="black")+
          coord_flip()+
          #geom_point(aes(color=p.adjust,size = Count)) +
          #scale_fill_gradient2(name="log10(FDR)",breaks=c(1,2,3), high=muted("dodgerblue4"), guide="colorbar")+
          scale_fill_gradient2(name="log10(FDR)",breaks=c(1,2,3), high=col2[2], guide="colorbar")+
          #scale_y_continuous(limits=c(0,23),breaks = c(0,20))+
          scale_y_continuous(limits = c(0,2.5), breaks=c(0,1,2))+
          theme_bw()+
          theme(
            plot.margin = unit(c(1, 1, 1, 1),"cm"),
            plot.title = element_text(color="black", size=15, hjust = 0.5),
            panel.border = element_rect(size=1, color="black"),
            panel.grid = element_blank(),
            strip.background = element_rect(size=1, color="black",fill="transparent"),
            strip.text = element_text(color="black", size=12),
            legend.position="right",
            legend.background = element_blank(),
            legend.title = element_text(color="black", size=8),
            legend.text= element_text(color="black", size=8),
            legend.key.size = unit(0.2,"cm"),
            axis.text.x = element_text(color="black", size=10),
            axis.text.y = element_text(color="black", size=10, margin=ggplot2::margin(1,0.5,0,0,"cm")),
            axis.title.x = element_text(color="black", size=10),
            axis.title.y = element_blank())+
          labs(title="CF-enriched",y="-NES")
        q
      }
      a=ggarrange(p,q,ncol=1, heights = c(1, 1),align="v")
      ggsave("result/paired/enrich/all-sample-plot/barplot-gsea.GO.EV-CF.pdf",a,width=5.7,height=6.5)
      
    }

    ### split cancer-GO
    {
      nc.gse.go <- read.table("result/paired/enrich/nc-sample-table/gse.GO.txt",sep="\t",header=T,quote="")
      nc.gse.go <- nc.gse.go %>% mutate(logp = -log10(p.adjust), Count = str_count(core_enrichment,"/")+1)
      nc.gse.go <- nc.gse.go[order(nc.gse.go$p.adjust, -abs(nc.gse.go$NES), -nc.gse.go$Count),]
      nc.gse.go.up <- nc.gse.go[nc.gse.go$NES>1,]
      nc.gse.go.down <- nc.gse.go[nc.gse.go$NES < -1,]
      
      crc.gse.go <- read.table("result/paired/enrich/crc-sample-table/gse.GO.txt",sep="\t",header=T,quote="")
      crc.gse.go <- crc.gse.go %>% mutate(logp = -log10(p.adjust), Count = str_count(core_enrichment,"/")+1)
      crc.gse.go <- crc.gse.go[order(crc.gse.go$p.adjust, -abs(crc.gse.go$NES), -crc.gse.go$Count),]
      crc.gse.go.up <- crc.gse.go[crc.gse.go$NES>1,]
      crc.gse.go.down <- crc.gse.go[crc.gse.go$NES < -1,]
      
      luad.gse.go <- read.table("result/paired/enrich/luad-sample-table/gse.GO.txt",sep="\t",header=T,quote="")
      luad.gse.go <- luad.gse.go %>% mutate(logp = -log10(p.adjust), Count = str_count(core_enrichment,"/")+1)
      luad.gse.go <- luad.gse.go[order(luad.gse.go$p.adjust, -abs(luad.gse.go$NES), -luad.gse.go$Count),]
      luad.gse.go.up <- luad.gse.go[luad.gse.go$NES>1,]
      luad.gse.go.down <- luad.gse.go[luad.gse.go$NES < -1,]
      
      up <- rbind(nc.gse.go.up %>% mutate(group="NC"),
                  crc.gse.go.up %>% mutate(group="CRC"),
                  luad.gse.go.up %>% mutate(group="LC"))
      up$logp <- as.numeric(up$logp)
      up <- up[order(up$p.adjust, abs(up$NES), up$Count),]
      select <- c(#MF
                  "DNA-binding transcription factor activity", 
                  "regulatory region nucleic acid binding",   
                  #BP
                  "defense response to other organism",
                  "regulation of immune response", 
                  "lymphocyte activation",
                  "B cell activation",
                  "immune response-regulating signaling pathway", 
                  "antigen receptor-mediated signaling pathway",
                  "activation of innate immune response",
                  "Fc receptor signaling pathway",
                  "T cell receptor signaling pathway"
                  #"nuclear transport"
                  
                  #CC
                  #"nuclear speck",
                  #"ribonucleoprotein complex",
                  #"spliceosomal complex",
                  #"anchoring junction",
                  #"cell-substrate junction",
                  #"focal adhesion"
                  )
      select <- c(#transcriprion activity
                  "DNA-binding transcription factor activity", 
                  "regulatory region nucleic acid binding",   
                  #immune response
                  "defense response to other organism",
                  "lymphocyte activation",
                  "immune response-regulating signaling pathway", 
                  "antigen receptor-mediated signaling pathway",
                  #"activation of innate immune response",
                  "Fc receptor signaling pathway",
                  "T cell receptor signaling pathway",
                  #"nuclear speck",
                  #"ribonucleoprotein complex",
                  #"spliceosomal complex",
                  
                  #cell migration related —— only LC???????
                  #"anchoring junction",
                  "focal adhesion",
                  "cell-substrate junction"
      )
      up <- up[up$Description %in% select,]
      up$Description <- factor(up$Description, levels=rev(select))
      up$group <- factor(up$group, levels=c("NC","CRC","LC"))
      up$signif <- ifelse(up$p.adjust<0.1, "*","")
      {
        q=ggplot(up, aes(x=Description,y=Count,fill=logp))+
          facet_grid(~group)+
          geom_bar(stat="identity",width=0.7,size=0.2,color="black")+
          #geom_text(aes(label=signif),size=10,vjust=0.8,hjust=1)+
          geom_text(aes(label=signif),size=8,vjust=0.8,hjust=0)+
          coord_flip()+
          #geom_point(aes(color=p.adjust,size = Count)) +
          #scale_fill_gradient2(name="log10(FDR)",breaks=c(2,4,6,8),low="dodgerblue3",high="dodgerblue4")+
          scale_fill_gradient2(name="log10(FDR)",breaks=c(2,4,6,8), high=muted("dodgerblue4"), guide="colorbar")+
          scale_y_continuous(limits = c(0,220), breaks=c(0,150))+
          theme_bw()+
          theme(
            plot.margin = unit(c(1, 1, 1, 1),"cm"),
            plot.title = element_text(color="black", size=15, hjust = 0.5),
            panel.border = element_rect(size=1, color="black"),
            panel.grid = element_blank(),
            strip.background = element_rect(size=1, color="black",fill="transparent"),
            strip.text = element_text(color="black", size=12),
            legend.position="right",
            legend.background = element_blank(),
            legend.title = element_text(color="black", size=8),
            legend.text= element_text(color="black", size=8),
            legend.key.size = unit(0.2,"cm"),
            axis.text.x = element_text(color="black", size=10),
            axis.text.y = element_text(color="black", size=10, margin=ggplot2::margin(1,0.5,0,0,"cm")),
            axis.title.x = element_text(color="black", size=10),
            axis.title.y = element_blank())+
          labs(title="EV-enriched",y="Gene count")
        q
      }
      down <- rbind(nc.gse.go.down %>% mutate(group="NC"),
                    crc.gse.go.down %>% mutate(group="CRC"),
                    luad.gse.go.down %>% mutate(group="LC"))
      down$logp <- as.numeric(down$logp)
      down <- down[order(down$p.adjust, abs(down$NES), down$Count),]
      select <- c(
                  #CC
                  "integral component of plasma membrane", # NC signif
                  "cytosolic ribosome",
                  "Sm-like protein family complex", 
                  "spliceosomal snRNP complex", 
                  "U1 snRNP",
                  "cytosolic large ribosomal subunit",
                  #BP
                  "antimicrobial humoral response",   # NC signif
                  "mRNA 5'-splice site recognition",
                  "organ or tissue specific immune response",
                  #"platelet activation",
                  "negative regulation of execution phase of apoptosis",
                  "formation of quadruple SL/U4/U5/U6 snRNP"
                  )
      select <- c(#splicing and RNP 
                  "Sm-like protein family complex", 
                  "spliceosomal snRNP complex", 
                  "mRNA 5'-splice site recognition",
                  "U1 snRNP",
                  "formation of quadruple SL/U4/U5/U6 snRNP",
                  #immune response
                  "antimicrobial humoral response",   # NC signif
                  "organ or tissue specific immune response",
                  #others
                  "negative regulation of execution phase of apoptosis",
                  "integral component of plasma membrane", # NC signif
                  "cytosolic ribosome"
      )
      down <- down[down$Description %in% select,]
      down$Description <- factor(down$Description, levels=rev(select))
      down$group <- factor(down$group, levels=c("NC","CRC","LC"))
      down$signif <- ifelse(down$p.adjust<0.1, "*","")
      {
        q=ggplot(down, aes(x=Description,y=Count,fill=logp))+
          facet_grid(~group)+
          geom_bar(stat="identity",width=0.7,size=0.2,color="black")+
          geom_text(aes(label=signif),size=8,vjust=0.8,hjust=0)+
          coord_flip()+
          #geom_point(aes(color=p.adjust,size = Count)) +
          #scale_fill_gradient2(name="log10(FDR)",breaks=c(2,4,6,8),low="dodgerblue3",high="dodgerblue4")+
          scale_fill_gradient2(name="log10(FDR)",breaks=c(1,2,3), high=muted("dodgerblue4"), guide="colorbar")+
          scale_y_continuous(limits=c(0,70),breaks = c(0,50))+
          theme_bw()+
          theme(
            plot.margin = unit(c(1, 1, 1, 1),"cm"),
            plot.title = element_text(color="black", size=15, hjust = 0.5),
            panel.border = element_rect(size=1, color="black"),
            panel.grid = element_blank(),
            strip.background = element_rect(size=1, color="black",fill="transparent"),
            strip.text = element_text(color="black", size=12),
            legend.position="right",
            legend.background = element_blank(),
            legend.title = element_text(color="black", size=8),
            legend.text= element_text(color="black", size=8),
            legend.key.size = unit(0.2,"cm"),
            axis.text.x = element_text(color="black", size=10),
            axis.text.y = element_text(color="black", size=10, margin=ggplot2::margin(1,0.5,0,0,"cm")),
            axis.title.x = element_text(color="black", size=10),
            axis.title.y = element_blank())+
          labs(title="CF-enriched",y="Gene count")
        q
      }
      #ggsave("result/paired/enrich/nc.crc.lc.gsea.CF.pdf",p,width=6.4,height=3.3)
      a=ggarrange(q,p,ncol=1, heights = c(1, 1),align="v")
      ggsave("result/paired/enrich/nc.crc.lc.gsea.GO.EV-CF-3.pdf",a,width=6.9,height=6.5)
      
      
    }
    ### split cancer-KEGG
    {
      nc.gse.kegg <- read.table("result/paired/enrich/nc-sample-table/gse.KEGG.txt",sep="\t",header=T,quote="")
      nc.gse.kegg <- nc.gse.kegg %>% mutate(logp = -log10(p.adjust), Count = str_count(core_enrichment,"/")+1)
      nc.gse.kegg <- nc.gse.kegg[order(nc.gse.kegg$p.adjust, -abs(nc.gse.kegg$NES), -nc.gse.kegg$Count),]
      nc.gse.kegg.up <- nc.gse.kegg[nc.gse.kegg$NES>1,]
      nc.gse.kegg.down <- nc.gse.kegg[nc.gse.kegg$NES < -1,]
      
      crc.gse.kegg <- read.table("result/paired/enrich/crc-sample-table/gse.KEGG.txt",sep="\t",header=T,quote="")
      crc.gse.kegg <- crc.gse.kegg %>% mutate(logp = -log10(p.adjust), Count = str_count(core_enrichment,"/")+1)
      crc.gse.kegg <- crc.gse.kegg[order(crc.gse.kegg$p.adjust, -abs(crc.gse.kegg$NES), -crc.gse.kegg$Count),]
      crc.gse.kegg.up <- crc.gse.kegg[crc.gse.kegg$NES>1,]
      crc.gse.kegg.down <- crc.gse.kegg[crc.gse.kegg$NES < -1,]
      
      luad.gse.kegg <- read.table("result/paired/enrich/luad-sample-table/gse.KEGG.txt",sep="\t",header=T,quote="")
      luad.gse.kegg <- luad.gse.kegg %>% mutate(logp = -log10(p.adjust), Count = str_count(core_enrichment,"/")+1)
      luad.gse.kegg <- luad.gse.kegg[order(luad.gse.kegg$p.adjust, -abs(luad.gse.kegg$NES), -luad.gse.kegg$Count),]
      luad.gse.kegg.up <- luad.gse.kegg[luad.gse.kegg$NES>1,]
      luad.gse.kegg.down <- luad.gse.kegg[luad.gse.kegg$NES < -1,]
      
      up <- rbind(nc.gse.kegg.up %>% mutate(group="NC"),
                  crc.gse.kegg.up %>% mutate(group="CRC"),
                  luad.gse.kegg.up %>% mutate(group="LC"))
      
      up$logp <- as.numeric(up$logp)
      up <- up[order(up$p.adjust, abs(up$NES), up$Count),]
      up.cc=filter(up,group=="CRC")
      
      select=c("Th1 and Th2 cell differentiation",
               "Th17 cell differentiation",
               "Hedgehog signaling pathway",
               "Systemic lupus erythematosus",
               "B cell receptor signaling pathway",
               "Proteasome",
               "NF-kappa B signaling pathway")
      up <- up[up$Description %in% select,]
      up$Description <- factor(up$Description, levels=rev(select))
      up$group <- factor(up$group, levels=c("NC","CRC","LC"))
      up$signif <- ifelse(up$p.adjust<0.1, "*","")
      {
        q=ggplot(up, aes(x=Description,y=Count,fill=logp))+
          facet_grid(~group)+
          geom_bar(stat="identity",width=0.7,size=0.2,color="black")+
          geom_text(aes(label=signif),size=8,vjust=0.8,hjust=0)+
          coord_flip()+
          scale_fill_gradient2(limits=c(0,1),name="log10(FDR)",breaks=c(0,1), high=muted("dodgerblue4"), guide="colorbar")+
          scale_y_continuous(limits = c(0,30), breaks=c(0,20))+
          theme_bw()+
          theme(
            plot.margin = unit(c(1, 1, 1, 1),"cm"),
            plot.title = element_text(color="black", size=15, hjust = 0.5),
            panel.border = element_rect(size=1, color="black"),
            panel.grid = element_blank(),
            strip.background = element_rect(size=1, color="black",fill="transparent"),
            strip.text = element_text(color="black", size=12),
            legend.position="right",
            legend.background = element_blank(),
            legend.title = element_text(color="black", size=8),
            legend.text= element_text(color="black", size=8),
            legend.key.size = unit(0.2,"cm"),
            axis.text.x = element_text(color="black", size=10),
            axis.text.y = element_text(color="black", size=10),
            axis.title.x = element_text(color="black", size=10),
            axis.title.y = element_blank())+
          labs(title="",y="Gene count")
        q
      }
      down <- rbind(nc.gse.kegg.down %>% mutate(group="NC"),
                    crc.gse.kegg.down %>% mutate(group="CRC"),
                    luad.gse.kegg.down %>% mutate(group="LC"))
      down$logp <- as.numeric(down$logp)
      down <- down[order(-down$p.adjust, abs(down$NES), down$Count),]
      down.cc <- filter(down, group=="LC")
      select=c("Systemic lupus erythematosus",
               "Alcoholism",
               "Retrograde endocannabinoid signaling",
               "N-Glycan biosynthesis",
               "Cardiac muscle contraction")
      down <- down[down$Description %in% select,]
      down$Description <- factor(down$Description, levels=rev(select))
      down$group <- factor(down$group, levels=c("NC","CRC","LC"))
      down$signif <- ifelse(down$p.adjust<0.1, "*","")
      {
        p=ggplot(down, aes(x=Description,y=Count,fill=logp))+
          facet_grid(~group)+
          geom_bar(stat="identity",width=0.7,size=0.2,color="black")+
          geom_text(aes(label=signif),size=8,vjust=0.8,hjust=0)+
          coord_flip()+
          scale_fill_gradient2(limits=c(0,1),name="log10(FDR)",breaks=c(0,1), high=muted("dodgerblue4"), guide="colorbar")+
          scale_y_continuous(limits=c(0,50),breaks = c(0,40))+
          theme_bw()+
          theme(
            plot.margin = unit(c(1, 1, 1, 1),"cm"),
            plot.title = element_text(color="black", size=15, hjust = 0.5),
            panel.border = element_rect(size=1, color="black"),
            panel.grid = element_blank(),
            strip.background = element_rect(size=1, color="black",fill="transparent"),
            strip.text = element_text(color="black", size=12),
            legend.position="right",
            legend.background = element_blank(),
            legend.title = element_text(color="black", size=8),
            legend.text= element_text(color="black", size=8),
            legend.key.size = unit(0.2,"cm"),
            axis.text.x = element_text(color="black", size=10),
            axis.text.y = element_text(color="black", size=10),
            axis.title.x = element_text(color="black", size=10),
            axis.title.y = element_blank())+
          labs(title="",y="Gene count")
        p
      }
      a=ggarrange(q,p,ncol=1, heights = c(1, 0.85),align="v")
      ggsave("result/paired/enrich/nc.crc.lc.gsea.KEGG.EV-CF.pdf",a,width=5.5,height=5.5)
      
    }
  }
  
  # GSEA plot
  {
    gse.GO
    plotdir="result/paired/enrich/nc-sample-plot/"
    # GSEA plot
    col=brewer.pal(8,"Set1")
    col1 = rev(brewer.pal(8,"Reds"))
    col2 = rev(brewer.pal(8,"Blues"))
    
    ###
    df = as.data.frame(gse.GO)
    df = filter(df, p.adjust<0.1)
    pdf(paste0("result/paired/enrich/all-sample-plot/gseaplot.GO.immune-pickup-2",".pdf"), width = 4.5, height = 3.5)
    id=c("GO:0098542", #defense response to other organism
         #"GO:0050851", #antigen receptor-mediated signaling pathway
         #"GO:0002764", #immune response-regulating signaling pathway
         "GO:0050852", #T cell receptor signaling pathway
         "GO:0038093", #Fc receptor signaling pathway
         #"GO:0002218", #activation of innate immune response
         "GO:0019730", #antimicrobial humoral response
         "GO:0002251" #organ or tissue specific immune response
    )
    gseaplot2(gse.GO, geneSetID = id, color = c(col2[1],col1[1],col1[3],col2[3],col1[6]), pvalue_table = F, 
              base_size = 12, rel_heights = c(1.5, 0.4, 0.6), subplots = 1:3)
    dev.off()
    ###### change plot function ######
    trace("gseaplot2", edit = TRUE) 
    # 26line:  theme(legend.position = "none") + scale_y_continuous(limits = c(-0.8, 0.5), breaks=c(-0.8,-0.4,0,0.4))
    untrace("gseaplot2")
    
    # KEGG
    pdf(paste0(plotdir,"gseaplot.KEGG",".pdf"), width = 8, height = 6)
    gseaplot2(gse.KEGG, geneSetID = 1:5, color = col[1:5], pvalue_table = F, 
              base_size = 10, rel_heights = c(1.5, 0.5, 0.8), subplots = 1:3)
    dev.off()
    # DO
    pdf(paste0(plotdir,"gseaplot.DO",".pdf"), width = 6, height = 6)
    gseaplot2(gse.DO, geneSetID = 1:5, color = col[1:5], pvalue_table = T, 
              base_size = 10, rel_heights = c(1.5, 0.5, 0.8), subplots = 1:3)
    dev.off()
  }
  
  ##### get gene list
  {
    
    #### GO term selection
    library(DOSE)
    library(GOSemSim)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    get_GO_data <- function(OrgDb, ont, keytype){
      GO_Env <- get_GO_Env()
      use_cached <- FALSE
      
      if (exists("organism", envir=GO_Env, inherits=FALSE) &&
          exists("keytype", envir=GO_Env, inherits=FALSE)) {
        
        org <- get("organism", envir=GO_Env)
        kt <- get("keytype", envir=GO_Env)
        
        if (org == DOSE:::get_organism(OrgDb) &&
            keytype == kt &&
            exists("goAnno", envir=GO_Env, inherits=FALSE)) {
          ## https://github.com/GuangchuangYu/clusterProfiler/issues/182
          ## && exists("GO2TERM", envir=GO_Env, inherits=FALSE)){
          
          use_cached <- TRUE
        }
      }
      
      if (use_cached) {
        goAnno <- get("goAnno", envir=GO_Env)
      } else {
        OrgDb <- GOSemSim:::load_OrgDb(OrgDb)
        kt <- keytypes(OrgDb)
        if (! keytype %in% kt) {
          stop("keytype is not supported...")
        }
        
        kk <- keys(OrgDb, keytype=keytype)
        goAnno <- suppressMessages(
          select(OrgDb, keys=kk, keytype=keytype,
                 columns=c("GOALL", "ONTOLOGYALL")))
        
        goAnno <- unique(goAnno[!is.na(goAnno$GOALL), ])
        
        assign("goAnno", goAnno, envir=GO_Env)
        assign("keytype", keytype, envir=GO_Env)
        assign("organism", DOSE:::get_organism(OrgDb), envir=GO_Env)
      }
      
      if (ont == "ALL") {
        GO2GENE <- unique(goAnno[, c(2,1)])
      } else {
        GO2GENE <- unique(goAnno[goAnno$ONTOLOGYALL == ont, c(2,1)])
      }
      
      GO_DATA <- DOSE:::build_Anno(GO2GENE, get_GO2TERM_table())
      
      goOnt.df <- goAnno[, c("GOALL", "ONTOLOGYALL")] %>% unique
      goOnt <- goOnt.df[,2]
      names(goOnt) <- goOnt.df[,1]
      assign("GO2ONT", goOnt, envir=GO_DATA)
      return(GO_DATA)
    }
    get_GO_Env <- function(){
      if (!exists(".GO_clusterProfiler_Env", envir = .GlobalEnv)) {
        pos <- 1
        envir <- as.environment(pos)
        assign(".GO_clusterProfiler_Env", new.env(), envir=envir)
      }
      get(".GO_clusterProfiler_Env", envir = .GlobalEnv)
    }
    get_GO2TERM_table <- function(){
      GOTERM.df <- get_GOTERM()
      GOTERM.df[, c("go_id", "Term")] %>% unique
    }
    get_GOTERM <- function(){
      pos <- 1
      envir <- as.environment(pos)
      if (!exists(".GOTERM_Env", envir=envir)) {
        assign(".GOTERM_Env", new.env(), envir)
      }
      GOTERM_Env <- get(".GOTERM_Env", envir = envir)
      if (exists("GOTERM.df", envir = GOTERM_Env)) {
        GOTERM.df <- get("GOTERM.df", envir=GOTERM_Env)
      } else {
        GOTERM.df <- toTable(GOTERM)
        assign("GOTERM.df", GOTERM.df, envir = GOTERM_Env)
      }
      return(GOTERM.df)
    }
    go_data <- get_GO_data("org.Hs.eg.db", "ALL", "SYMBOL") 　
    save(go_data, file="GO_data.RData")
    load("GO_data.RData")
    go_data$EXTID2PATHID[1]
    go_data$GO2ONT[1]
    go_data$PATHID2EXTID[1]
    go_data$PATHID2NAME[2]
    
    id <- unique(unlist(go_data$PATHID2EXTID["GO:0050852"]))
    id <- unique(unlist(go_data$PATHID2EXTID["GO:0038093"]))
    id.ensg <- bitr(id, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                        toType = "ENSEMBL", #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                        OrgDb = org.Hs.eg.db) #Orgdb是指对应的注释包是哪个
    unique(id.ensg$ENSEMBL)
    paste(unique(id.ensg$ENSEMBL),collapse = "\t")
    up[up$ID=="GO:0050852",]$core_enrichment

    
    
  }
  
  # 2.2 Cancer gene set Enrichment Analysis
  message("start gsea at cancer gene set")
  {
    CGC_1 <- read.table("/BioII/lulab_b/zhanqing/cancer_gene/v3/id/CGC_1.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
    CGC_2 <- read.table("/BioII/lulab_b/zhanqing/cancer_gene/v3/id/CGC_2.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
    CGC_all <- read.table("/BioII/lulab_b/zhanqing/cancer_gene/v3/id/CGC_all.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
    NCG_1 <- read.table("/BioII/lulab_b/zhanqing/cancer_gene/v3/id/NCG_1.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
    NCG_2 <- read.table("/BioII/lulab_b/zhanqing/cancer_gene/v3/id/NCG_2.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
    Uniprot <- read.table("/BioII/lulab_b/zhanqing/cancer_gene/v3/id/Uniprot.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
    Pansoft <- read.table("/BioII/lulab_b/zhanqing/cancer_gene/v3/id/PanSoftware.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
    IntOGene <- read.table("/BioII/lulab_b/zhanqing/cancer_gene/v3/id/IntOGen.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
    
    
    set1 <- read.table("/BioII/lulab_b/zhanqing/cancer_gene/v2/public/ensemblid/set/Set1.txt")[,1]
    set2 <- read.table("/BioII/lulab_b/zhanqing/cancer_gene/v2/public/ensemblid/set/Set2.txt")[,1]
    set3 <- read.table("/BioII/lulab_b/zhanqing/cancer_gene/v2/public/ensemblid/set/Set3.txt")[,1]
    set4 <- read.table("/BioII/lulab_b/zhanqing/cancer_gene/v2/public/ensemblid/set/Set4.txt")[,1]
    set5 <- read.table("/BioII/lulab_b/zhanqing/cancer_gene/v2/public/ensemblid/set/Set5.txt")[,1]
    set6 <- read.table("/BioII/lulab_b/zhanqing/cancer_gene/v2/public/ensemblid/set/Set6.txt")[,1]
    set7 <- read.table("/BioII/lulab_b/zhanqing/cancer_gene/v2/public/ensemblid/set/Set7.txt")[,1]
    
    geneset <- rbind(data.frame(Description = "CGC_1", Entrez = CGC_1$Entrez),
                     #data.frame(Description = "CGC_all", Entrez = CGC_all$Entrez),
                     data.frame(Description = "NCG_1", Entrez = NCG_1$ENTREZID),
                     #data.frame(Description = "NCG_all", Entrez = NCG_2$ENTREZID),
                     data.frame(Description = "Uniprot", Entrez = Uniprot$ENTREZID),
                     data.frame(Description = "Pansoft", Entrez = Pansoft$ENTREZID),
                     data.frame(Description = "IntOGene", Entrez = IntOGene$ENTREZID)
    )
    geneset2 <- rbind(data.frame(Description = "set1", ensg = set1),
                      data.frame(Description = "set2", ensg = set2),
                      data.frame(Description = "set3", ensg = set3),
                      data.frame(Description = "set4", ensg = set4),
                      data.frame(Description = "set5", ensg = set5),
                      data.frame(Description = "set6", ensg = set6),
                      data.frame(Description = "set7", ensg = set7))
    
    pValue = 1
    pAdjust = "BH"
    gsea <- clusterProfiler::GSEA(gene.list.ncbiid,  # ensg for geneset2, ncbi for geneset1
                                  TERM2GENE = geneset,
                                  nPerm  = 1000,
                                  minGSSize = 10,
                                  maxGSSize = 500,
                                  pvalueCutoff = pValue,
                                  pAdjustMethod = pAdjust)
    
    outdir="result/enrich/plot/"
    col = brewer.pal(8,"Set1")
    pdf(paste0(outdir,"gseaplot.Cancerdatabase",".pdf"), width = 8, height = 7)
    gseaplot2(gsea, geneSetID = 1:5, color = col[1:5], pvalue_table = T, 
              base_size = 10, rel_heights = c(1.5, 0.4, 0.7), subplots = 1:3)
    dev.off()
    
  }
  # fgsea
  {
    rank = gene.list.ensgid
    geneset = list(CGC_1 = CGC_1$Ensg,
                   CGC_all = CGC_2$Ensg,
                   NCG_1 = NCG_1$ENSEMBL,
                   NCG_2 = NCG_2$ENSEMBL,
                   Uniprot = Uniprot$ENSEMBL,
                   Pansoft = Pansoft$ENSEMBL,
                   IntOGene = IntOGene$ENSEMBL)
    library(fgsea)
    fgsea <- fgsea(geneset, 
                   rank, 
                   nperm = 10,
                   minSize = 0, 
                   maxSize = 500)
    summary(fgsea)
    
    plotGseaTable(fgsea$pathway,
                  rank,
                  fgsea, 
                  gseaParam=0.5)
    plotEnrichment(fgsea$pathway,rank)
    
  }
  {
    # emapplot
    enrichplot::emapplot(gse.GO)
    # cnetplot
    cnetplot(gsea)
    # 山峦图，展示每个geneset的基因logFC分布
    ridgeplot(gsea)
    # 选择单个gene set
    gsea3 <- data.frame(gsea)
    gseaplot2(gsea, geneSetID = 1, title=gsea3$Description[1])
    
  }
  
}
### Heterogenity
{
  RPM <- round(edgeR::cpm(y),3)
  geneid <- rownames(RPM)
  split <- strsplit(geneid,"|",fixed=T)
  length <- unlist(lapply(split,function(split) split[2]))
  y$genes <- data.frame(length=as.numeric(length))
  RPKM <- round(rpkm(y),3)
  TPM <- round((RPKM/colSums(RPKM))*10^6,3)
  cpm <- as.data.frame(RPM)
  tpm <- as.data.frame(TPM)
  rpkm <- as.data.frame(RPKM)
  
  actb = data.frame(sample_id=colnames(cpm),
                    actb=as.numeric(cpm["ENSG00000075624.17|4405|ACTB|protein_coding",]),
                    tmsb4x=as.numeric(cpm["ENSG00000205542.11|1703|TMSB4X|protein_coding",]),
                    cancer=cancer,
                    class=class)
  actb <- actb[order(class,cancer),]
  actb$sample_id <- factor(actb$sample_id, levels=actb$sample_id)
  actb <- melt(actb, id.vars = c("sample_id","class","cancer"))
  actb <- actb %>% mutate(actb, group=paste0(cancer,"-",class,"-",variable))
  {
    p=ggplot(actb,aes(x=sample_id,y=value,color=variable))+
      #facet_grid(~class)+
      geom_point(size=1)+
      geom_line(aes(group=group))+
      scale_color_aaas()+
      #scale_y_continuous(limits=c(0,40000))+
      theme_bw()+
      theme(
        plot.margin = unit(c(1, 1, 1, 0.5),"cm"),
        panel.border=element_rect(size=1, fill="transparent"),
        panel.grid=element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text= element_text(color="black", size=10),
        plot.title = element_text(hjust = 0.5, size=10, face="bold"),
        axis.text.x = element_text(color="black", size=7, angle=90),
        axis.text.y = element_text(color="black", size=10),
        axis.title.x = element_text(color="black", size=10),
        axis.title.y = element_text(color="black", size=10))+
      xlab("")+ylab("cpm")
    p
  }
  ggsave("result/actb_tmsb4x_variance.pdf",width=8,height=3.5)
  
}
### gene boxplot
{
  gene="hsa_circ_0013872"
  gene="ENSG00000271425" # NBPF10
  gene="NBPF10"
  
  gene="hsa_circ_0048555"
  gene="ENSG00000167658" # EEF2
  gene="EEF2"
  
  gene="hsa_circ_0138707"
  transcript="ENST00000602361.1"
  gene="ENSG00000269900" # RP11-331F9.10
  gene="RP11-331F9.10"
  gene="RMRP"
  
  gene="ENSG00000084693"
  gene="ENSG00000124942"
  
  ### top 4 expression genes
  gene="hsa_circ_0000722"
  gene="ENSG00000131149"
  gene="GSE1"
  gene="hsa_circ_0006940" 
  gene=""
  gene="hsa_circ_0001380" 
  gene="ENSG00000163960"
  gene="UBXN7"
  gene="hsa_circ_0001801" 
  gene="ENSG00000168300"
  gene="PCMTD1"
  
  aa=data.frame(cpm=edgeR::cpm(y)[substr(rownames(y),1,15)==gene,],
                sample=colnames(cpm(y)),
                patient,
                cancer=cancer,
                class=class)
  de.all$de[substr(rownames(de.all$de),1,15)==gene,]
  {
    p=ggplot(aa[aa$cancer=="LUAD",], aes(x=class,y=cpm,fill=class))+  # cpm or tpm
      #facet_wrap(~ cancer)+
      geom_boxplot(notch = FALSE, alpha = 0.9, size=0.4, outlier.shape = NA, position=position_dodge(0.9)) +
      geom_jitter(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9), alpha = 0.8, size=0.8)+
      geom_line(aes(group=patient), size=0.4, color="darkgrey", alpha=0.8,linetype="dashed")+
      #scale_y_continuous(limits = c(0,500))+
      #scale_color_brewer(palette="Set1")+
      scale_fill_brewer(palette="Set1")+
      #scale_color_aaas()+
      #scale_y_continuous(limits = c(0,850))+
      ##########修改图例内容
      theme_bw()+
      theme(
        plot.margin = unit(c(0.5, 1, 0.5, 0.5),"cm"),
        panel.grid = element_blank(),
        panel.border=element_rect(size=1,fill="transparent"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(color="black", size=12),
        plot.title = element_text(hjust = 0.5, size=12),
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.title.x = element_text(color="black", size=12),
        axis.title.y = element_text(color="black", size=12))+
      labs(x="",y="CPM",title=gene)  
    p
    p <- p+stat_compare_means(aes(label = paste0("p = ", ..p.format.., "\n    ", ..p.signif..)),
                              method = "wilcox.test",
                              label.x = 1.3,
                              label.y = 1680)+  
      geom_signif(annotations = c(""),
                  y_position = 1650,   # change coordinate everytime
                  xmin = 1, 
                  xmax = 2,
                  tip_length = c(0.02,0.02),
                  textsize=4)
    p
  } 
  ggsave("result/paired/PCMTD1_mRNA_boxplot_LUAD.pdf", p,height = 4,width = 3)
}
### circRNA
{
  data <- read.delim2("count_matrix/paired/count_matrix_data_filter.txt", sep="\t", check.names = F)
  circ <- as.data.frame(rownames(data)[substr(rownames(data),1,3)=="hsa"])
  #write.table(circ,"result/paired/circrna/circ_1065.txt",sep="\t",quote=F,col.names = F,row.names = F)
  #用gencode+circRNA比circRNA找到更多的diffgenes
  de.all$de.cf[substr(rownames(de.all$de.cf),1,3)=="hsa",]
  de.all$de.ev[substr(rownames(de.all$de.ev),1,3)=="hsa",]
  de.all.nc$de.cf[substr(rownames(de.all.nc$de.cf),1,3)=="hsa",]
  de.all.nc$de.ev[substr(rownames(de.all.nc$de.ev),1,3)=="hsa",]
  de.all.crc$de.cf[substr(rownames(de.all.crc$de.cf),1,3)=="hsa",]
  de.all.crc$de.ev[substr(rownames(de.all.crc$de.ev),1,3)=="hsa",]
  de.all.lc$de.cf[substr(rownames(de.all.lc$de.cf),1,3)=="hsa",]
  de.all.lc$de.ev[substr(rownames(de.all.lc$de.ev),1,3)=="hsa",]
  
  ### heatmap
  {
    if(FALSE){
      de.circ <- c("hsa_circ_0013872","hsa_circ_0048555") 
      de.gene <- c("ENSG00000271425.9|14362|NBPF10|protein_coding", "ENSG00000167658.16|4021|EEF2|protein_coding") 
      de.circ <- c("hsa_circ_0001801") 
      de.gene <- c("ENSG00000271425.9|14362|NBPF10|protein_coding") 
      de.circ <- c("hsa_circ_0001801","ENSG00000168300.14|7323|PCMTD1|protein_coding",
                   "hsa_circ_0001062","ENSG00000188177.14|11800|ZC3H6|protein_coding",
                   "hsa_circ_0011014","ENSG00000117713.20|15939|ARID1A|protein_coding")
    }
    #cutoff=1
    de.cf <- de.all$de.cf
    de.ev <- de.all$de.ev
    #cutoof=0
    #de.cf <- filter(de.all$de, logFC<0, FDR<0.1)
    #de.ev <- filter(de.all$de, logFC>0, FDR<0.1)
    
    de.circ <- rbind(de.cf[substr(rownames(de.cf),1,3)=="hsa",],
                     de.ev[substr(rownames(de.ev),1,3)=="hsa",])
    circ <- rownames(de.circ)
    
    circ
    #"hsa_circ_0000247" "hsa_circ_0007082" "hsa_circ_0138706" "hsa_circ_0002454" "hsa_circ_0112897"
    #"hsa_circ_0119456" "hsa_circ_0054674" "hsa_circ_0002743" "hsa_circ_0042961" "hsa_circ_0067678"
    #"hsa_circ_0069748" "hsa_circ_0031447" "hsa_circ_0008034" "hsa_circ_0109923" "hsa_circ_0007695"
    #"hsa_circ_0115643" "hsa_circ_0097631" "hsa_circ_0048555"
    gene=c("ENSG00000156026.14|4533|MCU|protein_coding","ENSG00000198648.11|3701|STK39|protein_coding","ENSG00000269900.3|268|RMRP|lncRNA","ENSG00000116675.16|6958|DNAJC6|protein_coding","ENSG00000117614.10|3773|SYF2|protein_coding",
           "ENSG00000119778.15|10501|ATAD2B|protein_coding","ENSG00000115392.12|2412|FANCL|protein_coding","ENSG00000167257.11|4407|RNF214|protein_coding","ENSG00000178691.11|6551|SUZ12|protein_coding","ENSG00000181804.15|4041|SLC9A9|protein_coding",
           "ENSG00000145216.16|5639|FIP1L1|protein_coding","ENSG00000196792.12|4965|STRN3|protein_coding","ENSG00000138592.14|20641|USP8|protein_coding","ENSG00000130254.12|6018|SAFB2|protein_coding","ENSG00000100485.12|6948|SOS2|protein_coding",
           "ENSG00000180530.11|8497|NRIP1|protein_coding","ENSG00000111011.18|8745|RSRC2|protein_coding","ENSG00000167658.16|4021|EEF2|protein_coding")
    

    z <- DGEList(counts=data, group = class, sample=patient, lib.size = libsize$lib.size)
    dim(z)
    z <- calcNormFactors(z, method="TMM")
    logRPM <- edgeR::cpm(z, log=T)
    
    logRPM.c <- logRPM[circ,]
    logRPM.c.scale <- scale(t(logRPM.c), center = T, scale = T)
    logRPM.c.scale <- t(logRPM.c.scale)
    logRPM.scale = logRPM.c.scale
    
    
    logRPM.g <- logRPM[gene,]
    logRPM.g.scale <- scale(t(logRPM.g), center = T, scale = T)
    logRPM.g.scale <- t(logRPM.g.scale)
    logRPM.scale = logRPM.g.scale
    
    class. = class[order(class)]
    cancer. = cancer[order(class,cancer)]
    logRPM.scale. = logRPM.scale[,order(class,cancer)]
    
    
    ## pheatmap
    {
      library(pheatmap)
      ann_col = data.frame(cancer = as.character(cancer.),
                           class = as.character(class.))
      rownames(ann_col) = colnames(logRPM.scale.)
      #labels_row = c(rep("",17),"hsa_circ_0048555")
      #set colors of each group
      ann_colors = list(cancer = brewer.pal(3,"Set2")[1:3], 
                        class = brewer.pal(3,"Set1")[c(2,1)])
      names(ann_colors$cancer) = c("NC","CRC","LC")
      names(ann_colors$class) = c("CF","EV")
      #GeneClass = c(Path1 = "#807DBA", Path2 = "#9E9AC8", Path3 = "#BCBDDC"))
      
      col = brewer.pal(9,"Set1")
      pdf("result/paired/circRNA/pheatmap_all_gene_FC1_pheatmap.pdf", width = 7, height = 3)
      par(mar=c(2,2,2,2))
      pheatmap(logRPM.scale., 
               #scale="row",
               #color = colorRampPalette(c(col[2],"white",col[1]))(1000),
               color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
               #cutree_col = 2, 
               #cutree_row = 3, #break up the heatmap by clusters you define
               cluster_rows = F,
               cluster_cols = F, #by default, pheatmap clusters by both row and col
               show_rownames = T,
               show_colnames = F,
               fontsize = 6,
               angle_col = 0,
               annotation_col = ann_col,
               #annotation_row = ann_row,
               annotation_colors = ann_colors,
               annotation_names_row = F,
               legend_breaks=-4:4,
               #labels_row = labels_row,
               border=F)
      dev.off()
    }
    ## ggplot
    {
      dd <- as.data.frame(logRPM.scale.) %>% mutate(feature=rownames(logRPM.scale.))
      dd <- melt(dd,id.vars = "feature")
      b=ggplot(dd, aes(variable, feature)) +
        geom_tile(aes(fill = value),colour = "white") + 
        scale_fill_gradient2(limits=c(-2.13,4.13),breaks=c(-2,0,2,4),low = "navy", mid="white",high = "firebrick3")+
        theme_grey(base_size = 10) + 
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(limits=rev(circ),expand = c(0, 0)) + 
        theme(plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
              legend.position = "right",
              panel.background = element_blank(),
              #panel.border = element_rect(),
              axis.ticks = element_blank(), 
              axis.title = element_blank(),
              #axis.text.x = element_text(size = 8, angle=90, hjust = 1),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 8),
              axis.text = element_text(color='black'),
              legend.key.height = unit(0.3,"cm"),
              legend.key.width = unit(0.15,"cm"))
      
      b
    }
    ggsave("result/paired/circrna/pheatmap_all_circ_FC1.pdf",b,width=5,height=3)
  }

  ### gene boxplot
  {
    ### top 4 expression genes
    {
      gene="hsa_circ_0000722"
      gene="ENSG00000131149"
      gene="GSE1"
      gene="hsa_circ_0006940" 
      gene=""
      gene="hsa_circ_0001380" 
      gene="ENSG00000163960"
      gene="UBXN7"
      gene="hsa_circ_0001801" 
      gene="ENSG00000168300"
      gene="PCMTD1"
    }
    
    ## old
    {
      bb=data.frame(PCMTD1=edgeR::cpm(z)["ENSG00000168300.14|7323|PCMTD1|protein_coding",],
                    hsa_circ_0001801=edgeR::cpm(z)["hsa_circ_0001801",],
                    sample=colnames(cpm(z)),
                    patient=patient,
                    cancer=cancer,
                    class=class)
      bb=data.frame(ZC3H6=edgeR::cpm(z)["ENSG00000188177.14|11800|ZC3H6|protein_coding",],
                    hsa_circ_0001062=edgeR::cpm(z)["hsa_circ_0001062",],
                    sample=colnames(cpm(z)),
                    patient=patient,
                    cancer=cancer,
                    class=class)
      bb=data.frame(ARID1A=edgeR::cpm(z)["ENSG00000117713.20|15939|ARID1A|protein_coding",],
                    hsa_circ_0011014=edgeR::cpm(z)["hsa_circ_0011014",],
                    sample=colnames(cpm(z)),
                    patient=patient,
                    cancer=cancer,
                    class=class)
      # 趋势相反
      bb=data.frame(ANKRD12=edgeR::cpm(z)["ENSG00000101745.17|12248|ANKRD12|protein_coding",],
                    hsa_circ_0046849=edgeR::cpm(z)["hsa_circ_0046849",],
                    sample=colnames(cpm(z)),
                    patient=patient,
                    cancer=cancer,
                    class=class)
      
      bb=data.frame(EEF2=edgeR::cpm(z)["ENSG00000167658.16|4021|EEF2|protein_coding",],
                    hsa_circ_0048555=edgeR::cpm(z)["hsa_circ_0048555",],
                    sample=colnames(cpm(z)),
                    patient=patient,
                    cancer=cancer,
                    class=class)
    }

     
    "hsa_circ_0112897"
    "ENSG00000117614.10|3773|SYF2|protein_coding"
    
    bb=data.frame(gene=edgeR::cpm(z)[gene[18],],
                  circ=edgeR::cpm(z)[circ[18],],
                  sample=colnames(cpm(z)),
                  patient=patient,
                  cancer=cancer,
                  class=class)
    
    bb <- melt(bb, id.vars = c("sample","patient","cancer","class"))
    {
      p=ggplot(bb[bb$variable=="gene",], aes(x=class,y=value,fill=class))+  # cpm or tpm
        facet_grid(~cancer)+
        geom_boxplot(notch = FALSE, alpha = 1, size=0.5, outlier.shape = 16, outlier.size = 0.8, position=position_dodge(0.9)) +
        geom_jitter(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9), alpha = 0.8, size=0.8)+
        scale_fill_manual(values=brewer_pal(4,"Set1")(4)[c(2,1)])+
        scale_y_continuous(limits = c(0,2100), breaks=c(0,500,1000,1500,2000))+
        #scale_y_continuous(limits = c(0,100), breaks=c(0,20,40))+
        ##########修改图例内容
        theme_bw()+
        theme(
          strip.text = element_text(size=12, color="black", hjust = 0.5),
          strip.background = element_blank(),
          plot.margin = unit(c(0.5, 1, 0.5, 0.5),"cm"),
          panel.background=element_rect(fill="white", colour="black", size=1),
          panel.grid = element_blank(),
          panel.border=element_blank(),
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(color="black", size=12),
          plot.title = element_text(hjust = 0.5, size=12),
          axis.line = element_blank(),
          axis.text.x = element_text(color="black", size=12),
          axis.text.y = element_text(color="black", size=12),
          axis.title.x = element_text(color="black", size=12),
          axis.title.y = element_text(color="black", size=12))+
        labs(x="",y="CPM",title=unlist(lapply(strsplit(gene[18],"|",fixed=T),function(x) x[3]))) 
        #labs(x="",y="CPM",title=circ[18])  
      
      p
      p <- p +stat_compare_means(aes(label = paste0(..p.signif..)),
                                #aes(label = paste0("p = ", ..p.format.., "\n    ", ..p.signif..)),
                                paired = TRUE,  ########## paired wilcox test
                                method = "wilcox.test",
                                label.x = 1.3,
                                label.y = 1900)+  
        geom_signif(annotations = c(""),
                    y_position = 1880,   # change coordinate everytime
                    xmin = 1, 
                    xmax = 2,
                    tip_length = c(0.02,0.02),
                    textsize=4)
      p
    }
    a=ggarrange(p,q,ncol=2,align="h")
    ggsave("result/paired/circrna/EEF2_gene_boxplot.pdf", a,height = 2.6,width = 6)
  }
  ### top circRNA
  {
    cpm.circ <- cpm(z,log=T)
    cpm.circ.ev <- cpm.circ[,class.=="EV"]
    cpm.circ.cf <- cpm.circ[,class.=="cf"]
    
    
    ################
    # top 10 genes
    #options(scipen=-1) # 999以内不使用科学技术法
    #options(digits=5) #有效数字两位
    genesum <- sort(rowSums(cpm.circ),decreasing = T)
    genesum[1:10]
    
    read.frac <- function(cpm.circ){
      mean <- sort(rowMeans(cpm.circ),decreasing = T)
      mean.top <- mean[names(genesum[c(1:10)])]
      reads.top <- data.frame(reads=mean.top)
      reads.top <- reads.top %>% mutate(fraction = reads / sum(reads))
      gene <- rownames(reads.top)
      reads.top$gene = gene
      reads.top$gene = factor(reads.top$gene, levels=gene[10:1])
      reads.top$reads = log10(reads.top$reads)
      return(reads.top)
    }
    df.cf = read.frac(cpm.circ.cf) %>% mutate(group="cf")
    df.ev = read.frac(cpm.circ.ev) %>% mutate(group="ev")
    
    col=brewer.pal(8,"Blues")
    {
      p=ggplot(df.cf, aes(x=gene, y=reads, fill=fraction))+
        #facet_wrap(~ cancer, ncol = 4)+
        #geom_bar(stat="identity",width=0.7, col="black")+
        geom_bar(stat="identity",position="dodge",width=0.7)+
        coord_flip()+
        scale_y_continuous(limits=c(0,1.5),breaks=c(0,1))+
        scale_fill_gradient(seq(0,0.1,0.2),low=col[2],high=col[7])+
        theme_bw()+
        theme(
          plot.margin = unit(c(1, 1, 1, 1),"cm"),
          panel.border=element_rect(size=1, fill="transparent"),
          panel.grid=element_blank(),
          #legend.position=c(0.25,0.85),
          legend.position="right",
          legend.background = element_blank(),
          legend.title = element_blank(),
          legend.text= element_text(face="bold", color="black", size=12),
          legend.key.height = unit(20, "pt"),
          legend.key.width = unit(10, "pt"),
          plot.title = element_text(face="bold",color="black", size=12),
          axis.text.x = element_text(face="bold",color="black", size=12),
          axis.text.y = element_text(face="bold",color="black", size=12),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          # 文字样式调整
          strip.text.x = element_text(size = 15, color = "black", face = "bold"),
          # 标签背景色调整，fill表示填充色，color表示边框颜色
          strip.background = element_blank())+
        labs(title="cell-free")
      p
    }     
    ggsave("result/paired/circRNA/top10gene-cf.pdf",width=4.2,height=4.5)
    
  }
  ### RCircos
  {
    #CMplot
    {
      #https://www.jianshu.com/p/50b3cf231b34
      #BiocManager::install("CMplot")
      library("CMplot")
      CMplot(df.cm, 
             plot.type = "c", 
             r = 0.5,    # 半径
             cex = 0.5,   # 点的大小
             #H = ,   # 设置每个圈的高度
             multracks = T,  # 设置是否需要绘制多个track
             LOG10 = T,
             threshold = c(0.01,0.05) / nrow(df.cm),  # 设置阈值并添加阈值线
             #threshold.col = c("red","orange"),   # 阈值线颜色
             #threshold.lwd = c(1,2),   # 阈值线的宽度
             #threshold.lty = c(1,2),   # 阈值线的类型
             amplify = T,   # 设置是否放大显著的点
             cir.chr.h = 2,    # 设置染色体边界的高度
             signal.cex = c(2,2),    # 设置显著点的大小
             signal.pch = c(19,20),  # 设置显著点的形状
             signal.col = c("red","green"),  # 设置显著点的颜色
             outward = TRUE  # 设置点的朝向是否向外
      )
      
      df.cm <- df %>% transmute(circRNA=rownames(df),
                                Chromosome=chromosome,
                                Position=chromStart,
                                cpm.cf=cpm.cf,
                                cpm.ev=cpm.ev)
      ### p-value not cpm
      head(df)
    }
    #RCircos
    {
      # https://blog.csdn.net/weixin_42360237/article/details/112581062
      #BiocManager::install("RCircos", ask=F, update=F)
      library("RCircos")
      data("UCSC.HG19.Human.CytoBandIdeogram")
      RCircos.Set.Core.Components(cyto.info = UCSC.HG19.Human.CytoBandIdeogram, 
                                  chr.exclude = NULL, 
                                  tracks.inside = 4, 
                                  tracks.outside = 0)
      #列出所有绘图参数
      params <- RCircos.Get.Plot.Parameters()
      #params$radius.len=3
      #params$plot.radius=4
      params$line.color="white"
      params$hist.color="DarkBlue"
      params$track.background="white"
      params$grid.line.color="white"
      params$hist.width=100
      params$track.height=0.4
      params$sub.tracks=4
      #绘制染色体图形，默认方法显示染色体名称
      RCircos.Reset.Plot.Parameters(params)
      #RCircos.List.Plot.Parameters()
      
      cpm.circ <- cpm(z,log=T)
      cpm.circ.ev.mean <- rowMeans(cpm.circ[,class=="EV"])
      cpm.circ.cf.mean <- rowMeans(cpm.circ[,class=="CF"])
      sort(cpm.circ.ev.mean,decreasing = T)[1:4]
      sort(cpm.circ.cf.mean,decreasing = T)[1:4]
      #top 4 expression circRNA
      #hsa_circ_0000722 
      #hsa_circ_0006940 
      #hsa_circ_0001380 
      #hsa_circ_0001801
      
      ## circbase - list search - export csv
      circgenes <- read.table("result/paired/circRNA/circBase_export_2144.csv",sep=",",header=F)
      colnames(circgenes) <- c("position","strand","circRNA ID","genomic length", "spliced length", "samples","scores","repeats","annotation","best transcript","gene symbol","circRNA study")
      rownames(circgenes) <- circgenes$`circRNA ID`
      df <- circgenes[rownames(cpm(z)),] %>% transmute(chromosome=unlist(lapply(strsplit(position, ":", fixed=T), function(x) x[1])),
                                                       chromStart=unlist(lapply(strsplit(position, "[:-]"), function(x) x[2])),
                                                       chromEnd=unlist(lapply(strsplit(position, "-", fixed=T), function(x) x[2])),
                                                       cpm.cf=cpm.circ.cf.mean,
                                                       cpm.ev=cpm.circ.ev.mean)
      df <- df[df$chromosome %in% c(paste0("chr",seq(1,22)),"chrX","chrY"),]
      unique(df$chromosome)
      df$chromStart <- as.numeric(df$chromStart)
      df$chromEnd <- as.numeric(df$chromEnd)
      
      #绘图
      pdf("result/paired/circRNA/circRNA_chr.pdf",width=5,height=5)
      RCircos.Set.Plot.Area()
      RCircos.Chromosome.Ideogram.Plot()
      params$hist.color="darkblue"
      RCircos.Reset.Plot.Parameters(params)
      RCircos.Histogram.Plot(df, data.col = 4, track.num = 2, side = "in")
      params$hist.color="brown"
      RCircos.Reset.Plot.Parameters(params)
      RCircos.Histogram.Plot(df, data.col = 5, track.num = 1, side = "in")
      dev.off()
    }
  }
  ### backsplicing ratio
  {
    cpm.gencode <- cpm(z)[1:9694,]
    cpm.circ <- cpm(z)[9695:length(rownames(z)),]
    
    # top expressed circRNA
    {
      dim(cpm.circ)
      cpm.circ.ev.mean <- rowMeans(cpm.circ[,class=="EV"])
      cpm.circ.cf.mean <- rowMeans(cpm.circ[,class=="CF"])
      sort(cpm.circ.ev.mean,decreasing = T)[1:4]
      sort(cpm.circ.cf.mean,decreasing = T)[1:4]
      #top 4 expression circRNA
      #hsa_circ_0000722 
      #hsa_circ_0006940 
      #hsa_circ_0001380 
      #hsa_circ_0001801
    }

    gene.list.ncbiid <- read.table("/Users/zhanqing/Desktop/lulab/reference/reference/geneset/gencode.v38.txt",sep = "\t",header = T)
    gene.list.ncbiid <- gene.list.ncbiid[,c("ensg","genename")]
    circgenes <- read.table("result/paired/circRNA/circBase_export.csv",sep=",",header=F)
    colnames(circgenes) <- c("position","strand","circRNA ID","genomic length", "spliced length", "samples","scores","repeats","annotation","best transcript","gene symbol","circRNA study")
    circ_gene_table <- data.frame(circRNA_ID=circgenes$`circRNA ID`,
                                  gene_symbol=circgenes$`gene symbol`,
                                  ensg=gene.list.ncbiid[match(circgenes$`gene symbol`, gene.list.ncbiid$genename),"ensg"])
    circ_gene_table <- circ_gene_table[order(circ_gene_table$circRNA_ID),]
    
    # Backsplicing ratio/Means
    {
      cpm.circ.mean <- data.frame(circ.ev.mean=rowMeans(cpm.circ[,class=="EV"]),
                                  circ.cf.mean=rowMeans(cpm.circ[,class=="CF"]))
      cpm.gencode.mean <- as.matrix(data.frame(gencode.ev.mean=rowMeans(cpm.gencode[,class=="EV"]),
                                               gencode.cf.mean=rowMeans(cpm.gencode[,class=="CF"])))
      rownames(cpm.gencode.mean) <- unlist(lapply(strsplit(rownames(cpm.gencode.mean),".",fixed=T),function(x) x[1]))
      
      circ_gene_table <- circ_gene_table[circ_gene_table$ensg %in% rownames(cpm.gencode.mean),] ### 只留下在gencode表达的circRNA
      circ_ensg <- circ_gene_table[match(rownames(cpm.circ.mean),circ_gene_table$circRNA_ID),"ensg"]
      
      df <- cbind(cpm.circ[!is.na(circ_ensg),],
                  cpm.gencode[circ_ensg[!is.na(circ_ensg)],])
      df <- df %>% mutate(ev.ratio=circ.ev.mean / gencode.ev.mean,
                          cf.ratio=circ.cf.mean / gencode.cf.mean)
      df. <- melt(df[,5:6])
      {
        q=ggplot(df., aes(x=variable,y=value,fill=variable))+  # cpm or tpm
          geom_boxplot(notch = FALSE, alpha = 0.9, size=0.4, outlier.shape = NA, position=position_dodge(0.9)) +
          geom_jitter(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9), alpha = 0.8, size=0.8)+
          scale_fill_brewer(palette="Set1")+
          scale_y_continuous(limits=c(0,0.5))+
          theme_bw()+
          theme(
            plot.margin = unit(c(0.5, 1, 0.5, 0.5),"cm"),
            panel.grid = element_blank(),
            panel.border=element_rect(size=1,fill="transparent"),
            legend.position = "none",
            legend.title = element_blank(),
            legend.text = element_text(color="black", size=12),
            plot.title = element_text(hjust = 0.5, size=12),
            axis.text.x = element_text(color="black", size=12),
            axis.text.y = element_text(color="black", size=12),
            axis.title.x = element_text(color="black", size=12),
            axis.title.y = element_text(color="black", size=12))+
          labs(x="",y="Backsplicing ratio")  
        q
        q <- q+stat_compare_means(aes(label = paste0("p = ", ..p.format.., "\n    ", ..p.signif..)),
                                  method = "wilcox.test",
                                  label.x = 1.28,
                                  label.y = 1700)+  
          geom_signif(annotations = c(""),
                      y_position = 1650,   # change coordinate everytime
                      xmin = 1, 
                      xmax = 2,
                      tip_length = c(0.02,0.02),
                      textsize=4)
        q
      }
      {
        df.2 <- df[,5:6]
        df.2$sum <- df.2$ev.ratio+df.2$cf.ratio
        df.2 <- df.2[order(df.2$sum, decreasing = F),]
        df.2$circRNA_number <- seq(1:nrow(df.2))
        df.2 <- melt(df.2, id.vars = "circRNA_number")
        #df.2 <- df.2[df.2$value<=1,]
        df.2 <- df.2[df.2$variable!="sum",]
        df.2 <- df.2[df.2$circRNA_number %in% c(seq(1,2262,20),2262),]
        q=ggplot(df.2, aes(x=circRNA_number,y=value,color=variable))+  # cpm or tpm
          geom_line()+
          #geom_smooth()+
          scale_color_aaas()+
          scale_y_continuous(limits = c(0,0.7))+
          theme_bw()+
          theme(
            plot.margin = unit(c(0.5, 1, 0.5, 0.5),"cm"),
            panel.grid = element_blank(),
            panel.border=element_rect(size=1,fill="transparent"),
            legend.position = "none",
            legend.title = element_blank(),
            legend.text = element_text(color="black", size=12),
            plot.title = element_text(hjust = 0.5, size=12),
            axis.text.x = element_text(color="black", size=12),
            axis.text.y = element_text(color="black", size=12),
            axis.title.x = element_text(color="black", size=12),
            axis.title.y = element_text(color="black", size=12))+
          labs(x="circRNA number",y="Backsplicing ratio")  
        q
        ggsave("result/paired/circRNA/Backsplicing_ratio.pdf",q,width=3.5,height=3.5)
      }
    }
    # Backsplicing ratio
    {
      colnames(cpm.circ)==colnames(cpm.gencode)
      rownames(cpm.gencode) <- unlist(lapply(strsplit(rownames(cpm.gencode),".",fixed=T),function(x) x[1]))
      circ_gene_table <- circ_gene_table[circ_gene_table$ensg %in% rownames(cpm.gencode),] ### 只留下在gencode表达的circRNA
      
      cpm.circ <- cpm.circ[circ_gene_table$circRNA_ID,]
      cpm.gencode <- cpm.gencode[circ_gene_table[match(rownames(cpm.circ),circ_gene_table$circRNA_ID),"ensg"],]
      dim(cpm.circ)
      dim(cpm.gencode)
      
      cpm.ratio <- cpm.circ/(cpm.gencode+cpm.circ)
      cpm.ratio <- as.data.frame(cpm.circ/(cpm.gencode+1))
      
      cpm.ratio.mean <- data.frame(all.cf=rowMeans(cpm.ratio[,class=="CF"]),
                                   all.ev=rowMeans(cpm.ratio[,class=="EV"]),
                                   nc.cf=rowMeans(cpm.ratio[,cancer=="NC"&class=="CF"]),
                                   nc.ev=rowMeans(cpm.ratio[,cancer=="NC"&class=="EV"]),
                                   crc.cf=rowMeans(cpm.ratio[,cancer=="CRC"&class=="CF"]),
                                   crc.ev=rowMeans(cpm.ratio[,cancer=="CRC"&class=="EV"]),
                                   lc.cf=rowMeans(cpm.ratio[,cancer=="LC"&class=="CF"]),
                                   lc.ev=rowMeans(cpm.ratio[,cancer=="LC"&class=="EV"]))
      a="all.cf"
      b="all.ev"
      cpm.ratio.mean.log <- log2(filter(cpm.ratio.mean[,c(a,b)], all.cf>0, all.ev>0))
      col = brewer.pal(9,"Set1")
      
      # ggscatter
      {
        t=cor.test(cpm.ratio.mean.log[,a],
                 cpm.ratio.mean.log[,b],
                 method = "pearson") # pearson比spearman算出的相关系数更高
        #label=paste("R = ",round(t$estimate,2),", p = ",signif(t$p.value,2),sep="")
        label=paste("R = ",round(t$estimate,2),", p < 2.2e-16",sep="")
        ggscatter(cpm.ratio.mean.log,
                  x = a,
                  y = b,
                  #add = "reg.line", #拟合曲线
                  add.params = list(color=col[2],size=1,alpha=0.6),
                  conf.int = TRUE, #置信区间阴影带
                  cor.coef = TRUE, #系数
                  cor.method = "pearson", #方法
                  xlab = paste("CPM ratio, ",sep=""),
                  ylab = paste("CPM ratio, ",sep=""),
                  shape=16,size=2,color=col[2])
      }

      # ggplot
      { 
        
        library(viridis)
        library(ggpointdensity)
        library(flextable)
        
        p=ggplot(cpm.ratio.mean.log, aes(all.cf, all.ev))+
          #geom_point(color=col[2])+
          geom_pointdensity(size = 1.3, shape = 16, alpha = 0.8)+
          geom_smooth(method = lm, size = 0.5, color = col[9])+
          geom_abline(aes(intercept=0,slope=1),linetype=5,size=0.5, color="black")+
          geom_hline(yintercept = 0,size=0.5,color=col[9], linetype="dashed")+
          geom_vline(xintercept = 0,size=0.5,color=col[9], linetype="dashed")+
          scale_x_continuous(limits = c(-12,4), breaks = c(-12,-8,-4,0,4))+
          scale_y_continuous(limits = c(-12,4), breaks = c(-12,-8,-4,0,4))+
          scale_color_viridis()+
          theme_bw()+
          theme(
            #plot.margin = unit(c(1, 1, 1, 1),"cm"),
            legend.position = "none",
            legend.title = element_text(color="black", size=15),
            legend.text= element_text(color="black", size=15),
            panel.grid = element_blank(),
            panel.border=element_blank(),
            axis.line = element_line(color="black", size=0.5),
            axis.text.x = element_text(color="black", size=15),
            axis.text.y = element_text(color="black", size=15),
            axis.title.x = element_text(color="black", size=15),
            axis.title.y = element_text(color="black", size=15))+
          labs(x="log2(CPM ratio), CF",y="log2(CPM ratio), EV",face="bold")+
        annotate("text",x=-6,y=4,label=label,size=5,angle=0,color="black")
        p
      }
      ggsave("result/paired/circrna/circratio-ALL.pdf", p, width=3.5, height=3.5)
    
      # ratio boxplot
      circ.="hsa_circ_0111150"
      gene.="ENSG00000197959.14|11329|DNM3|protein_coding"
      circ.="hsa_circ_0055926"
      gene.="ENSG00000115641.19|6626|FHL2|protein_coding"
      circ.="hsa_circ_0003553"
      gene.="ENSG00000117602.13|7138|RCAN3|protein_coding"
      circ.="hsa_circ_0011014"
      gene.="ENSG00000117713.20|15939|ARID1A|protein_coding"
      bb=data.frame(gene=edgeR::cpm(z)[gene.,],
                    circ=edgeR::cpm(z)[circ.,],
                    sample=colnames(cpm(z)),
                    patient=patient,
                    cancer=cancer,
                    class=class)
      bb$ratio <- bb$circ/(bb$gene+1)
      bb <- melt(bb, id.vars = c("sample","patient","cancer","class"))
      {
        q=ggplot(bb[bb$variable=="ratio",], aes(x=class,y=value,fill=class))+  # cpm or tpm
          facet_grid(~cancer)+
          geom_boxplot(notch = FALSE, alpha = 1, size=0.5, outlier.shape = NA, outlier.size = 0.8, position=position_dodge(0.9)) +
          geom_jitter(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9), alpha = 0.8, size=0.8)+
          #geom_line(aes(group=patient), size=0.4, color="darkgrey", alpha=0.8,linetype="dashed")+
          scale_fill_manual(values=brewer.pal(4,"Paired")[c(2,4)])+
          #scale_y_continuous(limits=c(0,1))+
          ##########修改图例内容
          theme_bw()+
          theme(
            strip.text = element_text(size=12, color="black", hjust = 0.5),
            strip.background = element_blank(),
            plot.margin = unit(c(0.5, 1, 0.5, 0.5),"cm"),
            panel.background=element_rect(fill="white", colour="black", size=1),
            panel.grid = element_blank(),
            panel.border=element_blank(),
            legend.position = "none",
            legend.title = element_blank(),
            legend.text = element_text(color="black", size=12),
            plot.title = element_text(hjust = 0.5, size=12),
            axis.line = element_blank(),
            axis.text.x = element_text(color="black", size=12),
            axis.text.y = element_text(color="black", size=12),
            axis.title.x = element_text(color="black", size=12),
            axis.title.y = element_text(color="black", size=12))+
          labs(x="",y="CPM ratio",title=circ.)  
        q
        q <- q+stat_compare_means(aes(label = paste0(..p.signif..)),
                                  #aes(label = paste0("p = ", ..p.format.., "\n    ", ..p.signif..)),
                                  paired = TRUE,  ########## paired wilcox test
                                  method = "wilcox.test",
                                  label.x = 2,
                                  label.y = 0.5)+  
          geom_signif(annotations = c(""),
                      y_position =  0.4,   # change coordinate everytime
                      xmin = 1, 
                      xmax = 2,
                      tip_length = c(0.02,0.02),
                      textsize=4)
        q
      }
      
      # ratio-test
      {
        y <- DGEList(counts=cpm.ratio, group=class)
        design <- model.matrix(~class)
        y <- estimateDisp(y, design)
        fit <- glmFit(y, design)
        lrt <- glmLRT(fit, coef=2)
        top <- topTags(lrt, n=nrow(y))$table

      
    }
  }
}




}


########### phaseI+phaseII
{
  annotate <- function(input){
    # input is library id, output is sample_id
    library_id <- input
    metadata <- read.table("/Users/zhanqing/Desktop/iPICO/output/PhaseII-multiplex/count_matrix/annotation.txt", sep = "\t", header = T)
    index=match(library_id, metadata$library_id)
    sample_id = c()
    for (i in 1:length(library_id)){
      sample_id[i] = metadata$sample_id[index[i]]
    }
    input <- sample_id
  }
  
  gencode.1.1 <- read.delim2("/Users/zhanqing/Desktop/iPICO/output/qhsky5/count_matrix/L1-L21/matrix/gencode.txt", header = T, row.names = 1, check.names = F) 
  gencode.1.2 <- read.delim2("/Users/zhanqing/Desktop/iPICO/output/qhsky5/count_matrix/reproducibility/gencode.txt", header = T, row.names = 1, check.names = F)[,-c(1:9)] 
  gencode.1 <- cbind(gencode.1.1, gencode.1.2)
  colnames(gencode.1) <- annotate(colnames(gencode.1))
  filter_sample_phaseI <- read.delim2("/Users/zhanqing/Desktop/iPICO/output/qhsky5/QC/L1-L21/filter_sample_phaseI.txt", header=T, check.names = F)[,1]
  libsize.1.1 <- read.delim2("/Users/zhanqing/Desktop/iPICO/output/qhsky5/count_matrix/L1-L21/matrix/libsize.txt",sep="\t",header=T,row.names = 1)
  libsize.1.2 <- read.delim2("/Users/zhanqing/Desktop/iPICO/output/qhsky5/QC/reproducibility/libsize.txt",sep="\t",header=T,row.names = 1)
  libsize.1 <- rbind(libsize.1.1, libsize.1.2)
  rownames(libsize.1) <- annotate(rownames(libsize.1))

  gencode.2 <- read.delim2("/Users/zhanqing/Desktop/iPICO/output/PhaseII-multiplex/count_matrix/gencode.txt", header = T, row.names = 1, check.names = F)
  colnames(gencode.2) <- annotate(colnames(gencode.2)) 
  filter_sample_phaseII <- read.delim2("/Users/zhanqing/Desktop/iPICO/output/PhaseII-multiplex/QC/filter_sample_phaseII.txt", header=T, check.names = F)[,1]
  libsize.2 <-  read.delim2("/Users/zhanqing/Desktop/iPICO/output/PhaseII-multiplex/count_matrix/libsize.txt",  header = T, row.names = 1, check.names = F)
  rownames(libsize.2) <- annotate(rownames(libsize.2))
  
  
  annotation <- read.delim2("/Users/zhanqing/Desktop/iPICO/output/PhaseII-multiplex/count_matrix/annotation.txt", header = T)
  
  data <- cbind(gencode.1[,filter_sample_phaseI], gencode.2[,filter_sample_phaseII])
  libsize <- c(libsize.1[filter_sample_phaseI,], libsize.2[filter_sample_phaseII,])
  names(libsize) <- c(filter_sample_phaseI, filter_sample_phaseII)
  
  cancer <- unlist(lapply(strsplit(colnames(data),"-",fixed=T),function(x) x[1]))
  cancer <- factor(cancer, levels = c("NC","CRC","LUAD","HCC","OSCC","GC"))
  class <- unlist(lapply(strsplit(colnames(data),"-",fixed=T),function(x) x[3]))
  class[is.na(class)] <- "cf"
  patient <- unlist(lapply(strsplit(colnames(data),"-EV",fixed=T),function(x) x[1]))
  group <- annotation[match(colnames(data),annotation$sample_id),"group"]
  
  
  y <- DGEList(counts=data)
  y <- DGEList(counts=data,lib.size = libsize)   ##### libsize
  RPM <- round(edgeR::cpm(y),3)
  geneid <- rownames(RPM)
  split <- strsplit(geneid,"|",fixed=T)
  length <- unlist(lapply(split,function(split) split[2]))
  y$genes <- data.frame(length=as.numeric(length))
  RPKM <- round(rpkm(y),3)
  TPM <- round((RPKM/colSums(RPKM))*10^6,3)
  cpm <- as.data.frame(RPM)
  tpm <- as.data.frame(TPM)
    
  ### edgeR: tpm cutoff to filter genes
  {
    # remove genes that are expressed less than 50% samples(tpm)
    y <- DGEList(counts=data, lib.size = libsize) 
    RPM <- round(edgeR::cpm(y),3)
    geneid <- rownames(RPM)
    split <- strsplit(geneid,"|",fixed=T)
    length <- unlist(lapply(split,function(split) split[2]))
    y$genes <- data.frame(length=as.numeric(length))
    RPKM <- round(rpkm(y),3)
    TPM <- round((RPKM/colSums(RPKM))*10^6,3)
    cpm <- as.data.frame(RPM)
    tpm <- as.data.frame(TPM)
    ###################
    ### group = class
    {
      # class
      cf = tpm[,class=="cf"]
      ev = tpm[,class=="EV"]
      n.cf = length(class[class=="cf"])/2
      n.ev = length(class[class=="EV"])/2
      
      detected_genes = data.frame(
        cellfree <- c(
          sum(rowSums(cf>1)>= n.cf),
          sum(rowSums(cf>2)>= n.cf),
          sum(rowSums(cf>5)>= n.cf)),
        EV <- c(
          sum(rowSums(ev>1)>= n.ev),
          sum(rowSums(ev>2)>= n.ev),
          sum(rowSums(ev>5)>= n.ev)))
      colnames(detected_genes) = c("cf","ev")
      rownames(detected_genes) = c("TPM>1(50%)","TPM>2(50%)","TPM>5(50%)")
      # stacked bar plot
      detected <- detected_genes %>% mutate(group=rownames(detected_genes))
      detected <- melt(detected)
      detected$group = factor(detected$group, levels=c("TPM>1(50%)","TPM>2(50%)","TPM>5(50%)"))
      {
        p=ggplot(detected, aes(x=variable, y=value, fill=group))+
          #geom_bar(stat="identity",width=0.7, col="black")+
          geom_bar(stat="identity",position="dodge",width=0.7, col="black")+
          #scale_y_continuous(limits=c(0,20000))+
          scale_fill_brewer()+
          theme_bw()+
          guides(fill=guide_legend(title=NULL, ncol=1))+
          theme(
            plot.margin = unit(c(1, 1, 0, 1),"cm"),
            panel.border=element_rect(size=1, fill="transparent"),
            panel.grid=element_blank(),
            #legend.position=c(0.25,0.85),
            legend.position="top",
            legend.background = element_blank(),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black", size=12),
            axis.text.x = element_text(face="bold",color="black", size=12),
            axis.text.y = element_text(face="bold",color="black", size=12),
            axis.title.x = element_text(face="bold",color="black", size=12),
            axis.title.y = element_text(face="bold",color="black", size=12))+
          ylab("Number of detected genes")+xlab("")
        #geom_vline(aes(xintercept=6.5))+
        p
      }             
      #ggsave("result/genenumber_cf-EV_TPMnormalized.pdf",p,width=4.5,height=4.5)               
      
      cf.genes <- rownames(cf)[rowSums(cf>1) >= n.cf]
      ev.genes <- rownames(ev)[rowSums(ev>1) >= n.ev]
      
      # venn plot
      {
        library(VennDiagram)
        color = brewer.pal(5, "Set1")
        venn <- venn.diagram(
          x = list(`cf`=cf.genes,
                   `EV`=ev.genes),
          filename = NULL,
          fill = color[1:2],
          alpha = 0.4,
          lwd = 2,
          lty = "longdash",
          col = "black", #线条色
          cex = 2.5,
          fontfamily = "sans",
          #fontface = "bold", #加粗
          cat.col = "black",
          cat.default.pos = "text", #位置, outer内 text外
          cat.cex = 2.5,
          cat.fontfamily = "sans",     
          #cat.fontface = "bold",     
          cat.dist = c(0.15, 0.15), #位置，用圆的距离
          cat.pos = c(-180,180), #位置，用圆的度数
          main = "",
          main.col = "black",
          main.cex = 1.5
          #main.fontface = "bold"
        );
        
        pdf("result/paired/gene_number_cf-EV_TPM>1_venn.pdf", width=5,height=5)
        grid.draw(venn);
        dev.off()
        ggsave("result/paired/gene_number_cf-EV_TPM>1_venn.pdf",venn,width = 5,height = 5)
      }  
      {
        venn=ggvenn(list(cf=cf.genes, EV=ev.genes),
                    fill_color=color[1:2],
                    stroke_linetype="longdash",
                    set_name_size=8,
                    text_size=6,
                    show_percentage = F)
        #ggsave("result/paired/gene_number_cf-EV_TPM>1_venn2.pdf",venn,width = 5,height = 5)
      }
    }
    filter = (rowSums(cf>1) >= n.cf) | (rowSums(ev>1) >= n.ev)
    table(filter)
    y <- y[filter,]
    # tpm>1: 7623
    ###################
    ### group = cancer
    {
      group = cancer
      crc = tpm[,group=="CRC"]
      luad = tpm[,group=="LUAD"]
      n.crc = length(group[group=="CRC"])/2
      n.luad = length(group[group=="LUAD"])/2
      
      detected_genes = data.frame(
        CRC <- c(sum(rowSums(crc>2)>= n.crc),
                 sum(rowSums(crc>5)>= n.crc),
                 sum(rowSums(crc>1)>= n.crc)),
        LUAD <- c(sum(rowSums(luad>2)>= n.luad),
                  sum(rowSums(luad>5)>= n.luad),
                  sum(rowSums(luad>1)>= n.luad)))
      colnames(detected_genes) = c("CRC","LUAD")
      rownames(detected_genes) = c("TPM>2(50%)","TPM>5(50%)","TPM>1(50%)")
      detected <- detected_genes %>% mutate(group=rownames(detected_genes))
      detected <- melt(detected)
      detected$group = factor(detected$group, levels=c("TPM>1(50%)","TPM>2(50%)","TPM>5(50%)"))
      {
        p=ggplot(detected, aes(x=variable, y=value, fill=group))+
          #geom_bar(stat="identity",width=0.7, col="black")+
          geom_bar(stat="identity",position="dodge",width=0.7, col="black")+
          #scale_y_continuous(limits=c(0,20000))+
          scale_fill_brewer()+
          theme_bw()+
          guides(fill=guide_legend(title=NULL, ncol=1))+
          theme(
            plot.margin = unit(c(1, 1, 0, 1),"cm"),
            panel.border=element_rect(size=1, fill="transparent"),
            panel.grid=element_blank(),
            #legend.position=c(0.25,0.85),
            legend.position="top",
            legend.background = element_blank(),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black", size=12),
            axis.text.x = element_text(face="bold",color="black", size=12),
            axis.text.y = element_text(face="bold",color="black", size=12),
            axis.title.x = element_text(face="bold",color="black", size=12),
            axis.title.y = element_text(face="bold",color="black", size=12))+
          ylab("Number of detected genes")+xlab("")
        #geom_vline(aes(xintercept=6.5))+
        p
      }   
      #ggsave("L1-L21/diff_exp_24v24/normalized/detected_barplot_(TPMnormLibsize).pdf",p,width=4.5,height=4.5) 
    }
    filter = (rowSums(nc>1) >= n.nc) | (rowSums(crc>1) >= n.crc) | (rowSums(luad>1) >= n.luad)
    table(filter)
    y <- y[filter,]
    
  }
  ### edgeR: filterbyExpr
  {
    y <- DGEList(counts=data,lib.size = libsize)
    keep <- filterByExpr(y, group = class, min.count=2, min.prop=0.1)
    y <- y[keep, ,keep.lib.size=TRUE]
    dim(y)   
  }
  ### PCA for all samples
  {
    library("preprocessCore")
    library("Rtsne")
    # normalize.quantiles: Using a normalization based upon quantiles, this function normalizes a matrix of probe level intensities.
    # normalize: Normalizes numeric data to a given scale.
    
    y <- DGEList(counts=data)
    y <- DGEList(counts=data,lib.size = libsize)
    y <- y[filter,]
    ###### X <- t(X)
    X <- t(edgeR::cpm(y))
    #X <- t(tpm)
    #X <- t(edgeR::cpm(y, log=T))

    
    # z-score scaling
    X_norm.4 <- scale(X, center = T, scale = T)

    ## Visualization
    {  
      Y <- cancer
      colors <- rainbow(length(unique(Y)))
      colors_plot <- as.character(Y)
      colors_plot[which(colors_plot == "NC")] <- colors[1]
      colors_plot[which(colors_plot == "CRC")] <- colors[2]
      colors_plot[which(colors_plot == "LUAD")] <- colors[3]
      colors_plot[which(colors_plot == "HCC")] <- colors[4]
    }
    {
      Y <- class
      colors <- rainbow(length(unique(Y)))
      colors_plot <- as.character(Y)
      colors_plot[which(colors_plot == "cf")] <- colors[1]
      colors_plot[which(colors_plot == "EV")] <- colors[2]
    }
    
    # PCA
    pca <- prcomp(X_norm.4, center = F, scale = F, rank. = 2)  
    pca.var <- pca$sdev^2  ## sdev是标准偏差，十个样本，就有十个标准偏差，平方是避免负数的干扰
    pca.var.per <- round(pca.var/sum(pca.var)*100, 1)  ##求每个样本的variation
    barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")  ##用柱状图可视化
    summary(pca)
    plot(pca$x, col = colors_plot, asp = 3, pch = 20, xlab = "component_1", ylab = "component_2", main = "PCA plot")
    ## ggplot
    pca.data <- data.frame(pca$x, class=class, cancer=cancer, group=group)
    percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
    percentage <- paste(colnames(pca.data)[1:9],"(", paste(as.character(percentage)), "%", ")", sep="")
    {
      p=ggplot(pca.data, aes(x=PC1,y=PC2,color=cancer)) + 
        geom_point(size=2,alpha=1) +
        scale_color_brewer(palette = "Set1") + 
        theme_bw()+
        theme(
          plot.margin = unit(c(1, 1, 1, 1),"cm"),
          panel.border=element_rect(size=1),
          panel.grid=element_blank(),
          legend.background = element_rect(color="grey"),
          legend.key.size = unit(0.1,"cm"),
          legend.position = c(0.8,0.85),
          legend.title = element_blank(),
          legend.text= element_text(color="black", size=10),
          axis.text.x = element_text(color="black", size=15),
          axis.text.y = element_text( color="black", size=15),
          axis.title.x = element_text(color="black", size=15),
          axis.title.y = element_text(color="black", size=15))+
        xlab(percentage[1]) + ylab(percentage[2])
      p
    }
    ggsave("result/phase12/pca.cpm.normlibsize.pdf",p,width = 4, height = 4)
    
    # tSNE
    tsne <- Rtsne(X_norm.4, dims=2, check_duplicates = FALSE, perplexity = 7) 
    # dims参数设置降维之后的维度，默认值为2
    # perplexity参数的取值必须小于(nrow(data) - 1 )/ 3
    # theta参数取值越大，结果的准确度越低，默认值为0.5
    # max_iter参数设置最大迭代次数。
    plot(tsne$Y, col = colors_plot, asp = 3, pch = 20,
         xlab = "tSNE_1", ylab = "tSNE_2", main = "tSNE plot")
    
  }
}


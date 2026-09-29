# RStudio工作目录设为本文件夹
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(RColorBrewer)
library(edgeR)
dir.create("generated/QC",recursive=TRUE,showWarnings=FALSE)
mt=read.table("inputs/1204/totalRNA.matrix.txt",sep="\t",header=TRUE,check.names=FALSE)
mt=mt[rowSums(mt)!=0,]
rna.type=read.table("inputs/1204/totalRNA.matrix.type.txt",sep="\t",header=TRUE)
rownames(rna.type)=rna.type[[1]]
cpm=as.data.frame(edgeR::cpm(mt))
test.type=rna.type[rownames(cpm),]$type2
group=sub("-.*","",colnames(cpm))
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
    
  }
#######  
  # 
  dfall <- dfall[order(match(dfall$celltype, unique(group)), tolower(as.character(dfall$type))),]
  dfall. <- dfall
  dfall.$type2 <- as.character(dfall.$type)
  dfall <- dfall[order(match(dfall$celltype, unique(group)), tolower(as.character(dfall$type))),]
  dfall. <- dfall.[dfall.$type!="IGTR",]
  dfall.$type2[dfall.$type2 %in% c("miRNA","piRNA","snoRNA","snRNA","scaRNA","sRNA","tRNA","Y_RNA")] = "sncRNA"
  dfall.$type2 <- factor(dfall.$type2, levels=c("protein_coding","lncRNA","sncRNA","pseudogene","srpRNA","tucpRNA"))
  dfall.$celltype <- factor(dfall.$celltype, levels = c("MDA5", "ARS", "HC"))
  
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
        legend.margin=margin(11,5.5,5.5,5.5),
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
        legend.margin=margin(11,5.5,5.5,5.5),
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
  ggsave("generated/QC/barplot-rnatype-count.pdf",p,width=9,height=4.6)
ggsave("generated/QC/RNA_biotype_total.pdf",a,width=4.5,height=4.6)
ggsave("generated/QC/RNA_biotype_small.pdf",b,width=4.5,height=4.6)

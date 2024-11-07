
## Fig6D: GSEA GO plot(pickup)
{
  tabledir="/Users/wangge/Documents/DM/cfRNA/enrich_all-sample-table/"
  plotdir="/Users/wangge/Documents/DM/cfRNA/plot/"
  dir.create(plotdir, showWarnings = T)
  col2 = pal_nejm()(8)
  all.gse.go <- read.table(paste0(tabledir,"gse.GO.txt"),sep="\t",header=T,quote="")
  all.gse.go <- all.gse.go %>% mutate(logp = -log10(p.adjust),
                                      Count = str_count(core_enrichment,"/")+1)
  all.gse.go <- all.gse.go[order(all.gse.go$p.adjust, -abs(all.gse.go$NES), -all.gse.go$Count),]
  # up
  up <- filter(all.gse.go, NES>1, p.adjust<0.1)
  up <- filter(up, logp<3)
  
  select <- c(#Muscle Activity
    "heart contraction",
    "cardiac muscle cell differentiation", 
    "regulation of muscle hypertrophy" ,
    "cardiac muscle contraction",  
    "actin binding",
    #ion Channels and Transport
    "calcium channel activity",
    "potassium ion transmembrane transport" ,        
    "ligand-gated monoatomic cation channel activity" ,
    #"lymphocyte activation",
    #"defense response to other organism",
    #"immune response-regulating signaling pathway", 
    #"Fc receptor signaling pathway",
    #"antigen receptor-mediated signaling pathway",
    #"activation of innate immune response",
    #"T cell receptor signaling pathway",
    #"nuclear speck",
    #"ribonucleoprotein complex",
    #"spliceosomal complex",
    
    #cell migration related —— only LC???????
    #"anchoring junction",
    "tRNA catabolic process",
    "ATP-dependent activity")
  up <- up[up$Description %in% select,]
  up$Description <- factor(up$Description, levels=rev(select))
  {
    p=ggplot(up, aes(x=Description,y=NES,fill=logp))+
      geom_bar(stat="identity",width=0.7,size=0.2,color="black")+
      #geom_text(aes(label=signif),size=10,vjust=0.8,hjust=1)+
      #geom_text(aes(label=signif),size=8,vjust=0.8,hjust=0)+
      coord_flip()+
      #geom_point(aes(color=p.adjust,size = Count)) +
      scale_fill_gradient2(name=bquote('-'~log[10]~'FDR'), breaks=c(1.5,2.5), high=col2[1], guide="colorbar")+
      #scale_y_continuous(limits = c(0,220), breaks=c(0,150))+
      scale_y_continuous(limits = c(0,2), breaks=c(0,1,2))+
      theme_bw()+
      theme(
        plot.margin = unit(c(1, 1, 1, 1),"cm"),
        plot.title = element_text(color="black", size=15, hjust = 0.5),
        panel.border = element_rect(size=1, color="black"),
        panel.grid = element_blank(),
        legend.position="right",
        legend.background = element_blank(),
        legend.title = element_text(color="black", size=10, angle=90, hjust=0.5),
        legend.text= element_text(color="black", size=10),
        legend.key.size = unit(0.2,"cm"),
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12, margin=ggplot2::margin(1,0.5,0,0,"cm")),
        axis.title.x = element_text(color="black", size=12),
        axis.title.y = element_blank())+
      labs(title="",y="NES")
    p
  }
  # down
  down <- filter(all.gse.go, NES < -1, p.adjust<0.1)
  down$logp <- as.numeric(down$logp)
  down <- down[order(down$p.adjust, abs(down$NES), down$Count),]
  select <- c(#splicing and RNPs 
    "pre-mRNA binding",
    "spliceosomal complex assembly", 
    "protein-RNA complex assembly",
    #"formation of quadruple SL/U4/U5/U6 snRNP",
    #"Sm-like protein family complex", 
    
    #immune response
    "innate immune response in mucosa",
    "killing by host of symbiont cells",# NC signif
    "organ or tissue specific immune response",
    "mucosal immune response",
    "antimicrobial humoral immune response mediated by antimicrobial peptide", 
    
    #others
    "negative regulation of megakaryocyte differentiation"
  )
  down <- down[down$Description %in% select,]
  down$Description <- factor(down$Description, levels=rev(select))
  {
    q=ggplot(down, aes(x=Description,y=abs(NES),fill=logp))+
      geom_bar(stat="identity",width=0.7,size=0.2,color="black")+
      coord_flip()+
      #geom_point(aes(color=p.adjust,size = Count)) +
      #scale_fill_gradient2(name="log10(FDR)",breaks=c(1,2,3), high=muted("dodgerblue4"), guide="colorbar")+
      scale_fill_gradient2(name=bquote('-'~log[10]~'FDR'), breaks=c(1.5,2.5), high=col2[2], guide="colorbar")+
      #scale_y_continuous(limits=c(0,23),breaks = c(0,20))+
      scale_y_continuous(limits = c(0,2.5), breaks=c(0,1,2))+
      theme_bw()+
      theme(
        plot.margin = unit(c(1, 1, 1, 1),"cm"),
        plot.title = element_text(color="black", size=15, hjust = 0.5),
        panel.border = element_rect(size=1, color="black"),
        panel.grid = element_blank(),
        legend.position="right",
        legend.background = element_blank(),
        legend.title = element_text(color="black", size=10, angle=90, hjust=0.5),
        legend.text= element_text(color="black", size=10),
        legend.key.size = unit(0.2,"cm"),
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12, margin=ggplot2::margin(1,0.5,0,0,"cm")),
        axis.title.x = element_text(color="black", size=12),
        axis.title.y = element_blank())+
      labs(title="",y="-NES")
    q
  }
  a=ggarrange(q,p,ncol=1, heights = c(1.2, 1.3),align="v")
  a
  ggsave("cfRNA/barplot-gsea.GO.EV-CF.pdf",a,width=8.5,height=6.5)
  #ggsave(paste0(figdir,"figure6/gsea-barplot-go-cfev-new.pdf"),a,width=5.7,height=6)
  
  
  
  
  #write.table(up,paste0(tabledir,"CFEV_GSE_GO_UP.txt"),sep="\t",quote = F,row.names = F)
  #write.table(down,paste0(tabledir,"CFEV_GSE_GO_DOWN.txt"),sep="\t",quote = F,row.names = F)
  
  
  {
    
    test <- strsplit(down$core_enrichment,"/",fixed=T)
    down1=test[[6]]
    down2=test[[4]]
    
    dev.off()
    venn <- venn.diagram(
      x = list(n1=down1,
               n2=down2),
      filename=NULL)
    grid.draw(venn);
    dev.off()
  }
  
  
  
}

library(ggplot)
library(tidyverse)

### miRPathDB
### the pathway of miRNA targets
path.miR.590.5p <- read.csv("miRPathDB_miR.590.5p.csv")
pathways <- c(
  "TGF-beta Receptor Signaling",
  "Signaling by TGF-beta Receptor Complex in Cancer",
  "Canonical and Non-Canonical TGF-B signaling",
  "TGF-beta receptor signaling activates SMADs",
  "Loss of Function of TGFBR1 in Cancer",
  "Loss of Function of TGFBR2 in Cancer",
  "MicroRNAs in cancer",
  "cardiac muscle cell proliferation",
  "muscle tissue development",
  "striated muscle tissue development",
  "muscle organ development",
  "TGF-beta signaling pathway",
  "blood vessel morphogenesis",
  "vasculature development",
  "blood vessel development",
  "cell migration",
  "anatomical structure morphogenesis",
  "cell morphogenesis involved in differentiation"
)
path.miR.590.5p <- path.miR.590.5p %>%
  filter(Pathway %in% pathways)
path.miR.590.5p <- path.miR.590.5p[c(2,6,7,10,15,16,17,18,19,21,22,23,28,29,30,34),]
path.miR.590.5p$Pathway <- factor(path.miR.590.5p$Pathway, levels = path.miR.590.5p$Pathway)
a1=ggplot(path.miR.590.5p, aes(x=Pathway, y=-log10(P.value))) + geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95),
        plot.margin = margin(0,0,0,15,unit = "pt")) 
a1
ggsave("Path_miR.590.5p.pdf",a1, width = 4.5, height = 6)

# 创建自定义颜色的x轴文字
custom_colors <- ifelse(1:nrow(path.miR.590) %in% c(2,3,10,11,12), "red", "black")

# 绘图
c1 <- ggplot(path.miR.590, aes(x=Pathway, y=-log10(P.value))) +
  geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95, color=custom_colors),
        plot.margin = margin(0,0,0,15,unit = "pt"))

c1
ggsave("Path_miR.590.5p_custom_text_colors.pdf", c1, width = 4.5, height = 6)



path.miR.590.3p <- read.csv("miRPathDB_targets/miRPathDB_miR.590.3p.csv")
pathways <- c(
  "MicroRNAs in cancer",
  "cellular response to steroid hormone stimulus",
  "cellular response to organic cyclic compound",
  "cellular response to DNA damage stimulus",
  "SMAD binding",
  "regulation of striated muscle tissue development",
  "Mesodermal Commitment Pathway",
  "heart development",
  "cell morphogenesis",
  "cellular response to lipid",
  "nucleocytoplasmic transport",
  "response to lipid",
  "response to steroid hormone",
  "urogenital system development",
  "eye development",
  "visual system development"
)
path.miR.590.3p <- path.miR.590.3p %>%
  filter(Pathway %in% pathways)
path.miR.590.3p <- path.miR.590.3p[-c(1,4,6,10,12,19),]
path.miR.590.3p$Pathway <- factor(path.miR.590.3p$Pathway, levels = path.miR.590.3p$Pathway)
a2=ggplot(path.miR.590.3p, aes(x=Pathway, y=-log10(P.value))) + geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95),
        plot.margin = margin(0,0,0,15,unit = "pt")) 
a2
ggsave("Path_miR.590.3p.pdf",a2, width = 4.5, height = 6)

custom_colors <- ifelse(1:nrow(path.miR.590.3p) %in% c(12), "red", "black")

# 绘图
c2 <- ggplot(path.miR.590.3p, aes(x=Pathway, y=-log10(P.value))) +
  geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95, color=custom_colors),
        plot.margin = margin(0,0,0,15,unit = "pt"))

c2
ggsave("Path_miR.590.3p_custom_text_colors.pdf", c2, width = 4.5, height = 5.5)



path.miR.409.5p <- read.csv("miRPathDB_targets/miRPathDB_miR.409.5p.csv")
pathways <- c(
  "regulation of RNA metabolic process",
  "transcription by RNA polymerase II",
  "homophilic cell adhesion via plasma membrane adhesion molecules",
  "intracellular transport",
  "gene expression",
  "cell-cell adhesion via plasma-membrane adhesion molecules",
  "regulation of RNA splicing",
  "glutamatergic synapse",
  "nuclear transport",
  "regulation of DNA-binding transcription factor activity",
  "Competing endogenous RNAs (ceRNAs) regulate PTEN translation",
  "PIP3 activates AKT signaling",
  "Post-transcriptional silencing by small RNAs",
  "Regulation of PTEN mRNA translation",
  "regulation of mitotic spindle organization",
  "Adherens junction",
  "nucleocytoplasmic transport"
)
path.miR.409.5p <- path.miR.409.5p %>%
  filter(Pathway %in% pathways)
path.miR.409.5p <- path.miR.409.5p[-c(1,3,8,11,13,18),]
path.miR.409.5p$Pathway <- factor(path.miR.409.5p$Pathway, levels = path.miR.409.5p$Pathway)
a3=ggplot(path.miR.409.5p, aes(x=Pathway, y=-log10(P.value))) + geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95),
        plot.margin = margin(0,0,0,60,unit = "pt")) 
a3
ggsave("Path_miR.409.5p.pdf",a3, width = 5, height = 7)


custom_colors <- ifelse(1:nrow(path.miR.409.5p) %in% c(1,3,4,5,8,9,10,11), "red", "black")

# 绘图
c3 <- ggplot(path.miR.409.5p, aes(x=Pathway, y=-log10(P.value))) +
  geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95, color=custom_colors),
        plot.margin = margin(0,0,0,60,unit = "pt"))

c3
ggsave("Path_miR.409.5p_custom_text_colors.pdf", c3, width = 5, height = 7)



path.miR.409.3p <- read.csv("miRPathDB_targets/miRPathDB_miR.409.3p.csv")
pathways <-c(
             "regulation of RNA splicing", 
             "cellular protein modification process", 
             "histone binding", 
             "covalent chromatin modification", 
             "cell morphogenesis involved in differentiation", 
             "regulation of cellular component organization", 
             "positive regulation of macromolecule biosynthetic process", 
             "phosphate-containing compound metabolic process", 
             "cell morphogenesis", 
             "cell development", 
             "regulation of phosphorus metabolic process", 
             "negative regulation of transcription by RNA polymerase II", 
             "anatomical structure morphogenesis", 
             "cell junction", 
             "adherens junction", 
             "positive regulation of cellular biosynthetic process", 
             "cytoplasmic stress granule")
pathways <- c(
  "IL-6 signaling pathway",
  "T-Cell antigen Receptor (TCR) pathway during Staphylococcus aureus infection",
  "B cell receptor signaling pathway",
  "cellular response to fibroblast growth factor stimulus",
  "Signaling by Interleukins",
  "muscle cell proliferation",
  "regulation of smooth muscle cell proliferation",
  "regulation of striated muscle tissue development",
  "VEGFA-VEGFR2 Signaling Pathway",
  "Wnt signaling pathway, calcium modulating pathway",
  "Adherens junction",
  "regulation of cell differentiation"
)
path.miR.409.3p <- path.miR.409.3p %>%
  filter(Pathway %in% pathways)
#path.miR.409.3p <- path.miR.409.3p[-c(3,6,7,8,9,11,20),]
path.miR.409.3p <- path.miR.409.3p[-c(1,8,4,2),]
path.miR.409.3p$Pathway <- factor(path.miR.409.3p$Pathway, levels = path.miR.409.3p$Pathway)
a3=ggplot(path.miR.409.3p, aes(x=Pathway, y=-log10(P.value))) + geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95),
        plot.margin = margin(0,0,0,30,unit = "pt")) 
a3
ggsave("Path_miR.409.3p_2.pdf",a3, width = 4.5, height = 7)

custom_colors <- ifelse(1:nrow(path.miR.409.3p) %in% c(1,3,4,5,8,9,10,11), "red", "black")

# 绘图
c3 <- ggplot(path.miR.409.3p, aes(x=Pathway, y=-log10(P.value))) +
  geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95, color=custom_colors),
        plot.margin = margin(0,0,0,30,unit = "pt"))

c3
ggsave("Path_miR.409.3p_custom_text_colors.pdf", c3, width = 4, height = 7.5)

path.miR.28.3p <- read.csv("miRPathDB_targets/miRPathDB_miR.28.3p.csv")
pathways <-c(
  "ErbB Signaling Pathway", 
  "Non-small cell lung cancer", 
  "Oncostatin M Signaling Pathway", 
  "RAC1/PAK1/p38/MMP2 Pathway", 
  "MECP2 regulates neuronal receptors and channels", 
  "VEGFA-VEGFR2 Signaling Pathway", 
  "mTOR signaling pathway", 
  "Regulation of TP53 Activity", 
  "SUMOylation", 
  "FoxO signaling pathway", 
  "Thyroid hormone signaling pathway", 
  "Phosphatidylinositol signaling system", 
  "Muscle structure development", 
  "VEGFR2 mediated cell proliferation", 
  "Hippo signaling pathway")
path.miR.28.3p <- path.miR.28.3p %>%
  filter(Pathway %in% pathways)
path.miR.28.3p <- path.miR.28.3p[-c(6,8,12,16,19,20,21),]
path.miR.28.3p$Pathway <- factor(path.miR.28.3p$Pathway, levels = path.miR.28.3p$Pathway)
a4=ggplot(path.miR.28.3p, aes(x=Pathway, y=-log10(P.value))) + geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95),
        plot.margin = margin(0,0,0,15,unit = "pt")) 
a4
ggsave("Path_miR.28.3p.pdf",a4, width = 4.5, height = 6)

custom_colors <- ifelse(1:nrow(path.miR.28.3p) %in% c(1,3,4,5,8,9,10,11), "red", "black")

# 绘图 #没用到
c4 <- ggplot(path.miR.28.3p, aes(x=Pathway, y=-log10(P.value))) +
  geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95, color=custom_colors),
        plot.margin = margin(0,0,0,15,unit = "pt"))

c4
ggsave("Path_miR.28.3p_custom_text_colors.pdf", c4, width = 4.5, height = 6)


path.miR.28.5p <- read.csv("miRPathDB_targets/miRPathDB_miR.28.5p.csv")
pathways <-c("Cytokine signaling in immune system",
             "IL-9 Signaling Pathway",
             "IL-7 Signaling Pathway",
             "Chemokine signaling pathway",
             "Signaling by interleukins",
             "positive regulation of lymphocyte proliferation",#
             "positive regulation of mononuclear cell proliferation",#
             "leukocyte proliferation",#
             "lymphocyte proliferation",#
             "mononuclear cell proliferation",#
             "Positive regulation of activated t cell proliferation",
             "muscle structure development",#
             "muscle organ development",#
             "regulation of muscle organ development",#
             "positive regulation of muscle cell differentiation",
             "Muscle hypertrophy",
             "positive regulation of cytoskeleton organization",#
             "Non-small cell lung cancer",
             "HIF-1 signaling pathway",#
             "ErbB signaling pathway")#
path.miR.28.5p <- path.miR.28.5p %>%
  filter(Pathway %in% pathways)
path.miR.28.5p <- path.miR.28.5p[-c(16,27,19,15,21,25,26,20,5,14,18,23),]
path.miR.28.5p$Pathway <- factor(path.miR.28.5p$Pathway, levels = path.miR.28.5p$Pathway)
a5=ggplot(path.miR.28.5p, aes(x=Pathway, y=-log10(P.value))) + geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95),
        plot.margin = margin(0,0,0,15,unit = "pt")) 
a5
ggsave("Path_miR.28.5p.pdf",a5, width = 4.5, height = 6)

custom_colors <- ifelse(1:nrow(path.miR.28.5p) %in% c(1,3,4,7,8,9,10,11,12,13,14,15), "red", "black")

# 绘图 #没用到
c5 <- ggplot(path.miR.28.5p, aes(x=Pathway, y=-log10(P.value))) +
  geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95, color=custom_colors),
        plot.margin = margin(0,0,0,15,unit = "pt"))

c5
ggsave("Path_miR.28.5p_custom_text_colors.pdf", c5, width = 4.5, height = 6)


path.miR.146b.3p <- read.csv("miRPathDB_targets/miRPathDB_miR.146b.3p.csv")
pathways <-c(
  "TGF-beta Signaling Pathway",#
  "regulation of small GTPase mediated signal transduction",#
  "Melatonin metabolism and effects",#
  "cellular response to organic cyclic compound",#
  "regulation of muscle organ development",#
  "positive regulation of cardiac muscle tissue development",#
  "positive regulation of heart growth",#
  "Toll-like Receptor Signaling",#
  "MicroRNAs in cancer",
  "Hippo signaling pathway",
  "cellular response to steroid hormone stimulus",#
  "regulation of cellular localization"#
)
path.miR.146b.3p <- path.miR.146b.3p %>%
  filter(Pathway %in% pathways)
path.miR.146b.3p <- path.miR.146b.3p[-c(1,11),]
path.miR.146b.3p$Pathway <- factor(path.miR.146b.3p$Pathway, levels = path.miR.146b.3p$Pathway)
a6=ggplot(path.miR.146b.3p, aes(x=Pathway, y=-log10(P.value))) + geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95),
        plot.margin = margin(0,0,0,15,unit = "pt")) 
a6
ggsave("Path_miR.146b.3p.pdf",a6, width = 4.5, height = 6)

custom_colors <- ifelse(1:nrow(path.miR.146b.3p) %in% c(5,6,10), "red", "black")

# 绘图 
c6 <- ggplot(path.miR.146b.3p, aes(x=Pathway, y=-log10(P.value))) +
  geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95, color=custom_colors),
        plot.margin = margin(0,0,0,15,unit = "pt"))

c6
ggsave("Path_miR.146b.3p_custom_text_colors.pdf", c6, width = 4.5, height = 6)


path.miR.146b.5p <- read.csv("miRPathDB_targets/miRPathDB_miR.146b.5p.csv")
pathways <-c(
  "IL-1 signaling pathway",
  "Toll-like receptor signaling pathway",
  "Interleukin-1 Induced Activation of NF-kappa-B",
  "TLR4 Signaling and Tolerance",
  "Signal transduction through IL1R",
  "miRNAs involvement in the immune response in sepsis",
  "NF-kB is activated and signals survival",
  "ApoE and miR-146 in inflammation and atherosclerosis",
  "Toll-like receptor signaling pathway",
  "Interleukin-1 signaling",
  "Interleukin-1 family signaling",
  "Signaling by Interleukins",
  "Cytokine Signaling in Immune system",
  "NOD-like receptor signaling pathway",
  "T-Cell antigen Receptor (TCR) Signaling Pathway",
  "p75 NTR receptor-mediated signalling",
  "TRAF6 mediated induction of NFkB and MAP kinases upon TLR7/8 or 9 activation",
  "IRAK1 recruits IKK complex upon TLR7/8 or 9 stimulation"
)
path.miR.146b.5p <- path.miR.146b.5p %>%
  filter(Pathway %in% pathways)
path.miR.146b.5p <- path.miR.146b.5p[-c(12,31,14,19,3,20,23,30,8,10,32,7,26,4,18,27),]
path.miR.146b.5p$Pathway <- factor(path.miR.146b.5p$Pathway, levels = path.miR.146b.5p$Pathway)
a7=ggplot(path.miR.146b.5p, aes(x=Pathway, y=-log10(P.value))) + geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95),
        plot.margin = margin(0,0,0,25,unit = "pt")) 
a7
ggsave("Path_miR.146b.5p.pdf",a7, width = 5, height = 7.5)

custom_colors <- ifelse(1:nrow(path.miR.146b.5p) %in% c(1:11,13:16), "red", "black")

# 绘图 
c7 <- ggplot(path.miR.146b.5p, aes(x=Pathway, y=-log10(P.value))) +
  geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95, color=custom_colors),
        plot.margin = margin(0,0,0,25,unit = "pt"))

c7
ggsave("Path_miR.146b.5p_custom_text_colors.pdf", c7, width = 5, height = 7.5)


path.miR.181a.5p <- read.csv("miRPathDB_targets/miRPathDB_miR.181a.5p.csv")
pathways <-c(
  "MicroRNAs in cancer",#
  "IL-5 Signaling Pathway",#
  "IL-6 signaling pathway",#
  "Cellular response to organic substance",
  "Signaling Pathways in Glioblastoma",#
  "Negative regulation of MAPK pathway",#
  "Non-small cell lung cancer",#
  "Gene expression (Transcription)",#
  "RNA Polymerase II Transcription",#
  "positive regulation of cell population proliferation",#
  "cell population proliferation",#
  "cell development",#
  "regulation of cellular component movement",#
  "cell morphogenesis involved in differentiation",#
  "regulation of cell death"#
)
path.miR.181a.5p <- path.miR.181a.5p %>%
  filter(Pathway %in% pathways)
path.miR.181a.5p <- path.miR.181a.5p[c(6,12,13,14,15,16,23,24,27,28,29,31,32,33),]
path.miR.181a.5p$Pathway <- factor(path.miR.181a.5p$Pathway, levels = path.miR.181a.5p$Pathway)
a8=ggplot(path.miR.181a.5p, aes(x=Pathway, y=-log10(P.value))) + geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95),
        plot.margin = margin(0,0,0,0,unit = "pt")) 
a8
ggsave("Path_miR.181a.5p.pdf",a8, width = 4.5, height = 6)

custom_colors <- ifelse(1:nrow(path.miR.181a.5p) %in% c(3,4), "red", "black")

# 绘图 
c8 <- ggplot(path.miR.181a.5p, aes(x=Pathway, y=-log10(P.value))) +
  geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95, color=custom_colors),
        plot.margin = margin(0,0,0,0,unit = "pt"))

c8
ggsave("Path_miR.181a.5p_custom_text_colors.pdf", c8, width = 4.2, height = 6)


path.miR.181a.3p <- read.csv("miRPathDB_targets/miRPathDB_miR.181a.3p.csv")
pathways <-c("regulation of gene expression", 
             "regulation of transcription, DNA-templated", 
             "regulation of RNA metabolic process", 
             "Wnt signaling pathway", 
             "regulation of metabolic process", 
             "cell morphogenesis involved in differentiation", 
             "cellular metabolic process", 
             "multicellular organism development", 
             "vasculature development", 
             "tube development", 
             "blood vessel morphogenesis", 
             "cell surface receptor signaling pathway involved in cell-cell signaling", 
             "cell differentiation", 
             "animal organ morphogenesis")
path.miR.181a.3p <- path.miR.181a.3p %>%
  filter(Pathway %in% pathways)
path.miR.181a.3p$Pathway < factor(path.miR.181a.3p$Pathway, levels = path.miR.181a.3p$Pathway)
a9=ggplot(path.miR.181a.3p, aes(x=Pathway, y=-log10(P.value))) + geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95),
        plot.margin = margin(0,0,0,0,unit = "pt")) 
a9
ggsave("Path_miR.181a.3p.pdf",a9, width = 5, height = 7)

custom_colors <- ifelse(1:nrow(path.miR.181a.3p) %in% c(7), "red", "black")

# 绘图 
c9 <- ggplot(path.miR.181a.3p, aes(x=Pathway, y=-log10(P.value))) +
  geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95, color=custom_colors),
        plot.margin = margin(0,0,0,0,unit = "pt"))

c9
ggsave("Path_miR.181a.3p_custom_text_colors.pdf", c9, width = 5, height = 7)




### miRNA expression in human tissue
### PD-specifc, iRBD-specific
biomarks3 <- c("miR-590-3p", "miR-590-5p", "miR-28-3p", "miR-28-5p", "miR-409-3p", "miR-409-5p",
               "miR-146b-3p", "miR-146b-5p", "miR-181a-3p", "miR-181a-5p")

human.tissue.miRNA <- read.csv("miRNA_tissue_matrix_quantile.csv")
human.tissue.miRNA <- human.tissue.miRNA[,c(4,5,7,8)]
#human.tissue.miRNA$acc <- gsub("hsa-", "", human.tissue.miRNA$acc)
avg_expression <- human.tissue.miRNA %>%
  group_by(organ, acc) %>%
  summarise(mean_expression = mean(expression, na.rm = TRUE)) %>%
  arrange(organ, acc)
expression_matrix <- avg_expression %>%
  pivot_wider(names_from = acc, values_from = mean_expression)
normalized_matrix <- expression_matrix %>%
  mutate(across(-organ, ~ . / sum(.)))
write.table(normalized_matrix, file = "human.tissue.miRNA.normalized.expression.txt", sep="\t", quote = F)

normalized_matrix <- normalized_matrix[-c(7,18,23,25),]
human.tissue.miRNA.marker.melt <- reshape2::melt(normalized_matrix, ids.var="organ")

human.tissue.miRNA.marker2 <-  human.tissue.miRNA.marker.melt %>% 
  as_tibble() %>% 
  group_by(organ, variable) %>% 
  dplyr::summarise(value=mean(value))


#human.tissue.miRNA.mark <- human.tissue.miRNA[rownames(human.tissue.miRNA)%in% biomarks3 , ] # lack miR-3615, miR-4433b-3p
#human.tissue.miRNA.mark <- as.data.frame(t(human.tissue.miRNA.mark))
#human.tissue.miRNA.mark$group <- rownames(human.tissue.miRNA.mark)

x <- strsplit(rownames(human.tissue.miRNA.mark), "[.]")
tissueName <- NULL
for (i in x) {
 tissueName <- c(tissueName, i[1])  
}

human.tissue.miRNA.mark$group <- tissueName
#write.table(human.tissue.miRNA.mark, file = "human.tissue.miRNA.markers2.txt", sep="\t", quote = F)

#human.tissue.miRNA.marker <- read.csv("human.tissue.miRNA.markers2.txt", sep = "\t", header = T)
#human.tissue.miRNA.marker.melt <- reshape2::melt(human.tissue.miRNA.mark, ids.var="group")
#human.tissue.miRNA.marker2 <-  human.tissue.miRNA.marker.melt %>% as_tibble()  %>% 
#  group_by(group, variable) %>% dplyr::summarise(value=mean(value)) 

pdf("human.tissue.miRNA.biomarkers222.pdf", width = 10, height = 3.8)
ggplot(human.tissue.miRNA.marker2, aes(x=factor(organ), y=variable, size=value)) +
  #geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", )
  geom_point(color="#A50F15") + #
  scale_color_gradient(low = "#FEE5D9", high = "#A50F15") + 
  scale_fill_gradient(low = "#FEE5D9", high = "#A50F15") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 0.95)) +
  labs(y=NULL, x="Human tissue")

dev.off()



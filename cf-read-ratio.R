
################################
##### RNA type reads ratio #####
#####        wangge        #####
#####       20240910       #####
################################



setwd("/Users/wangge/Documents/DM/cfRNA/")

mt <- read.table("totalRNA.matrix.txt",sep="\t",header=T,check.names = F)
libsize <- colSums(mt)
rtype <- data.frame(libsize=libsize)
df <- read.table("totalRNA.matrix.type.txt",sep="\t",header=T,check.names = F)
identical(rownames(mt), df$ncrnaid)
mt$type2 <- df$type2
mt_sums <- aggregate(. ~ type2, data = mt, FUN = sum)
rownames(mt_sums) <- mt_sums$type2
mt_sums <- mt_sums[, -which(names(mt_sums) == "type2")]
mt_sums <- t(mt_sums)
mt_sums <- cbind(mt_sums,rtype)
mt_ratio <- sweep(mt_sums, 1, mt_sums$libsize, FUN = "/")
mt_ratio <- mt_ratio[, -which(names(mt_ratio) == "libsize")]
mt_ratio$type <- str_split(rownames(mt_ratio),"\\-",simplify=T)[,3]
mt_ratio$type[mt_ratio$type == "CF"] = "Plasma"
mt_ratio <- mt_ratio %>%
  rename(mRNA = protein_coding)

library(ggplot2)
library(ggsignif)
library(dplyr)

rna_type = "lncRNA"
compare_mt_ratio <- function(rna_type, mt_ratio, y_limits, y_breaks, y_labels) {
  rna_data <- mt_ratio[, c(rna_type, "type")]
  colnames(rna_data) <- c("value", "group")
  #rna_data$group <- gsub("CF", "Plasma", rna_data$group)
  rna_data$group <- factor(rna_data$group, levels = c("Plasma", "EV"))

  summary_stats <- rna_data %>%
    group_by(group) %>%
    summarize(mean_y = mean(value),
              q25 = quantile(value, 0.25),
              q75 = quantile(value, 0.75))
  

  bb = pal_nejm()(8)
  p = ggplot(rna_data, aes(x = group, y = value, fill = group)) +
    #stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(notch = FALSE, size = 0.4, outlier.shape = 16, outlier.size = 0.8, position = position_dodge(0.9)) +
    theme_bw() +
    scale_fill_manual(values = c(bb[2], bb[1])) +
    scale_x_discrete(labels = c("Plasma", "EV")) +
    theme_classic() +
    theme(panel.background = element_rect(fill = "white", colour = "black", size = 1),
          strip.background = element_blank(),
          strip.placement = "inside",
          strip.text = element_text(size = 12, color = "black"),
          strip.switch.pad.grid = unit(0, "cm"),
          strip.switch.pad.wrap = unit(0, "cm"),
          axis.line = element_blank(),
          axis.title = element_text(color = "black", size = 20),
          axis.text.y = element_text(color = "black", size = 20),
          axis.text.x = element_text(color = "black", size = 20),
          legend.position = "none",
          legend.title = element_text(color = "black", size = 20),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          legend.text = element_text(color = "black", size = 20)) +
    labs(title = rna_type, y = "Reads Ratio", x = "") +
    geom_signif(comparisons = list(c("EV", "Plasma")),
                map_signif_level = TRUE,
                vjust = 0.25)
  
  # 保存图像
  ggsave(paste0("plot/", rna_type, "-ratio.pdf"), p, width = 3.3, height = 4)
  
  return(p)
}

rna_types <- colnames(mt_ratio[,-15])

for (rna_type in rna_types) {
  print(paste("Plotting for:", rna_type))
  compare_mt_ratio(rna_type, mt_ratio, NULL, NULL, NULL)
}
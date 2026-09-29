# RStudio工作目录设为本文件夹
# 先运行severity_RNA.R、severity_pathway.R和severity_clinical_plot.R
library(data.table); library(ggplot2)
root <- file.path("generated", "figure6E")
dir.create(root,recursive=TRUE,showWarnings=FALSE)
out <- "generated"
gene <- fread(file.path(out,"transcript/results/Pooled_severity_participant_mean_predictions.tsv"))[
 comparison=="Pooled_severity" & validation=="repeated_nested_5fold"]
path <- fread(file.path(out,"clinical_context/participant_mean_oof_scores.tsv"))[
 model=="EV_RNA_pathway_only" & penalty_rule=="lambda_min"]
gm <- fread(file.path(out,"transcript/results/Fig2_Fig3_model_performance.tsv"))[
 comparison=="Pooled_severity" & validation=="repeated_nested_5fold"]
pm <- fread(file.path(out,"DM_severity_reviewer_completion_20260830/output/severity_pr_auc_auc_brier.tsv"))[
 model=="EV_RNA_pathway_only"]
stopifnot(nrow(gene)==87, nrow(path)==87, !anyDuplicated(gene$sample_id),
  !anyDuplicated(path$sample_id), nrow(gm)==1,nrow(pm)==1)
joined <- merge(gene[,.(sample_id,truth,transcript=score)],
  path[,.(sample_id,path_truth=truth,pathway=mean_oof_score)], by='sample_id')
stopifnot(nrow(joined)==87, all(joined$truth==joined$path_truth),sum(joined$truth)==39)
auc <- function(y,p) (sum(rank(p)[y==1])-sum(y==1)*(sum(y==1)+1)/2)/(sum(y==1)*sum(y==0))
stopifnot(abs(auc(gene$truth,gene$score)-gm$AUC)<1e-12,
          abs(auc(path$truth,path$mean_oof_score)-pm$AUC)<1e-12)
# Preserve the submitted ROC coordinate and line convention.
roc <- function(y,p,model) {
  y<-y[order(p,decreasing=TRUE)]
  data.table(FPR=c(0,cumsum(y==0)/sum(y==0),1),
             TPR=c(0,cumsum(y==1)/sum(y==1),1),Model=model)
}
d<-rbind(roc(gene$truth,gene$score,'Transcript'),roc(path$truth,path$mean_oof_score,'Pathway'))
labs <- c(sprintf('Transcript-level: %.3f (%.3f-%.3f)',gm$AUC,gm$AUC_low,gm$AUC_high),
          sprintf('Pathway-level: %.3f (%.3f-%.3f)',pm$AUC,pm$AUC_low,pm$AUC_high))
d[,Model:=factor(Model,levels=c('Transcript','Pathway'),labels=labs)]
p<-ggplot(d,aes(FPR,TPR,color=Model))+
  geom_abline(intercept=0,slope=1,linewidth=0.55,color='black')+
  geom_line(linewidth=1.05)+
  scale_color_manual(values=c('#777777','#377EB8'))+
  scale_x_continuous(breaks=seq(0,1,.25),labels=function(x)sprintf('%.2f',x))+
  scale_y_continuous(breaks=seq(0,1,.25),labels=function(x)sprintf('%.2f',x))+
  coord_equal(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+
  labs(x='False positive rate',y='True positive rate',color=NULL)+
  theme_bw(base_family='Helvetica',base_size=13)+
  theme(panel.grid=element_blank(),panel.border=element_rect(color='black',fill=NA,linewidth=.8),
    axis.text=element_text(color='black',size=16),axis.title=element_text(color='black',size=18),
    legend.position='top',legend.background=element_blank(),legend.key.width=grid::unit(6,'mm'),
    legend.text=element_text(size=12.5),legend.margin=margin(0,0,3,0),
    legend.location='plot',legend.justification='left',legend.box.just='left',
    plot.margin=margin(7,26,7,9))+
  guides(color=guide_legend(ncol=1,byrow=TRUE))
# Explicit physical dimensions; no changes to the other Figure 6 panels.
ggsave(file.path(root,'Figure6E_corrected_large_text_v4.pdf'),p,width=4.0,height=4.25,units='in',device=cairo_pdf,bg='white')
png(file.path(root,'Figure6E_corrected_large_text_v4.png'),width=4.0,height=4.25,units='in',res=300,bg='white',type='cairo')
print(p)
dev.off()
fwrite(d,file.path(root,'Figure6E_ROC_coordinates.tsv'),sep='\t')
fwrite(as.data.table(ggplot_build(p)$data[[2]])[,.(x,y,colour,group)],
       file.path(root,'Figure6E_rendered_line_coordinates.tsv'),sep='\t')
fwrite(rbind(data.table(model='Transcript-level',AUC=gm$AUC,CI_low=gm$AUC_low,CI_high=gm$AUC_high),
             data.table(model='Pathway-level',AUC=pm$AUC,CI_low=pm$AUC_low,CI_high=pm$AUC_high)),
       file.path(root,'Figure6E_statistics.tsv'),sep='\t')
writeLines(c('Figure 6E. Participant-level mean out-of-fold ROC curves for transcript-level and pathway-level severity models.',
 'Both models used 87 participants (39 moderate/severe and 48 mild). Legend values are AUC (95% CI).',
 'Pathway confidence limits match the final manuscript and repaired Supplementary Data 14.',
 'Helvetica; 4.0 x 4.25 inches; vector PDF and 300-dpi PNG.'),file.path(root,'README.txt'))
writeLines(capture.output(sessionInfo()),file.path(root,'sessionInfo.txt'))
print(labs)

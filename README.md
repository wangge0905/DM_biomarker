在RStudio中把工作目录设为本文件夹，分别打开R文件运行。
inputs是分析输入，folds是交叉验证分配，results是R2最终结果，供核对。
运行产生的文件放在generated；这里不附带运行中间结果。

三组biomarker分析
biomarker.R → biomarker_panel.R
最终RNA组合为4、2、6个。biomarker_permutation.R用于两患者亚型的置换检验，先完成上述两步。

Figure 3–5
ComBat_seq.R → differential_expression.R：A–C，描述性差异分析。
selected_gene_expression.R：D–F，先运行biomarker.R和biomarker_panel.R。
ComBat-seq调整计数只用于描述性差异图，不作为分类器输入。

Figure 2
cf-biotype.R：A，RNA类型构成。
cf-read-ratio.R：B，各类RNA的配对样本箱线图。
TE.R：D，Plasma/EV的TE比例箱线图，同时输出三组TE比例附图。
cf-differential.R → cf-differential-plot.R：E，配对差异RNA数量。

测序质控附图
QC_region.R：A，外显子、内含子和基因间区比例。
QC.R：B、D、E，clean/usable reads、rRNA比例及miRNA种类数。
RNA_biotype.R：C，三组RNA类型构成。

TE和miRNA附图
TE_heatmap.R：TE热图。
miRNA_organ.R：miRNA组织表达点图。

Figure 6和临床背景比较附图
ev-pathway-score.R：Hallmark评分，先运行。
severity_RNA.R：转录本模型，重复嵌套交叉验证及留一测序批次验证。
severity_pathway.R：抗体、临床、通路和联合模型。
severity_clinical_plot.R：ROC、PR、校准和决策曲线，接着severity_pathway.R运行。
severity_plot.R：Figure 6A、B、D和通路模型混淆矩阵，先完成以上四步。
severity_confusion.R：Figure 6转录本模型混淆矩阵，先运行severity_RNA.R。
severity_ROC.R：Figure 6E，先完成severity_RNA.R、severity_pathway.R和severity_clinical_plot.R。
severity_permutation.R → severity_permutation_plot.R：严重性通路模型置换检验及图，先完成通路评分。
severity_rank.R：fold-wise秩评分敏感性分析及图，先完成通路模型。

固定panel的协变量敏感性附图
age_sex_sensitivity.R：年龄/性别调整，先完成biomarker.R和biomarker_panel.R。
clinical_sensitivity.R → clinical_sensitivity_plot.R：患者亚型比较的临床协变量调整。

各小图分别保存。QC和RNA类型图另外保留原先的双面板导出，方便按原PPT裁剪。
原PPT中的裁剪、覆盖文字及拼接不在脚本内处理。
不包含Figure 2C细胞/组织构成饼图和富集分析图。
富集图的排除不包括Hallmark严重性评分及模型。

R 4.3.2。
主要包版本：data.table 1.15.2，edgeR 4.0.16，e1071 1.7-14，glmnet 4.1-9，
sva 3.35.2，ggplot2 3.5.0，pheatmap 1.0.12，ggpubr 0.6.0，ggsignif 0.6.4，
ggsci 3.0.1，ragg 1.2.7，patchwork 1.2.0，ggrepel 0.9.5。
作图使用Helvetica。跨系统字体替代可能改变文字排版。

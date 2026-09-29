# RStudio工作目录设为本文件夹
library(data.table);library(sva)
dm <- "inputs"
out <- file.path("generated", "descriptive_counts")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
raw <- fread(file.path(dm,"1204/gencode.txt"),data.table=FALSE,check.names=FALSE)
feature <- raw[[1]]
raw <- raw[,-1,drop=FALSE]
raw <- raw[,-c(40:44,50:54),drop=FALSE]
colnames(raw)[35:44] <- paste0(substr(colnames(raw)[35:44],1,12),substr(colnames(raw)[35:44],16,20))
rownames(raw) <- feature
meta <- read.delim(file.path(dm,"annotation.txt"),check.names=FALSE)
names(meta)[2] <- "library_id"
meta <- meta[-c(117,155:158),]
meta <- meta[order(factor(meta$group,levels=c("MDA5","ARS","HC"))),]
stopifnot(nrow(meta)==155L,all(meta$library_id %in% colnames(raw)))
counts <- as.matrix(raw[,meta$library_id,drop=FALSE])
counts <- counts[rowSums(counts)>1,,drop=FALSE]
adjusted <- ComBat_seq(counts,batch=substr(meta$library_id,1,6))
write.table(adjusted,file.path(out,"gencode_rmbatch.all.txt"),sep="\t",quote=FALSE,row.names=TRUE)
reference <- read.delim(file.path(dm,"gencode_rmbatch.all.txt"),row.names=1,check.names=FALSE)
stopifnot(identical(rownames(adjusted),rownames(reference)),setequal(colnames(adjusted),colnames(reference)))
delta <- max(abs(adjusted-as.matrix(reference[,colnames(adjusted)])))
fwrite(data.table(check="ComBat-seq adjusted counts",max_absolute_difference=delta,
 status=ifelse(delta==0,"PASS","FAIL")),file.path(out,"verification.tsv"),sep="\t")
if(delta!=0) stop("ComBat-seq count reproduction differs; see verification.tsv")

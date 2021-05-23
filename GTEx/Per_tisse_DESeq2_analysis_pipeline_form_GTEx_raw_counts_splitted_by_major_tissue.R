# set working directory to folder with input raw counts files
setwd("/home/danko_0001/projects/gb446/ANTXR2/reads_by_tissues/")
# get a list of files this directory
input_files <- list.files("/home/danko_0001/projects/gb446/ANTXR2/reads_by_tissues/", pattern = "[.]gct.gz")

library("DESeq2")
# loop for reading the files at the work directory and outputing normalized counts and fold-change data
for(i in 1:length(input_files)){
  cts <- as.matrix(read.csv(input_files[i],sep="\t",row.names="Name"))
  tis <- paste0(strsplit(input_files[i], "_")[[1]][1])
  print(tis)
  coldata <- read.csv("/home/danko_0001/projects/gb446/ANTXR2/DEseq2_AFR_EUR_meta_match_counts.csv", row.names=1)
  coldata <- coldata[colnames(cts),]
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ maternal)
  dds$maternal <- relevel(dds$maternal, ref="AFR")
  dds <- DESeq(dds)
  normalized_counts <- counts(dds, normalized=TRUE)
  write.table(normalized_counts, file = file.path(paste("/home/danko_0001/projects/gb446/ANTXR2/DESeq2_NC_per_tissue/",tis,".csv",sep="")), sep="\t", quote=F, col.names=NA)
  res <- results(dds)
  write.csv(res, file = file.path(paste("/home/danko_0001/projects/gb446/ANTXR2/DESeq2_DE_per_tissue/",tis,".csv",sep="")))}
 ################################################################################
require(DESeq2)

PVAL <- 0.01
FOLD <- 1

## Read count data.
source("/local/storage/projects/NHP/annotation/readData.R")


counts <- ca[,c(12:20)]
countsPI <- ca[,c(22:30)]
pc_genes  <- ca[,c(1:10)]

## Build experimental design matrix
sampleID <- c("Human 1", "Human 2", "Human 4", "Chimp 3", "Chimp 4", "Chimp 5", "R. Macaque 2", "R. Macaque 3", "R. Macaque 4")
prepDay  <- c(1, 3, 5, 2, 3, 4, 2, 3, 5)
subject  <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
species  <- c("H", "H", "H", "C", "C", "C", "M", "M", "M")
Design <- data.frame(sampleID, prepDay, subject, species)

## Fit the model 
fitModel <- function(species, cnts) {
  Design <- data.frame(sampleID, prepDay, subject, species)

  ## Estimate dispersions.
  dds <- DESeqDataSetFromMatrix(countData= cnts, colData= Design, design= ~ species)
  dds <- DESeq(dds, betaPrior=FALSE)
  dds <- replaceOutliersWithTrimmedMean(dds) ## NAs otherwise.
  resultsNames(dds)
  res <- results(dds, name="species_SS_vs_SO", independentFiltering = FALSE)
  ss <- data.frame(genes, res)

  return(ss)
}

## Treating alternatie primates as a single 'group'.
species <- c("SS", "SS", "SS", "SO", "SO", "SO", "SO", "SO", "SO")
hs <- fitModel(species, counts)
head(hs[order(hs$padj), ])

##GO analysis
library("biomaRt")

#get go terms for all protein coding genes
go_terms1=do.call("rbind", mclapply(1:NROW(pc_genes), function(i){
  getBM(attributes=c('hgnc_symbol', 'go_id', 'name_1006'), filters = 'hgnc_symbol', values=  pc_genes$mgi[i], mart= ensembl)
}, mc.cores=1))

#find genes labeled as 'integral component of membrane'
icom=which(go_terms1$name_1006=="integral component of membrane")
icom_genes=go_terms1$hgnc_symbol[icom] 


#see how icom genes overlap with significantly different genes
membrane=which(hs_pc$mgi %in% icom_genes)

significant=which(hs_pc$padj<PVAL)

sig_mem=intersect(membrane, significant)



#plot ma plot with icom genes highlighted
pdf("ma_icom_genes_sig_antxr2.pdf")
plot(hs_pc$baseMean, hs_pc$log2FoldChange, pch=19, log="x", col=rgb(red=0.66, green=0.66, blue=0.66, alpha=0.5))
points(hs_pc$baseMean[hs_pc$padj<PVAL], hs_pc$log2FoldChange[hs_pc$padj<PVAL], col=rgb(red=0,green=0,blue=0, alpha=0.5), pch=19)
points(hs_pc$baseMean[sig_mem], hs_pc$log2FoldChange[sig_mem], col=rgb(red=0.43, green=0.91, blue=1, alpha=0.6), pch=19)
points(hs_pc$baseMean[which(hs_pc$mgi=="ANTXR2")], hs_pc$log2FoldChange[which(hs_pc$mgi =="ANTXR2")], col=rgb(red=0, green=0, blue=1, alpha=0.6), pch=19)
dev.off()


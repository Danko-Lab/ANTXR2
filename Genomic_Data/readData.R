## 
## Read non-human primate count data, for dendrogram and expression analysis.
rowMax <- function(x) { sapply(1:NROW(x), function(i) {return(max(x[i,], na.rm=TRUE))}) }
rowMin <- function(x) { sapply(1:NROW(x), function(i) {return(min(x[i,], na.rm=TRUE))}) }

ca <- read.table("countall_rnaseq.tsv")
ca <- cbind(ca[,1:9], "gc18", ca[,10:NCOL(ca)])
gap <- read.table("genes.inGap")[!is.na(ca[,11]),] ## Read gap data.
ca <- ca[!is.na(ca[,11]),]

## Get pause counts.
ps <- read.table("countpause_rnaseq.tsv")
#ps[,7] <- paste(ps[,7],"_PauseSite", sep="")
ps <- cbind(ps[,1:6], "PauseSite", paste(ps[,4],"_PauseSite",sep=""), ps[,7], "ps", ps[,8:NCOL(ps)])
colnames(ps) <- colnames(ca)

## Get dREG counts.
ts <- read.table("counttss_rnaseq.tsv")
ts[,5] <- rowMax(ts[7:12])

## Rough classes...
stab <- rowMax(ts[,17:18])
dist <- ts[,13]

class <- rep("tss", NROW(ts)) ## tss is then unclassified as a promoter or enhancer
class[stab < 0.1 & dist < 500]  <- "Prox_Stab" ## Clearly protein coding promoter
class[stab > 0.1  & dist > 10000] <- "Dist_UnSt" ## Clearly distal enhancer
summary(as.factor(class))

ts <- cbind(ts[,1:6], class, ts[,c(4,20,19,21:NCOL(ts))])
#ts <- cbind(ts[,1:9], "tss", ts[,10:NCOL(ts)])
colnames(ts) <- colnames(ca)

## Join em
#ca <- rbind(ca, ps, ts)
ca <- rbind(ca, ts)  ## Don't use pause sites here!! Analyze these separately.

## Rename.
names(ca) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "type", "mgi", "mapSize", "annot_type",
                        "HF1", "HM1", "CF1", "CM1", "RF1", "RM1")
Condition <- as.factor(c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", 
                                                "U", "U", "U", "U", "U", "U", "U"))
Species  <- as.factor(c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", 
                                                "Human", "Human", "Chimp", "Chimp", "RMacaque", "RMacaque"))
Labels <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA",
                                                "Human Female1", "Human Male1", "Chimp Female1", "Chimp Male1", "Rhesus Female1", "Rhesus Male1")

## Cleanup...
ca <- ca[!is.na(ca[,11]),] ## Removes those which are not mappable/ orthologues in at least one species.
ca <- ca[grep("random", ca$chrom, invert=TRUE),]


## Used for getting useful subsets of the data.
#indx.all <- c(11:35)## ALL
#indx.unt <- c(11:20,31:35)## ONLY UNTREATED.
#indx.good <- c(11:35) ## "Good?!"  Remove M1-PI

# Get RPKM
rpkm_df <- as.matrix(ca) ## "Good?!"  Remove H2-U, H3-PI, C2-U+PI, M1-PI
for(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- 1000*(as.numeric(rpkm_df[,i])+0)/sum(as.numeric(rpkm_df[,i])) *1000/(as.numeric(ca[,"mapSize"]))


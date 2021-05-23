#!/bin/bash

# To get chromInfo, run the command below.  NOTE: MUST REMOVE HEADER LINE.
#twoBitInfo /gbdb/hg19/hg19.2bit chromInfo.hg19
#twoBitInfo /gbdb/rheMac3/rheMac3.2bit chromInfo.rheMac3
#twoBitInfo /gbdb/panTro4/panTro4.2bit chromInfo.panTro4


#files
chainfile=/home/danko_0001/projects/comparaReg/external_resources/rhesus/rhesus.human.rbest.chain.gz
bwaindex=/home/danko_0001/projects/comparaReg/external_resources/rhesus/bwa/index
chromInfo=/home/danko_0001/projects/comparaReg/external_resources/rhesus/chrominfo.txt

#proseq files


#rhesus macaque 1: SRX2008420 (SRR4012415, SRR4012416)
cat SRR4012415.fastq SRR4012416.fastq > RM1_proseq.fastq

#rhesus macaque 2: SRX2008421( SRR4012419, SRR4012418, SRR4012417)
fasterq-dump SRR4012419
fasterq-dump SRR4012418
fasterq-dump SRR4012417
cat SRR4012419.fastq SRR4012418.fastq SRR4012417.fastq > RM2_proseq.fastq

#rhesus macaque 3: SRX2008423 (SRR4012423,SRR4012424, SRR4012425)
fasterq-dump SRR4012423
fasterq-dump SRR4012424
fasterq-dump SRR4012425
cat SRR4012423.fastq SRR4012424.fastq SRR4012425.fastq > RM3_proseq.fastq


gzip *.fastq

export human_genome=/home/danko_0001/projects/comparaReg/external_resources/rhesus/bwa/index/genome
export human_chinfo=/home/danko_0001/projects/lac334/rhemac10/chromInfo.rhemac10

mkdir My_output_dir
bash proseq_mapper.sh -I RM1_proseq.fastq.gz -i $human_genome -c $human_chinfo -T /workdir/ -O ./My_output_dir/


export human_genome=/home/danko_0001/projects/comparaReg/external_resources/rhesus/bwa/index/genome
export human_chinfo=/home/danko_0001/projects/lac334/rhemac10/chromInfo.rhemac10
mkdir My_output_dir1
bash proseq_mapper.sh -I RM2_proseq.fastq.gz -i $human_genome -c $human_chinfo -T /workdir/ -O ./My_output_dir1/


#rhesus macaque 4 already on server


export human_genome=/home/danko_0001/projects/comparaReg/external_resources/rhesus/bwa/index/genome
export human_chinfo=/home/danko_0001/projects/lac334/rhemac10/chromInfo.rhemac10
mkdir My_output_dir2
bash proseq_mapper.sh -I /home/danko_0001/projects/RawSequenceFiles/2017-04-11-CD4_nhp4/M4-U.fastq.gz -i $human_genome -c $human_chinfo -T /workdir/ -O ./My_output_dir2/

/home/danko_0001/projects/RawSequenceFiles/2017-04-11-CD4_nhp4/M4-U.fastq.gz

#make bwa index for rhemac10

TMPDIR=./
BWAIDX=/home/danko_0001/projects/lac334/rhemac10/proseq/bwaindex/rhemac10
for fastq in `ls ./*.fastq.gz`
 do
  name=`echo $fastq | awk -F"/" '{print $NF}' | cut -d \. -f 1` ## | cut -d \/ -f 2`

  ## Align using BWA.
  bwa aln -t 32 ${BWAIDX} $fastq > ${TMPDIR}/$name.sai
  bwa samse -n 1 -f ${TMPDIR}/$name.sam ${BWAIDX} ${TMPDIR}/$name.sai $fastq
  samtools view -b -S ${TMPDIR}/$name.sam > ${TMPDIR}/$name.bam
  samtools sort -@ 32 ${TMPDIR}/$name.bam -o ${TMPDIR}/$name.sort.bam ## Writes the sorted BAM to the output folder. # add ".bam" for new version of samtools
  rm ${TMPDIR}/$name.bam ${TMPDIR}/$name.sam ${TMPDIR}/$name.sai
  cp ${TMPDIR}/$name.sort.bam ${OUTPUT} ## Saves the sorted BAM in the output file.  Make this optional?
 done


###CURRENTLY RUNNING ON ALL FILES####### #have bam files in rhemac10. need to convert to hg19.
CHINFO=/home/danko_0001/projects/lac334/rhemac10/chromInfo.rhemac10

FILES=`ls *.bam`
function converttohg19 {
for f in $FILES
 do 
   echo $f
   g=`echo $f | cut -d \. -f 1`
   echo $g
#make bed
   bedtools bamtobed -i $f > $g.sortedByCoord.bed
   gzip $g.sortedByCoord.bed

#make bedgraph
   ## Remove rRNA and reverse the strand (PRO-seq).
   zcat $f | grep "rRNA" -v | \
         awk 'BEGIN{OFS="\t"} {print $1,$2,($2+1),$4,$5,$6=="+"?"-":"+"}' | sort-bed - | gzip > $g.nr.rs.bed.gz
 
   ## Convert to bedGraph ... Can't gzip these, unfortunately.
   bedtools genomecov -bg -i $g.nr.rs.bed.gz -g $CHINFO | sort-bed - | gzip > $g.bedGraph.gz

#crossmap
   /home/danko_0001/projects/lac334/rhemac10/CrossMap-0.4.2/bin/CrossMap.py bed $MAPCHAIN $g.bedGraph.gz $g.hg19.bedGraph
   gzip $g.hg19.bedGraph

#convert to stranded bedgraph
   zcat $g.hg19.bedGraph.gz | grep -v "random" | grep "\+$" | awk 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{OFS="\t"} {print $1,$2,$3,abs($5)}' | sort-bed - > $g.plus.hg19.bedGraph
   zcat $g.hg19.bedGraph.gz | grep -v "random" | grep "\-$" | awk 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{OFS="\t"} {print $1,$2,$3,-1*abs($5)}' | sort-bed - > $g.plus.hg19.bedGraph

 done
}


FILES='RM1_proseq.sort.bam'
CHINFO=/home/danko_0001/projects/lac334/rhemac10/chromInfo.rhemac10
function converttohg19 {
for f in $FILES
 do 
   echo $f
   g=`echo $f | cut -d \. -f 1`
   echo $g
#make bed
   bedtools bamtobed -i $f > $g.sortedByCoord.bed
   gzip $g.sortedByCoord.bed

#make bedgraph
   ## Convert to bedGraph ... Can't gzip these, unfortunately.
   bedtools genomecov -bg -i $g.sortedByCoord.bed.gz -g $CHINFO | sort-bed - | gzip > $g.bedGraph.gz

#crossmap
   /home/danko_0001/projects/lac334/rhemac10/CrossMap-0.4.2/bin/CrossMap.py bed $MAPCHAIN $g.bedGraph.gz $g.hg19.bedGraph
   gzip $g.hg19.bedGraph

#convert to stranded bedgraph
   zcat $g.hg19.bedGraph.gz | grep -v "random" | grep "\+$" | awk 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{OFS="\t"} {print $1,$2,$3,abs($5)}' | sort-bed - > $g.plus.hg19.bedGraph
   zcat $g.hg19.bedGraph.gz | grep -v "random" | grep "\-$" | awk 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{OFS="\t"} {print $1,$2,$3,-1*abs($5)}' | sort-bed - > $g.plus.hg19.bedGraph

 done
}

FILES=`ls *.bam`
for f in $FILES
 do 
   echo $f
   g=`echo $f | cut -d \. -f 1`
   echo $g
bedtools bamtobed -i $f > $g.sortedByCoord.bed
done

FILES='RM2_proseq.sort.bam'

FILES='RM3_proseq.sort.bam'

FILES='RM4_proseq.sort.bam'
CHINFO=/home/danko_0001/projects/lac334/rhemac10/chromInfo.rhemac10
function converttohg19 {
for f in $FILES
 do 
   echo $f
   g=`echo $f | cut -d \. -f 1`
   echo $g
#make bed
   gzip $g.sortedByCoord.bed

#make bedgraph
   ## Convert to bedGraph ... Can't gzip these, unfortunately.
   bedtools genomecov -bg -i $g.sortedByCoord.bed.gz -g $CHINFO | sort-bed - | gzip > $g.bedGraph.gz

#crossmap
   /home/danko_0001/projects/lac334/rhemac10/CrossMap-0.4.2/bin/CrossMap.py bed $MAPCHAIN $g.bedGraph.gz $g.hg19.bedGraph
   gzip $g.hg19.bedGraph

#convert to stranded bedgraph
   zcat $g.hg19.bedGraph.gz | grep -v "random" | grep "\+$" | awk 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{OFS="\t"} {print $1,$2,$3,abs($5)}' | sort-bed - > $g.plus.hg19.bedGraph
   zcat $g.hg19.bedGraph.gz | grep -v "random" | grep "\-$" | awk 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{OFS="\t"} {print $1,$2,$3,-1*abs($5)}' | sort-bed - > $g.plus.hg19.bedGraph

 done
}






###REDO RM4 with barcode
#use proseq2.0.bsh
export PATH=$PATH:/programs/prinseq-lite-0.20.2:/programs:/home/zw355/lib/bin:/home/zw355/lib/ucsc

bash proseq2.0.bsh -i /home/danko_0001/projects/lac334/rhemac10/proseq/bwaindex/rhemac10 -c /home/danko_0001/projects/lac334/rhemac10/chromInfo.rhemac10 -SE -G -T myOutput1 -O myOutput1 --UMI1=6


function converttohg19 {
for f in $FILES
 do 
   echo $f
   g=`echo $f | cut -d \. -f 1`
   echo $g
#make bed

#make bedgraph
   ## Convert to bedGraph ... Can't gzip these, unfortunately.
   bedtools genomecov -bg -i $g.sortedByCoord.bed.gz -g $CHINFO | sort-bed - | gzip > $g.bedGraph.gz

#crossmap
g="RM4_proseq"
MAPCHAIN=/home/danko_0001/projects/lac334/rhemac10/rheMac10.hg19.rbest.chain.gz

   /home/danko_0001/projects/lac334/rhemac10/CrossMap-0.4.2/bin/CrossMap.py bed $MAPCHAIN $g.bedGraph.gz $g.hg19.bedGraph
   gzip $g.hg19.bedGraph

#convert to stranded bedgraph
   zcat $g.hg19.bedGraph.gz | grep -v "random" | grep "\+$" | awk 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{OFS="\t"} {print $1,$2,$3,abs($5)}' | sort-bed - > $g.plus.hg19.bedGraph
   zcat $g.hg19.bedGraph.gz | grep -v "random" | grep "\-$" | awk 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{OFS="\t"} {print $1,$2,$3,-1*abs($5)}' | sort-bed - > $g.minus.hg19.bedGraph

cat $g.plus.hg19.bedGraph | sed '/fix/d' > $g.plus.hg19.updated.bedGraph
cat $g.minus.hg19.bedGraph | sed '/fix/d' > $g.minus.hg19.updated.bedGraph


 done
}

#final steps after liftover


FILES=`ls *.hg19.bedGraph.gz`

FILES='RM2_proseq.hg19.bedGraph.gz'

function convertbw {
for f in $FILES
 do 
   echo $f
   g=`echo $f | cut -d \. -f 1`
   echo $g

   zcat $g.hg19.bedGraph.gz | grep -v "random" | grep "\+$" | awk 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{OFS="\t"} {print $1,$2,$3,abs($5)}' | sort-bed - > $g.plus.hg19.bedGraph
   zcat $g.hg19.bedGraph.gz | grep -v "random" | grep "\-$" | awk 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{OFS="\t"} {print $1,$2,$3,-1*abs($5)}' | sort-bed - > $g.minus.hg19.bedGraph

cat $g.plus.hg19.bedGraph | sed '/fix/d' > $g.plus.hg19.updated.bedGraph
cat $g.minus.hg19.bedGraph | sed '/fix/d' > $g.minus.hg19.updated.bedGraph


   bedGraphToBigWig $g.plus.hg19.updated.bedGraph /fs/cbsudanko/storage/projects/NHP/hg19.chromInfo $g.plus.hg19.updated.bw
   bedGraphToBigWig $g.minus.hg19.updated.bedGraph /fs/cbsudanko/storage/projects/NHP/hg19.chromInfo $g.minus.hg19.updated.bw

 done
}



#############rna-seq files
/home/danko_0001/projects/RawSequenceFiles/CD4_nhp_rnaseq/M4_U_.fastq.gz
/home/danko_0001/projects/RawSequenceFiles/CD4_nhp_rnaseq/M5_U_.fastq.gz

twoBitInfo rheMac10.2bit chromInfo.rhemac10

#align rna-seq

export PATH=/programs/STAR:$PATH

STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --genomeDir /home/danko_0001/projects/lac334/rhesusRNAseq/star/STAR-align --outFileNamePrefix M4_rhemac10 --readFilesIn /home/danko_0001/projects/RawSequenceFiles/CD4_nhp_rnaseq/M4_U_.fastq.gz
STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --genomeDir /home/danko_0001/projects/lac334/rhesusRNAseq/star/STAR-align --outFileNamePrefix M5_rhemac10 --readFilesIn /home/danko_0001/projects/RawSequenceFiles/CD4_nhp_rnaseq/M5_U_.fastq.gz

bedtools bamtobed -i M5_rhemac10Aligned.sortedByCoord.out.bam > M5_rhemac10Aligned.sortedByCoord.bed
bedtools bamtobed -i M4_rhemac10Aligned.sortedByCoord.out.bam > M4_rhemac10Aligned.sortedByCoord.bed

gzip M4_rhemac10Aligned.sortedByCoord.bed
gzip M5_rhemac10Aligned.sortedByCoord.bed


FILES=`ls M*.bed.gz`
CHINFO=/home/danko_0001/projects/lac334/rhemac10/chromInfo.rhemac10
makeBigWig


#download rhemac10/hg19 chain files
wget https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/vsHg19/reciprocalBest/hg19.rheMac10.rbest.chain.gz

MAPCHAIN=/home/danko_0001/projects/lac334/rhemac10/rheMac10.hg19.rbest.chain.gz
# Use rbest
FILES=`ls M*.bed.gz`
cat M4_rhemac10Aligned_plus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"n",$4,"+"}' > M4.tmp.bedGraph
cat M4_rhemac10Aligned_minus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"n",$4,"-"}' >> M4.tmp.bedGraph
cat M4.tmp.bedGraph | sort-bed - | gzip > M4.bedGraph.gz
/home/danko_0001/projects/lac334/rhemac10/CrossMap-0.4.2/bin/CrossMap.py bed $MAPCHAIN M4.bedGraph.gz M4.hg19.bedGraph


cat M5_rhemac10Aligned_plus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"n",$4,"+"}' > M5.tmp.bedGraph
cat M5_rhemac10Aligned_minus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"n",$4,"-"}' >> M5.tmp.bedGraph
cat M5.tmp.bedGraph | sort-bed - | gzip > M5.bedGraph.gz
/home/danko_0001/projects/lac334/rhemac10/CrossMap-0.4.2/bin/CrossMap.py bed $MAPCHAIN M5.bedGraph.gz M5.hg19.bedGraph


   zcat M4.hg19.bedGraph.gz | grep -v "random" | grep "\+$" | awk 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{OFS="\t"} {print $1,$2,$3,abs($5)}' | sort-bed - > M4_plus.hg19.bedGraph
   zcat M4.hg19.bedGraph.gz | grep -v "random" | grep "\-$" | awk 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{OFS="\t"} {print $1,$2,$3,-1*abs($5)}' | sort-bed - > M4_minus.hg19.bedGraph



   zcat M5.hg19.bedGraph.gz | grep -v "random" | grep "\+$" | awk 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{OFS="\t"} {print $1,$2,$3,abs($5)}' | sort-bed - > M5_plus.hg19.bedGraph
   zcat M5.hg19.bedGraph.gz | grep -v "random" | grep "\-$" | awk 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{OFS="\t"} {print $1,$2,$3,-1*abs($5)}' | sort-bed - > M5_minus.hg19.bedGraph

#remove chroms with fix in name. not found in hg19. won't work in next step
cat M4_plus.hg19.bedGraph | sed '/fix/d' > M4_plus.hg19.updated.bedGraph
cat M4_minus.hg19.bedGraph | sed '/fix/d' > M4_minus.hg19.updated.bedGraph



cat M5_plus.hg19.bedGraph | sed '/fix/d' > M5_plus.hg19.updated.bedGraph
cat M5_minus.hg19.bedGraph | sed '/fix/d' > M5_minus.hg19.updated.bedGraph




   bedGraphToBigWig M4_plus.hg19.updated.bedGraph /fs/cbsudanko/storage/projects/NHP/hg19.chromInfo M4_plus.updated.hg19.bw
   bedGraphToBigWig M4_minus.hg19.updated.bedGraph /fs/cbsudanko/storage/projects/NHP/hg19.chromInfo M4_minus.updated.hg19.bw



   bedGraphToBigWig M5_plus.hg19.updated.bedGraph /fs/cbsudanko/storage/projects/NHP/hg19.chromInfo M5_plus.updated.hg19.bw
   bedGraphToBigWig M5_minus.hg19.updated.bedGraph /fs/cbsudanko/storage/projects/NHP/hg19.chromInfo M5_minus.updated.hg19.bw




##########done. have bedgraph and bigwigs for m4 and m5 rna-seq based in rhemac10 and hg19.

#get counts
#soft link from //local/storage/projects/NHP/annotations
cp /fs/cbsudanko/storage/projects/NHP/annotations/pause.increments.bed ./
ln -s /fs/cbsudanko/storage/projects/NHP/annotations/pausesites.tsv
ln -s /fs/cbsudanko/storage/projects/NHP/annotations/genebodies.tsv
ln -s /fs/cbsudanko/storage/projects/NHP/annotations/gene.increments.bed
ln -s /fs/cbsudanko/storage/projects/NHP/annotations/gene.inGap
cp /home/danko_0001/projects/lac334/naive/sorted/getCounts.R ./
cp /home/danko_0001/projects/lac334/naive/sorted/bedmap-strand.bsh ./

export GENEBODIES=genebodies.tsv
export INCREMENTS=gene.increments.bed
export COUNTALL=countall.tsv

export PAUSESITES=pausesites.tsv
export PAUSEINCR=pause.increments.bed
export COUNTPAUSE=countpause.tsv

export TSSSITES=tsssites.tsv
export TSSINCR=tss.increments.bed
export COUNTTSS=counttss.tsv

function initCounts {
 bedmap --bases-uniq $1 $2 > olSize.sum
 paste $1 olSize.sum > $3
 rm olSize.sum
}

initCounts $GENEBODIES $INCREMENTS $COUNTALL
initCounts $PAUSESITES $PAUSEINCR  $COUNTPAUSE
initCounts $TSSSITES   $TSSINCR    $COUNTTSS

function getCounts {
 R --no-save --args $2 $COUNTFILE $BWPLUS $BWMINUS < getCounts.R
 bash bedmap-strand.bsh $1 $COUNTFILE | sed s/^.*\|//g | sed "s/NAN/NA/g" > $COUNTFILE.sum 
 paste $3 $COUNTFILE.sum > tmp
 mv tmp $3
 rm $COUNTFILE $COUNTFILE.sum
}

function getPauseCounts {
 R --no-save --args $1 $COUNTFILE $BWPLUS $BWMINUS < getCounts.R
 cat $COUNTFILE | awk '{print $5}' | sed "s/NAN/NA/g" > $COUNTFILE.tmp
 paste $2 $COUNTFILE.tmp > tmp
 mv tmp $2
 rm $COUNTFILE $COUNTFILE.tmp
}

#function getTSSCounts {
# ## Get plus counts.
# R --no-save --args $2 $COUNTFILE $BWPLUS $BWMINUS < getCounts.R
# bash bedmap-strand.bsh $1 $COUNTFILE | sed s/^.*\|//g | sed "s/NAN/NA/g" > $COUNTFILE.sum
#  
# ## Get minus counts.
#
#}

BWPATH=./

## Get counts, and place this in 'score'.
export COUNTFILE=M4.count.bed
export BWPLUS=M4_plus.updated.hg19.bw
export BWMINUS=M4_minus.updated.hg19.bw
getCounts $GENEBODIES $INCREMENTS $COUNTALL


export COUNTFILE=M5.count.bed
export BWPLUS=M5_plus.updated.hg19.bw
export BWMINUS=M5_minus.updated.hg19.bw
getCounts $GENEBODIES $INCREMENTS $COUNTALL


getCounts $PAUSESITES $PAUSEINCR $COUNTPAUSE
#getPauseCounts $PAUSESITES $COUNTPAUSE
getCounts $TSSSITES $TSSINCR $COUNTTSS

function makeBigWig {
for f in $FILES
 do 
   echo $f
   g=`echo $f | cut -d \. -f 1`

   ## Remove rRNA and reverse the strand (PRO-seq).
   zcat $f | grep "rRNA" -v | \
         awk 'BEGIN{OFS="\t"} {print $1,$2,($2+1),$4,$5,$6=="+"?"-":"+"}' | sort-bed - | gzip > $g.nr.rs.bed.gz
 
   ## Convert to bedGraph ... Can't gzip these, unfortunately.
   bedtools genomecov -bg -i $g.nr.rs.bed.gz -g $CHINFO -strand + | sort-bed - > $g\_plus.bedGraph
   bedtools genomecov -scale -1.0 -bg -i $g.nr.rs.bed.gz -g $CHINFO -strand - | sort-bed - > $g\_minus.bedGraph
 
   ## Then to bigWig
   bedGraphToBigWig $g\_plus.bedGraph $CHINFO $g\_plus.bw
   bedGraphToBigWig $g\_minus.bedGraph $CHINFO $g\_minus.bw

 done
}


function liftToHg19 {
 for i in $FILES
 do
   f=`echo $i | cut -d \. -f 1`
   echo $f

   ## Force each entry to represent precisely one base.
#   cat $f\_plus.bedGraph | perl ~/NHP/lib/bedGraph2Wig.pl | perl ~/NHP/lib/wig2bedGraph.pl + | less
#   cat $f\_minus.bedGraph | perl ~/NHP/lib/bedGraph2Wig.pl | perl ~/NHP/lib/wig2bedGraph.pl -  

   ## Convert to bedGraph
   cat $f\_plus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"n",$4,"+"}' > $f.tmp.bedGraph
   cat $f\_minus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"n",$4,"-"}' >> $f.tmp.bedGraph
   cat $f.tmp.bedGraph | sort-bed - | gzip > $f.bedGraph.gz

   ## Use liftOver.  Other options: -tab -minMatch=0.1
   #python2.7 /home/lac334/lib/CrossMap/usr/bin/CrossMap.py bed $MAPCHAIN $f.bedGraph.gz $f.hg19.bedGraph
   #python2.7 /home/lac334/crossmap/bin/CrossMap.py bed $MAPCHAIN $f.bedGraph.gz $f.hg19.bedGraph
   /home/lac334/lib/CrossMap/usr/bin/CrossMap.py bed $MAPCHAIN $f.bedGraph.gz $f.hg19.bedGraph
   #liftOver $f.bedGraph $MAPCHAIN $f.hg19.bedGraph unmap
   #rm $f.tmp.bedGraph.gz $f.bedGraph.gz
   gzip $f.hg19.bedGraph

   ## Split back into separate files.
   zcat $f.hg19.bedGraph.gz | grep -v "random" | grep "\+$" | awk 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{OFS="\t"} {print $1,$2,$3,abs($5)}' | sort-bed - > $f\_plus.hg19.bedGraph
   zcat $f.hg19.bedGraph.gz | grep -v "random" | grep "\-$" | awk 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{OFS="\t"} {print $1,$2,$3,-1*abs($5)}' | sort-bed - > $f\_minus.hg19.bedGraph

   ## Convert bedGarph to bigWig
   bedGraphToBigWig $f\_plus.hg19.bedGraph $CHINFOhg $f\_plus.hg19.bw
   bedGraphToBigWig $f\_minus.hg19.bedGraph $CHINFOhg $f\_minus.hg19.bw
 done
}



## Preclean
#rm *.wig.gz *.bw *.bedGraph *.bedGraph.gz *.nr.rs.* *tmp* tmp

UNMAPPED=tmp
CHINFOhg=/local/storage/projects/NHP/hg19.chromInfo

## Human
FILES='ls /home/danko_0001/projects/lac334/naive/sorted/H*_dedup*.bed.gz'
CHINFO=/local/storage/projects/NHP/hg19.chromInfo
makeBigWig

## Chimpanzee
FILES=`ls C*_dedup*.bam.bed.gz`
CHINFO=/local/storage/projects/NHP/panTro4.chromInfo
makeBigWig

MAPCHAIN=/local/storage/projects/NHP/makeRecipBest/hg19.panTro4/panTro4.hg19.rbest.chain.gz # Use rbest
FILES=`ls C*_dedup*.bam.bed.gz`
liftToHg19
#FILES=`ls C-PI.bed.gz`
#liftToHg19 &


## Rhesus Macaque
FILES=`ls R*_dedup*.bam.bed.gz`
CHINFO=/local/storage/projects/NHP/rheMac3.chromInfo
makeBigWig

MAPCHAIN=/local/storage/projects/NHP/makeRecipBest/hg19.rheMac3/rheMac3.hg19.rbest.chain.gz # Use rbest
FILES=`ls R*_dedup*.bam.bed.gz`
liftToHg19


## Cleanup.
#rm *.nr.rs.* tmp *.bedGraph *.bedGraph.gz

## Softlinks are required in just a few places.
#ln -s H-U_plus.bw H-U_plus.hg19.bw
#n -s H-U_minus.bw H-U_minus.hg19.bw
#ln -s H-PI_plus.bw H-PI_plus.hg19.bw
#sln -s H-PI_minus.bw H-PI_minus.hg19.bw



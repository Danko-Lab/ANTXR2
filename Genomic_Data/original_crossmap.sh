#!/bin/bash

# To get chromInfo, run the command below.  NOTE: MUST REMOVE HEADER LINE.
#twoBitInfo /gbdb/hg19/hg19.2bit chromInfo.hg19
#twoBitInfo /gbdb/rheMac3/rheMac3.2bit chromInfo.rheMac3
#twoBitInfo /gbdb/panTro4/panTro4.2bit chromInfo.panTro4


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



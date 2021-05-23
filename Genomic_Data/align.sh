#!/bin/sh
#SBATCH -J hg19.bwa.align # Job Name.
#SBATCH -o align.o%j      # Name of the output file (eg. myMPI.oJobID).
#SBATCH -e align.err%j    # Direct error to the error file.
#SBATCH -p normal         # Queue name.
#SBATCH -t 24:00:00       # Run time (hh:mm:ss) - 24 hours.
#SBATCH -N 1              # Requests 1 MPI node.
#SBATCH -n 16             # 3 tasks total. ... up to 16?!
#SBATCH -A TG-MCB130131   # Siepel lab account ID to bill.
#SBATCH --mail-user=dankoc@gmail.com
#SBATCH --mail-type=all

## Set some enviroment variables.
STARTDIR=$WORK/alignments/CD4/rnaseq/

## Set path.
export PATH=$PATH:$HOME/src/bwa/bwa-0.7.5a/:$HOME/src/samtools/samtools-0.1.19:$HOME/src/sratk/sratoolkit.2.3.2-5-ubuntu64/bin

## Move to the TMP directory.
cd /tmp

## Copy data files to the working directory.
cp $WORK/genomes/hg19/bwa/hg19* /tmp
cp $STARTDIR/*.sra /tmp

## Align reads.
for i in `ls *.sra`
do
  fastq-dump --gzip -Z -B ./$i > $i.fastq.gz   ## DESIGNED TO CONVERT FROM COLORSPACE!
  bwa aln -t 16 hg19 $i.fastq.gz > $SCRATCH/$i.sai              ## Align.
  bwa samse -n 1 -f $SCRATCH/$i.sam hg19 $SCRATCH/$i.sai $i.fastq.gz    ## Move to sam format.
  samtools view -b -S $SCRATCH/$i.sam > $SCRATCH/$i.bam
  samtools sort $SCRATCH/$i.bam $i.sort.bam
  rm $SCRATCH/$i.* ## Cleanup.
done

## Move result files from the working directory to my home directory.
cp *sort.bam.bam $STARTDIR


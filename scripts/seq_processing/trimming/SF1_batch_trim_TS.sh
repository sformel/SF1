#!/bin/bash

#SBATCH --job-name=SF1_batch_trim_TS
#SBATCH --output=job_reports/SF1_batch_trim_TS.output
#SBATCH --error=job_reports/SF1_batch_trim_TS.error
#SBATCH --qos=normal
#SBATCH --time=0-24:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

#batch trim all pairs of seqs

SOFTPATH='/lustre/project/svanbael/steve/software'
DATAPATH='/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data'
OUTPATH='/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out_TS'


cd $DATAPATH

for R1 in *_1.fq.gz
do
   R2=${R1//_1.fq.gz/_2.fq.gz}
   R1paired=${R1//.fq.gz/_paired.fq.gz}
   R1unpaired=${R1//.fq.gz/_unpaired.fq.gz}	
   R2paired=${R2//.fq.gz/_paired.fq.gz}
   R2unpaired=${R2//.fq.gz/_unpaired.fq.gz}	
   java -jar $SOFTPATH/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 20 -phred33 $R1 $R2 $OUTPATH/$R1paired $OUTPATH/$R1unpaired $OUTPATH/$R2paired $OUTPATH/$R2unpaired ILLUMINACLIP:$SOFTPATH/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25

done

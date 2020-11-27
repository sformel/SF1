#!/bin/bash

#SBATCH --job-name=SF1_test_QA
#SBATCH --output=job_reports/SF1_test_QA.output
#SBATCH --error=job_reports/SF1_test_QA.error
#SBATCH --qos=normal
#SBATCH --time=0-24:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

module load trinity/2.4.0
module load bowtie2/2.3.3 
module load samtools/1.5

F='/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R1/SF1_01_1_paired.fq.gz'
R='/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R2/SF1_01_2_paired.fq.gz'
O='/lustre/project/svanbael/steve/SF1/Trinity_test_out_SF1_26Nov2018'

bowtie2-build $O/Trinity.fasta Trinity.fasta
bowtie2 -p 10 -q --no-unal -k 20 -x Trinity.fasta -1 $F -2 $R  \
     2>align_stats.txt| samtools view -@10 -Sb -o bowtie2.bam
cat 2>&1 align_stats.txt

# sort the alignments by coordinate
samtools sort bowtie2.bam -o bowtie2.coordSorted.bam

# index the coordinate-sorted bam file
samtools index bowtie2.coordSorted.bam

# index the Trinity.fasta file
samtools faidx Trinity.fasta


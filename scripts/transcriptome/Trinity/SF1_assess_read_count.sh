#!/bin/bash
	
#SBATCH --qos=long
#SBATCH --job-name SF1_assess_read_count_Feb20
#SBATCH --error job_reports/SF1_assess_read_count_Feb20.error          
#SBATCH --output job_reports/SF1_assess_read_count_Feb20.output  
#SBATCH --time=5-23:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

module load trinity/2.8.4

#trinity output
TO='/lustre/project/svanbael/steve/SF1/Trinity_test_out_SF1_20Feb2019_test'

bowtie2-build $TO/Trinity.fasta $TO/Trinity.fasta

bowtie2 -p 20 -q --no-unal -k 20 -x /lustre/project/svanbael/steve/SF1/Trinity_test_out_SF1_20Feb2019_test/Trinity.fasta -1 /lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R1/SF1_01_1_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R1/SF1_04_1_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R1/SF1_10_1_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R1/SF1_13_1_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R1/SF1_17_1_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R1/SF1_19_1_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R1/SF1_26_1_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R1/SF1_30_1_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R1/SF1_34_1_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R1/SF1_37_1_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R1/SF1_40_1_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R1/SF1_44_1_paired.fq.gz -2 /lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R2/SF1_01_2_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R2/SF1_04_2_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R2/SF1_10_2_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R2/SF1_13_2_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R2/SF1_17_2_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R2/SF1_19_2_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R2/SF1_26_2_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R2/SF1_30_2_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R2/SF1_34_2_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R2/SF1_37_2_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R2/SF1_40_2_paired.fq.gz,/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R2/SF1_44_2_paired.fq.gz \
2>$TO/align_stats.txt| samtools view -@10 -Sb -o $TO/bowtie2.bam

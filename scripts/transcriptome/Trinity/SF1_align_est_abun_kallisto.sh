#!/bin/bash
	
#SBATCH --qos=long
#SBATCH --job-name SF1_align_est_abun_kallisto
#SBATCH --error job_reports/SF1_align_est_abun_kallisto.error          
#SBATCH --output job_reports/SF1_align_est_abun_kallisto.output  
#SBATCH --time=5-23:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

module load trinity/2.8.4

#trinity output
TO='/lustre/project/svanbael/steve/SF1/Trinity_test_out_SF1_20Feb2019_test'

#sample list
S='/lustre/project/svanbael/steve/SF1/sample_lists/SF1_Trinity_samples_list.txt'


$TRINITY_HOME/util/align_and_estimate_abundance.pl \
--transcripts $TO/Trinity.fasta \
--est_method kallisto \
--seqType fq \
--samples_file $S \
--output_dir $TO/est_abund_kallisto \
--SS_lib_type RF \
--thread_count 20 \
--gene_trans_map $TO/Trinity.fasta.gene_trans_map \
--prep_reference




#!/bin/bash
	
#SBATCH --qos=long
#SBATCH --job-name SF1_kallisto_abun_matrix
#SBATCH --error job_reports/SF1_kallisto_abun_matrix.error
#SBATCH --output job_reports/SF1_kallisto_abun_matrix.output
#SBATCH --time=5-23:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

module load trinity/2.8.4

#trinity output
TO='/lustre/project/svanbael/steve/SF1/Trinity_test_out_SF1_20Feb2019_test'

$TRINITY_HOME/util/abundance_estimates_to_matrix.pl \
--est_method kallisto \
--gene_trans_map $TO/Trinity.fasta.gene_trans_map \
--name_sample_by_basedir \
--quant_files est_abund_kallisto/kallisto_quant_files.txt \

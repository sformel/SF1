#!/bin/bash
	
#SBATCH --qos=long
#SBATCH --job-name SF1_trembl_BLAST
#SBATCH --error job_reports/SF1_trembl_BLAST.error          
#SBATCH --output job_reports/SF1_trembl_BLAST.output  
#SBATCH --time=6-23:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu
	
module load ncbi-blast/2.5.0+

#trinity output
TO='/lustre/project/svanbael/steve/SF1/Trinity_test_out_SF1_20Feb2019_test'

export BLASTDB='/lustre/project/svanbael/steve/SF1/BLAST_db'

blastx -query $TO/Trinity.fasta -db uniprot_trembl.fasta -out $TO/uniprot_trembl_blastx.outfmt6 -evalue 1e-20 -num_threads 20 -max_target_seqs 1 -outfmt 6

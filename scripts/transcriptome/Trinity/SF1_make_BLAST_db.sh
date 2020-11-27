#!/bin/bash
	
#SBATCH --qos=long
#SBATCH --job-name SF1_make_BLAST_db
#SBATCH --error job_reports/SF1_make_BLAST_db.error          
#SBATCH --output job_reports/SF1_make_BLAST_db.output  
#SBATCH --time=23:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu


module load ncbi-blast/2.5.0+
	
cd /lustre/project/svanbael/steve/SF1/BLAST_db

makeblastdb -in uniprot_trembl.fasta -dbtype prot 

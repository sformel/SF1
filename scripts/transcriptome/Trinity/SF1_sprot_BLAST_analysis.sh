#!/bin/bash
	
#SBATCH --qos=long
#SBATCH --job-name SF1_sprot_BLAST_analysis
#SBATCH --error job_reports/SF1_sprot_BLAST_analysis.error          
#SBATCH --output job_reports/SF1_sprot_BLAST_analysis.output  
#SBATCH --time=1-23:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu


module load ncbi-blast/2.5.0+
module load trinity/2.8.4
	
#trinity output
TO='/lustre/project/svanbael/steve/SF1/Trinity_test_out_SF1_20Feb2019_test'

export BLASTDB='/lustre/project/svanbael/steve/SF1/BLAST_db'

$TRINITY_HOME/util/analyze_blastPlus_topHit_coverage.pl $TO/uniprot_sprot_blastx.outfmt6 $TO/Trinity.fasta $BLASTDB/uniprot_sprot.fasta

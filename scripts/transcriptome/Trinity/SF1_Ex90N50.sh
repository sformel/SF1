#!/bin/bash
	
#SBATCH --qos=long
#SBATCH --job-name SF1_Ex90N50
#SBATCH --error job_reports/SF1_Ex90N50.error          
#SBATCH --output job_reports/SF1_Ex90N50.output  
#SBATCH --time=1-23:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu


module load trinity/2.8.4

#trinity output
TO='/lustre/project/svanbael/steve/SF1/Trinity_test_out_SF1_20Feb2019_test'

export BLASTDB='/lustre/project/svanbael/steve/SF1/BLAST_db'

$TRINITY_HOME/util/misc/contig_ExN50_statistic.pl \
     salmon.isoform.TMM.EXPR.matrix $TO/Trinity.fasta | tee $TO/ExN50.stats

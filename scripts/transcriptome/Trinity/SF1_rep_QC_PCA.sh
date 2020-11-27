#!/bin/bash
	
#SBATCH --qos=long
#SBATCH --job-name SF1_rep_QC_PCA
#SBATCH --error job_reports/SF1_rep_QC_PCA.error          
#SBATCH --output job_reports/SF1_rep_QC_PCA.output  
#SBATCH --time=1-23:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu


module load trinity/2.8.4

#trinity output
TO='/lustre/project/svanbael/steve/SF1/Trinity_test_out_SF1_20Feb2019_test'

$TRINITY_HOME/Analysis/DifferentialExpression/PtR \
    --matrix salmon.isoform.counts.matrix \
    -s sample_lists/SF1_samples_for_QC.txt --min_rowSums 10 --log2 \
    --CPM --center_rows \
    --prin_comp 3

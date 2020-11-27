#!/bin/bash
	
#SBATCH --qos=long
#SBATCH --job-name SF1_rep_QC_gene
#SBATCH --error job_reports/SF1_rep_QC_gene.error          
#SBATCH --output job_reports/SF1_rep_QC_gene.output  
#SBATCH --time=1-23:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu


module load trinity/2.8.4

#trinity output
TO='/lustre/project/svanbael/steve/SF1/Trinity_test_out_SF1_20Feb2019_test'

$TRINITY_HOME/Analysis/DifferentialExpression/PtR --matrix  salmon.gene.counts.matrix \
                  --samples sample_lists/SF1_samples_for_QC.txt --log2 \
                  --min_rowSums 10 \
                  --compare_replicates

$TRINITY_HOME/Analysis/DifferentialExpression/PtR \
          --matrix salmon.gene.counts.matrix \
          --min_rowSums 10 \
          -s sample_lists/SF1_samples_for_QC.txt --log2 --CPM --sample_cor_matrix

$TRINITY_HOME/Analysis/DifferentialExpression/PtR \
    --matrix salmon.gene.counts.matrix \
    -s sample_lists/SF1_samples_for_QC.txt --min_rowSums 10 --log2 \
    --CPM --center_rows \
    --prin_comp 3 
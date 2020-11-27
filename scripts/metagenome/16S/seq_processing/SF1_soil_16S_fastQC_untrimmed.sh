#!/bin/bash

#SBATCH --job-name=SF1_soil_16S_fastQC_untrimmed
#SBATCH --output=job_reports/SF1_soil_16S_fastQC_untrimmed.output
#SBATCH --error=job_reports/SF1_soil_16S_fastQC_untrimmed.error
#SBATCH --qos=normal
#SBATCH --time=0-24:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

#fastqc -help for options

DATAPATH='/lustre/project/svanbael/steve/SF1/soil_16S/H202SC19080294/raw_data'
mkdir $DATAPATH/fastqc_untrimmed

ls $DATAPATH/*.fq.gz | xargs -P 2 -n 1 fastqc -o $DATAPATH/fastqc_untrimmed

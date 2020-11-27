#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --qos=normal
#SBATCH -e job_reports/SF1_Trinity_26Nov2018.error            # File to which STDERR will be written
#SBATCH -o job_report/SF1_Trinity_26Nov2018.output           # File to which STDOUT will be written
#SBATCH -J SF1_Trinity_26Nov2018               # Job name
#SBATCH --mem=64000
#SBATCH --cpus-per-task=20
#SBATCH --time=1:00:00              # Runtime in D-HH:MM:SS
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=sformel@tulane.edu # Email to send notifications to

echo 'it worked'

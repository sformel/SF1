#Script for SF1 project

#Purpose: Generate descriptive stats of sequencing results
#Last updated by Steve FOrmel
#Jan 15, 2020

#script name
subtitle <- "Exploratory Analysis, Script: Seq_descriptive_stats.R"

#Import and clean data-----
source("scripts/metagenome/ITS/import_and_clean_SF1_ITS.R")

library(Rmisc)

df.seqs <- data.frame(sample_data(SF1), "Total Reads" = sample_sums(SF1))

summarySE(data = df.seqs , measurevar = "Total Reads", groupvars = c(""))
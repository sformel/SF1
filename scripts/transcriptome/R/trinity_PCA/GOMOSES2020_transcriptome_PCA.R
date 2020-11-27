#Make NMDS of Root Transcriptome for SF1 project
#script modified from Trinity script: salmon.isoform.counts.matrix.R
#modified for figure on poster for GOMOSES2020 conference

#by Steve Formel
#last updated: Jan 29, 2020

#metadata for figures
subtitle = "Exploratory Analysis, transcriptome_PCA.R"

library(cluster)
library(Biobase)
library(qvalue)
library(fastcluster)
library(plyr)
library(tidyverse)
library(phyloseq)

options(stringsAsFactors = FALSE)
NO_REUSE = F

#load in data----
primary_data = read.table("data/transcriptome/SF1_Trinity_24Mar2019_output/salmon_output/salmon.isoform.counts.matrix", header=T, com='', row.names=1, check.names=F, sep='\t')

primary_data = as.matrix(primary_data)

#load in Trinity R scripts
source("scripts/transcriptome/R/trinity_PCA/heatmap.3.R")
source("scripts/transcriptome/R/trinity_PCA/misc_rnaseq_funcs.R")
source("scripts/transcriptome/R/trinity_PCA/pairs3.R")
source("scripts/transcriptome/R/trinity_PCA/vioplot2.R")


data = primary_data

samples_data = read.csv("data/transcriptome/SF1_Trinity_24Mar2019_output/salmon_output/SF1_done_sample_sheet.csv", fileEncoding="UTF-8-BOM")

colnames(data) = samples_data$sample_ID


data = data[rowSums(data)>=10,]
initial_matrix = data # store before doing various data transformations
cs = colSums(data)
data = t( t(data)/cs) * 1e6;
data = log2(data+1)
data = as.matrix(data) # convert to matrix
 
# Centering rows
data = t(scale(t(data), scale=F))

prin_comp_data = data
pca = prcomp(prin_comp_data, center = FALSE, scale. = FALSE)
pc_pct_variance = (pca$sdev^2)/sum(pca$sdev^2)
def.par <- par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1));

PCA.loadings=pca$x
PCA.scores = pca$rotation

#reorder treatments
samples_data$Treatment <- interaction(samples_data$soil, samples_data$oil)

#rename

samples_data$Treatment <- revalue(samples_data$Treatment, c("L.N"="Live Soil , No Oil", "S.N"="Autoclaved Soil , No Oil" , "L.Y" = "Live Soil , Oiled", "S.Y" = "Autoclaved Soil , Oiled"))

#reorder levels
samples_data$Treatment <- factor(samples_data$Treatment, levels = c("Live Soil , No Oil", "Live Soil , Oiled" , "Autoclaved Soil , No Oil", "Autoclaved Soil , Oiled"))

#get summary data
pca.sum <- as.data.frame(summary(pca)$importance)
df <- cbind(PCA.scores,samples_data)

ggplot(data = df, 
                aes(x = PC1, y = PC2, 
                    fill = Treatment)) +
  geom_point(size = 6,
             shape = 21) +
  geom_polygon(alpha = 0.5) + 
  scale_fill_manual(values = c("#D55E00", "#56B4E9", "#F0E442","black")) +
  labs(title = "PCA on Root Transcriptome",
       x = paste("PC1 (", round(pca.sum$PC1[2]*100, digits = 1), "%)"),
       y = paste("PC1 (", round(pca.sum$PC2[2]*100, digits = 1), "%)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("figs/GOMOSES_2020/PCA_post.png", width = 10, height = 7, units = "in")

#PERMANOVA on transcriptome----
t.mat <- t(initial_matrix)
genes.df <- as.data.frame(t.mat)

library(vegan)

mod1 <- adonis(genes.df ~ soil*oil, data = samples_data, method = "bray") #took about 1 min at home

mod1

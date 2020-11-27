#Make NMDS of Root Transcriptome for SF1 project
#script modified from Trinity script: salmon.isoform.counts.matrix.R

#by Steve Formel
#last updated: Jan 21, 2020

#metadata for figures
subtitle = "Exploratory Analysis, transcriptome_PCA.R"

#library(parallelDist)
library(cluster)
library(Biobase)
library(qvalue)
library(fastcluster)
library(tidyverse)
library(phyloseq)
options(stringsAsFactors = FALSE)
NO_REUSE = F

#load in data----
primary_data = read.table("data/transcriptome/SF1_Trinity_24Mar2019_output/salmon_output/salmon.isoform.counts.matrix", header=T, com='', row.names=1, check.names=F, sep='\t')

primary_data = as.matrix(primary_data)

#load in Trinity R scripts
source("scripts/transcriptome/R/heatmap.3.R")
source("scripts/transcriptome/R/misc_rnaseq_funcs.R")
source("scripts/transcriptome/R/pairs3.R")
source("scripts/transcriptome/R/vioplot2.R")


data = primary_data

#myheatcol = colorpanel(75, 'purple','black','yellow')

samples_data = read.csv("data/transcriptome/SF1_Trinity_24Mar2019_output/salmon_output/SF1_done_sample_sheet.csv")

#samples_data = samples_data[samples_data[,2] != '',]
colnames(data) = colnames(samples_data)

#drop trt names from rep names
#samples_data[,2] <- substring(samples_data[,2], 14)
#make trt columns
#samples_data <- as.data.frame(samples_data)
#samples_data$soil <- str_sub(samples_data$sample_name, 6,6)
#samples_data$oil <- str_sub(samples_data$sample_name, 12,12)

#this still doesn't work:  error memory exhausted
#df.dist <- parDist(x = data,method = "bray", threads = 16)

#Back to PCA!

#sample_types = as.character(unique(samples_data[,1]))
#rep_names = as.character(samples_data[,2])
#data = data[, colnames(data) %in% rep_names, drop=F ]
#nsamples = length(sample_types)
#sample_colors = rainbow(nsamples)
#names(sample_colors) = sample_types
#sample_type_list = list()
#for (i in 1:nsamples) {
# samples_want = samples_data[samples_data[,1]==sample_types[i], 2]
#sample_type_list[[sample_types[i]]] = as.vector(samples_want)
#}
#sample_factoring = colnames(data)
#for (i in 1:nsamples) {
#  sample_type = sample_types[i]
#  replicates_want = sample_type_list[[sample_type]]
#  sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
#}

data = data[rowSums(data)>=10,]
initial_matrix = data # store before doing various data transformations
cs = colSums(data)
data = t( t(data)/cs) * 1e6;
data = log2(data+1)
#sample_factoring = colnames(data)
#for (i in 1:nsamples) {
#sample_type = sample_types[i]
#replicates_want = sample_type_list[[sample_type]]
#sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
#}
#sampleAnnotations = matrix(ncol=ncol(data),nrow=nsamples)
# for (i in 1:nsamples) {
#   sampleAnnotations[i,] = colnames(data) %in% sample_type_list[[sample_types[i]]]
# }
# sampleAnnotations = apply(sampleAnnotations, 1:2, function(x) as.logical(x))
# sampleAnnotations = sample_matrix_to_color_assignments(sampleAnnotations, col=sample_colors)
# rownames(sampleAnnotations) = as.vector(sample_types)
# colnames(sampleAnnotations) = colnames(data)
data = as.matrix(data) # convert to matrix
 
# Centering rows
data = t(scale(t(data), scale=F))

#write.table(data, file="salmon.isoform.counts.matrix.minRow10.CPM.log2.centered.dat", quote=F, sep='	');
#pdf("salmon.isoform.counts.matrix.minRow10.CPM.log2.centered.prcomp.principal_components_SF_mod.pdf")
prin_comp_data = data
pca = prcomp(prin_comp_data, center = FALSE, scale. = FALSE)
pc_pct_variance = (pca$sdev^2)/sum(pca$sdev^2)
def.par <- par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1));
#write.table(pca$rotation, file="salmon.isoform.counts.matrix.minRow10.CPM.log2.centered.PCA.prcomp.scores", quote=F, sep="	")
#write.table(pca$x, file="salmon.isoform.counts.matrix.minRow10.CPM.log2.centered.PCA.prcomp.loadings", quote=F, sep="	")
PCA.loadings=pca$x
PCA.scores = pca$rotation
#for (i in 1:(max(3,2)-1)) {
 # xrange = range(PCA.scores[,i])
  #yrange = range(PCA.scores[,i+1])
  #samples_want = rownames(PCA.scores) %in% sample_type_list[[sample_types[1]]]
  #pc_i_pct_var = sprintf("(%.2f%%)", pc_pct_variance[i]*100)
  #pc_i_1_pct_var = sprintf("(%.2f%%)", pc_pct_variance[i+1]*100)
  
#  plot(PCA.scores[samples_want,i], PCA.scores[samples_want,i+1], xlab=paste('PC',i, pc_i_pct_var), ylab=paste('PC',i+1, pc_i_1_pct_var), xlim=xrange, ylim=yrange, col=sample_colors[1])
#  text(PCA.scores[samples_want,i], PCA.scores[samples_want,i+1], labels=gsub(".*r","r",rownames(PCA.scores)[samples_want]), pos=1)
 # for (j in 2:nsamples) {
  #  samples_want = rownames(PCA.scores) %in% sample_type_list[[sample_types[j]]]
   # points(PCA.scores[samples_want,i], PCA.scores[samples_want,i+1], col=sample_colors[j], pch=j)
  #  text(PCA.scores[samples_want,i], PCA.scores[samples_want,i+1], labels=gsub(".*r","r",rownames(PCA.scores)[samples_want]), pos=1)
  #}
  #plot.new()
  #legend('topleft', as.vector(sample_types), col=sample_colors, pch=1:nsamples, ncol=2)
#}

#my version
#reorder treatments
samples_data$Treatment <- interaction(samples_data$soil, samples_data$oil)

#rename
library(plyr)
samples_data$Treatment <- revalue(samples_data$Treatment, c("L.N"="Live, No Oil", "S.N"="Sterile, No Oil", "L.Y" = "Live, Oiled", "S.Y" = "Sterile, Oiled"))

#reorder levels
samples_data$Treatment <- factor(samples_data$Treatment, levels = c("Live, No Oil", "Live, Oiled", "Sterile, No Oil", "Sterile, Oiled"))

df <- cbind(PCA.scores,samples_data)

ggplot(data = df, 
                aes(x = PC1, y = PC2, 
                    fill = Treatment, 
                    label = rep)) +
  geom_point(color = "black", 
             shape = 21,
             size = 4) +
  geom_text(hjust = 0, nudge_x = 0.005) +
  theme_bw() +
  scale_fill_manual(values = c("#D55E00", "#56B4E9", "#F0E442","black")) +
  labs(title = "SF1 Transcriptome PCA - post experiment", 
       subtitle = subtitle,
       caption = "seq counts were normalized to colsums and mult. by 1e6, then log2 transformed (+1 added for zeros.")

ggsave("figs/transcriptome/exploratory/PCA_post.png", width = 10, height = 7, units = "in")


#par(def.par)
#pcloadings_mat_vals = PCA.loadings[,1:3]
#print(dim(pcloadings_mat_vals))
#pcloadings_mat = matrix_to_color_assignments(pcloadings_mat_vals, col=colorpanel(256,'purple','black','yellow'), by='col')
#print(dim(pcloadings_mat))
#colnames(pcloadings_mat) = paste('PC', 1:ncol(pcloadings_mat))
#dev.off()
#gene_cor = NULL

#PERMANOVA on transcriptome----
t.mat <- t(initial_matrix)
genes.df <- as.data.frame(t.mat)

library(vegan)

mod1 <- adonis(genes.df ~ soil*oil, data = samples_data, method = "bray") #took about 1 min at home

mod1

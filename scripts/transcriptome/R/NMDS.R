#Make NMDS of Root Transcriptome for SF1 project
#script modified from Trinity script: salmon.isoform.counts.matrix.R

#by Steve Formel
#last updated: Sept 13, 2019

#library(parallelDist)
library(cluster)
library(Biobase)
library(qvalue)
library(fastcluster)
library(tidyverse)
library(phyloseq)
options(stringsAsFactors = FALSE)
NO_REUSE = F

#load in data
primary_data = read.table("Data/SF1_Trinity_24Mar2019_output/salmon_output/salmon.isoform.counts.matrix", header=T, com='', row.names=1, check.names=F, sep='\t')
primary_data = as.matrix(primary_data)

#load in Trinity R scripts
source("Data/SF1_Trinity_24Mar2019_output/salmon_output/heatmap.3.R")
source("Data/SF1_Trinity_24Mar2019_output/salmon_output/misc_rnaseq_funcs.R")
source("Data/SF1_Trinity_24Mar2019_output/salmon_output/pairs3.R")
source("Data/SF1_Trinity_24Mar2019_output/salmon_output/vioplot2.R")


data = primary_data
#myheatcol = colorpanel(75, 'purple','black','yellow')
samples_data = read.csv("Data/SF1_Trinity_24Mar2019_output/salmon_output/SF1_done_sample_sheet.csv", header=T, check.names=F, fill=T)
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
# 
# # Centering rows
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

trans <- ggplot(data = df, aes(x = PC1, y = PC2, fill = Treatment, label = rep)) +
  geom_polygon(color = "black", alpha = 0.5) +
  geom_text(hjust = 0, nudge_x = 0.005) +
  ggtitle("Root Transcriptome PCA") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9", "#F0E442","black" ))
  
#ggsave(file = "trans_PCA.png", width = 8, height = 6, unit = "in")


#par(def.par)
#pcloadings_mat_vals = PCA.loadings[,1:3]
#print(dim(pcloadings_mat_vals))
#pcloadings_mat = matrix_to_color_assignments(pcloadings_mat_vals, col=colorpanel(256,'purple','black','yellow'), by='col')
#print(dim(pcloadings_mat))
#colnames(pcloadings_mat) = paste('PC', 1:ncol(pcloadings_mat))
#dev.off()
#gene_cor = NULL

#add 16S soil NMDS
#make phyloseq object

#read in data
library(data.table)
seqtab.nochim <- readRDS("Data/Soil_DNA/soil_16S/16S_soil_otu.rds")
samdf <- readRDS("Data/Soil_DNA/soil_16S/16S_soil_env.rds")
taxa <- readRDS("Data/Soil_DNA/soil_16S/16S_soil_taxa.rds")

#add rep num
samdf$rep <- c(rep(c(1,2),4),rep(c(1,2,3),4))
#reorder time factor
samdf$Time <- factor(samdf$Time, levels = c("pre-exp", "post-exp"))
#rename
library(plyr)
samdf$Treatment <- interaction(samdf$Soil, samdf$Oil)
samdf$Treatment <- revalue(samdf$Treatment, c("L.N"="Live, No Oil", "S.N"="Sterile, No Oil", "L.Y" = "Live, Oiled", "S.Y" = "Sterile, Oiled"))

#reorder levels
samdf$Treatment <- factor(samdf$Treatment, levels = c("Live, No Oil", "Live, Oiled", "Sterile, No Oil", "Sterile, Oiled"))
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

#abbreviate ASV names
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#visualize diversity
rich <- plot_richness(ps, x="Oil", measures=c("Shannon"))

ggplot(data = rich$data, aes(x = Oil, y = value, fill = Treatment)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(~Time) +
  ggtitle("16S Soil Shannon Diversity") +theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9", "#F0E442","#0072B2" )) + 
  theme(plot.title = element_text(hjust = 0.5))

#ggsave("16S_soil_Shannon.png",width = 8,height = 8,units = "in")

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

#plot
soil_NMDS <- plot_ordination(ps.prop, ord.nmds.bray)

soil_NMDS_post <- ggplot(data = subset(soil_NMDS$data, soil_NMDS$data$Time!="pre-exp"), aes(x = NMDS1, y = NMDS2, fill = Treatment,label = rep)) +
  geom_polygon(color = "black", alpha = 0.5) + 
  #facet_grid(~Time) +
  geom_text(hjust = 0, nudge_x = 0.01, color = "black") +
  ggtitle("Soil 16S Community NMDS") +theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9", "#F0E442","#0072B2" )) + 
  scale_color_manual(values = c("#D55E00", "#56B4E9", "#F0E442","#0072B2" )) + 
  theme(plot.title = element_text(hjust = 0.5))

soil_NMDS_pre <- ggplot(data = subset(soil_NMDS$data, soil_NMDS$data$Time!="post-exp"), aes(x = NMDS1, y = NMDS2, fill = Treatment,label = rep)) +
  geom_point(aes(color = Treatment), size = 4) +
  geom_polygon(color = "black", alpha = 0.5) + 
  #facet_grid(~Time) +
  geom_text(hjust = 0, nudge_x = 0.01, color = "black") +
  ggtitle("Soil 16S Community at time = 0") +theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9", "#F0E442","#0072B2" )) + 
  scale_color_manual(values = c("#D55E00", "#56B4E9", "#F0E442","#0072B2" )) + 
  theme(plot.title = element_text(hjust = 0.5))

#plot with transcriptome
library(cowplot)

A <- trans + theme(legend.position = "none")
plot_grid(A, soil_NMDS_post)
#ggsave("trans_and_soil.png", width = 10, height = 8, units = "in")

B <- soil_NMDS_pre + theme(legend.position = "none")
plot_grid(B , soil_NMDS_post)

#ggsave("pre_and_post_soil.png", width = 10, height = 8, units = "in")

#import and combine all DE genes
filenames <- list.files(path = "Data/SF1_Trinity_24Mar2019_output/DEseq_output/DE_results/",full.names = TRUE, pattern = "*results")

DE_genes <- do.call("rbind", lapply(filenames, fread, data.table = FALSE))

#Group by main effects of treatments------
#DE_genes.Live <- DE_genes %>%
#  filter(sampleA=="soil_L_oil_N" & sampleB=="soil_L_oil_Y" & baseMeanA > baseMeanB & padj < 0.05)

DE_genes.Live <- DE_genes %>%
  filter(sampleA=="soil_L_oil_N" & sampleB=="soil_L_oil_Y" & padj<0.05)

DE_genes.Sterile <- DE_genes %>%
  filter(sampleA=="soil_S_oil_N" & sampleB=="soil_S_oil_Y" & padj<0.05)

Oil_genes <- intersect(DE_genes.Live$gene_name, DE_genes.Sterile$gene_name)

DE_genes.oil <- DE_genes %>%
  filter(sampleA=="soil_L_oil_Y" & sampleB=="soil_S_oil_Y" & padj<0.05)

DE_genes.no_oil <- DE_genes %>%
  filter(sampleA=="soil_L_oil_N" & sampleB=="soil_S_oil_N" & padj<0.05)

Soil_genes <- intersect(DE_genes.oil$gene_name, DE_genes.no_oil$gene_name)

#Make plot
#diff_genes <- c(Oil_genes,Soil_genes)
int_gene <- intersect(Oil_genes, Soil_genes)
DE_plot <- DE_genes[(DE_genes$sampleA=="soil_L_oil_Y" & DE_genes$sampleB=="soil_S_oil_Y") | (DE_genes$sampleA=="soil_L_oil_N" & DE_genes$sampleB=="soil_S_oil_N") | (DE_genes$sampleA=="soil_L_oil_N" & DE_genes$sampleB=="soil_L_oil_Y") | (DE_genes$sampleA=="soil_S_oil_N" & DE_genes$sampleB=="soil_S_oil_Y"),]

DE_plot <- DE_plot[DE_plot$gene_name %in% int_gene,]

DE_plot$Trt <- ifelse((DE_plot$sampleA=="soil_L_oil_Y" & DE_plot$sampleB=="soil_S_oil_Y") | (DE_plot$sampleA=="soil_L_oil_N" & DE_plot$sampleB=="soil_S_oil_N") , "Soil", "Oil")

#view by gene
DE_plot %>%
  select(gene_name, sampleA, sampleB, log2FoldChange, Trt) %>%
  group_split(gene_name)

DE_plot$COMP <- as.factor(paste(DE_plot$sampleA, DE_plot$sampleB, sep = "_"))

DE_plot$log2FoldChange <- round(DE_plot$log2FoldChange,2)
#This is silly but rught now this is the only way my brain can see to do this
plot.data <- DE_plot %>%
  select(gene_name,COMP, log2FoldChange) %>%
  spread(key = gene_name, value = log2FoldChange) %>%
  t()

colnames(plot.data) <- plot.data[1,]
plot.data <- as.data.frame(plot.data[-1,])
plot.data$gene_name <- row.names(plot.data)

#bar plots

ggplot(data = DE_plot, aes(x = gene_name, y = log2FoldChange, fill = COMP)) +
  geom_bar(stat = "identity", position = position_dodge())

#point plots

A <- ggplot(data = plot.data, aes(x = as.numeric(soil_L_oil_N_soil_L_oil_Y), y = as.numeric(soil_S_oil_N_soil_S_oil_Y), color = gene_name)) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  xlab(expression(atop('No Oil' %<->% Oil, paste("Live Soil")))) +
  ylab(expression(atop('No Oil' %<->% Oil, paste("Sterile Soil"))))


ggsave("Soil_Trt.png", width = 10, height = 8, units = "in")

B <- ggplot(data = plot.data, aes(x = as.numeric(soil_L_oil_N_soil_S_oil_N), y = as.numeric(soil_L_oil_Y_soil_S_oil_Y), color = gene_name)) +
  geom_point(size = 4) +
geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  xlab(expression(atop(Sterile %<->% Live, paste("No Oil")))) +
  ylab(expression(atop(Sterile %<->% Live, paste("Oil Added"))))

ggsave("Oil_Trt.png", width = 10, height = 8, units = "in")

#C <- get_legend(A)

plot_grid(A,B, nrow = 2)

#D <- plot_grid(A + theme(legend.position = "none"), B + theme(legend.position = "none"), nrow = 2)

#plot_grid(D, C, ncol = 2, rel_widths = c(1.5,0.5))

#filter to DE genes between LN and LY
DE.filtered <- DE_genes[DE_genes$padj<0.01,]

#make sure to refresh data at the top so it's not transformed
data1 = primary_data
#data1 = subset(data1, rownames(data1) %in% DE.filtered$gene_name)

# # Centering rows
data1 = t(scale(t(data1), scale=F))

prin_comp_data = data1
pca = prcomp(prin_comp_data, center = FALSE, scale. = FALSE)
pc_pct_variance = (pca$sdev^2)/sum(pca$sdev^2)
def.par <- par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1));
PCA.loadings=pca$x
PCA.scores = pca$rotation
#my version
#reorder treatments
samples_data$Treatment <- interaction(samples_data$soil, samples_data$oil)
#rename
samples_data$Treatment <- revalue(samples_data$Treatment, c("L.N"="Live, No Oil", "S.N"="Sterile, No Oil", "L.Y" = "Live, Oiled", "S.Y" = "Sterile, Oiled"))

#reorder levels
samples_data$Treatment <- factor(samples_data$Treatment, levels = c("Live, No Oil", "Live, Oiled", "Sterile, No Oil", "Sterile, Oiled"))
df <- cbind(PCA.scores,samples_data)

trans <- ggplot(data = df, aes(x = PC1, y = PC2, fill = Treatment, label = rep)) +
  geom_polygon(color = "black", alpha = 0.5) +
  geom_text(hjust = 0, nudge_x = 0.005) +
  ggtitle("Root Transcriptome PCA: filtered to DE genes, padj < 0.01") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9", "#F0E442","#0072B2" ))

ggsave("trans_DE_filt.png", width = 8, height = 7, units = "in")

#NMDS on bray curtis
data1 = subset(data, rownames(data) %in% DE.filtered$gene_name)

library(vegan)
DE.filt_NMDS <- metaMDS(comm = t(data1), distance = "bray")

NMDS.scores <- scores(DE.filt_NMDS)

NMDS.DE<- cbind(NMDS.scores, samples_data)

ggplot(data = NMDS.DE, aes(x = NMDS1, y = NMDS2, fill = Treatment, label = rep)) +
  geom_polygon(color = "black", alpha = 0.5) +
  geom_text(hjust = 0, nudge_x = 0.005) +
  ggtitle("Root Transcriptome NMDS: filtered to DE genes, padj < 0.01") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9", "#F0E442","#0072B2" ))

ggsave("trans_DE_filt_NMDS.png", width = 8, height = 7, units = "in")

#mantel test

#subset soil to remove pre-exp
library(vegan)
ps.post <- prune_samples(sample_data(ps)$Time=="post-exp", ps)
DE_filt_dist <- vegdist(x = t(data1), method = "bray")
MG_dist <- vegdist(otu_table(ps.post), method = "bray")

#on full genes
adonis(t(data1) ~ Oil*Soil, data = bac.env, method='eu')

bac.env <- data.frame(sample_data(ps.post))

adonis(MG_dist~Oil*Soil, data = bac.env)
adonis(DE_filt_dist~Oil*Soil, data = bac.env)

mantel(xdis = DE_filt_dist, ydis = MG_dist,method = "pearson")

# Call:
#   mantel(xdis = DE_filt_dist, ydis = MG_dist, method = "pearson") 
# 
# Mantel statistic r: 0.2836 
# Significance: 0.029 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
#   0.182 0.239 0.286 0.336 
# Permutation: free
# Number of permutations: 999

mantel(xdis = DE_filt_dist, ydis = MG_dist,method = method = "spearman")

#Mantel statistic based on Spearman's rank correlation rho 

# Call:
# mantel(xdis = DE_filt_dist, ydis = MG_dist, method = "spearman") 
# 
# Mantel statistic r: 0.2196 
#       Significance: 0.065 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.176 0.241 0.287 0.315 
# Permutation: free
# Number of permutations: 999

mantel(xdis = DE_filt_dist, ydis = MG_dist, method = "kendall")

# Mantel statistic based on Kendall's rank correlation tau 
# 
# Call:
# mantel(xdis = DE_filt_dist, ydis = MG_dist, method = "kendall") 
# 
# Mantel statistic r: 0.1534 
#       Significance: 0.057 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.115 0.157 0.191 0.234 
# Permutation: free
# Number of permutations: 999
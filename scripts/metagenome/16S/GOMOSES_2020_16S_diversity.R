#Begin to explore Fungal diversity for SF1 project
#Last updated by Steve FOrmel
#Jan 29, 2020

#modified for figure for GOMOSES 2020 poster

library(ggplot2)

#Import and clean data-----
source("scripts/metagenome/16S/import_and_clean_SF1_16S.R")

#make subtitle for metadata on all figures-----
#subtitle <- "Exploratory Analysis, Script: 16S_diversity.R"

#Alpha div-----
library(ggplot2)

B <- plot_richness(SF1,
              measures = c("Observed"),
              x = "time",
              color = "soil")

ggplot(data = B$data,
       aes(x = time,
           y = value,
           color = oil
           )) +
  geom_boxplot(aes(fill = oil), alpha = 0.2) +
  facet_grid(~soil) +
  theme_bw() +
  labs(title = "SF1 16S Richness", subtitle = subtitle)

ggsave("figs/metagenome/16S/exploratory/SF1_richness_all.png", width = 10, height = 7, units = "in")

#Shannon Div
B <- plot_richness(SF1,
                   measures = c("Shannon"),
                   x = "time",
                   color = "soil")

ggplot(data = B$data,
       aes(x = time,
           y = value,
           color = oil
       )) +
  geom_boxplot(aes(fill = oil), alpha = 0.2) +
  facet_grid(~soil) +
  theme_bw() +
  labs(title = "SF1 16S Shannon Diversity", subtitle = subtitle)

ggsave("figs/metagenome/16S/exploratory/SF1_shannon_all.png", width = 10, height = 7, units = "in")

#Beta div----

SF1.ord <- ordinate(SF1, method = "NMDS", distance = "bray")

B <- plot_ordination(SF1, SF1.ord, type="samples")

ggplot(data = B$data,
       aes(x = NMDS1*-1,
           y = NMDS2,
           fill = Treatment,
           shape = time)) +
  geom_point(size = 6) +
  geom_polygon(alpha = 0.5) +
  scale_shape_manual(values = c(21,24)) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9", "#F0E442","black")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "NMDS of Bacterial Communities (Bray-Curtis)",
       caption = paste("stress = ", round(SF1.ord$stress, digits = 2), sep = " "),
       shape = "Time")

ggsave("figs/GOMOSES_2020/16S_NMDS_all.png", width = 10, height = 7, units = "in")

#NMDS post-exp only

SF1.post <- prune_samples(sample_data(SF1)$time=="Post-Experiment", SF1)
SF1.ord <- ordinate(SF1.post, method = "NMDS", distance = "bray")

B <- plot_ordination(SF1.post, SF1.ord, type="samples", color="oil")

ggplot(data = B$data,
       aes(x = NMDS1,
           y = NMDS2,
           color = oil,
           shape = time)) +
  geom_point(aes(color = interaction(oil,soil)),
             size = 4) +
  theme_bw() +
  labs(title = "SF1 16S NMDS on Bray-Curtis - post experiment", 
       subtitle = subtitle,
       caption = paste("stress = ", round(SF1.ord$stress, digits = 2), sep = " "))

ggsave("figs/metagenome/16S/exploratory/NMDS_post_only.png", width = 10, height = 7, units = "in")

#PERMANOVA-----

#pre-experiment samples

ps <- prune_samples(sample_data(SF1)$time=="Pre-Experiment", SF1)

#Extract Data
library(vegan)

#extract otu table
otu <- otu_table(ps)
if (taxa_are_rows(otu)) {
  otu <- t(otu)
}
mat <- as(otu, "matrix")
df <- as.data.frame(mat)

#extract sample data
sam.sub <- data.frame(sample_data(ps))

#replace spaces in colnames
colnames(sam.sub) <- str_replace_all(colnames(sam.sub), " ", ".")

#can't add in oil data until I get the rest... (Jan 20, 2020)

mod1 <- adonis(df ~ soil*oil + pH + Calcium + Copper + Magnesium + Phosphorus + Potassium + Sodium + Sulfur + Zinc + Percent.Carbon + Percent.Nitrogen, 
               data = sam.sub,
               distance = "bray")

mod1  #nothing significant, but it's a ridiculous model with too many part

#hypothesis based model (i.e. things I think are sginificant)

mod2 <- adonis(df ~ soil*oil, 
               data = sam.sub,
               distance = "bray")

mod2  #soil (autoclaving) becomes significant if I get rid of all soil chemistry from the model.

#check for dispersion problems
dist.mat <- vegdist(df, method = "bray")

BD <- betadisper(dist.mat, group = sam.sub$soil)
anova(BD) #not significant


#post samples
ps <- prune_samples(sample_data(SF1)$time=="Post-Experiment", SF1)

#Extract Data
library(vegan)

#extract otu table
otu <- otu_table(ps)
if (taxa_are_rows(otu)) {
  otu <- t(otu)
}
mat <- as(otu, "matrix")
df <- as.data.frame(mat)

#extract sample data
sam.sub <- data.frame(sample_data(ps))

#replace spaces in colnames
colnames(sam.sub) <- str_replace_all(colnames(sam.sub), " ", ".")

#can't add in oil data until I get the rest... (Jan 20, 2020)

mod1 <- adonis(df ~ soil*oil + pH + Calcium + Copper + Magnesium + Phosphorus + Potassium + Sodium + Sulfur + Zinc + Percent.Carbon + Percent.Nitrogen, 
               data = sam.sub,
               distance = "bray")

mod1  #nothing significant, but it's a ridiculous model with too many part

#hypothesis based model (i.e. things I think are sginificant)

mod2 <- adonis(df ~ soil*oil + pH + Percent.Carbon + Percent.Nitrogen, 
               data = sam.sub,
               distance = "bray")

mod2  #soil, oil, and % nitrogen are sig, but not interaction of soil and oil

#check for dispersion problems
dist.mat <- vegdist(df, method = "bray")

BD <- betadisper(dist.mat, group = sam.sub$oil)
anova(BD) #not significant

BD <- betadisper(dist.mat, group = sam.sub$soil)
anova(BD)                                          ###significant#### - need to address

#post samples live only

ps <- prune_samples(sample_data(SF1)$time=="Post-Experiment" &
                      sample_data(SF1)$soil=="Live Soil", SF1)

#Extract Data
library(vegan)

#extract otu table
otu <- otu_table(ps)
if (taxa_are_rows(otu)) {
  otu <- t(otu)
}
mat <- as(otu, "matrix")
df <- as.data.frame(mat)

#extract sample data
sam.sub <- data.frame(sample_data(ps))

#replace spaces in colnames
colnames(sam.sub) <- str_replace_all(colnames(sam.sub), " ", ".")

#can't add in oil data until I get the rest... (Jan 20, 2020)

#hypothesis based model (i.e. things I think are sginificant)

mod1 <- adonis(df ~ oil + pH + Percent.Carbon + Percent.Nitrogen, 
               data = sam.sub,
               distance = "bray")

mod1  #soil, oil, and % nitrogen are sig, but not interaction of soil and oil

#check for dispersion problems
dist.mat <- vegdist(df, method = "bray")

BD <- betadisper(dist.mat, group = sam.sub$oil)
anova(BD) #not significant

#post samples autoclaved only

ps <- prune_samples(sample_data(SF1)$time=="Post-Experiment" &
                      sample_data(SF1)$soil=="Autoclaved Soil", SF1)

#Extract Data
library(vegan)

#extract otu table
otu <- otu_table(ps)
if (taxa_are_rows(otu)) {
  otu <- t(otu)
}
mat <- as(otu, "matrix")
df <- as.data.frame(mat)

#extract sample data
sam.sub <- data.frame(sample_data(ps))

#replace spaces in colnames
colnames(sam.sub) <- str_replace_all(colnames(sam.sub), " ", ".")

#can't add in oil data until I get the rest... (Jan 20, 2020)

#hypothesis based model (i.e. things I think are sginificant)

mod1 <- adonis(df ~ oil, 
               data = sam.sub,
               distance = "bray")

mod1  #soil, oil, and % nitrogen are sig, but not interaction of soil and oil

#check for dispersion problems
dist.mat <- vegdist(df, method = "bray")

BD <- betadisper(dist.mat, group = sam.sub$oil)
anova(BD) #not significant

#Specific taxa for oil----

oil.tax.list <- c("Alteromonadales", 
                  "Chromatiales", 
                  "Desulfobacterales",
                  "Desulfovibrionales", 
                  "Desulfuromonadales", 
                  "Enterobacteriales", 
                  "Methylococcales", 
                  "Oceanospirillales",
                  "Pseudomonadales", 
                  "Rhizobiales", 
                  "Rhodobacterales", 
                  "Thiotrichales")

#subset to post samples only
ps <- prune_samples(sample_data(SF1)$time=="Post-Experiment", SF1)
                    
#Subset to live_oiled and autoclaved_oiled samples
SF1.oiled <- prune_samples(sample_data(ps)$oil=="Oil Added",
                        ps)

#Subset to taxa that are in both live-unoiled and autoclave_unoiled
SF1.not_oiled <- prune_samples(sample_data(ps)$oil=="No Oil",
                        ps)
#prune empty taxa
SF1.oiled <- prune_taxa(taxa_sums(SF1.oiled)>0, SF1.oiled)
SF1.not_oiled <- prune_taxa(taxa_sums(SF1.not_oiled)>0, SF1.not_oiled)

#Get taxa that are in oiled but not in unoiled
SF1.oil.tax <- prune_taxa(setdiff(taxa_names(SF1.oiled), taxa_names(SF1.not_oiled)),
                              SF1.oiled)

#look at ordination of oil taxa

SF1.ord <- ordinate(SF1.oil.tax, method = "NMDS", distance = "bray")

B <- plot_ordination(SF1.oil.tax, SF1.ord, type="samples", color="oil")

ggplot(data = B$data,
       aes(x = NMDS1,
           y = NMDS2,
           color = oil)) +
  geom_point(aes(color = interaction(oil,soil)),
             size = 4) +
  theme_bw() +
  labs(title = "SF1 16S NMDS on Bray-Curtis - all samples", 
       subtitle = subtitle,
       caption = paste("stress = ", round(SF1.ord$stress, digits = 2), sep = " "))

ggsave("figs/metagenome/16S/NMDS_common_oil_taxa.png", width = 10, height = 7, units = "in")

#Look at specific taxa
SF1.oil.indicators <- subset_taxa(SF1.oil.tax, Order %in% oil.tax.list)

#agglomerate at order level
SF1.order <- tax_glom(SF1.oil.indicators, "Order")

#NMDS
SF1.ord <- ordinate(SF1.oil.indicators, method = "NMDS", distance = "bray")

B <- plot_ordination(SF1.oil.indicators, SF1.ord, type="samples", color="oil")

ggplot(data = B$data,
       aes(x = NMDS1,
           y = NMDS2,
           color = oil)) +
  geom_point(aes(color = interaction(oil,soil)),
             size = 4) +
  theme_bw() +
  labs(title = "SF1 16S NMDS on Bray-Curtis - all samples", 
       subtitle = subtitle,
       caption = paste("stress = ", round(SF1.ord$stress, digits = 2), sep = " "))

ggsave("figs/metagenome/16S/NMDS_oil_tax.indicators.png", width = 10, height = 7, units = "in")

#relative abundance of indicator taxa-----

#per https://joey711.github.io/phyloseq/preprocess.html

#transform
SF1.order.RA <- transform_sample_counts(SF1.order, function(OTU) OTU/sum(OTU) )

#Standardize abundances to the median sequencing depth

#total = median(sample_sums(SF1.order))
#standf = function(x, t=total) round(t * (x / sum(x)))
#SF1.order.std <-  transform_sample_counts(SF1.order, standf)

#Remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean & trivially large C.V.
#gpsf = filter_taxa(gps, function(x) sd(x)/mean(x) > 3.0, TRUE)

#plot 
plot_bar(SF1.order.RA,  "soil", "Abundance", "Order") +
  facet_wrap(~ Order)

ggsave("figs/metagenome/16S/barplot_oil_tax.indicators_in_oiled_samples.png", width = 10, height = 7, units = "in")





#But how do those orders look in unoiled samples?----

#Look at specific taxa in unoiled samples
SF1.oil.indicators <- subset_taxa(SF1.not_oiled, Order %in% oil.tax.list)

#agglomerate at order level
SF1.order <- tax_glom(SF1.oil.indicators, "Order")

#NMDS
SF1.ord <- ordinate(SF1.oil.indicators, method = "NMDS", distance = "bray")

B <- plot_ordination(SF1.oil.indicators, SF1.ord, type="samples", color="oil")

ggplot(data = B$data,
       aes(x = NMDS1,
           y = NMDS2,
           color = oil)) +
  geom_point(aes(color = interaction(oil,soil)),
             size = 4) +
  theme_bw() +
  labs(title = "SF1 16S NMDS on Bray-Curtis - all samples", 
       subtitle = subtitle,
       caption = paste("stress = ", round(SF1.ord$stress, digits = 2), sep = " "))

ggsave("figs/metagenome/16S/NMDS_oil_tax.indicators_in_unoiled_samples.png", width = 10, height = 7, units = "in")

#relative abundance of indicator taxa-----

#per https://joey711.github.io/phyloseq/preprocess.html

#transform
SF1.order.RA <- transform_sample_counts(SF1.order, function(OTU) OTU/sum(OTU) )

#Standardize abundances to the median sequencing depth

#total = median(sample_sums(SF1.order))
#standf = function(x, t=total) round(t * (x / sum(x)))
#SF1.order.std <-  transform_sample_counts(SF1.order, standf)

#Remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean & trivially large C.V.
#gpsf = filter_taxa(gps, function(x) sd(x)/mean(x) > 3.0, TRUE)

#plot 
plot_bar(SF1.order.RA,  "soil", "Abundance", "Order") +
  facet_wrap(~ oil*sampleID)

ggsave("figs/metagenome/16S/barplot_oil_tax.indicators_in_unoiled_samples.png", width = 10, height = 7, units = "in")

#But how do those orders look in all samples?----

#Look at specific taxa in unoiled samples
SF1.oil.indicators <- subset_taxa(ps, Order %in% oil.tax.list)

#agglomerate at order level
SF1.order <- tax_glom(SF1.oil.indicators, "Order")

#NMDS
SF1.ord <- ordinate(SF1.oil.indicators, method = "NMDS", distance = "bray")

B <- plot_ordination(SF1.oil.indicators, SF1.ord, type="samples", color="oil")

ggplot(data = B$data,
       aes(x = NMDS1,
           y = NMDS2,
           color = oil)) +
  geom_point(aes(color = interaction(oil,soil)),
             size = 4) +
  theme_bw() +
  labs(title = "SF1 16S NMDS on Bray-Curtis - all samples", 
       subtitle = subtitle,
       caption = paste("stress = ", round(SF1.ord$stress, digits = 2), sep = " "))

ggsave("figs/metagenome/16S/NMDS_oil_tax.indicators_in_all_samples.png", width = 10, height = 7, units = "in")

#relative abundance of indicator taxa-----

#per https://joey711.github.io/phyloseq/preprocess.html

#transform
SF1.order.RA <- transform_sample_counts(SF1.order, function(OTU) OTU/sum(OTU) )

#Standardize abundances to the median sequencing depth

#total = median(sample_sums(SF1.order))
#standf = function(x, t=total) round(t * (x / sum(x)))
#SF1.order.std <-  transform_sample_counts(SF1.order, standf)

#Remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean & trivially large C.V.
#gpsf = filter_taxa(gps, function(x) sd(x)/mean(x) > 3.0, TRUE)

#plot 
plot_bar(SF1.order.RA,  "soil", "Abundance", "Order") +
  facet_wrap(~ oil)


ggsave("figs/metagenome/16S/barplot_oil_tax.indicators_in_all_samples.png", width = 10, height = 7, units = "in")

#heatmap of taxa of interest
library(pheatmap)

#extract otu table
otu <- otu_table(ps)
if (taxa_are_rows(otu)) {
  otu <- t(otu)
}
mat <- as(otu, "matrix")
df <- as.data.frame(mat)

#extract sample data
sam.sub <- data.frame(sample_data(ps))

#replace spaces in colnames
colnames(sam.sub) <- str_replace_all(colnames(sam.sub), " ", ".")

#make color scheme for treatments
annotation_colors = list(
  oil = c(Oiled="black", Not_Oiled="gray"),
  soil = c(Autoclaved="burlywood4", Live="darkgreen"))

pheatmap(counts_scale.genes_that_go_up_with_autoclave,
         show_rownames=F,
         annotation_col= samdf[,2:3],
         annotation_colors = annotation_colors,
         main = "SF1 DEgenes when autoclave is present \n exploratory_heatmap.R",
         filename = "figs/transcriptome/exploratory/heatmap_autoclave_added_genes.png", width = 9, height = 7)

#look at oil taxa orders compared to all other taxa-----

#Look at specific taxa
SF1.oil.indicators <- subset_taxa(ps, Order %in% oil.tax.list)
SF1.others <- subset_taxa(ps, !(Order %in% oil.tax.list))

#agglomerate at order level
SF1.order.oil <- tax_glom(SF1.oil.indicators, "Order")
SF1.order.others <- tax_glom(SF1.others, "Kingdom")

#remove non-bacteria
SF1.order.others <- subset_taxa(SF1.order.others, Kingdom=="Bacteria")

#rename Order to "Others"
tt <- tax_table(SF1.order.others)

mat <- as(tt, "matrix")
df <- as.data.frame(mat)
df$Order <- factor("All Other Orders")

#put table back into ps object
tax_table(SF1.order.others) <- tax_table(as.matrix(df))

#merge with oil indicator taxa
SF1.all.orders <- merge_phyloseq(SF1.order.oil, SF1.order.others)

#NMDS
SF1.ord <- ordinate(SF1.all.orders, method = "NMDS", distance = "bray")

B <- plot_ordination(SF1.all.orders, SF1.ord, type="samples", color="oil")

ggplot(data = B$data,
       aes(x = NMDS1,
           y = NMDS2,
           color = oil)) +
  geom_point(aes(color = interaction(oil,soil)),
             size = 4) +
  theme_bw() +
  labs(title = "SF1 16S NMDS on Bray-Curtis - collapsed to order",
       caption = paste("stress = ", round(SF1.ord$stress, digits = 2), sep = " "))

ggsave("figs/GOMOSES_2020/NMDS_oil_tax.indicators_by_order.png", width = 10, height = 7, units = "in")

#relative abundance of indicator taxa-----

#per https://joey711.github.io/phyloseq/preprocess.html

#transform
SF1.order.RA <- transform_sample_counts(SF1.all.orders, function(OTU) OTU/sum(OTU) )

#Standardize abundances to the median sequencing depth

#total = median(sample_sums(SF1.order))
#standf = function(x, t=total) round(t * (x / sum(x)))
#SF1.order.std <-  transform_sample_counts(SF1.order, standf)

#Remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean & trivially large C.V.
#gpsf = filter_taxa(gps, function(x) sd(x)/mean(x) > 3.0, TRUE)

#add rep number for facetting
sample_data(SF1.order.RA)$rep <- rep(c("rep 1", "rep 2", "rep 3"), 4)

#plot 
plot_bar(SF1.order.RA,  "soil", "Abundance", "Order") +
  facet_wrap(~ oil, ncol = 3) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme_bw()

ggsave("figs/GOMOSES_2020/barplot_oil_tax.indicators_in_oiled_samples.png", width = 10, height = 7, units = "in")

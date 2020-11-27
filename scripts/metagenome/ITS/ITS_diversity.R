#Begin to explore Fungal diversity for SF1 project
#Last updated by Steve FOrmel
#Jan 20, 2020

#Import and clean data-----
source("scripts/metagenome/ITS/import_and_clean_SF1_ITS.R")

#make subtitle for metadata on all figures-----
subtitle <- "Exploratory Analysis, Script: ITS_diversity.R"

#load libraries----
library(vegan)
library(ggplot2)
library(tidyverse)

#Alpha div-----

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
  labs(title = "SF1 ITS Richness", subtitle = subtitle)

ggsave("figs/metagenome/ITS/exploratory/SF1_richness_all.png", width = 10, height = 7, units = "in")

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
  labs(title = "SF1 ITS Richness", subtitle = subtitle)

ggsave("figs/metagenome/ITS/exploratory/SF1_shannon_all.png", width = 10, height = 7, units = "in")

#Beta div----

SF1.ord <- ordinate(SF1, method = "NMDS", distance = "bray")

B <- plot_ordination(SF1, SF1.ord, type="samples", color="oil")

ggplot(data = B$data,
       aes(x = NMDS1,
           y = NMDS2,
           color = oil,
           shape = time)) +
  geom_point(aes(color = interaction(oil,soil)),
                 size = 4) +
  theme_bw() +
  labs(title = "SF1 ITS NMDS on Bray-Curtis - all samples", 
       subtitle = subtitle,
       caption = paste("stress = ", round(SF1.ord$stress, digits = 2), sep = " "))

ggsave("figs/metagenome/ITS/exploratory/NMDS_all.png", width = 10, height = 7, units = "in")

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
  labs(title = "SF1 ITS NMDS on Bray-Curtis - post experiment", 
       subtitle = subtitle,
       caption = paste("stress = ", round(SF1.ord$stress, digits = 2), sep = " "))

ggsave("figs/metagenome/ITS/exploratory/NMDS_post_only.png", width = 10, height = 7, units = "in")

#PERMANOVA-----

#pre-experiment samples

#post samples
ps <- prune_samples(sample_data(SF1)$time=="Pre-Experiment", SF1)

#Extract Data


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


#post samples-----
ps <- prune_samples(sample_data(SF1)$time=="Post-Experiment", SF1)

#Extract Data

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
anova(BD) #not significant


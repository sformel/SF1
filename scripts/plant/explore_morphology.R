#Script to explore morphology data collected for SF1 Project
#by Steve Formel Oct 19, 2018
#Last updated by Steve Formel Oct 19,2018

#Read in data and packages-------------
library(readxl)
library(dplyr)
library(ggplot2)

SF1_GC_data <- read_excel("Data/SF1_GC_data.xlsx", 
                                   sheet = "data", col_types = c("numeric", 
                                                                 "numeric", "text", "text", "date", 
                                                                 "text", "text", "numeric", "numeric", 
                                                                 "numeric", "numeric", "text"), na = "NA")

#Sum total of stemplant trait grouped by plant ID, date and total stems (for normalization)------------
stem.ht <- SF1_GC_data %>% 
  group_by(plantID, date, soil, oil, number_of_stems) %>% 
  summarise(tot_stem_ht = sum(stem_ht))

stem.diam <- SF1_GC_data %>% 
  group_by(plantID, date, soil, oil, number_of_stems) %>% 
  summarise(tot_stem_diam = sum(stem_diam))

num.leaves <- SF1_GC_data %>% 
  group_by(plantID, date, soil, oil, number_of_stems) %>% 
  summarise(tot_num_leaves = sum(num_leaves))

#plot stem trait over time by treatment, normalized to number of stems-----

#note that stem diam and height were problematic as sometimes they were unmeasurable if the were below the bucket lip.  As of Oct 19, 2018 I'm still deciding how to deal with these NA values so as to not throw out the other measurements.

#One possibility: #SF1_GC_data[is.na(SF1_GC_data)] <- 1

ggplot(data = stem.ht, aes(x = date, y = tot_stem_ht/number_of_stems, color = interaction(oil,soil))) +
  geom_jitter() +
  geom_smooth(method = lm, aes(linetype = soil)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#D55E00", "#56B4E9", "#F0E442","black" ))

ggplot(data = stem.diam, aes(x = date, y = tot_stem_diam/number_of_stems, color = interaction(oil,soil))) +
  geom_jitter() +
  geom_smooth(method = lm, aes(linetype = soil))

ggplot(data = num.leaves, aes(x = date, y = tot_num_leaves/number_of_stems, color = oil, shape = soil)) +
  geom_jitter() +
  geom_smooth(method = lm, aes(linetype = soil))

ggplot(data = num.leaves, aes(x = date, y = number_of_stems, color = oil, shape = soil)) +
  geom_jitter() +
  geom_smooth(method = lm, aes(linetype = soil))


#plot another way of thinking about biomass:  multiplying total stems by each of the traits-----
biomass_proxy <- SF1_GC_data %>% 
  group_by(plantID, date, soil, oil, number_of_stems, stem_diam, stem_ht, num_leaves) %>% 
  summarise(biomass_proxy = sum(num_leaves*stem_diam*stem_ht*number_of_stems))

ggplot(data = biomass_proxy, aes(x = date, y = biomass_proxy, color = oil, shape = soil)) +
  geom_jitter() +
  geom_smooth(method = lm, aes(linetype = soil))

#mixed effects model-----

library(nlme)

A <- lme(biomass_proxy ~ oil*soil, random = ~ 1|plantID, data = na.omit(biomass_proxy), method = "ML")

summary(A)
anova(A)


#post-hoc
library(multcomp)
summary(glht(A, linfct=mcp(month ="Tukey")))

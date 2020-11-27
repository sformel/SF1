#SF1 by Steve Formel
#Last updated: Oct 30, 2020

# Post-experiment plant Biomass to traits relationship for SF1 project

#load libraries----
library(readxl)
library(tidyverse)
library(ggplot2)
library(brms)

#load and clean data----
df.mass <- read_excel("data/plant/SF1_GC_data.xlsx", 
                          sheet = "biomass")

df.morph <- read_excel("data/plant/SF1_GC_data.xlsx", 
                          sheet = "morphology")

#convert to numeric
df.morph[,9:12] <- lapply(df.morph[,9:12], as.numeric)

#convert stem diam to cm
df.morph$stem_diam <- df.morph$stem_diam/10

#Make total mass by plantID
df.mass.total <- df.mass %>%
  group_by(plantID) %>%
  mutate(total_mass = sum(dry_mass)) %>%
  select(plantID,soil,oil,total_mass) %>%
  unique.data.frame()

#filter to final week
df.morph.lastweek <- df.morph %>%
  filter(week==13) %>%
  select(-notes) %>%
  na.omit() %>%
  mutate(total_traits = stem_ht*stem_diam*num_leaves) %>%
  group_by(plantID) %>%
  mutate(sum_total_traits = sum(total_traits)) %>%
  select(plantID,soil,oil,sum_total_traits) %>%
  unique.data.frame()
  

#compare traits from final week with total mass
plot(df.mass.total$total_mass, df.morph.lastweek$sum_total_traits)
m <- lm(df.mass.total$total_mass ~ df.morph.lastweek$sum_total_traits)
summary(m)

#compare traits from final week with above and belowground mass
plot(df.mass[df.mass$Tissue=="AG",]$dry_mass, df.morph.lastweek$sum_total_traits)
m <- lm(df.mass[df.mass$Tissue=="AG",]$dry_mass ~ df.morph.lastweek$sum_total_traits)
summary(m)

plot(df.mass[df.mass$Tissue=="IR",]$dry_mass, df.morph.lastweek$sum_total_traits)
m <- lm(df.mass[df.mass$Tissue=="IR",]$dry_mass ~ df.morph.lastweek$sum_total_traits)
summary(m)

#stems over time----
#filter to unique num stems
df.morph.stems <- df.morph %>%
  select(-notes) %>%
  na.omit() %>%
  group_by(plantID) %>%
  select(plantID,soil,oil,week,number_of_stems) %>%
  unique.data.frame()

M2<- brm(number_of_stems ~ oil + soil + unclass(week) + (1|plantID), 
         data = df.morph.stems, 
         family = "poisson",
         chains = 4, cores = 4,
         seed = 1)

M <- M2
pp_check(object = M, type = "dens_overlay", nsamples = 100)
pp_check(M, type = "stat", stat = 'median', nsamples = 100)
pp_check(M, type = "stat", stat = 'mean', nsamples = 100)
pp_check(M,type = 'intervals')

plot(M)
as.data.frame(posterior_summary(M))

#R2
as.data.frame(bayes_R2(object = M,resp = NULL,summary = TRUE,robust = FALSE,probs = c(0.025, 0.975)))

#This isn't a great model, but it doesn't lie.  Sterile and Oil increased stem counts, although that didn't translate to biomass necessarily.  This is also rate that I'm modeling, not total stem counts.
#SN has the most biomass, but that isn't evident from stem counts
df.morph.stems %>%
ggplot(aes(x = week,
           y = number_of_stems,
           color = interaction(soil,oil))) +
  geom_jitter() +
  facet_wrap(~soil*oil)

#What does it look like when it's stems*traits?

df.morph %>%
  select(-notes) %>%
  na.omit() %>%
  mutate(total_traits = stem_ht*stem_diam*num_leaves) %>%
  ggplot(aes(x = week,
           y = total_traits,
           color = interaction(soil,oil))) +
  geom_jitter() +
  geom_smooth(method = "lm", color = "black", formula=y~x-1 ) +
  facet_wrap(~soil*oil)

#What about stem diam?
df.morph %>%
  ggplot(aes(x = oil, 
             y = stem_diam, 
             fill = soil)) +
  geom_boxplot() +
  xlab("Oil Added")

#What about num leaves?
df.morph %>%
  ggplot(aes(x = oil, 
             y = num_leaves, 
             fill = soil)) +
  geom_boxplot() +
  xlab("Oil Added")

#What about stem height?
df.morph %>%
  ggplot(aes(x = oil, 
             y = stem_ht, 
             fill = soil)) +
  geom_boxplot() +
  xlab("Oil Added")
#Looks like a winner, let's see how it corresponds to biomass

df.morph.ht <- df.morph %>%
  filter(week==13) %>%
  select(-notes) %>%
  na.omit() %>%
  group_by(plantID) %>%
  mutate(max_ht = max(stem_ht)) %>%
  select(plantID,soil,oil,week,max_ht) %>%
  unique.data.frame()

#compare traits from final week with total mass
plot(df.mass.total$total_mass, df.morph.ht$max_ht)
m <- lm(df.mass.total$total_mass ~ df.morph.ht$max_ht)
summary(m)

#no dice.  But what about max height over time?

#filter to stem heights greater than 0
df.morph.ht <- df.morph %>%
  filter(stem_ht > 0) %>%
  select(-notes) %>%
  na.omit() %>%
  group_by(plantID) %>%
  select(plantID,soil,oil,week,stem_ht) %>%
  unique.data.frame()

M2 <- brm(stem_ht ~ oil + soil + unclass(week) + (1|plantID), 
         data = df.morph.ht, 
         family = "gamma",
         chains = 4, cores = 4,
         seed = 1)

M <- M2
pp_check(object = M, type = "dens_overlay", nsamples = 100)
pp_check(M, type = "stat", stat = 'median', nsamples = 100)
pp_check(M, type = "stat", stat = 'mean', nsamples = 100)
pp_check(M,type = 'intervals')

plot(M)
as.data.frame(posterior_summary(M))

#R2
as.data.frame(bayes_R2(object = M,resp = NULL,summary = TRUE,robust = FALSE,probs = c(0.025, 0.975)))


#stems*traits

#What does it look like when it's stems*traits?

df.morph %>%
  select(-notes) %>%
  na.omit() %>%
  mutate(total_traits = stem_ht*stem_diam*num_leaves) %>%
  ggplot(aes(x = week,
           y = total_traits,
           color = interaction(soil,oil))) +
  geom_jitter() +
  geom_smooth(method = "lm", color = "black", formula=y~x-1 ) +
  facet_wrap(~soil*oil)

df.alltraits <- df.morph %>%
  select(-notes) %>%
  na.omit() %>%
  mutate(total_traits = stem_ht*stem_diam*num_leaves) %>%
  select(plantID,soil,oil,week,total_traits)

#replace zeros with 1 (because the stem was there even if height or diam = 0)
df.alltraits$total_traits[df.alltraits$total_traits==0] <- 1

M2 <- brm(total_traits ~ oil + soil + unclass(week) + (1|plantID), 
         data = df.alltraits, 
         family = "gamma",
         chains = 4, cores = 4,
         seed = 1)

M <- M2
pp_check(object = M, type = "dens_overlay", nsamples = 1000)
pp_check(M, type = "stat", stat = 'median', nsamples = 1000)
pp_check(M, type = "stat", stat = 'mean', nsamples = 1000)
pp_check(M,type = 'intervals')

plot(M)
as.data.frame(posterior_summary(M))

#R2
as.data.frame(bayes_R2(object = M,resp = NULL,summary = TRUE,robust = FALSE,probs = c(0.025, 0.975)))

#But the above doesn't really account for the repeated measures apsect of the individual stems, so it's probalby best to sum the info by plant

# df.alltraits <- df.morph %>%
#   select(-notes) %>%
#   na.omit() %>%
#   mutate(total_traits = stem_ht*stem_diam*num_leaves) %>%
#   select(plantID,soil,oil,week,total_traits) %>%
#   group_by(plantID) %>%
#   mutate(sum_total_traits = sum(total_traits)) %>%
#   select(plantID,soil,oil,week,sum_total_traits) %>%
#   unique.data.frame()
# 
# M2 <- brm(sum_total_traits ~ oil + soil + unclass(week) + (1|plantID), 
#          data = df.alltraits, 
#          family = "gaussian",
#          chains = 4, cores = 4,
#          seed = 1, 
#          iter = 10000, 
#          control = list(max_treedepth = 15))
# 
# M <- M2
# pp_check(object = M, type = "dens_overlay", nsamples = 100)
# pp_check(M, type = "stat", stat = 'median', nsamples = 100)
# pp_check(M, type = "stat", stat = 'mean', nsamples = 100)
# pp_check(M,type = 'intervals')
# 
# plot(M)
# as.data.frame(posterior_summary(M))
# 
# #R2
# as.data.frame(bayes_R2(object = M,resp = NULL,summary = TRUE,robust = FALSE,probs = c(0.025, 0.975)))

#I couldn't get these models to converge with reasonable parameters
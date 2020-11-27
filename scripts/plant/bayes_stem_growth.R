#SF1 by Steve Formel
#Last updated: Oct 24, 2019

#Mixed effects model to look at relationship between oil, and morphology.  Built off S5 model
library(brms)
library(readxl)
library(tidyverse)
library(tidybayes)
library(ggplot2)


plant.morph <- read_excel("Data/SF1_GC_data.xlsx", 
                          sheet = "data")
plant.morph$week <- as.factor(plant.morph$week)

plant.morph$number_of_stems <- as.integer(plant.morph$number_of_stems)
plant.morph$num_leaves <- as.integer(plant.morph$num_leaves)

hist(plant.morph$num_leaves)
hist(plant.morph$number_of_stems)

#both are probably negative binomial distribution

#beginning Full model for selection 

#number of stems ~ oil*soil + (1|plantID) + (1|block)
#anova on tot_stems ~ block is barely significant - but seems to be making the model struggle, if I remove it the models run successfully.  With block, I could not get the model to converge, no mtter how high the adapt_delta parameter.

M1<- brm(number_of_stems ~ unclass(week) + (1|plantID), data = na.omit(stems), save_all_pars = TRUE, control = list(adapt_delta = 0.99, max_treedepth = 15), family = "negbinomial")

M1.poisson <- brm(number_of_stems ~ unclass(week) + (1|plantID), data = na.omit(stems), save_all_pars = TRUE, control = list(adapt_delta = 0.99, max_treedepth = 15), family = "poisson")

M2<- brm(number_of_stems ~ oil + soil + unclass(week) + (1|plantID), data = na.omit(stems), save_all_pars = TRUE, control = list(adapt_delta = 0.99, max_treedepth = 15), family = "poisson")

M3<- brm(number_of_stems ~ oil*soil + unclass(week) + (1|plantID), data = na.omit(stems), save_all_pars = TRUE, control = list(adapt_delta = 0.99, max_treedepth = 15), family = "poisson")

saveRDS(M1, "SF1_M1.rds")
saveRDS(M1.poisson, "SF1_M1poisson.rds")
saveRDS(M2, "SF1_M2.rds")
saveRDS(M3, "SF1_M3.rds")

loo_all <- loo(M1,M2,M3, reloo = TRUE) #null without treatments is the best option :-(
loo_poisson <- loo(M1, M1.poisson) #oops, maybe a poisson is a better fit?

loo_all_poisson <- loo(M1.poisson, M2, M3) #With a poisson fit, M2 is no worse than null.

plot(M2, pars = c("oil", "soil"))

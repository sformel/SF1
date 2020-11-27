#SF1 by Steve Formel
#Last updated: Oct 30, 2020

# Post-experiment plant Biomass for SF1 project

#load libraries----
library(readxl)
library(tidyverse)
library(ggplot2)

#load and clean data----
df <- read_excel("data/plant/SF1_GC_data.xlsx", 
                          sheet = "biomass")

#Make total mass
df <- df %>%
  group_by(plantID) %>%
  mutate(total_mass = sum(dry_mass))

#Exploratory plots----

#All tissue types
df %>%
  ggplot(aes(x = oil, 
             y = dry_mass, 
             fill = Tissue)) +
  geom_boxplot() +
  facet_grid(~ ifelse(soil=="L", "Live Soil", "Sterile Soil")) +
  xlab("Oil Added")

#Total Biomass
df %>%
  ggplot(aes(x = oil, 
             y = total_mass, 
             fill = soil)) +
  geom_boxplot() +
  xlab("Oil Added")

#Oil might have affected sterile soil, but did not effect live.

m <- lm(total_mass ~ soil*oil, data = df)
summary(m)
hist(resid(m))
shapiro.test(resid(m))

#Another not quite normal data set.  Awesome.  Let's see if we can fit it with Bayes.

library(brms)
library(modelr)

b0 <- brm(total_mass ~ 1 + soil*oil, 
          data = df,
          family = "skew_normal",
          chains = 4, cores = 4,
          seed = 1)

#https://bayesed-and-confused.netlify.app/post/model-fit-checks/

M <- b0
pp_check(object = M, type = "dens_overlay", nsamples = 100)
pp_check(M, type = "stat", stat = 'median', nsamples = 100)
pp_check(M, type = "stat", stat = 'mean', nsamples = 100)
pp_check(M,type = 'intervals')

plot(M)
as.data.frame(posterior_summary(M))

#R2
as.data.frame(bayes_R2(object = M,resp = NULL,summary = TRUE,robust = FALSE,probs = c(0.025, 0.975)))

#loo_R2
loo_R2(object = M)

#Plot with tidy_bayes-----

#still working on this
library(tidybayes)

M %>%
  spread_draws(b_Intercept, b_soilS) %>%
  mutate(condition_mean = b_Intercept + b_soilS) %>%
  ggplot(aes(y = condition, x = condition_mean, fill = stat(abs(x) < .8))) +
  stat_halfeye() +
  geom_vline(xintercept = c(-.8, .8), linetype = "dashed") +
  scale_fill_manual(values = c("gray80", "skyblue"))

M %>%
  spread_draws(r_condition[condition,]) %>%
  compare_levels(r_condition, by = condition) %>%
  ungroup() %>%
  mutate(condition = reorder(condition, r_condition)) %>%
  ggplot(aes(y = condition, x = r_condition)) +
  stat_halfeye() +
  geom_vline(xintercept = 0, linetype = "dashed")

m_mpg = brm(
  mpg ~ hp * cyl,
  data = mtcars)

mtcars %>%
  group_by(cyl) %>%
  data_grid(hp = seq_range(hp, n = 101)) %>%
  add_fitted_draws(m_mpg, n = 100) %>%
  ggplot(aes(x = hp, y = mpg, color = ordered(cyl))) +
  geom_line(aes(y = .value, group = paste(cyl, .draw)), alpha = .1) +
  geom_point(data = mtcars) +
  scale_color_brewer(palette = "Dark2")

#this script is for analyzing the effect of playback 

#load packages
library(tidyverse)
library(ggplot2)
library(emmeans)
library(glmmTMB)
library(DHARMa)
library(performance)
library(lme4)
library(dplyr)
library(broom.mixed)


#read in clean data
sjdf_clean <- read_csv("clean_data.csv")

#subset playback data
sjdf_clean <- sjdf_clean |>
  filter(SEASON == "WINTER")

#set variables classes
sjdf_clean <- sjdf_clean |> #fix class of variables
  mutate(across(c(ALARM, MOB, FLEE, INTEREST, HYPOTENUSE, LATENCY.ALARM, 
                  GROUP.SIZE, NUMBER.MOB, NUMBER.ALARM), as.numeric)) |>
  mutate(across(c(SPECIES, SEASON, PLAYBACK, SITE, SUBSITE, TREATMENT), 
                as.factor))

#Data analysis begins here


sjdf_clean <- sjdf_clean |> #mutate the group size and hypotenuse variables to 
  #standardize
  mutate(GROUP.SIZE = as.numeric(scale(GROUP.SIZE))) |>
  mutate(HYPOTENUSE = as.numeric(scale(HYPOTENUSE)))


# ALARM occurrence model
m.alarm <- glmer(ALARM ~ PLAYBACK + GROUP.SIZE + HYPOTENUSE +
                   (1 | SUBSITE),
                 data = sjdf_clean,
                 family = binomial)  

check_collinearity(m.alarm) #check for colinearity between predictors

sim.m.alarm <- simulateResiduals(fittedModel = m.alarm) #check model fit
testDispersion(sim.m.alarm)
testZeroInflation(sim.m.alarm)
testOutliers(sim.m.alarm) 

summary(m.alarm) # get model output

#rerun model without HYPOTENUSE, it is not significant and removes 7 hawk trials


# MOB occurrence model
m.mob <- glmer(MOB ~ PLAYBACK + GROUP.SIZE + HYPOTENUSE +
                 (1 | SUBSITE),
               data = sjdf_clean,
               family = binomial)

check_collinearity(m.mob)

sim.m.mob <- simulateResiduals(fittedModel = m.mob)
testDispersion(sim.m.mob)
testZeroInflation(sim.m.mob)
testOutliers(sim.m.mob)

summary(m.mob)


# INTEREST occurrence model
m.interest <- glmer(INTEREST ~ PLAYBACK + GROUP.SIZE + HYPOTENUSE +
                      (1 | SUBSITE),
                    data = sjdf_clean,
                    family = binomial)

check_collinearity(m.interest)

sim.m.interest <- simulateResiduals(fittedModel = m.interest)
testDispersion(sim.m.interest)
testZeroInflation(sim.m.interest)
testOutliers(sim.m.interest)

summary(m.interest)





#playback has no effect on our behavior probability repsonses.

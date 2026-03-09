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



# now check if playback affects alarm call counts or latency

#check alarm counts
m.alarm.call <- glmmTMB(NUMBER.ALARM ~ PLAYBACK + GROUP.SIZE + HYPOTENUSE +
                          (1 | SUBSITE),
                        ziformula = ~1,
                        data = sjdf_clean, 
                        family = nbinom2()
)

check_collinearity(m.alarm.call) #check for colinearity between predictors

sim.m.alarm.call <- simulateResiduals(fittedModel = m.alarm.call) #check model fit
testDispersion(sim.m.alarm.call)
testZeroInflation(sim.m.alarm.call)
testOutliers(sim.m.alarm.call) 

summary(m.alarm.call)
#no

#check mob counts
m.mob.call <- glmmTMB(NUMBER.MOB ~ PLAYBACK + GROUP.SIZE + HYPOTENUSE +
                        (1 | SUBSITE),
                      ziformula = ~0,
                      data = sjdf_mob, 
                      family = nbinom2()
)

check_collinearity(m.mob.call) # check for collinearity between predictors

sim.m.mob.call <- simulateResiduals(fittedModel = m.mob.call) # check model fit
testDispersion(sim.m.mob.call)
testZeroInflation(sim.m.mob.call)
testOutliers(sim.m.mob.call) 

summary(m.mob.call)

#no


#last check latency

m.latency <- glmmTMB(LATENCY.ALARM ~ PLAYBACK + GROUP.SIZE + HYPOTENUSE +
                       (1 | SUBSITE),
                     family = gaussian(link = "log"),
                     data = sjdf_clean)

check_collinearity(m.latency) # check for collinearity between predictors

sim.m.latency <- simulateResiduals(fittedModel = m.latency) # check model fit
testDispersion(sim.m.latency)
testZeroInflation(sim.m.latency)
testOutliers(sim.m.latency) 


summary(m.latency)

#no


#we can conclude that playback has no effect on our measured behaviors






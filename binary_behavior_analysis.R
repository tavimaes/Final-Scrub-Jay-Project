#this script is for analyzing the effect of species on probability 
#to observe behaviors


#load packages
library(tidyverse)
library(ggplot2)
library(ggeffects)
library(glmmTMB)
library(DHARMa)
library(performance)
library(lme4)
library(dplyr)
library(broom.mixed)


#read in clean data
sjdf_clean <- read_csv("clean_data.csv")

#set variables classes
sjdf_clean <- sjdf_clean |> #fix class of variables
  mutate(across(c(ALARM, MOB, FLEE, INTEREST, HYPOTENUSE, LATENCY.ALARM, 
                  GROUP.SIZE, NUMBER.MOB, NUMBER.ALARM), as.numeric)) |>
  mutate(across(c(SPECIES, SEASON, PLAYBACK, SITE, SUBSITE, TREATMENT), 
                as.factor))

#Data analysis begins here

#make hawk the reference level for model estimates
sjdf_clean$TREATMENT <- relevel(sjdf_clean$TREATMENT, ref = "HAWK")

sjdf_clean <- sjdf_clean |> #mutate the group size and hypotenuse variables to 
  #standardize
  mutate(GROUP.SIZE = as.numeric(scale(GROUP.SIZE))) |>
  mutate(HYPOTENUSE = as.numeric(scale(HYPOTENUSE)))


# ALARM occurrence model
m.alarm <- glmer(ALARM ~  SPECIES + TREATMENT + GROUP.SIZE + HYPOTENUSE +
                   SEASON +
                   (1 | SUBSITE),
                 data = sjdf_clean,
                 family = binomial)  

check_collinearity(m.alarm) #check for colinearity between predictors

sim.m.alarm <- simulateResiduals(fittedModel = m.alarm) #check model fit
testDispersion(sim.m.alarm)
testZeroInflation(sim.m.alarm)
testOutliers(sim.m.alarm) 

summary(m.alarm) # get model output
alarm_hypotenuse_est <- tidy(m.alarm) |> filter(term == "HYPOTENUSE")

#rerun model without HYPOTENUSE, it is not significant and removes 7 hawk trials

m.alarm <- glmer(ALARM ~  SPECIES + TREATMENT + GROUP.SIZE +
                   SEASON +
                   (1 | SUBSITE),
                 data = sjdf_clean,
                 family = binomial)

alarm_est <- tidy(m.alarm) |> mutate(model = "ALARM")


model_output_tibble <- bind_rows(alarm_est, alarm_hypotenuse_est)


# MOB occurrence model
m.mob <- glmer(MOB ~ SPECIES + TREATMENT + GROUP.SIZE + HYPOTENUSE +
                 SEASON +
                 (1 | SUBSITE),
               data = sjdf_clean,
               family = binomial)

check_collinearity(m.mob)

sim.m.mob <- simulateResiduals(fittedModel = m.mob)
testDispersion(sim.m.mob)
testZeroInflation(sim.m.mob)
testOutliers(sim.m.mob)

summary(m.mob)
mob_hypotenuse_est <- tidy(m.mob) |> filter(term == "HYPOTENUSE")

m.mob <- glmer(MOB ~ SPECIES + TREATMENT + GROUP.SIZE +
                 SEASON +
                 (1 | SUBSITE),
               data = sjdf_clean,
               family = binomial)
summary((m.mob))
mob_est <- tidy(m.mob) |> mutate(model = "MOB")

model_output_tibble <- bind_rows(model_output_tibble, mob_est, 
                                 mob_hypotenuse_est)


# INTEREST occurrence model
m.interest <- glmer(INTEREST ~ SPECIES + TREATMENT + GROUP.SIZE + HYPOTENUSE +
                      SEASON +
                      (1 | SUBSITE),
                    data = sjdf_clean,
                    family = binomial)

check_collinearity(m.interest)

sim.m.interest <- simulateResiduals(fittedModel = m.interest)
testDispersion(sim.m.interest)
testZeroInflation(sim.m.interest)
testOutliers(sim.m.interest)

summary(m.interest)
interest_hypotenuse_est <- tidy(m.interest) |> filter(term == "HYPOTENUSE")

m.interest <- glmer(INTEREST ~ SPECIES + TREATMENT + GROUP.SIZE +
                      SEASON +
                      (1 | SUBSITE),
                    data = sjdf_clean,
                    family = binomial)

interest_est <- tidy(m.interest) |> mutate(model = "INTEREST")

model_output_tibble <- bind_rows(model_output_tibble, interest_est, 
                                 interest_hypotenuse_est)



# FLEE occurrence model
m.flee <- glmer(FLEE ~ SPECIES + TREATMENT + GROUP.SIZE + HYPOTENUSE +
                  SEASON +
                  (1 | SUBSITE),
                data = sjdf_clean,
                family = binomial)

check_collinearity(m.flee)

sim.m.flee <- simulateResiduals(fittedModel = m.flee)
testDispersion(sim.m.flee)
testZeroInflation(sim.m.flee)
testOutliers(sim.m.flee)

summary(m.flee)
flee_hypotenuse_est <- tidy(m.flee) |> filter(term == "HYPOTENUSE")

m.flee <- glmer(FLEE ~ SPECIES + TREATMENT + GROUP.SIZE +
                  SEASON +
                  (1 | SUBSITE),
                data = sjdf_clean,
                family = binomial)

flee_est <- tidy(m.flee) |> mutate(model = "FLEE")

model_output_tibble <- bind_rows(model_output_tibble, flee_est, 
                                 flee_hypotenuse_est)


# tidy tibble up
model_output_tibble$model[7] <- "ALARM"
model_output_tibble$model[14] <- "MOB"
model_output_tibble$model[21] <- "INTEREST"
model_output_tibble$model[28] <- "FLEE"

model_output_tibble <- model_output_tibble |> 
  select(-effect, - group) |>
  filter(term != "(Intercept)") |> 
  filter(term != "sd__(Intercept)")
  
model_output_tibble <- model_output_tibble |> select(model, everything())

# save as a file 
write.csv(model_output_tibble, "binary_behav_stats.csv")


















all_pred <- merge(all_pred, labels, by = "MODEL")
ggplot(all_pred, aes(x = SPECIES, y = pred_prob, color = SPECIES)) +
  geom_point(size = 2) +
  geom_errorbar(data = subset(all_pred, MODEL != "Model 4"),  # exclude model 4
                aes(ymin = lower, ymax = upper),
                width = 0.2) +
  geom_text(
    aes(label = PANEL_LABEL),
    x = -Inf,  # far left
    y = 1.0,  # just above top of y-axis (adjust as needed)
    hjust = -.5, # nudges left
    vjust = 1,
    inherit.aes = FALSE) +
  facet_wrap(~ MODEL) +
  ylim(0,1) +
  labs(
    x = "Species",
    y = "Predicted Probability of Behaviour"
  ) +
  theme_classic() +
  theme(strip.text = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = 1))



#this script is for analyzing the effect of species on probability 
#to observe behaviors


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
library(ggdist)


#read in clean data
sjdf_clean <- read_csv("clean_data.csv")

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
m.alarm <- glmer(ALARM ~  SPECIES*TREATMENT + GROUP.SIZE + HYPOTENUSE + PLAYBACK +
                   (1 | SUBSITE),
                 data = sjdf_clean,
                 family = binomial) 

check_collinearity(m.alarm) #check for colinearity between predictors

sim.m.alarm <- simulateResiduals(fittedModel = m.alarm) #check model fit
testDispersion(sim.m.alarm)
testZeroInflation(sim.m.alarm)
testOutliers(sim.m.alarm) 

summary(m.alarm) # get model output


#save playback summary stats for alarm
alarm_playback_est <- tidy(m.alarm) |> filter(term == "PLAYBACKYES") |> 
  mutate(model = "ALARM")

#rerun model without PLAYBACK, it is not significant. Now we can include SEASON

m.alarm <- glmer(ALARM ~  SPECIES*TREATMENT + GROUP.SIZE + HYPOTENUSE + SEASON +
                   (1 | SUBSITE),
                 data = sjdf_clean,
                 family = binomial) 

check_collinearity(m.alarm) #check for colinearity between predictors

sim.m.alarm <- simulateResiduals(fittedModel = m.alarm) #check model fit
testDispersion(sim.m.alarm)
testZeroInflation(sim.m.alarm)
testOutliers(sim.m.alarm) 

summary(m.alarm) # get model output

alarm_hypotenuse_est <- tidy(m.alarm) |> filter(term == "HYPOTENUSE") |> 
  mutate(model = "ALARM")

#now rerun once more without hypotenuse, it has no significance and removing 
#allows 7 extra data points to be used

m.alarm <- glmer(ALARM ~  SPECIES*TREATMENT + GROUP.SIZE + SEASON +
                   (1 | SUBSITE),
                 data = sjdf_clean,
                 family = binomial) 

check_collinearity(m.alarm) #check for colinearity between predictors

sim.m.alarm <- simulateResiduals(fittedModel = m.alarm) #check model fit
testDispersion(sim.m.alarm)
testZeroInflation(sim.m.alarm)
testOutliers(sim.m.alarm) 


summary(m.alarm) # get model output


#get all other stats from final model
alarm_est <- tidy(m.alarm) |> mutate(model = "ALARM")


model_output_tibble <- bind_rows(alarm_est, alarm_hypotenuse_est, 
                                 alarm_playback_est)


# MOB occurrence model
m.mob <- glmer(MOB ~ SPECIES*TREATMENT + GROUP.SIZE + HYPOTENUSE + PLAYBACK +
                 (1 | SUBSITE),
               data = sjdf_clean,
               family = binomial)

check_collinearity(m.mob) # check collinearity between predictors

sim.m.mob <- simulateResiduals(fittedModel = m.mob) # check model fit
testDispersion(sim.m.mob)
testZeroInflation(sim.m.mob)
testOutliers(sim.m.mob)

summary(m.mob) # get model output

# save playback summary stats for MOB
mob_playback_est <- tidy(m.mob) |> 
  filter(term == "PLAYBACKYES") |> 
  mutate(model = "MOB")

# rerun model without PLAYBACK (if not significant), now include SEASON
m.mob <- glmer(MOB ~ SPECIES*TREATMENT + GROUP.SIZE + HYPOTENUSE + SEASON +
                 (1 | SUBSITE),
               data = sjdf_clean,
               family = binomial)

check_collinearity(m.mob)
sim.m.mob <- simulateResiduals(fittedModel = m.mob)
testDispersion(sim.m.mob)
testZeroInflation(sim.m.mob)
testOutliers(sim.m.mob)

summary(m.mob)

mob_hypotenuse_est <- tidy(m.mob) |> 
  filter(term == "HYPOTENUSE") |> 
  mutate(model = "MOB")

# rerun final model without HYPOTENUSE if not significant
m.mob <- glmer(MOB ~ SPECIES*TREATMENT + GROUP.SIZE + SEASON +
                 (1 | SUBSITE),
               data = sjdf_clean,
               family = binomial)

check_collinearity(m.mob)
sim.m.mob <- simulateResiduals(fittedModel = m.mob)
testDispersion(sim.m.mob)
testZeroInflation(sim.m.mob)
testOutliers(sim.m.mob)

summary(m.mob)

# get all other stats from final model
mob_est <- tidy(m.mob) |> mutate(model = "MOB")

model_output_tibble <- bind_rows(model_output_tibble, mob_hypotenuse_est, 
                               mob_est, mob_playback_est)



# INTEREST occurrence model
m.interest <- glmer(INTEREST ~ SPECIES*TREATMENT + GROUP.SIZE + HYPOTENUSE + PLAYBACK +
                      (1 | SUBSITE),
                    data = sjdf_clean,
                    family = binomial)

check_collinearity(m.interest)  # check collinearity between predictors

sim.m.interest <- simulateResiduals(fittedModel = m.interest) # check model fit
testDispersion(sim.m.interest)
testZeroInflation(sim.m.interest)
testOutliers(sim.m.interest)

summary(m.interest)  # get model output

# save playback summary stats for INTEREST
interest_playback_est <- tidy(m.interest) |> 
  filter(term == "PLAYBACKYES") |> 
  mutate(model = "INTEREST")

# rerun model without PLAYBACK if not significant, include SEASON
m.interest <- glmer(INTEREST ~ SPECIES*TREATMENT + GROUP.SIZE + HYPOTENUSE + SEASON +
                      (1 | SUBSITE),
                    data = sjdf_clean,
                    family = binomial)

check_collinearity(m.interest)
sim.m.interest <- simulateResiduals(fittedModel = m.interest)
testDispersion(sim.m.interest)
testZeroInflation(sim.m.interest)
testOutliers(sim.m.interest)

summary(m.interest)

interest_hypotenuse_est <- tidy(m.interest) |> 
  filter(term == "HYPOTENUSE") |> 
  mutate(model = "INTEREST")

# rerun final model without HYPOTENUSE if not significant
m.interest <- glmer(INTEREST ~ SPECIES*TREATMENT + GROUP.SIZE + SEASON +
                      (1 | SUBSITE),
                    data = sjdf_clean,
                    family = binomial)

check_collinearity(m.interest)
sim.m.interest <- simulateResiduals(fittedModel = m.interest)
testDispersion(sim.m.interest)
testZeroInflation(sim.m.interest)
testOutliers(sim.m.interest)

summary(m.interest)

# get all other stats from final model
interest_est <- tidy(m.interest) |> mutate(model = "INTEREST")

model_output_tibble <- bind_rows(model_output_tibble, interest_hypotenuse_est, 
                                interest_est, interest_playback_est)

# FLEE occurrence model
m.flee <- glmer(FLEE ~ SPECIES*TREATMENT + GROUP.SIZE + HYPOTENUSE + PLAYBACK +
                  (1 | SUBSITE),
                data = sjdf_clean,
                family = binomial)

check_collinearity(m.flee)  # check collinearity between predictors

sim.m.flee <- simulateResiduals(fittedModel = m.flee) # check model fit
testDispersion(sim.m.flee)
testZeroInflation(sim.m.flee)
testOutliers(sim.m.flee)

summary(m.flee)  # get model output

# save playback summary stats for FLEE
flee_playback_est <- tidy(m.flee) |> 
  filter(term == "PLAYBACKYES") |> 
  mutate(model = "FLEE")

# rerun model without PLAYBACK if not significant, include SEASON
m.flee <- glmer(FLEE ~ SPECIES*TREATMENT + GROUP.SIZE + HYPOTENUSE + SEASON +
                  (1 | SUBSITE),
                data = sjdf_clean,
                family = binomial)

check_collinearity(m.flee)
sim.m.flee <- simulateResiduals(fittedModel = m.flee)
testDispersion(sim.m.flee)
testZeroInflation(sim.m.flee)
testOutliers(sim.m.flee)

summary(m.flee)

flee_hypotenuse_est <- tidy(m.flee) |> 
  filter(term == "HYPOTENUSE") |> 
  mutate(model = "FLEE")

# rerun final model without HYPOTENUSE if not significant
m.flee <- glmer(FLEE ~ SPECIES*TREATMENT + GROUP.SIZE + SEASON +
                  (1 | SUBSITE),
                data = sjdf_clean,
                family = binomial)

check_collinearity(m.flee)
sim.m.flee <- simulateResiduals(fittedModel = m.flee)
testDispersion(sim.m.flee)
testZeroInflation(sim.m.flee)
testOutliers(sim.m.flee)

summary(m.flee)

# get all other stats from final model
flee_est <- tidy(m.flee) |> mutate(model = "FLEE")

model_output_tibble <- bind_rows(model_output_tibble, flee_hypotenuse_est, 
                                 flee_est, flee_playback_est)


# GRAPH building

library(ggbeeswarm)

ggplot() +
  
  # Replace stat_slab with quasirandom points
  geom_quasirandom(
    data = sj_long,
    aes(
      x = SPECIES,
      y = VALUE,
      color = TREATMENT,
      group = TREATMENT
    ),
    dodge.width = 0.6,  
    width = 0.2,
    size = 1,
    alpha = 0.1
  ) +
  
  # Predicted points
  geom_point(
    data = pred_df,
    aes(
      x = SPECIES,
      y = prob,
      color = TREATMENT,
      shape = TREATMENT
    ),
    size = 3,
    position = position_dodge(width = 0.6)
  ) +
  
  # Predicted error bars
  geom_errorbar(
    data = pred_df, show.legend = FALSE,
    aes(
      x = SPECIES,
      ymin = asymp.LCL,
      ymax = asymp.UCL,
      color = TREATMENT
    ),
    width = 0.3,
    position = position_dodge(width = 0.6)
  ) +
  
  # Facet labels
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0,
    size = 5
  ) +
  
  facet_wrap(~ factor(BEHAVIOR, levels = c("ALARM", "MOB", "INTEREST", "FLEE")),
             ncol = 2) +
  
  scale_color_manual(values = c("CONTROL" = "gray60",
                                "HAWK" = "darkblue")) +
  scale_fill_manual(values = c("CONTROL" = "gray60",
                               "HAWK" = "darkblue")) +
  scale_shape_manual(values = c("CONTROL" = 15,
                                "HAWK" = 17)) +
  
  scale_x_discrete(labels = c(
    CASJ = expression(italic("A. californica")),
    ISSJ = expression(italic("A. insularis"))
  )) +
  
  scale_y_continuous(limits = c(0,1), expand = c(0.02,0.05)) +
  
  theme_classic(base_size = 14) +
  labs(x = "Species", y = "Probability") +
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))

ggsave("fig1.png", last_plot(), width = 7, height = 5, dpi = 300)







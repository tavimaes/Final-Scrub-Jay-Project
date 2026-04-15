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
library(brms)

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

check_collinearity(m.alarm) #check for collinearity between predictors

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

model_output_tibble <- bind_rows(model_output_tibble,  mob_est, mob_hypotenuse_est, 
                              mob_playback_est)



# INTEREST occurrence model
m.interest <- glmer(INTEREST ~ SPECIES + TREATMENT + GROUP.SIZE + HYPOTENUSE + PLAYBACK +
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
m.interest <- glmer(INTEREST ~ SPECIES + TREATMENT + GROUP.SIZE + HYPOTENUSE + SEASON +
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
m.interest <- glmer(INTEREST ~ SPECIES + TREATMENT + GROUP.SIZE + SEASON +
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

model_output_tibble <- bind_rows(model_output_tibble, interest_est, interest_hypotenuse_est, 
                                interest_playback_est)

# FLEE occurrence model
m.flee <- glmer(FLEE ~ SPECIES + TREATMENT + GROUP.SIZE + HYPOTENUSE + PLAYBACK +
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
m.flee <- glmer(FLEE ~ SPECIES + TREATMENT + GROUP.SIZE + HYPOTENUSE + SEASON +
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
m.flee <- glmer(FLEE ~ SPECIES + TREATMENT + GROUP.SIZE + SEASON +
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

model_output_tibble <- bind_rows(model_output_tibble, flee_est, 
                                 flee_hypotenuse_est, flee_playback_est)

model_output_tibble_final <- model_output_tibble |>
  filter(effect != "ran_pars") |>
  filter(term != "(Intercept)")


write_csv(model_output_tibble_final, "binary_behav_stats.csv")



# GRAPH building

alarm_df <- emmeans(m.alarm, ~ SPECIES * TREATMENT, type = "response") |>
  as.data.frame() |>
  mutate(BEHAVIOR = "ALARM")

mob_df <- emmeans(m.mob, ~ SPECIES * TREATMENT, type = "response") |>
  as.data.frame() |>
  mutate(BEHAVIOR = "MOB")

interest_df <- emmeans(m.interest, ~ SPECIES * TREATMENT, type = "response") |>
  as.data.frame() |>
  mutate(BEHAVIOR = "INTEREST")

flee_df <- emmeans(m.flee, ~ SPECIES * TREATMENT, type = "response") |>
  as.data.frame() |>
  mutate(BEHAVIOR = "FLEE")

pred_df <- bind_rows(alarm_df, mob_df, interest_df, flee_df)

#create facet labels
label_df <- data.frame(
  BEHAVIOR = c("ALARM", "MOB", "INTEREST", "FLEE"),
  label = c("a)", "b)", "c)", "d)"),
  x = -Inf,
  y = Inf
)

#fix error bars where there was no variance

pred_df$asymp.LCL[13] <- 1
pred_df$asymp.UCL[13] <- 1
pred_df$asymp.LCL[15] <- 1
pred_df$asymp.UCL[15] <- 1
pred_df$asymp.LCL[16] <- 1
pred_df$asymp.UCL[16] <- 1

#raw data
raw_df <- sjdf_clean |>
  pivot_longer(
    cols = c(ALARM, MOB, INTEREST, FLEE),
    names_to = "BEHAVIOR",
    values_to = "RESPONSE"
  ) |>
  group_by(SPECIES, TREATMENT, BEHAVIOR) |>
  summarise(
    prop = mean(RESPONSE),
    n = n(),
    .groups = "drop"
  )


ggplot() +
  # Predicted points
  geom_col(
    data = raw_df,
    aes(
      x = SPECIES,
      y = prop,
      fill = TREATMENT
    ),
    position = position_dodge(width = 0.6),
    width = 0.5,
    alpha = 0.2, show.legend = FALSE   # transparency so predictions show on top
  ) +
  
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
  
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.3,   # pushes slightly inside from left
    vjust = 1.5,    # pushes slightly down from top
    size = 5
  ) +
  
  facet_wrap(~ factor(BEHAVIOR, levels = c("ALARM", "MOB", "INTEREST", "FLEE")),
             ncol = 2) +
  
  scale_color_manual(
    name = "Treatment",
    values = c("CONTROL" = "red", "HAWK" = "darkblue"),
    labels = c("CONTROL" = "Control", "HAWK" = "Hawk")
  ) +
  
  scale_shape_manual(
    name = "Treatment",
    values = c("CONTROL" = 15, "HAWK" = 17),
    labels = c("CONTROL" = "Control", "HAWK" = "Hawk")
  ) +
  
  scale_x_discrete(labels = c(
    CASJ = expression(italic("A. californica")),
    ISSJ = expression(italic("A. insularis"))
  )) +
  
  scale_y_continuous(limits = c(0,1), expand = c(0.02,0.05)) +
  
  theme_classic(base_size = 14) +
  labs(x = "Species", y = "Probability of behavior occurring",
       color = "Treatment", shape = "Treatment") +
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.2))

ggsave("fig1.png", last_plot(), width = 7, height = 5, dpi = 300)



contrast(
  emmeans(m.alarm, ~ SPECIES * TREATMENT),
  interaction = "pairwise", type = "response"
)

#----- Post-hoc tests -----------

#test within treatment means

emm.alarm <- emmeans(m.alarm, ~ SPECIES | TREATMENT)

pairs(emm.alarm, type = "response")

emm.mob <- emmeans(m.mob, ~ SPECIES | TREATMENT)

pairs(emm.mob, type = "response")

emm.interest <- emmeans(m.interest, ~ SPECIES | TREATMENT)
pairs(emm.interest, type = "response")

emm.flee <- emmeans(m.flee, ~ SPECIES | TREATMENT)
pairs(emm.flee, type = "response")

#test within species means

emm.alarm <- emmeans(m.alarm, ~ TREATMENT | SPECIES)
pairs(emm.alarm, type = "response")

emm.mob <- emmeans(m.mob, ~ TREATMENT | SPECIES)
pairs(emm.mob, type = "response")

emm.interest <- emmeans(m.interest, ~ TREATMENT | SPECIES)
pairs(emm.interest, type = "response")

emm.flee <- emmeans(m.flee, ~ TREATMENT | SPECIES)
pairs(emm.flee, type = "response")

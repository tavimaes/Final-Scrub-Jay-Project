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


# GRAPH building


# Generate predictions for each model

# Create a function to get predicted probabilities for a model
get_pred <- function(model, response_name) {
  # SPECIES varies, TREATMENT fixed at "hawk", other covariates averaged
  #(including weighted mean of SEASON)
  emm <- emmeans(model, ~ SPECIES, at = list(TREATMENT = "HAWK"))
  
  # Convert to data frame with predicted probabilities and 95% CI
  emm_df <- summary(emm, type = "response", infer = TRUE) |> 
    as.data.frame() |> 
    mutate(response = response_name)
  
  return(emm_df)
}

# Apply to all four models
pred_alarm    <- get_pred(m.alarm, "ALARM")
pred_mob      <- get_pred(m.mob, "MOB")
pred_interest <- get_pred(m.interest, "INTEREST")
pred_flee     <- get_pred(m.flee, "FLEE")

# Combine all into one table
predictions_table <- bind_rows(pred_alarm, pred_mob, pred_interest, pred_flee)

# View
predictions_table


# graphing time
predictions_table <- predictions_table |> 
  mutate(panel_label = factor(response, 
                              levels = c("Alarm", "Mob", "Interest", "FLEE"),
                              labels = c("a)", "b)", "c)", "d)")))
predictions_table <- predictions_table |>
  filter(response != "FLEE")




predictions_table <- predictions_table %>%
  mutate(response = recode(response,
                           "ALARM" = "Alarm",
                           "MOB" = "Mob",
                           "INTEREST" = "Interest"))

sjdf_issj <- sjdf_clean |>
  filter(TREATMENT == "HAWK") |>
  filter(SPECIES == "ISSJ")

sjdf_casj <- sjdf_clean |>
  filter(TREATMENT == "HAWK") |>
  filter(SPECIES == "CASJ")

issj_alarm_mean <- mean(sjdf_issj$ALARM)
casj_alarm_mean <- mean(sjdf_casj$ALARM)
issj_mob_mean <- mean(sjdf_issj$MOB)
casj_mob_mean <- mean(sjdf_casj$MOB)
issj_interest_mean <- mean(sjdf_issj$INTEREST)
casj_interest_mean <- mean(sjdf_casj$INTEREST)

raw_data <- data.frame(
  response = c("Alarm", "Alarm", "Mob","Mob", "Interest", "Interest"), 
  mean = c(issj_alarm_mean, casj_alarm_mean, issj_mob_mean,
           casj_mob_mean, issj_interest_mean, casj_interest_mean),
  SPECIES = c("ISSJ", "CASJ", "ISSJ", "CASJ", "ISSJ", "CASJ")
)

sig_df <- data.frame(
  response = c("Alarm", "Mob", "Interest"),  
  SPECIES = c("ISSJ", "ISSJ", "ISSJ"),
  y = c(0.68, 0.03, 0.28),
  label = c("*", "*", "*")
)



ggplot(predictions_table, aes(x = response, y = prob, shape = SPECIES, 
                              group = SPECIES, color = SPECIES)) +
  geom_point(data = raw_data, aes(x = response, y = mean), 
             show.legend = FALSE,
             color = "darkgray",
             position = position_nudge(x = 0.2),
             size = 4) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, 
                linewidth = 0.75, show.legend = FALSE) +
  geom_text(data = sig_df, aes(x = response, y = y, label = label),
            inherit.aes = FALSE, size = 6) +
  labs(x = "Behaviour", y = "Predicted probability", color = NULL, shape = NULL) +
  scale_shape_manual(values = c("ISSJ" = 16, "CASJ" = 1)) +
  theme_classic(base_size = 14) +
  scale_color_manual(values = c("red", "darkblue")) +
  scale_x_discrete(limits = c("Alarm", "Mob", "Interest"))
  


ggsave("fig1.png", last_plot(), width = 8, height = 5, dpi = 300)



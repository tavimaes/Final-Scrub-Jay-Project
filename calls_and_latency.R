#this script is for analyzing call and latency data

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

#read in data
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

# ALARM call count model

m.alarm.call <- glmmTMB(NUMBER.ALARM ~ SPECIES + TREATMENT + GROUP.SIZE + SEASON + 
                          HYPOTENUSE + (1 | SUBSITE),
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

alarm.call_hypotenuse_est <- tidy(m.alarm.call) |> filter(term == "HYPOTENUSE")

# rerun without hypotenuse
m.alarm.call <- glmmTMB(NUMBER.ALARM ~ SPECIES + TREATMENT + GROUP.SIZE + SEASON
                        + (1 | SUBSITE),
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

alarm.call_est <- tidy(m.alarm.call) |> mutate(model = "ALARM.CALL")

model_output_tibble <- bind_rows(alarm.call_est, alarm.call_hypotenuse_est)


# MOB call count model

sjdf_mob <- sjdf_clean |>
  filter(NUMBER.MOB != 179) #remove outlier, shot flew directly at nesting area


m.mob.call <- glmmTMB(NUMBER.MOB ~ SPECIES + TREATMENT + GROUP.SIZE + SEASON + 
                        HYPOTENUSE + (1 | SUBSITE),
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

mob.call_hypotenuse_est <- tidy(m.mob.call) |> 
  filter(term == "HYPOTENUSE")

# rerun without hypotenuse
m.mob.call <- glmmTMB(NUMBER.MOB ~ SPECIES + TREATMENT + GROUP.SIZE + SEASON
                      + (1 | SUBSITE),
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

mob.call_est <- tidy(m.mob.call) |> 
  mutate(model = "MOB.CALL")

model_output_tibble <- bind_rows(model_output_tibble, mob.call_est, mob.call_hypotenuse_est)


#latency analysis

m.latency <- glmmTMB(LATENCY.ALARM ~ SPECIES + TREATMENT + GROUP.SIZE + 
                       SEASON + HYPOTENUSE + (1|SUBSITE),
        family = gaussian(link = "log"),
        data = sjdf_clean)

check_collinearity(m.latency) # check for collinearity between predictors

sim.m.latency <- simulateResiduals(fittedModel = m.latency) # check model fit
testDispersion(sim.m.latency)
testZeroInflation(sim.m.latency)
testOutliers(sim.m.latency) 


summary(m.latency)

latency_est <- tidy(m.latency) |> 
  mutate(model = "LATENCY")

model_output_tibble <- bind_rows(model_output_tibble, latency_est)

model_output_tibble <- model_output_tibble |>
  select(-group, -component, -effect) |>
  filter(term != "(Intercept)") |>
  filter(term != "sd__(Intercept)") |>
  filter(term != "sd__Observation")

model_output_tibble$model[5] <- "ALARM.CALL"
model_output_tibble$model[10] <- "MOB.CALL"

model_output_tibble <- model_output_tibble |> select(model, everything())

# save as a file 
write.csv(model_output_tibble, "call_and_latency_stats.csv")

#graphing time

get_pred <- function(model, response_name) {
  
  emm <- emmeans(model, ~ SPECIES, at = list(TREATMENT = "HAWK"))
  
  emm_df <- summary(emm, type = "response") |>
    as.data.frame() |>
    mutate(model = response_name)
  
  return(emm_df)
}

# Apply to all four models
pred_alarm.call    <- get_pred(m.alarm.call, "NUMBER.ALARM")
pred_mob.call     <- get_pred(m.mob.call, "NUMBER.MOB")
pred_latency <- get_pred(m.latency, "LATENCY.ALARM")

# Combine all into one table
predictions_table <- bind_rows(pred_alarm.call, pred_mob.call, pred_latency)

# View
predictions_table

#graph building

# Prepare the data for faceting
predictions_table_facet <- predictions_table %>%
  mutate(
    facet_panel = case_when(
      model == "LATENCY.ALARM"      ~ "Latency",
      model %in% c("NUMBER.ALARM","NUMBER.MOB") ~ "Call counts"
    ),
    behavior = case_when(
      model == "LATENCY.ALARM"      ~ "Alarming latency (s)",
      model == "NUMBER.ALARM"       ~ "Alarm calls",
      model == "NUMBER.MOB"         ~ "Mobbing calls"
    )
  )

# plot
sig_df <- data.frame(
  behavior = c("Alarming latency (s)", "Mobbing calls"),       
  facet_panel = c("Latency", "Call counts"), 
  SPECIES = c("ISSJ", "ISSJ"),          
  y = c(1.75, 13),                  
  label = c("*", "*")
)



raw_data_long <- sjdf_clean |>
  pivot_longer(
    cols = c(NUMBER.ALARM, NUMBER.MOB, LATENCY.ALARM),
    names_to = "behavior",
    values_to = "value"
  ) |>
  filter(TREATMENT != "CONTROL") |>
  mutate(
    facet_panel = case_when(
      behavior == "LATENCY.ALARM"      ~ "Latency",
      behavior %in% c("NUMBER.ALARM","NUMBER.MOB") ~ "Call counts"
    ),
    behavior = recode(behavior,
                      "NUMBER.ALARM" = "Alarm calls",
                      "NUMBER.MOB" = "Mobbing calls",
                      "LATENCY.ALARM" = "Alarming latency (s)")
      
  ) |>
  filter(value != 179)


ggplot(predictions_table_facet, aes(x = behavior, y = response, color = SPECIES, shape = SPECIES)) +
  geom_point(data = raw_data_long,
             aes(x = behavior, y = value, shape = SPECIES),
             size = 2,
             position = position_jitter(height = 0, width = 0.25),
             color = "darkgray",
             alpha = 0.6,
             show.legend = FALSE) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = response - SE, ymax = response + SE),
                linewidth = 0.75, width = 0.25, show.legend = FALSE) +
  geom_text(data = sig_df, aes(x = behavior, y = y, label = label),
            inherit.aes = FALSE, size = 8) +
  facet_wrap(~ facet_panel, scales = "free") +
  scale_shape_manual(values = c("ISSJ" = 16, "CASJ" = 1)) +
  labs(x = "Vocal response",
       y = "Predicted value",
       color = NULL,
       shape = NULL) +
  theme_classic(base_size = 14) +
  theme(strip.text = element_blank()) +
  scale_color_manual(values = c("red", "darkblue")) +
  scale_y_continuous(
    limits = c(0, NA)
  )
  
  

ggsave("fig2.png", last_plot(), width = 8, height = 5, dpi = 300)







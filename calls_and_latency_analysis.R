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
library(ggeffects)

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

m.alarm.call <- glmmTMB(NUMBER.ALARM ~ SPECIES*TREATMENT + GROUP.SIZE + 
                          HYPOTENUSE + PLAYBACK + (1 | SUBSITE),
  data = sjdf_clean, 
  family = nbinom2()
  )

check_collinearity(m.alarm.call) #check for colinearity between predictors

sim.m.alarm.call <- simulateResiduals(fittedModel = m.alarm.call) #check model fit
testDispersion(sim.m.alarm.call)
testZeroInflation(sim.m.alarm.call)
testOutliers(sim.m.alarm.call) 

summary(m.alarm.call)

alarm.call_playback_est <- tidy(m.alarm.call) |> filter(term == "PLAYBACKYES")

# rerun without playback
m.alarm.call1 <- glmmTMB(NUMBER.ALARM ~ SPECIES*TREATMENT + GROUP.SIZE + SEASON 
                        + HYPOTENUSE
                        + (1 | SUBSITE),
                        data = sjdf_clean, 
                        family = nbinom2()
)

m.alarm.call <- glmmTMB(NUMBER.ALARM ~ SPECIES*TREATMENT + GROUP.SIZE + SEASON 
                         + HYPOTENUSE
                         + (1 | SUBSITE),
                         ziformula = ~1,
                         data = sjdf_clean, 
                         family = nbinom2()
)

AIC(m.alarm.call1, m.alarm.call) #zero inflated model has lower AIC

check_collinearity(m.alarm.call) #check for colinearity between predictors

sim.m.alarm.call <- simulateResiduals(fittedModel = m.alarm.call) #check model fit
testDispersion(sim.m.alarm.call)
testZeroInflation(sim.m.alarm.call)
testOutliers(sim.m.alarm.call) 

summary(m.alarm.call)

alarm.call_hypotenuse_est <- tidy(m.alarm.call) |> filter(term == "HYPOTENUSE")

# rerun without hypotenuse, final model
m.alarm.call <- glmmTMB(NUMBER.ALARM ~ SPECIES*TREATMENT + GROUP.SIZE + SEASON 
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

model_output_tibble <- bind_rows(alarm.call_est, alarm.call_hypotenuse_est, 
                                 alarm.call_playback_est)



# MOB call count model

sjdf_clean <- sjdf_clean |> #remove outlier 
  filter(NUMBER.MOB != 179)


m.mob.call <- glmmTMB(
  NUMBER.MOB ~ SPECIES*TREATMENT + GROUP.SIZE + HYPOTENUSE + PLAYBACK + (1 | SUBSITE),
  data = sjdf_clean,
  family = nbinom2()
)

# Check collinearity
check_collinearity(m.mob.call)

# Check model fit
sim.m.mob <- simulateResiduals(fittedModel = m.mob.call)
testDispersion(sim.m.mob)
testZeroInflation(sim.m.mob)
testOutliers(sim.m.mob)

# Inspect summary
summary(m.mob.call)

# Extract estimate for PLAYBACK
mob.call_playback_est <- tidy(m.mob.call) |> filter(term == "PLAYBACKYES")


m.mob.call <- glmmTMB(
  NUMBER.MOB ~ SPECIES*TREATMENT + GROUP.SIZE + SEASON + HYPOTENUSE + (1 | SUBSITE),
  data = sjdf_clean,
  family = nbinom2()
)

check_collinearity(m.mob.call)

sim.m.mob <- simulateResiduals(fittedModel = m.mob.call)
testDispersion(sim.m.mob)
testZeroInflation(sim.m.mob)
testOutliers(sim.m.mob)

summary(m.mob.call)

# Extract estimate for HYPOTENUSE
mob.call_hypotenuse_est <- tidy(m.mob.call) |> filter(term == "HYPOTENUSE")


m.mob.call <- glmmTMB(
  NUMBER.MOB ~ SPECIES*TREATMENT + GROUP.SIZE + SEASON + (1 | SUBSITE),
  data = sjdf_clean,
  family = nbinom2()
)

check_collinearity(m.mob.call)

sim.m.mob <- simulateResiduals(fittedModel = m.mob.call)
testDispersion(sim.m.mob)
testZeroInflation(sim.m.mob)
testOutliers(sim.m.mob)

summary(m.mob.call)

mob.call_final_est <- tidy(m.mob.call) |> mutate(model = "MOB.CALL")


mob_model_output <- bind_rows(
  mob.call_final_est,
  mob.call_hypotenuse_est,
  mob.call_playback_est
)

model_output_tibble <- bind_rows(model_output_tibble, mob_model_output)

model_output_tibble <- model_output_tibble |>
  filter(effect != "ran_pars") |>
  filter(term != "(Intercept)")

write_csv(model_output_tibble, "call_count_stats.csv")


#latency analysis

# LATENCY model



# --- FULL MODEL ---
m.latency <- glmmTMB(
  LATENCY.ALARM ~ SPECIES*TREATMENT + GROUP.SIZE + HYPOTENUSE + PLAYBACK + (1 | SUBSITE),
  data = sjdf_clean,
  family = gaussian(link = "log")
)

# Check collinearity
check_collinearity(m.latency)

# Check model fit
sim.m.latency <- simulateResiduals(fittedModel = m.latency)
testDispersion(sim.m.latency)
testZeroInflation(sim.m.latency)
testOutliers(sim.m.latency)

# Inspect summary
summary(m.latency)

# Extract estimate for PLAYBACK
latency_playback_est <- tidy(m.latency) |> 
  filter(term == "PLAYBACKYES")


# --- REMOVE PLAYBACK ---
m.latency <- glmmTMB(
  LATENCY.ALARM ~ SPECIES*TREATMENT + GROUP.SIZE + SEASON + HYPOTENUSE + (1 | SUBSITE),
  data = sjdf_clean,
  family = gaussian(link = "log")
)

check_collinearity(m.latency)

sim.m.latency <- simulateResiduals(fittedModel = m.latency)
testDispersion(sim.m.latency)
testZeroInflation(sim.m.latency)
testOutliers(sim.m.latency)

summary(m.latency)

# Extract estimate for HYPOTENUSE
latency_hypotenuse_est <- tidy(m.latency) |> 
  filter(term == "HYPOTENUSE")


# --- FINAL MODEL ---
m.latency <- glmmTMB(
  LATENCY.ALARM ~ SPECIES*TREATMENT + GROUP.SIZE + SEASON + (1 | SUBSITE),
  data = sjdf_clean,
  family = gaussian(link = "log")
)

check_collinearity(m.latency)

sim.m.latency <- simulateResiduals(fittedModel = m.latency)
testDispersion(sim.m.latency)
testZeroInflation(sim.m.latency)
testOutliers(sim.m.latency)

summary(m.latency)

latency_final_est <- tidy(m.latency) |> 
  mutate(model = "LATENCY.ALARM")


# --- COMBINE OUTPUTS ---
latency_model_output <- bind_rows(
  latency_final_est,
  latency_hypotenuse_est,
  latency_playback_est
)




# --- PREDICTIONS ---
pred_latency <- emmeans(m.latency, ~ SPECIES * TREATMENT, type = "response") |>
  summary() |>
  as.data.frame()


# --- PLOT ---
ggplot() +
  
  # raw data
  geom_quasirandom(
    data = sjdf_clean,
    aes(x = SPECIES,
        y = LATENCY.ALARM,
        color = TREATMENT),
    dodge.width = 0.75,
    width = 0.35,
    alpha = 0.3,
    size = 1
  ) +
  
  # predicted means
  geom_point(
    data = pred_latency,
    aes(x = SPECIES,
        y = response,
        color = TREATMENT,
        shape = TREATMENT),
    size = 5,
    position = position_dodge(width = 0.75)
  ) +
  
  # confidence intervals
  geom_errorbar(
    data = pred_latency,
    aes(x = SPECIES,
        ymin = asymp.LCL,
        ymax = asymp.UCL,
        color = TREATMENT),
    width = 0.25,
    linewidth = 1,
    position = position_dodge(width = 0.75)
  ) +
  
  scale_color_manual(values = c(
    CONTROL = "gray60",
    HAWK = "darkblue"
  )) +
  
  scale_shape_manual(values = c(
    CONTROL = 15,
    HAWK = 17
  )) +
  
  scale_x_discrete(labels = c(
    CASJ = expression(italic("A. californica")),
    ISSJ = expression(italic("A. insularis"))
  )) +
  
  theme_classic(base_size = 14) +
  labs(
    x = "Species",
    y = "Latency to alarm call")





# ------- final graph -------

# Combine predictions for both models
alarm_call_df <- emmeans(m.alarm.call, ~ SPECIES * TREATMENT, type = "response") |>
  as.data.frame() |>
  mutate(BEHAVIOR = "ALARM")

mob_call_df <- emmeans(m.mob.call, ~ SPECIES * TREATMENT, type = "response") |>
  as.data.frame() |>
  mutate(BEHAVIOR = "MOB")

pred_df <- bind_rows(alarm_call_df, mob_call_df)

# Raw data in long format
raw_df <- sjdf_clean |>
  pivot_longer(
    cols = c(NUMBER.ALARM, NUMBER.MOB),
    names_to = "BEHAVIOR",
    values_to = "COUNT"
  ) |>
  mutate(
    BEHAVIOR = recode(BEHAVIOR,
                      "NUMBER.ALARM" = "ALARM",
                      "NUMBER.MOB" = "MOB")
  )

# Facet labels
label_df <- data.frame(
  BEHAVIOR = c("ALARM", "MOB"),
  label = c("a)", "b)"),
  x = -Inf,
  y = Inf
)


# Make TREATMENT a factor so Control comes first
pred_df$TREATMENT <- factor(pred_df$TREATMENT, levels = c("CONTROL", "HAWK"))
raw_df$TREATMENT <- factor(raw_df$TREATMENT, levels = c("CONTROL", "HAWK"))

library("scales")

ggplot() +
  
  # Raw data points (jittered)
  geom_jitter(
    data = raw_df,
    aes(
      x = SPECIES,
      y = COUNT,
      color = TREATMENT
    ),
 
    alpha = 0.1,
    shape = 16,
    size = 2,
    position = position_jitterdodge(
      dodge.width = 0.6,  
      jitter.width = 0.5, 
      jitter.height = 0    
    ),
    show.legend = FALSE
  ) +
  
  geom_violin(
    data = raw_df,
    aes(
      x = SPECIES,
      y = COUNT,
      fill = TREATMENT
    ),
    position = position_dodge(width = 0.6),  
    alpha = 0.1, linetype = 0,                             
    show.legend = FALSE,
    scale = "width"
  ) + 
  
  # Predicted points
  geom_point(
    data = pred_df,
    aes(
      x = SPECIES,
      y = response,
      color = TREATMENT,
      shape = TREATMENT
    ),
    size = 3,
    position = position_dodge(width = 0.6)
  ) +
  
  # Predicted error bars
  geom_errorbar(
    data = pred_df,
    aes(
      x = SPECIES,
      ymin = asymp.LCL,  # emmeans CI names
      ymax = asymp.UCL,
      color = TREATMENT
    ),
    width = 0.3,
    position = position_dodge(width = 0.6),
    show.legend = FALSE
  ) +
  
  # Panel labels
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.3,
    vjust = 1.5,
    size = 5
  ) +
  
  # Facets
  facet_wrap(~ BEHAVIOR, ncol = 2) +
  
  # Colors and shapes
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
  
  # X axis
  scale_x_discrete(labels = c(
    CASJ = expression(italic("A. californica")),
    ISSJ = expression(italic("A. insularis"))
  )) +
  
  # Log scale with space below zero
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1),  # smooth linear near 0, log for higher values
    expand = expansion(mult = c(0.05, 0.1)),  # adds some space below 0 and above max
    labels = comma_format()
  ) +
  
  # Theme
  theme_classic(base_size = 14) +
  labs(
    x = "Species",
    y = "Number of calls"
  ) +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2)
  )



ggsave("fig2.png", last_plot(), width = 8, height = 5, dpi = 300)


# ----- alarm post-hoc ------
em_alarm <- emmeans(m.alarm.call, ~ SPECIES * TREATMENT)

#differences within species
alarm_treatment_diff <- contrast(em_alarm, "pairwise", by = "SPECIES", type = "response")  
alarm_treatment_diff

#difference of the differences
diff_of_diffs_alarm <- contrast(em_alarm, interaction = "trt.vs.ctrl", type = "response" )
diff_of_diffs_alarm

# ----- MOB post-hoc ------

# Estimated marginal means for SPECIES × TREATMENT
em_mob <- emmeans(m.mob.call, ~ SPECIES * TREATMENT)

# Differences within each species: (Hawk - Control)
mob_treatment_diff <- contrast(em_mob, "pairwise", by = "SPECIES", type = "response")
mob_treatment_diff

# Difference of the differences: (Hawk-Control for ISSJ) - (Hawk-Control for CASJ)
mob_diff_of_diffs <- contrast(em_mob, interaction = "trt.vs.ctrl", type = "response")
mob_diff_of_diffs




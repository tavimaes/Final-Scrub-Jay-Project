#this script is for analyzing the effect of species on probability 
#to observe alarming using Bayesian models


####----Load packages----

library(tidyverse)
library(ggplot2)
library(dplyr)
library(brms)
library(tidybayes)
library(marginaleffects)

#read in clean data
sjdf_clean <- read_csv("clean_data.csv")

#set variables classes
sjdf_clean <- sjdf_clean |>
  mutate(across(c(ALARM, MOB, FLEE, INTEREST, HYPOTENUSE, LATENCY.ALARM, 
                  GROUP.SIZE, NUMBER.MOB, NUMBER.ALARM), as.numeric)) |>
  mutate(across(c(SPECIES, SEASON, PLAYBACK, SITE, SUBSITE, TREATMENT), 
                as.factor)) |>
  #scale continuous predictors 
  mutate(GROUP.SIZE = as.numeric(scale(GROUP.SIZE))) |>
  mutate(HYPOTENUSE = as.numeric(scale(HYPOTENUSE)))

####----Data analysis begins here----

####----ALARM----

#set priors 
priors <- c(
  prior(normal(0, 1.5), class = "b"),
  prior(normal(0, 1.5), class = "Intercept"),
  prior(exponential(1), class = "sd")
)

#set seed for reproducibility
set.seed(123)

m_alarm_v1 <- brm(
  ALARM ~ SPECIES * TREATMENT + GROUP.SIZE + PLAYBACK + HYPOTENUSE + 
    (1 | SUBSITE),
  data = sjdf_clean,
  family = bernoulli(link = "logit"),
  chains = 4,
  cores = 4,
  prior = priors,
  iter = 4000
)

saveRDS(m_alarm_v1, "models/m_alarm_v1.rds")

#check model convergence
summary(m_alarm_v1)
plot(m_alarm_v1)

#check model fit
pp_check(m_alarm_v1, type = "bars")
check_collinearity(m_alarm_v1)
mcmc_acf(m_alarm_v1)

#save playback stats
playback_effect <- as_tibble(fixef(m_alarm_v1), rownames = "term") |>
  filter(grepl("PLAYBACK", term)) |>
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

#rerun without playback and with season
set.seed(123)

m_alarm_v2 <- brm(
  ALARM ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + HYPOTENUSE + 
    (1 | SUBSITE),
  data = sjdf_clean,
  family = bernoulli(link = "logit"),
  chains = 4,
  cores = 4,
  prior = priors,
  iter = 4000
)

saveRDS(m_alarm_v2, "models/m_alarm_v2.rds")

#check model convergence
summary(m_alarm_v2)
plot(m_alarm_v2)

#check model fit
pp_check(m_alarm_v2, type = "bars")

hypotenuse_effect <- as_tibble(fixef(m_alarm_v2), rownames = "term") |>
  filter(grepl("HYPOTENUSE", term)) |>
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

#rerun final model without hypotenuse
set.seed(123)

m_alarm_final <- brm(
  ALARM ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + 
    (1 | SUBSITE),
  data = sjdf_clean,
  family = bernoulli(link = "logit"),
  chains = 4,
  cores = 4,
  prior = priors,
  iter = 4000
)

saveRDS(m_alarm_final, "models/m_alarm_final.rds")

#check model convergence
summary(m_alarm_final)
plot(m_alarm_final)

#check model fit
pp_check(m_alarm_final, type = "bars")
check_collinearity(m_alarm_final)

#save all other stats
coef_df <- as_tibble(fixef(m_alarm_final), rownames = "term") |>
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

alarm_coefs <- bind_rows(coef_df, hypotenuse_effect, playback_effect)

#coef plot for alarm
ggplot(alarm_coefs, aes(x = estimate, y = reorder(term, estimate))) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(
    x = "Effect size (log-odds)",
    y = "Predictor",
    title = "Effects on Alarm Probability"
  )

fitted_draws <- sjdf_clean %>%
  select(SEASON, GROUP.SIZE, SUBSITE) %>%
  expand_grid(
    SPECIES = c("CASJ", "ISSJ"),
    TREATMENT = c("CONTROL", "HAWK")
  ) %>%
  add_epred_draws(m_alarm_final, re_formula = NA, ndraws = 1000) %>%
  group_by(SPECIES, TREATMENT, .draw) %>%
  summarise(.epred = mean(.epred), .groups = "drop")

plot1 <- fitted_draws %>%
  group_by(SPECIES, TREATMENT) %>%
  median_qi(.epred, .width = 0.95) %>%
  ggplot(aes(x = SPECIES, y = .epred, color = TREATMENT, group = TREATMENT)) +
  geom_point(position = position_dodge(0.3), size = 3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper),
                position = position_dodge(0.3), width = 0.15) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  labs(
    y = "Predicted probability of alarm",
    x = "Species",
    color = "Treatment"
  ) +
  theme_classic()

plot1

plot2 <- fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(SPECIES) %>%
  median_qi(difference, .width = 0.95) %>%
  ggplot(aes(x = SPECIES, y = difference, color = SPECIES)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    y = "Difference in alarm probability (Hawk − Control)",
    x = "Species",
    color = "Species"
  ) +
  theme_classic() +
  theme(legend.position = "none")

plot2

plot3 <- fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES = as.character(SPECIES)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(treatment_effect = HAWK - CONTROL) %>%
  select(SPECIES, .draw, treatment_effect) %>%
  pivot_wider(names_from = SPECIES, values_from = treatment_effect) %>%
  mutate(interaction = ISSJ - CASJ) %>%
  median_qi(interaction, .width = 0.95) %>%
  ggplot(aes(x = "ISSJ − CASJ", y = interaction)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    y = "Difference in treatment effect (ISSJ − CASJ)",
    x = "",
    title = "Contrast of contrasts"
  ) +
  theme_classic()

plot3

#compare difference of differences for treatment effect by species
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES = as.character(SPECIES)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(treatment_effect = HAWK - CONTROL) %>%
  select(SPECIES, .draw, treatment_effect) %>%
  pivot_wider(names_from = SPECIES, values_from = treatment_effect) %>%
  mutate(interaction = ISSJ - CASJ) %>%
  median_qi(interaction, .width = 0.95)

fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES = as.character(SPECIES)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(treatment_effect = HAWK - CONTROL) %>%
  select(SPECIES, .draw, treatment_effect) %>%
  pivot_wider(names_from = SPECIES, values_from = treatment_effect) %>%
  mutate(interaction = ISSJ - CASJ) %>%
  summarise(p_greater = mean(interaction > 0))

#for treatments by species
# Median and CI
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = ISSJ - CASJ) %>%
  group_by(TREATMENT) %>%
  median_qi(difference, .width = 0.95)

# Posterior probability ISSJ > CASJ per treatment
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = ISSJ - CASJ) %>%
  group_by(TREATMENT) %>%
  summarise(p_greater = mean(difference > 0))
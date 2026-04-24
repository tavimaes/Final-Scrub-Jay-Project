#this script is for analyzing the effect of species on probability 
#to observe alarming using Bayesian models


####----Load packages----

library(tidyverse)
library(ggplot2)
library(dplyr)
library(brms)
library(tidybayes)
library(marginaleffects)
library(performance)
library(bayesplot)

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

####----ALARM PROBABILITY ANALYSIS ----

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
  iter = 4000,
  file = "models/m_alarm_v1"
)

#check model convergence
summary(m_alarm_v1)
plot(m_alarm_v1)

#check model fit
pp_check(m_alarm_v1, type = "bars")
check_collinearity(m_alarm_v1)


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
  iter = 4000,
  file = "models/m_alarm_v2"
)

#check model convergence
summary(m_alarm_v2)
plot(m_alarm_v2)

#check model fit
pp_check(m_alarm_v2, type = "bars")
check_collinearity(m_alarm_v2)

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
  iter = 4000,
  file = "models/m_alarm_final"
)

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

#save all coefficients
alarm_coefs <- bind_rows(coef_df, hypotenuse_effect, playback_effect) |>
  filter(term != "Intercept") |>
  mutate(model = "alarm_prob")

#save coefficient data
write_csv(alarm_coefs, file = "alarm_prob_coefs.csv")


#generate data to plot
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
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  median_qi(difference, .width = 0.95) %>%
  ggplot(aes(x = TREATMENT, y = difference, color = TREATMENT)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    y = "Difference in alarm probability (CASJ − ISSJ)",
    x = "Treatment",
    color = "Treatment"
  ) +
  theme_classic() +
  theme(legend.position = "none")

plot3

# ── Plot 2 numerical summaries: treatment effect (HAWK - CONTROL) per species ──

# Median and CI
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(SPECIES) %>%
  median_qi(difference, .width = 0.95)

# Posterior probability treatment effect > 0 per species
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(SPECIES) %>%
  summarise(p_greater = mean(difference > 0))


# ── Plot 3 numerical summaries: species difference (CASJ - ISSJ) per treatment ──

# Median and CI
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  median_qi(difference, .width = 0.95)

# Posterior probability CASJ > ISSJ per treatment
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  summarise(p_greater = mean(difference > 0))














###----MOB PROBABILITY ANALYSIS -----

m_mob_v1 <- brm(
  MOB ~ SPECIES * TREATMENT + GROUP.SIZE + PLAYBACK + HYPOTENUSE + 
    (1 | SUBSITE),
  data = sjdf_clean,
  family = bernoulli(link = "logit"),
  chains = 4,
  cores = 4,
  prior = priors,
  iter = 4000,
  file = "models/m_mob_v1"
)

#check model convergence
summary(m_mob_v1)
plot(m_mob_v1)

#check model fit
pp_check(m_mob_v1, type = "bars")
check_collinearity(m_mob_v1)


#save playback stats
playback_effect <- as_tibble(fixef(m_mob_v1), rownames = "term") |>
  filter(grepl("PLAYBACK", term)) |>
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

#rerun without playback and with season
set.seed(123)

m_mob_v2 <- brm(
  MOB ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + HYPOTENUSE + 
    (1 | SUBSITE),
  data = sjdf_clean,
  family = bernoulli(link = "logit"),
  chains = 4,
  cores = 4,
  prior = priors,
  iter = 4000,
  file = "models/m_mob_v2"
)

#check model convergence
summary(m_mob_v2)
plot(m_mob_v2)

#check model fit
pp_check(m_mob_v2, type = "bars")
check_collinearity(m_mob_v2)

hypotenuse_effect <- as_tibble(fixef(m_mob_v2), rownames = "term") |>
  filter(grepl("HYPOTENUSE", term)) |>
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

#rerun final model without hypotenuse
set.seed(123)

m_mob_final <- brm(
  MOB ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + 
    (1 | SUBSITE),
  data = sjdf_clean,
  family = bernoulli(link = "logit"),
  chains = 4,
  cores = 4,
  prior = priors,
  iter = 4000,
  file = "models/m_mob_final"
)

#check model convergence
summary(m_mob_final)
plot(m_mob_final)

#check model fit
pp_check(m_mob_final, type = "bars")
check_collinearity(m_mob_final)

#save all other stats
coef_df <- as_tibble(fixef(m_mob_final), rownames = "term") |>
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

#save all coefficients
mob_coefs <- bind_rows(coef_df, hypotenuse_effect, playback_effect) |>
  filter(term != "Intercept") |>
  mutate(model = "mob_prob")

#save coefficient data
write_csv(mob_coefs, file = "mob_prob_coefs.csv")


#generate data to plot
fitted_draws <- sjdf_clean %>%
  select(SEASON, GROUP.SIZE, SUBSITE) %>%
  expand_grid(
    SPECIES = c("CASJ", "ISSJ"),
    TREATMENT = c("CONTROL", "HAWK")
  ) %>%
  add_epred_draws(m_mob_final, re_formula = NA, ndraws = 1000) %>%
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
    y = "Predicted probability of mobbing",
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
    y = "Difference in mob probability (Hawk − Control)",
    x = "Species",
    color = "Species"
  ) +
  theme_classic() +
  theme(legend.position = "none")

plot2

plot3 <- fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  median_qi(difference, .width = 0.95) %>%
  ggplot(aes(x = TREATMENT, y = difference, color = TREATMENT)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    y = "Difference in mob probability (CASJ − ISSJ)",
    x = "Treatment",
    color = "Treatment"
  ) +
  theme_classic() +
  theme(legend.position = "none")

plot3

# ── Plot 2 numerical summaries: treatment effect (HAWK - CONTROL) per species ──

# Median and CI
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(SPECIES) %>%
  median_qi(difference, .width = 0.95)

# Posterior probability treatment effect > 0 per species
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(SPECIES) %>%
  summarise(p_greater = mean(difference > 0))


# ── Plot 3 numerical summaries: species difference (CASJ - ISSJ) per treatment ──

# Median and CI
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  median_qi(difference, .width = 0.95)

# Posterior probability CASJ > ISSJ per treatment
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  summarise(p_greater = mean(difference > 0))








###----INTEREST PROBABILITY ANALYSIS -----

m_interest_v1 <- brm(
  INTEREST ~ SPECIES * TREATMENT + GROUP.SIZE + PLAYBACK + HYPOTENUSE + 
    (1 | SUBSITE),
  data = sjdf_clean,
  family = bernoulli(link = "logit"),
  chains = 4,
  cores = 4,
  prior = priors,
  iter = 4000,
  file = "models/m_interest_v1"
)

#check model convergence
summary(m_interest_v1)
plot(m_interest_v1)

#check model fit
pp_check(m_interest_v1, type = "bars")
check_collinearity(m_interest_v1)

#save playback stats
playback_effect <- as_tibble(fixef(m_interest_v1), rownames = "term") |>
  filter(grepl("PLAYBACK", term)) |>
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

#rerun without playback and with season
set.seed(123)

m_interest_v2 <- brm(
  INTEREST ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + HYPOTENUSE + 
    (1 | SUBSITE),
  data = sjdf_clean,
  family = bernoulli(link = "logit"),
  chains = 4,
  cores = 4,
  prior = priors,
  iter = 4000,
  file = "models/m_interest_v2"
)

#check model convergence
summary(m_interest_v2)
plot(m_interest_v2)

#check model fit
pp_check(m_interest_v2, type = "bars")
check_collinearity(m_interest_v2)

hypotenuse_effect <- as_tibble(fixef(m_interest_v2), rownames = "term") |>
  filter(grepl("HYPOTENUSE", term)) |>
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

#rerun final model without hypotenuse
set.seed(123)

m_interest_final <- brm(
  INTEREST ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + 
    (1 | SUBSITE),
  data = sjdf_clean,
  family = bernoulli(link = "logit"),
  chains = 4,
  cores = 4,
  prior = priors,
  iter = 4000,
  file = "models/m_interest_final"
)

#check model convergence
summary(m_interest_final)
plot(m_interest_final)

#check model fit
pp_check(m_interest_final, type = "bars")
check_collinearity(m_interest_final)

#save all other stats
coef_df <- as_tibble(fixef(m_interest_final), rownames = "term") |>
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

#save all coefficients
interest_coefs <- bind_rows(coef_df, hypotenuse_effect, playback_effect) |>
  filter(term != "Intercept") |>
  mutate(model = "interest_prob")

#save coefficient data
write_csv(interest_coefs, file = "interest_prob_coefs.csv")

#generate data to plot
fitted_draws <- sjdf_clean %>%
  select(SEASON, GROUP.SIZE, SUBSITE) %>%
  expand_grid(
    SPECIES = c("CASJ", "ISSJ"),
    TREATMENT = c("CONTROL", "HAWK")
  ) %>%
  add_epred_draws(m_interest_final, re_formula = NA, ndraws = 1000) %>%
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
    y = "Predicted probability of interest",
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
    y = "Difference in interest probability (Hawk − Control)",
    x = "Species",
    color = "Species"
  ) +
  theme_classic() +
  theme(legend.position = "none")

plot2

plot3 <- fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  median_qi(difference, .width = 0.95) %>%
  ggplot(aes(x = TREATMENT, y = difference, color = TREATMENT)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    y = "Difference in interest probability (CASJ − ISSJ)",
    x = "Treatment",
    color = "Treatment"
  ) +
  theme_classic() +
  theme(legend.position = "none")

plot3

# ── Plot 2 numerical summaries: treatment effect (HAWK - CONTROL) per species ──

# Median and CI
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(SPECIES) %>%
  median_qi(difference, .width = 0.95)

# Posterior probability treatment effect > 0 per species
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(SPECIES) %>%
  summarise(p_greater = mean(difference > 0))

# ── Plot 3 numerical summaries: species difference (CASJ - ISSJ) per treatment ──

# Median and CI
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  median_qi(difference, .width = 0.95)

# Posterior probability CASJ > ISSJ per treatment
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  summarise(p_greater = mean(difference > 0))









###----FLEE PROBABILITY ANALYSIS -----

m_flee_v1 <- brm(
  FLEE ~ SPECIES * TREATMENT + GROUP.SIZE + PLAYBACK + HYPOTENUSE + 
    (1 | SUBSITE),
  data = sjdf_clean,
  family = bernoulli(link = "logit"),
  chains = 4,
  cores = 4,
  prior = priors,
  iter = 4000,
  file = "models/m_flee_v1"
)

#check model convergence
summary(m_flee_v1)
plot(m_flee_v1)

#check model fit
pp_check(m_flee_v1, type = "bars")
check_collinearity(m_flee_v1)

#save playback stats
playback_effect <- as_tibble(fixef(m_flee_v1), rownames = "term") |>
  filter(grepl("PLAYBACK", term)) |>
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

#rerun without playback and with season
set.seed(123)

m_flee_v2 <- brm(
  FLEE ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + HYPOTENUSE + 
    (1 | SUBSITE),
  data = sjdf_clean,
  family = bernoulli(link = "logit"),
  chains = 4,
  cores = 4,
  prior = priors,
  iter = 4000,
  file = "models/m_flee_v2"
)

#check model convergence
summary(m_flee_v2)
plot(m_flee_v2)

#check model fit
pp_check(m_flee_v2, type = "bars")
check_collinearity(m_flee_v2)

hypotenuse_effect <- as_tibble(fixef(m_flee_v2), rownames = "term") |>
  filter(grepl("HYPOTENUSE", term)) |>
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

#rerun final model without hypotenuse
set.seed(123)

m_flee_final <- brm(
  FLEE ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + 
    (1 | SUBSITE),
  data = sjdf_clean,
  family = bernoulli(link = "logit"),
  chains = 4,
  cores = 4,
  prior = priors,
  iter = 4000,
  file = "models/m_flee_final"
)

#check model convergence
summary(m_flee_final)
plot(m_flee_final)

#check model fit
pp_check(m_flee_final, type = "bars")
check_collinearity(m_flee_final)

#save all other stats
coef_df <- as_tibble(fixef(m_flee_final), rownames = "term") |>
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

#save all coefficients
flee_coefs <- bind_rows(coef_df, hypotenuse_effect, playback_effect) |>
  filter(term != "Intercept") |>
  mutate(model = "flee_prob")

#save coefficient data
write_csv(flee_coefs, file = "flee_prob_coefs.csv")

#generate data to plot
fitted_draws <- sjdf_clean %>%
  select(SEASON, GROUP.SIZE, SUBSITE) %>%
  expand_grid(
    SPECIES = c("CASJ", "ISSJ"),
    TREATMENT = c("CONTROL", "HAWK")
  ) %>%
  add_epred_draws(m_flee_final, re_formula = NA, ndraws = 1000) %>%
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
    y = "Predicted probability of fleeing",
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
    y = "Difference in flee probability (Hawk − Control)",
    x = "Species",
    color = "Species"
  ) +
  theme_classic() +
  theme(legend.position = "none")

plot2

plot3 <- fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  median_qi(difference, .width = 0.95) %>%
  ggplot(aes(x = TREATMENT, y = difference, color = TREATMENT)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    y = "Difference in flee probability (CASJ − ISSJ)",
    x = "Treatment",
    color = "Treatment"
  ) +
  theme_classic() +
  theme(legend.position = "none")

plot3

# ── Plot 2 numerical summaries: treatment effect (HAWK - CONTROL) per species ──

# Median and CI
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(SPECIES) %>%
  median_qi(difference, .width = 0.95)

# Posterior probability treatment effect > 0 per species
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(SPECIES) %>%
  summarise(p_greater = mean(difference > 0))

# ── Plot 3 numerical summaries: species difference (CASJ - ISSJ) per treatment ──

# Median and CI
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  median_qi(difference, .width = 0.95)

# Posterior probability CASJ > ISSJ per treatment
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  summarise(p_greater = mean(difference > 0))

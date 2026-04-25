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
library(patchwork)

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

#save for final plot
alarm_draws <- fitted_draws %>% mutate(model = "Alarm")












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


#save for final plot
mob_draws <- fitted_draws %>% mutate(model = "Mob")





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


#save for final plot
interest_draws <- fitted_draws %>% mutate(model = "Interest")






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

#save for final plot
flee_draws <- fitted_draws %>% mutate(model = "Flee")




#----FINAL PLOTTING----




# ---- SETUP ----
all_draws <- bind_rows(alarm_draws, mob_draws, interest_draws, flee_draws) %>%
  mutate(
    model = fct_relevel(model, "Alarm", "Mob", "Flee", "Interest"),
    SPECIES = factor(SPECIES, levels = c("CASJ", "ISSJ")),
    TREATMENT = factor(TREATMENT, levels = c("CONTROL", "HAWK"))
  )

raw_proportions <- bind_rows(
  sjdf_clean %>% group_by(SPECIES, TREATMENT) %>%
    summarise(prop = mean(ALARM), n = n(), .groups = "drop") %>% mutate(model = "Alarm"),
  sjdf_clean %>% group_by(SPECIES, TREATMENT) %>%
    summarise(prop = mean(MOB), n = n(), .groups = "drop") %>% mutate(model = "Mob"),
  sjdf_clean %>% group_by(SPECIES, TREATMENT) %>%
    summarise(prop = mean(FLEE), n = n(), .groups = "drop") %>% mutate(model = "Flee"),
  sjdf_clean %>% group_by(SPECIES, TREATMENT) %>%
    summarise(prop = mean(INTEREST), n = n(), .groups = "drop") %>% mutate(model = "Interest")
) %>%
  mutate(
    model = fct_relevel(model, "Alarm", "Mob", "Flee", "Interest"),
    SPECIES = factor(SPECIES, levels = c("CASJ", "ISSJ")),
    TREATMENT = factor(TREATMENT, levels = c("CONTROL", "HAWK"))
  ) |>
  mutate(group = interaction(SPECIES, TREATMENT, sep = "."))

dodge <- position_dodge(width = 0.8)

tag_theme <- theme(
  plot.tag = element_text(face = "bold", size = 15),
  plot.tag.position = c(0, 0.98)
)

common_theme <- theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(0.9, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "none"
  )

diff_scale <- scale_y_continuous(
  limits = c(-0.75, 0.75),
  labels = scales::percent_format()
)

sig_color <- scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray60"))

# ---- PANEL A ----
p1 <- all_draws %>%
  mutate(group = interaction(SPECIES, TREATMENT)) %>%
  group_by(model, SPECIES, TREATMENT, group) %>%
  median_qi(.epred, .width = 0.95) %>%
  ggplot(aes(x = SPECIES, y = .epred, group = group, color = group)) +
  geom_col(
    data = raw_proportions,
    aes(x = SPECIES, y = prop, fill = group, group = group),
    position = dodge,
    width = 0.7,
    alpha = 0.2,   
    inherit.aes = FALSE
  ) +
  geom_point(position = dodge, size = 2.5) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper),
                position = dodge, width = 0.15, linewidth = 0.6) +
  geom_text(
    data = raw_proportions %>% mutate(group = interaction(SPECIES, TREATMENT)),
    aes(x = SPECIES, y = Inf, label = n, group = group),
    position = dodge, vjust = -0.5, size = 3, inherit.aes = FALSE
  ) +
  scale_color_manual(values = c(
    "CASJ.CONTROL" = "plum3", "CASJ.HAWK" = "darkorchid4",
    "ISSJ.CONTROL" = "lightblue3", "ISSJ.HAWK" = "dodgerblue3"
  )) +
  scale_fill_manual(values = c(
    "CASJ.CONTROL" = "plum3",
    "CASJ.HAWK" = "darkorchid4",
    "ISSJ.CONTROL" = "lightblue3",
    "ISSJ.HAWK" = "dodgerblue3"
  )) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.08)),
    labels = scales::percent_format(accuracy = 1)
  ) +
  facet_wrap(~ model, ncol = 1) +
  labs(y = "Response probability", x = "Species", tag = "a)") +
  coord_cartesian(clip = "off") +
  common_theme + tag_theme

# ---- PANEL B ----
p2 <- all_draws %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(model, SPECIES) %>%
  median_qi(difference, .width = 0.95) %>%
  mutate(significant = .lower > 0 | .upper < 0) %>%
  ggplot(aes(x = SPECIES, y = difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = significant), size = 2.5) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper, color = significant), width = 0.15) +
  sig_color +
  diff_scale +
  facet_wrap(~ model, ncol = 1) +
  labs(y = "Percentage difference by treatment (Hawk − Control)", x = "Species", tag = "b)") +
  common_theme + tag_theme

# ---- PANEL C ----
p3 <- all_draws %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(model, TREATMENT) %>%
  median_qi(difference, .width = 0.95) %>%
  mutate(significant = .lower > 0 | .upper < 0) %>%
  ggplot(aes(x = TREATMENT, y = difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = significant), size = 2.5) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper, color = significant), width = 0.15) +
  sig_color +
  diff_scale +
  facet_wrap(~ model, ncol = 1) +
  labs(y = "Percentage difference by species (CASJ − ISSJ)", x = "Treatment", tag = "c)") +
  common_theme + tag_theme

p3 <- p3 +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# ---- COMBINE ----
final_plot <- (p1 | p2 | p3) +
  plot_layout(widths = c(2.2, 1.2, 1.2))

final_plot

ggsave("fig5.png", final_plot, width = 10, height = 10, dpi = 300)


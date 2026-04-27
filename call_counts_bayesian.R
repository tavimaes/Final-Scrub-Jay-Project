# this script is for running byesian models on call count data

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


###----NUMBER OF ALARM CALLS ANALYSIS -----

m_nalarm_v1 <- brm(
  NUMBER.ALARM ~ SPECIES * TREATMENT + GROUP.SIZE + PLAYBACK + HYPOTENUSE + 
    (1 | SUBSITE),
  data = sjdf_clean,
  family = negbinomial(),
  chains = 4,
  cores = 4,
  prior = priors,
  iter = 4000,
  file = "models/m_nalarm_v1"
)

#check model convergence
summary(m_nalarm_v1)
plot(m_nalarm_v1)

#check model fit
pp_check(m_nalarm_v1)
check_collinearity(m_nalarm_v1)

#save playback stats
playback_effect <- as_tibble(fixef(m_nalarm_v1), rownames = "term") |>
  filter(grepl("PLAYBACK", term)) |>
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

#rerun without playback and with season
set.seed(123)

m_nalarm_v2 <- brm(
  NUMBER.ALARM ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + HYPOTENUSE + 
    (1 | SUBSITE),
  data = sjdf_clean,
  family = negbinomial(),
  chains = 4,
  cores = 4,
  prior = priors,
  iter = 4000,
  file = "models/m_nalarm_v2"
)

#check model convergence
summary(m_nalarm_v2)
plot(m_nalarm_v2)

#check model fit
pp_check(m_nalarm_v2)
check_collinearity(m_nalarm_v2)

hypotenuse_effect <- as_tibble(fixef(m_nalarm_v2), rownames = "term") |>
  filter(grepl("HYPOTENUSE", term)) |>
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

#rerun final model without hypotenuse
set.seed(123)

m_nalarm_final <- brm(
  NUMBER.ALARM ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + 
    (1 | SUBSITE),
  data = sjdf_clean,
  family = negbinomial(),
  chains = 4,
  cores = 4,
  prior = priors,
  iter = 4000,
  file = "models/m_nalarm_final"
)

#check model convergence
summary(m_nalarm_final)
plot(m_nalarm_final)

#check model fit
pp_check(m_nalarm_final)
check_collinearity(m_nalarm_final)

#save all other stats
coef_df <- as_tibble(fixef(m_nalarm_final), rownames = "term") |>
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

#save all coefficients
nalarm_coefs <- bind_rows(coef_df, hypotenuse_effect, playback_effect) |>
  filter(term != "Intercept") |>
  mutate(model = "number_alarm")

#save coefficient data
write_csv(nalarm_coefs, file = "number_alarm_coefs.csv")

#generate data to plot
fitted_draws <- sjdf_clean %>%
  select(SEASON, GROUP.SIZE, SUBSITE) %>%
  expand_grid(
    SPECIES = c("CASJ", "ISSJ"),
    TREATMENT = c("CONTROL", "HAWK")
  ) %>%
  add_epred_draws(m_nalarm_final, re_formula = NA, ndraws = 1000) %>%
  group_by(SPECIES, TREATMENT, .draw) %>%
  summarise(.epred = mean(.epred), .groups = "drop")

plot1 <- fitted_draws %>%
  group_by(SPECIES, TREATMENT) %>%
  median_qi(.epred, .width = 0.95) %>%
  ggplot(aes(x = SPECIES, y = .epred, color = TREATMENT, group = TREATMENT)) +
  geom_point(position = position_dodge(0.3), size = 3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper),
                position = position_dodge(0.3), width = 0.15) +
  labs(
    y = "Predicted number of alarm calls",
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
  labs(
    y = "Difference in alarm calls (Hawk − Control)",
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
  labs(
    y = "Difference in alarm calls (CASJ − ISSJ)",
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






###----NUMBER OF MOB CALLS ANALYSIS -----

m_nmob_v1 <- brm(
  NUMBER.MOB ~ SPECIES * TREATMENT + GROUP.SIZE + PLAYBACK + HYPOTENUSE + 
    (1 | SUBSITE),
  data = sjdf_clean,
  family = negbinomial(),
  chains = 4,
  cores = 4,
  prior = priors,
  iter = 4000,
  file = "models/m_nmob_v1"
)

#check model convergence
summary(m_nmob_v1)
plot(m_nmob_v1)

#check model fit
pp_check(m_nmob_v1)
check_collinearity(m_nmob_v1)

#save playback stats
playback_effect <- as_tibble(fixef(m_nmob_v1), rownames = "term") |>
  filter(grepl("PLAYBACK", term)) |>
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

#rerun without playback and with season
set.seed(123)

m_nmob_v2 <- brm(
  NUMBER.MOB ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + HYPOTENUSE + 
    (1 | SUBSITE),
  data = sjdf_clean,
  family = negbinomial(),
  chains = 4,
  cores = 4,
  prior = priors,
  iter = 4000,
  file = "models/m_nmob_v2"
)

#check model convergence
summary(m_nmob_v2)
plot(m_nmob_v2)

#check model fit
pp_check(m_nmob_v2)
check_collinearity(m_nmob_v2)

hypotenuse_effect <- as_tibble(fixef(m_nmob_v2), rownames = "term") |>
  filter(grepl("HYPOTENUSE", term)) |>
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

#rerun final model without hypotenuse
set.seed(123)

m_nmob_final <- brm(
  NUMBER.MOB ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + 
    (1 | SUBSITE),
  data = sjdf_clean,
  family = negbinomial(),
  chains = 4,
  cores = 4,
  prior = priors,
  iter = 4000,
  file = "models/m_nmob_final"
)

#check model convergence
summary(m_nmob_final)
plot(m_nmob_final)

#check model fit
pp_check(m_nmob_final)
check_collinearity(m_nmob_final)

#save all other stats
coef_df <- as_tibble(fixef(m_nmob_final), rownames = "term") |>
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

#save all coefficients
nmob_coefs <- bind_rows(coef_df, hypotenuse_effect, playback_effect) |>
  filter(term != "Intercept") |>
  mutate(model = "number_mob")

#save coefficient data
write_csv(nmob_coefs, file = "number_mob_coefs.csv")

#generate data to plot
fitted_draws <- sjdf_clean %>%
  select(SEASON, GROUP.SIZE, SUBSITE) %>%
  expand_grid(
    SPECIES = c("CASJ", "ISSJ"),
    TREATMENT = c("CONTROL", "HAWK")
  ) %>%
  add_epred_draws(m_nmob_final, re_formula = NA, ndraws = 1000) %>%
  group_by(SPECIES, TREATMENT, .draw) %>%
  summarise(.epred = mean(.epred), .groups = "drop")

plot1 <- fitted_draws %>%
  group_by(SPECIES, TREATMENT) %>%
  median_qi(.epred, .width = 0.95) %>%
  ggplot(aes(x = SPECIES, y = .epred, color = TREATMENT, group = TREATMENT)) +
  geom_point(position = position_dodge(0.3), size = 3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper),
                position = position_dodge(0.3), width = 0.15) +
  labs(
    y = "Predicted number of mob calls",
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
  labs(
    y = "Difference in mob calls (Hawk − Control)",
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
  labs(
    y = "Difference in mob calls (CASJ − ISSJ)",
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




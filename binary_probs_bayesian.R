#this script analyzes the effect of species and treatment on antipredator
#behaviour probability using Bayesian mixed effects models

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
library(ggtext)

#read in clean data
sjdf_clean <- read_csv("clean_data.csv")

#set variable classes and scale continuous predictors
sjdf_clean <- sjdf_clean |>
  mutate(across(c(ALARM, MOB, FLEE, INTEREST, HYPOTENUSE, LATENCY.ALARM,
                  GROUP.SIZE, NUMBER.MOB, NUMBER.ALARM), as.numeric)) |>
  mutate(across(c(SPECIES, SEASON, PLAYBACK, SITE, SUBSITE, TREATMENT), as.factor)) |>
  mutate(GROUP.SIZE = as.numeric(scale(GROUP.SIZE))) |>
  mutate(HYPOTENUSE = as.numeric(scale(HYPOTENUSE)))

####----Data analysis begins here----

#set priors — weakly informative
priors <- c(
  prior(normal(0, 1.5), class = "b"),
  prior(normal(0, 1.5), class = "Intercept"),
  prior(exponential(1), class = "sd")
)

####----ALARM PROBABILITY----

#v1: includes playback and hypotenuse
set.seed(123)
m_alarm_v1 <- brm(
  ALARM ~ SPECIES * TREATMENT + GROUP.SIZE + PLAYBACK + HYPOTENUSE + (1 | SUBSITE),
  data = sjdf_clean, family = bernoulli(link = "logit"),
  chains = 4, cores = 4, prior = priors, iter = 4000,
  file = "models/m_alarm_v1"
)

summary(m_alarm_v1)
plot(m_alarm_v1)
pp_check(m_alarm_v1, type = "bars")
check_collinearity(m_alarm_v1)

#save playback effect before dropping it
playback_effect <- as_tibble(fixef(m_alarm_v1), rownames = "term") |>
  filter(grepl("PLAYBACK", term)) |>
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5)

#v2: drop playback, add season, keep hypotenuse
set.seed(123)
m_alarm_v2 <- brm(
  ALARM ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + HYPOTENUSE + (1 | SUBSITE),
  data = sjdf_clean, family = bernoulli(link = "logit"),
  chains = 4, cores = 4, prior = priors, iter = 4000,
  file = "models/m_alarm_v2"
)

summary(m_alarm_v2)
plot(m_alarm_v2)
pp_check(m_alarm_v2, type = "bars")
check_collinearity(m_alarm_v2)

#save hypotenuse effect before dropping it
hypotenuse_effect <- as_tibble(fixef(m_alarm_v2), rownames = "term") |>
  filter(grepl("HYPOTENUSE", term)) |>
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5)

#final model: drop hypotenuse
set.seed(123)
m_alarm_final <- brm(
  ALARM ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + (1 | SUBSITE),
  data = sjdf_clean, family = bernoulli(link = "logit"),
  chains = 4, cores = 4, prior = priors, iter = 4000,
  file = "models/m_alarm_final"
)

summary(m_alarm_final)
plot(m_alarm_final)
pp_check(m_alarm_final, type = "bars")
check_collinearity(m_alarm_final)

#save all coefficients and write to csv
alarm_coefs <- as_tibble(fixef(m_alarm_final), rownames = "term") |>
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5) |>
  bind_rows(hypotenuse_effect, playback_effect) |>
  filter(term != "Intercept") |>
  mutate(model = "alarm_prob")

write_csv(alarm_coefs, file = "alarm_prob_coefs.csv")

#generate posterior predicted draws for plotting
fitted_draws <- sjdf_clean %>%
  select(SEASON, GROUP.SIZE, SUBSITE) %>%
  expand_grid(SPECIES = c("CASJ", "ISSJ"), TREATMENT = c("CONTROL", "HAWK")) %>%
  add_epred_draws(m_alarm_final, re_formula = NA, ndraws = 1000) %>%
  group_by(SPECIES, TREATMENT, .draw) %>%
  summarise(.epred = mean(.epred), .groups = "drop")

#numerical summaries: treatment effect (HAWK - CONTROL) per species
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(SPECIES) %>%
  median_qi(difference, .width = 0.95)

fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(SPECIES) %>%
  summarise(p_greater = mean(difference > 0))

#difference of differences
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES   = as.character(SPECIES)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  select(SPECIES, .draw, difference) %>%
  pivot_wider(names_from = SPECIES, values_from = difference) %>%
  mutate(diff_of_diffs = CASJ - ISSJ) %>%
  median_qi(diff_of_diffs, .width = 0.95)

fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES   = as.character(SPECIES)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  select(SPECIES, .draw, difference) %>%
  pivot_wider(names_from = SPECIES, values_from = difference) %>%
  mutate(diff_of_diffs = CASJ - ISSJ) %>%
  summarise(p_greater = mean(diff_of_diffs > 0))

#numerical summaries: species difference (CASJ - ISSJ) per treatment
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT), SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  median_qi(difference, .width = 0.95)

fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT), SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  summarise(p_greater = mean(difference > 0))

alarm_draws <- fitted_draws %>% mutate(model = "Alarm")

####----MOB PROBABILITY----

#v1: includes playback and hypotenuse
set.seed(123)
m_mob_v1 <- brm(
  MOB ~ SPECIES * TREATMENT + GROUP.SIZE + PLAYBACK + HYPOTENUSE + (1 | SUBSITE),
  data = sjdf_clean, family = bernoulli(link = "logit"),
  chains = 4, cores = 4, prior = priors, iter = 4000,
  file = "models/m_mob_v1"
)

summary(m_mob_v1)
plot(m_mob_v1)
pp_check(m_mob_v1, type = "bars")
check_collinearity(m_mob_v1)

playback_effect <- as_tibble(fixef(m_mob_v1), rownames = "term") |>
  filter(grepl("PLAYBACK", term)) |>
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5)

#v2: drop playback, add season, keep hypotenuse
set.seed(123)
m_mob_v2 <- brm(
  MOB ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + HYPOTENUSE + (1 | SUBSITE),
  data = sjdf_clean, family = bernoulli(link = "logit"),
  chains = 4, cores = 4, prior = priors, iter = 4000,
  file = "models/m_mob_v2"
)

summary(m_mob_v2)
plot(m_mob_v2)
pp_check(m_mob_v2, type = "bars")
check_collinearity(m_mob_v2)

hypotenuse_effect <- as_tibble(fixef(m_mob_v2), rownames = "term") |>
  filter(grepl("HYPOTENUSE", term)) |>
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5)

#final model: drop hypotenuse
set.seed(123)
m_mob_final <- brm(
  MOB ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + (1 | SUBSITE),
  data = sjdf_clean, family = bernoulli(link = "logit"),
  chains = 4, cores = 4, prior = priors, iter = 4000,
  file = "models/m_mob_final"
)

summary(m_mob_final)
plot(m_mob_final)
pp_check(m_mob_final, type = "bars")
check_collinearity(m_mob_final)

mob_coefs <- as_tibble(fixef(m_mob_final), rownames = "term") |>
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5) |>
  bind_rows(hypotenuse_effect, playback_effect) |>
  filter(term != "Intercept") |>
  mutate(model = "mob_prob")

write_csv(mob_coefs, file = "mob_prob_coefs.csv")

fitted_draws <- sjdf_clean %>%
  select(SEASON, GROUP.SIZE, SUBSITE) %>%
  expand_grid(SPECIES = c("CASJ", "ISSJ"), TREATMENT = c("CONTROL", "HAWK")) %>%
  add_epred_draws(m_mob_final, re_formula = NA, ndraws = 1000) %>%
  group_by(SPECIES, TREATMENT, .draw) %>%
  summarise(.epred = mean(.epred), .groups = "drop")

#numerical summaries: treatment effect per species
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(SPECIES) %>%
  median_qi(difference, .width = 0.95)

fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(SPECIES) %>%
  summarise(p_greater = mean(difference > 0))

#difference of differences
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES   = as.character(SPECIES)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  select(SPECIES, .draw, difference) %>%
  pivot_wider(names_from = SPECIES, values_from = difference) %>%
  mutate(diff_of_diffs = CASJ - ISSJ) %>%
  median_qi(diff_of_diffs, .width = 0.95)

fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES   = as.character(SPECIES)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  select(SPECIES, .draw, difference) %>%
  pivot_wider(names_from = SPECIES, values_from = difference) %>%
  mutate(diff_of_diffs = CASJ - ISSJ) %>%
  summarise(p_greater = mean(diff_of_diffs > 0))

#numerical summaries: species difference per treatment
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT), SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  median_qi(difference, .width = 0.95)

fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT), SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  summarise(p_greater = mean(difference > 0))

mob_draws <- fitted_draws %>% mutate(model = "Mob")

####----INTEREST PROBABILITY----

#v1: includes playback and hypotenuse
set.seed(123)
m_interest_v1 <- brm(
  INTEREST ~ SPECIES * TREATMENT + GROUP.SIZE + PLAYBACK + HYPOTENUSE + (1 | SUBSITE),
  data = sjdf_clean, family = bernoulli(link = "logit"),
  chains = 4, cores = 4, prior = priors, iter = 4000,
  file = "models/m_interest_v1"
)

summary(m_interest_v1)
plot(m_interest_v1)
pp_check(m_interest_v1, type = "bars")

check_collinearity(m_interest_v1)

playback_effect <- as_tibble(fixef(m_interest_v1), rownames = "term") |>
  filter(grepl("PLAYBACK", term)) |>
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5)

#v2: drop playback, add season, keep hypotenuse
set.seed(123)
m_interest_v2 <- brm(
  INTEREST ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + HYPOTENUSE + (1 | SUBSITE),
  data = sjdf_clean, family = bernoulli(link = "logit"),
  chains = 4, cores = 4, prior = priors, iter = 4000,
  file = "models/m_interest_v2"
)

summary(m_interest_v2)
plot(m_interest_v2)
pp_check(m_interest_v2, type = "bars")
check_collinearity(m_interest_v2)

hypotenuse_effect <- as_tibble(fixef(m_interest_v2), rownames = "term") |>
  filter(grepl("HYPOTENUSE", term)) |>
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5)

#final model: drop hypotenuse
set.seed(123)
m_interest_final <- brm(
  INTEREST ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + (1 | SUBSITE),
  data = sjdf_clean, family = bernoulli(link = "logit"),
  chains = 4, cores = 4, prior = priors, iter = 4000,
  file = "models/m_interest_final"
)

summary(m_interest_final)
plot(m_interest_final)
pp_check(m_interest_final, type = "bars")
check_collinearity(m_interest_final)

interest_coefs <- as_tibble(fixef(m_interest_final), rownames = "term") |>
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5) |>
  bind_rows(hypotenuse_effect, playback_effect) |>
  filter(term != "Intercept") |>
  mutate(model = "interest_prob")

write_csv(interest_coefs, file = "interest_prob_coefs.csv")

fitted_draws <- sjdf_clean %>%
  select(SEASON, GROUP.SIZE, SUBSITE) %>%
  expand_grid(SPECIES = c("CASJ", "ISSJ"), TREATMENT = c("CONTROL", "HAWK")) %>%
  add_epred_draws(m_interest_final, re_formula = NA, ndraws = 1000) %>%
  group_by(SPECIES, TREATMENT, .draw) %>%
  summarise(.epred = mean(.epred), .groups = "drop")

#numerical summaries: treatment effect per species
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(SPECIES) %>%
  median_qi(difference, .width = 0.95)

fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(SPECIES) %>%
  summarise(p_greater = mean(difference > 0))

#difference of differences
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES   = as.character(SPECIES)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  select(SPECIES, .draw, difference) %>%
  pivot_wider(names_from = SPECIES, values_from = difference) %>%
  mutate(diff_of_diffs = CASJ - ISSJ) %>%
  median_qi(diff_of_diffs, .width = 0.95)

fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES   = as.character(SPECIES)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  select(SPECIES, .draw, difference) %>%
  pivot_wider(names_from = SPECIES, values_from = difference) %>%
  mutate(diff_of_diffs = CASJ - ISSJ) %>%
  summarise(p_greater = mean(diff_of_diffs > 0))

#numerical summaries: species difference per treatment
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT), SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  median_qi(difference, .width = 0.95)

fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT), SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  summarise(p_greater = mean(difference > 0))

interest_draws <- fitted_draws %>% mutate(model = "Interest")

####----FLEE PROBABILITY----

#v1: includes playback and hypotenuse
set.seed(123)
m_flee_v1 <- brm(
  FLEE ~ SPECIES * TREATMENT + GROUP.SIZE + PLAYBACK + HYPOTENUSE + (1 | SUBSITE),
  data = sjdf_clean, family = bernoulli(link = "logit"),
  chains = 4, cores = 4, prior = priors, iter = 4000,
  file = "models/m_flee_v1"
)

summary(m_flee_v1)
plot(m_flee_v1)
pp_check(m_flee_v1, type = "bars")
check_collinearity(m_flee_v1)

playback_effect <- as_tibble(fixef(m_flee_v1), rownames = "term") |>
  filter(grepl("PLAYBACK", term)) |>
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5)

#v2: drop playback, add season, keep hypotenuse
set.seed(123)
m_flee_v2 <- brm(
  FLEE ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + HYPOTENUSE + (1 | SUBSITE),
  data = sjdf_clean, family = bernoulli(link = "logit"),
  chains = 4, cores = 4, prior = priors, iter = 4000,
  file = "models/m_flee_v2"
)

summary(m_flee_v2)
plot(m_flee_v2)
pp_check(m_flee_v2, type = "bars")
check_collinearity(m_flee_v2)

hypotenuse_effect <- as_tibble(fixef(m_flee_v2), rownames = "term") |>
  filter(grepl("HYPOTENUSE", term)) |>
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5)

#final model: drop hypotenuse
set.seed(123)
m_flee_final <- brm(
  FLEE ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + (1 | SUBSITE),
  data = sjdf_clean, family = bernoulli(link = "logit"),
  chains = 4, cores = 4, prior = priors, iter = 4000,
  file = "models/m_flee_final"
)

summary(m_flee_final)
plot(m_flee_final)
pp_check(m_flee_final, type = "bars")
check_collinearity(m_flee_final)

flee_coefs <- as_tibble(fixef(m_flee_final), rownames = "term") |>
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5) |>
  bind_rows(hypotenuse_effect, playback_effect) |>
  filter(term != "Intercept") |>
  mutate(model = "flee_prob")

write_csv(flee_coefs, file = "flee_prob_coefs.csv")

fitted_draws <- sjdf_clean %>%
  select(SEASON, GROUP.SIZE, SUBSITE) %>%
  expand_grid(SPECIES = c("CASJ", "ISSJ"), TREATMENT = c("CONTROL", "HAWK")) %>%
  add_epred_draws(m_flee_final, re_formula = NA, ndraws = 1000) %>%
  group_by(SPECIES, TREATMENT, .draw) %>%
  summarise(.epred = mean(.epred), .groups = "drop")

#numerical summaries: treatment effect per species
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(SPECIES) %>%
  median_qi(difference, .width = 0.95)

fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(SPECIES) %>%
  summarise(p_greater = mean(difference > 0))

#difference of differences
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES   = as.character(SPECIES)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  select(SPECIES, .draw, difference) %>%
  pivot_wider(names_from = SPECIES, values_from = difference) %>%
  mutate(diff_of_diffs = CASJ - ISSJ) %>%
  median_qi(diff_of_diffs, .width = 0.95)

fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT),
         SPECIES   = as.character(SPECIES)) %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  select(SPECIES, .draw, difference) %>%
  pivot_wider(names_from = SPECIES, values_from = difference) %>%
  mutate(diff_of_diffs = CASJ - ISSJ) %>%
  summarise(p_greater = mean(diff_of_diffs > 0))

#numerical summaries: species difference per treatment
fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT), SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  median_qi(difference, .width = 0.95)

fitted_draws %>%
  mutate(TREATMENT = as.character(TREATMENT), SPECIES = as.character(SPECIES)) %>%
  select(SPECIES, TREATMENT, .draw, .epred) %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(TREATMENT) %>%
  summarise(p_greater = mean(difference > 0))

flee_draws <- fitted_draws %>% mutate(model = "Flee")

####----FINAL PLOTTING----

#combine all draws and set factor levels
all_draws <- bind_rows(alarm_draws, mob_draws, interest_draws, flee_draws) %>%
  mutate(
    model     = fct_relevel(model, "Alarm", "Mob", "Flee", "Interest"),
    SPECIES   = factor(SPECIES,   levels = c("CASJ", "ISSJ")),
    TREATMENT = factor(TREATMENT, levels = c("CONTROL", "HAWK"))
  )

#raw observed proportions for background bars
raw_proportions <- bind_rows(
  sjdf_clean %>% group_by(SPECIES, TREATMENT) %>%
    summarise(prop = mean(ALARM),    n = n(), .groups = "drop") %>% mutate(model = "Alarm"),
  sjdf_clean %>% group_by(SPECIES, TREATMENT) %>%
    summarise(prop = mean(MOB),      n = n(), .groups = "drop") %>% mutate(model = "Mob"),
  sjdf_clean %>% group_by(SPECIES, TREATMENT) %>%
    summarise(prop = mean(FLEE),     n = n(), .groups = "drop") %>% mutate(model = "Flee"),
  sjdf_clean %>% group_by(SPECIES, TREATMENT) %>%
    summarise(prop = mean(INTEREST), n = n(), .groups = "drop") %>% mutate(model = "Interest")
) %>%
  mutate(
    model     = fct_relevel(model, "Alarm", "Mob", "Flee", "Interest"),
    SPECIES   = factor(SPECIES,   levels = c("CASJ", "ISSJ")),
    TREATMENT = factor(TREATMENT, levels = c("CONTROL", "HAWK")),
    group     = interaction(SPECIES, TREATMENT, sep = ".")
  )

#shared plot elements
dodge      <- position_dodge(width = 0.8)
diff_scale <- scale_y_continuous(limits = c(-0.75, 0.75), labels = scales::percent_format())
sig_color  <- scale_color_manual(values = c("TRUE" = "#E15759", "FALSE" = "gray60"))


common_theme <- theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_blank(),
    panel.spacing    = unit(2, "cm"),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position  = "none"
  )

#panel a: predicted probabilities with raw data bars and sample sizes
p1 <- all_draws %>%
  mutate(group = interaction(SPECIES, TREATMENT)) %>%
  group_by(model, SPECIES, TREATMENT, group) %>%
  median_qi(.epred, .width = 0.95) %>%
  ggplot(aes(x = SPECIES, y = .epred, group = group, color = group)) +
  geom_col(
    data = raw_proportions,
    aes(x = SPECIES, y = prop, fill = group, group = group),
    position = dodge, width = 0.7, alpha = 0.25, inherit.aes = FALSE
  ) +
  geom_point(position = dodge, size = 5) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper),
                position = dodge, width = 0, linewidth = 1.2) +
  geom_text(
    data = raw_proportions %>%
      mutate(group = interaction(SPECIES, TREATMENT)) %>%
      filter(model == "Alarm"),
    aes(x = SPECIES, y = Inf, label = n, group = group),
    position = dodge, vjust = -0.5, size = 6, inherit.aes = FALSE
  ) +
  scale_color_manual(values = c(
    "CASJ.CONTROL" = "#59A14F", "CASJ.HAWK" = "#4E79A7",
    "ISSJ.CONTROL" = "#59A14F", "ISSJ.HAWK" = "#4E79A7"
  )) +
  scale_fill_manual(values = c(
    "CASJ.CONTROL" = "#59A14F", "CASJ.HAWK" = "#4E79A7",
    "ISSJ.CONTROL" = "#59A14F", "ISSJ.HAWK" = "#4E79A7"
  )) +
  scale_x_discrete(
    labels = c(
      "CASJ" = expression("A. californica"),
      "ISSJ" = expression("A. insularis")
    ),
    expand = expansion(add = 0.1)  # 
  ) +
  scale_y_continuous(
    breaks = c(0, 0.50, 1.00),
    expand = expansion(mult = c(0, 0.04))
  ) +
  facet_wrap(~ model, ncol = 1) +
  labs(x = NULL, y = NULL) +  # y label removed, x set manually via scale
  coord_cartesian(clip = "off") +
  common_theme +
  theme(
    axis.text.x  = element_text(size = 20, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title = element_blank(),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin = margin(20, 2, 5, 2)  # extra top margin
  )

#panel b: treatment effect (hawk - control) per species
diff_scale <- scale_y_continuous(
  limits = c(-0.8, 0.8),
  breaks = c(-0.8, -0.4, 0, 0.4, 0.8)
)
p2 <- all_draws %>%
  pivot_wider(names_from = TREATMENT, values_from = .epred) %>%
  mutate(difference = HAWK - CONTROL) %>%
  group_by(model, SPECIES) %>%
  median_qi(difference, .width = 0.95) %>%
  mutate(significant = .lower > 0 | .upper < 0) %>%
  ggplot(aes(x = SPECIES, y = difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = significant), size = 5) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper, color = significant), width = 0, linewidth = 1.2) +
  sig_color + diff_scale +
  scale_x_discrete(labels = c(
    "CASJ" = expression("A. californica"),
    "ISSJ" = expression("A. insularis")
  )) +
  facet_wrap(~ model, ncol = 1) +
  labs(x = NULL, y = NULL) +
  common_theme +
  coord_cartesian(clip = "off") +
  theme(
    axis.text.x  = element_text(size = 20, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title  = element_blank(),
    plot.margin = margin(20, 2, 5, 2)
  )
#panel c: species difference (CASJ - ISSJ) per treatment
p3 <- all_draws %>%
  pivot_wider(names_from = SPECIES, values_from = .epred) %>%
  mutate(difference = CASJ - ISSJ) %>%
  group_by(model, TREATMENT) %>%
  median_qi(difference, .width = 0.95) %>%
  mutate(significant = .lower > 0 | .upper < 0) %>%
  ggplot(aes(x = TREATMENT, y = difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = significant), size = 5) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper, color = significant), width = 0, linewidth = 1.2) +
  sig_color + diff_scale +
  scale_x_discrete(labels = c("CONTROL" = "A. californica", "HAWK" = "Hawk")) +
  facet_wrap(~ model, ncol = 1) +
  labs(x = NULL, y = NULL) +
  common_theme +
  coord_cartesian(clip = "off") +
  theme(
    axis.text.x  = element_text(size = 20, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title   = element_blank(),
    plot.margin  = margin(20, 2, 5, 2)
  )


#save each panel separately for assembly 
ggsave("fig4_a.png", p1, width = 2.2, height = 10, dpi = 300)
ggsave("fig4_b.png", p2, width = 2.2, height = 10, dpi = 300)
ggsave("fig4_c.png", p3, width = 2.2, height = 10, dpi = 300)



# this script runs Bayesian negative binomial models on alarm and mob call count data

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

####----NUMBER OF ALARM CALLS----

#v1: includes playback and hypotenuse
set.seed(123)
m_nalarm_v1 <- brm(
  NUMBER.ALARM ~ SPECIES * TREATMENT + GROUP.SIZE + PLAYBACK + HYPOTENUSE + (1 | SUBSITE),
  data = sjdf_clean, family = negbinomial(),
  chains = 4, cores = 4, prior = priors, iter = 4000,
  file = "models/m_nalarm_v1"
)

summary(m_nalarm_v1)
plot(m_nalarm_v1)
pp_check(m_nalarm_v1)
check_collinearity(m_nalarm_v1)

#save playback effect before dropping it
playback_effect <- as_tibble(fixef(m_nalarm_v1), rownames = "term") |>
  filter(grepl("PLAYBACK", term)) |>
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5)

#v2: drop playback, add season, keep hypotenuse
set.seed(123)
m_nalarm_v2 <- brm(
  NUMBER.ALARM ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + HYPOTENUSE + (1 | SUBSITE),
  data = sjdf_clean, family = negbinomial(),
  chains = 4, cores = 4, prior = priors, iter = 4000,
  file = "models/m_nalarm_v2"
)

summary(m_nalarm_v2)
plot(m_nalarm_v2)
pp_check(m_nalarm_v2)
check_collinearity(m_nalarm_v2)

#save hypotenuse effect before dropping it
hypotenuse_effect <- as_tibble(fixef(m_nalarm_v2), rownames = "term") |>
  filter(grepl("HYPOTENUSE", term)) |>
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5)

#final model: drop hypotenuse
set.seed(123)
m_nalarm_final <- brm(
  NUMBER.ALARM ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + (1 | SUBSITE),
  data = sjdf_clean, family = negbinomial(),
  chains = 4, cores = 4, prior = priors, iter = 4000,
  file = "models/m_nalarm_final"
)

summary(m_nalarm_final)
plot(m_nalarm_final)
pp_check(m_nalarm_final)
check_collinearity(m_nalarm_final)

#save all coefficients and write to csv
nalarm_coefs <- as_tibble(fixef(m_nalarm_final), rownames = "term") |>
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5) |>
  bind_rows(hypotenuse_effect, playback_effect) |>
  filter(term != "Intercept") |>
  mutate(model = "number_alarm")

write_csv(nalarm_coefs, file = "number_alarm_coefs.csv")

#generate posterior predicted draws for plotting and summaries
fitted_draws <- sjdf_clean %>%
  select(SEASON, GROUP.SIZE, SUBSITE) %>%
  expand_grid(SPECIES = c("CASJ", "ISSJ"), TREATMENT = c("CONTROL", "HAWK")) %>%
  add_epred_draws(m_nalarm_final, re_formula = NA, ndraws = 1000) %>%
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

nalarm_draws <- fitted_draws %>% mutate(model = "Alarm calls")

####----NUMBER OF MOB CALLS----

#v1: includes playback and hypotenuse
set.seed(123)
m_nmob_v1 <- brm(
  NUMBER.MOB ~ SPECIES * TREATMENT + GROUP.SIZE + PLAYBACK + HYPOTENUSE + (1 | SUBSITE),
  data = sjdf_clean, family = negbinomial(),
  chains = 4, cores = 4, prior = priors, iter = 4000,
  file = "models/m_nmob_v1"
)

summary(m_nmob_v1)
plot(m_nmob_v1)
pp_check(m_nmob_v1)
check_collinearity(m_nmob_v1)

#save playback effect before dropping it
playback_effect <- as_tibble(fixef(m_nmob_v1), rownames = "term") |>
  filter(grepl("PLAYBACK", term)) |>
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5)

#v2: drop playback, add season, keep hypotenuse
set.seed(123)
m_nmob_v2 <- brm(
  NUMBER.MOB ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + HYPOTENUSE + (1 | SUBSITE),
  data = sjdf_clean, family = negbinomial(),
  chains = 4, cores = 4, prior = priors, iter = 4000,
  file = "models/m_nmob_v2"
)

summary(m_nmob_v2)
plot(m_nmob_v2)
pp_check(m_nmob_v2)
check_collinearity(m_nmob_v2)

#save hypotenuse effect before dropping it
hypotenuse_effect <- as_tibble(fixef(m_nmob_v2), rownames = "term") |>
  filter(grepl("HYPOTENUSE", term)) |>
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5)

#final model: drop hypotenuse
set.seed(123)
m_nmob_final <- brm(
  NUMBER.MOB ~ SPECIES * TREATMENT + GROUP.SIZE + SEASON + (1 | SUBSITE),
  data = sjdf_clean, family = negbinomial(),
  chains = 4, cores = 4, prior = priors, iter = 4000,
  file = "models/m_nmob_final"
)

summary(m_nmob_final)
plot(m_nmob_final)
pp_check(m_nmob_final)
check_collinearity(m_nmob_final)

#save all coefficients and write to csv
nmob_coefs <- as_tibble(fixef(m_nmob_final), rownames = "term") |>
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5) |>
  bind_rows(hypotenuse_effect, playback_effect) |>
  filter(term != "Intercept") |>
  mutate(model = "number_mob")

write_csv(nmob_coefs, file = "number_mob_coefs.csv")

#generate posterior predicted draws for plotting and summaries
fitted_draws <- sjdf_clean %>%
  select(SEASON, GROUP.SIZE, SUBSITE) %>%
  expand_grid(SPECIES = c("CASJ", "ISSJ"), TREATMENT = c("CONTROL", "HAWK")) %>%
  add_epred_draws(m_nmob_final, re_formula = NA, ndraws = 1000) %>%
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

nmob_draws <- fitted_draws %>% mutate(model = "Mob calls")

####----FINAL PLOTTING----

#combine draws and set factor levels
all_draws <- bind_rows(nalarm_draws, nmob_draws) %>%
  mutate(
    model     = fct_relevel(model, "Alarm calls", "Mob calls"),
    SPECIES   = factor(SPECIES,   levels = c("CASJ", "ISSJ")),
    TREATMENT = factor(TREATMENT, levels = c("CONTROL", "HAWK"))
  )

#raw count data for violin and jitter layers
raw_points <- bind_rows(
  sjdf_clean %>% select(SPECIES, TREATMENT, value = NUMBER.ALARM) %>% mutate(model = "Alarm calls"),
  sjdf_clean %>% select(SPECIES, TREATMENT, value = NUMBER.MOB)   %>% mutate(model = "Mob calls")
) %>%
  mutate(
    model     = fct_relevel(model, "Alarm calls", "Mob calls"),
    SPECIES   = factor(SPECIES,   levels = c("CASJ", "ISSJ")),
    TREATMENT = factor(TREATMENT, levels = c("CONTROL", "HAWK")),
    group     = as.character(interaction(SPECIES, TREATMENT, sep = "."))
  )

#sample counts for panel a labels — only shown on first facet
raw_counts <- raw_points %>%
  group_by(model, SPECIES, TREATMENT) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(group = as.character(interaction(SPECIES, TREATMENT, sep = ".")))

#pre-compute model summary for panel a
model_summary <- all_draws %>%
  mutate(group = as.character(interaction(SPECIES, TREATMENT, sep = "."))) %>%
  group_by(model, SPECIES, TREATMENT, group) %>%
  median_qi(.epred, .width = 0.95) %>%
  ungroup()

#shared plot elements
dodge     <- position_dodge(width = 0.8)
sig_color <- scale_color_manual(values = c("TRUE" = "#E15759", "FALSE" = "gray60"))

group_colors <- c(
  "CASJ.CONTROL" = "#59A14F", "CASJ.HAWK" = "#4E79A7",
  "ISSJ.CONTROL" = "#59A14F", "ISSJ.HAWK" = "#4E79A7"
)

common_theme <- theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_blank(),
    panel.spacing    = unit(0.9, "cm"),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position  = "none"
  )

#shared difference scale for panels b and c
count_diff_scale <- scale_y_continuous(
  limits = c(-5, 35),
  labels = scales::label_number()
)

#panel a: raw data violin + jitter with model estimates overlaid
p1 <- ggplot() +
  geom_violin(
    data = raw_points,
    aes(x = SPECIES, y = value, group = group, fill = group),
    position = position_dodge(width = 0.8),
    color = NA, alpha = 0.25, scale = "width"
  ) +
  geom_jitter(
    data = raw_points,
    aes(x = SPECIES, y = value, group = group, color = group),
    position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8),
    size = 1.0, alpha = 0.15
  ) +
  geom_point(
    data = model_summary,
    aes(x = SPECIES, y = .epred, group = group, color = group),
    position = position_dodge(width = 0.8), size = 5
  ) +
  geom_errorbar(
    data = model_summary,
    aes(x = SPECIES, ymin = .lower, ymax = .upper, group = group, color = group),
    position = position_dodge(width = 0.8), width = 0, linewidth = 1.2
  ) +
  geom_text(
    data = raw_counts %>% filter(model == "Alarm calls"),
    aes(x = SPECIES, y = Inf, label = n, group = group),
    position = position_dodge(width = 0.8), vjust = -0.5, size = 6,
    color = "black"
  ) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values  = group_colors) +
  scale_x_discrete(
    labels = c(
      "CASJ" = expression(italic("A. c.")),
      "ISSJ" = expression(italic("A. i."))
    ),
    expand = expansion(add = 0.1)
  ) +
  scale_y_continuous(
    trans   = "pseudo_log",
    labels  = scales::label_number(),
    expand  = expansion(mult = c(0.02, 0.04))
  ) +
  facet_wrap(~ model, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  common_theme +
  theme(
    axis.text.x       = element_text(size = 20),
    axis.text.y       = element_text(size = 18),
    axis.title        = element_blank(),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin       = margin(20, 2, 5, 2)
  )

#panel b: treatment effect (hawk - control) per species
#red = CI excludes zero
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
  sig_color + count_diff_scale +
  scale_x_discrete(labels = c(
    "CASJ" = expression(italic("A. c.")),
    "ISSJ" = expression(italic("A. i."))
  )) +
  facet_wrap(~ model, ncol = 1) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  common_theme +
  theme(
    axis.text.x       = element_text(size = 20),
    axis.text.y       = element_text(size = 18),
    axis.title        = element_blank(),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin       = margin(20, 2, 5, 2)
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
  sig_color + count_diff_scale +
  scale_x_discrete(labels = c("CONTROL" = "C", "HAWK" = "H")) +
  facet_wrap(~ model, ncol = 1) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  common_theme +
  theme(
    axis.text.x       = element_text(size = 20),
    axis.text.y       = element_text(size = 18),
    axis.title        = element_blank(),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin       = margin(20, 2, 5, 2)
  )

#save each panel separately
ggsave("fig5_a.png", p1, width = 2.2, height = 7, dpi = 300)
ggsave("fig5_b.png", p2, width = 2.2, height = 7, dpi = 300)
ggsave("fig5_c.png", p3, width = 2.2, height = 7, dpi = 300)


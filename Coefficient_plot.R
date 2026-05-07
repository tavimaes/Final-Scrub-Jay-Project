####This script is for making a coefficient plot for all models

#load packages
library(tidyverse)
library(ggplot2)
library(dplyr)

#read in data for all models

alarm_prob_data <- read.csv("alarm_prob_coefs.csv")
mob_prob_data <- read.csv("mob_prob_coefs.csv")
interest_prob_data <- read.csv("interest_prob_coefs.csv")
flee_prob_data <- read.csv("flee_prob_coefs.csv")
alarm_call_data <- read.csv("number_alarm_coefs.csv")
mob_call_data <- read.csv("number_mob_coefs.csv")

#bind all dataframes

all_data <- bind_rows(alarm_prob_data, mob_prob_data, interest_prob_data, 
                      flee_prob_data, alarm_call_data, mob_call_data)


#clean data titles



#plot

panel_labels <- data.frame(
  model = factor(
    c("alarm_prob", "mob_prob", "flee_prob",
      "interest_prob", "number_alarm", "number_mob"),
    levels = c("alarm_prob", "mob_prob", "flee_prob",
               "interest_prob", "number_alarm", "number_mob")
  ),
  label = c("a)", "b)", "c)", "d)", "e)", "f)")
)

all_data <- all_data %>%
  mutate(term = as.factor(term))

n_terms <- nlevels(all_data$term)

band_data <- data.frame(
  ymin = seq(0.5, n_terms - 0.5, by = 1),
  ymax = seq(1.5, n_terms + 0.5, by = 1),
  fill = rep(c("gray95", "white"), length.out = n_terms)
)


all_data %>%
  mutate(
    significant = lower > 0 | upper < 0,
    model = factor(model, levels = c(
      "alarm_prob", "mob_prob", "flee_prob",
      "interest_prob", "number_alarm", "number_mob"
    )),
    term = factor(term,
                  levels = c(
                    "PLAYBACKYES",
                    "HYPOTENUSE",
                    "SEASONWINTER",
                    "GROUP.SIZE",
                    "SPECIESISSJ:TREATMENTHAWK",
                    "TREATMENTHAWK",
                    "SPECIESISSJ"
                  ),
                  labels = c(
                    "Playback*",
                    "Hypotenuse*",
                    "Season (Winter)",
                    "Group Size",
                    "Species × Treatment",
                    "Treatment (Hawk)",
                    "Species (A. insularis)"
                  )
    )
  ) %>%
  ggplot(aes(x = estimate, y = term, color = significant)) +
  geom_rect(
    data = band_data,
    aes(ymin = ymin, ymax = ymax, fill = fill),
    xmin = -Inf, xmax = Inf,
    inherit.aes = FALSE
  ) +
  scale_fill_identity(guide = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3) +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red")) +
  facet_wrap(~ model, scales = "fixed") +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  scale_color_manual(
    values = c("FALSE" = "gray60", "TRUE" = "#E15759"),
    labels = c("FALSE" = "Not credible", "TRUE" = "Credible"),
    name = NULL
  ) +
  geom_text(
    data = panel_labels,
    aes(label = label, x = -Inf, y = Inf),
    hjust = -0.3,
    vjust = 1.5,
    color = "black",
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  scale_y_discrete(expand = expansion(add = c(0.5, 1))) +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = "Effect size (log-odds or log-count)",
    y = "Predictor"
  )


ggsave("fig3.png", last_plot(), width = 7, height = 5, dpi = 300)









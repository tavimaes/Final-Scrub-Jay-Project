####----Coefficient plot for all models----

library(tidyverse)
library(ggplot2)

#read in and combine all model coefficients
all_data <- bind_rows(
  read.csv("alarm_prob_coefs.csv"),
  read.csv("mob_prob_coefs.csv"),
  read.csv("interest_prob_coefs.csv"),
  read.csv("flee_prob_coefs.csv"),
  read.csv("number_alarm_coefs.csv"),
  read.csv("number_mob_coefs.csv")
) %>%
  mutate(
    significant = lower > 0 | upper < 0,
    model = factor(model, levels = c(
      "alarm_prob", "mob_prob", "flee_prob",
      "interest_prob", "number_alarm", "number_mob"
    )),
    term = factor(as.character(term),
                  levels = c(
                    "PLAYBACKYES", "HYPOTENUSE", "SEASONWINTER",
                    "GROUP.SIZE", "SPECIESISSJ:TREATMENTHAWK",
                    "TREATMENTHAWK", "SPECIESISSJ"
                  ),
                  labels = c(
                    "Playback*", "Hypotenuse*", "Season (Winter)",
                    "Group Size", "Species × Treatment",
                    "Treatment (Hawk)", "Species (A. insularis)"
                  )
    )
  )

#panel labels
panel_labels <- data.frame(
  model = factor(
    c("alarm_prob", "mob_prob", "flee_prob",
      "interest_prob", "number_alarm", "number_mob"),
    levels = c("alarm_prob", "mob_prob", "flee_prob",
               "interest_prob", "number_alarm", "number_mob")
  ),
  label = c("a)", "b)", "c)", "d)", "e)", "f)")
)

#plot
ggplot(all_data, aes(x = estimate, y = term, color = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), width = 0, linewidth = 1) +
  geom_text(
    data = panel_labels,
    aes(label = label, x = -Inf, y = Inf, size = 5),
    hjust = -0.3, vjust = 1.5,
    color = "black", fontface = "bold",
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "#E15759")) +
  scale_x_continuous(
    breaks = c(-3, 0, 3)) +
  facet_wrap(~ model, scales = "fixed") +
  theme_classic() +
  scale_y_discrete(expand = expansion(add = c(0.5, 1.5))) +
  theme(
    legend.position  = "none",
    strip.text       = element_blank(),
    strip.background = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text.x      = element_text(size = 14),
    axis.text.y      = element_text(size = 14),
    axis.title       = element_blank()
  )

ggsave("fig3.png", last_plot(), width = 8, height = 5, dpi = 300)

####----Coefficient plot — flipped layout, 6 stacked panels----

library(tidyverse)
library(ggplot2)

#read in and combine all model coefficients
all_data <- bind_rows(
  read.csv("alarm_prob_coefs.csv"),
  read.csv("mob_prob_coefs.csv"),
  read.csv("interest_prob_coefs.csv"),
  read.csv("flee_prob_coefs.csv"),
  read.csv("number_alarm_coefs.csv"),
  read.csv("number_mob_coefs.csv")
) %>%
  mutate(
    significant = lower > 0 | upper < 0,
    model = factor(model, levels = c(
      "alarm_prob", "mob_prob", "flee_prob",
      "interest_prob", "number_alarm", "number_mob"
    )),
    # reversed order so Species (A. insularis) appears at top
    term = factor(as.character(term),
                  levels = c(
                    "SPECIESISSJ", "TREATMENTHAWK", "SPECIESISSJ:TREATMENTHAWK",
                    "GROUP.SIZE", "SEASONWINTER", "HYPOTENUSE", "PLAYBACKYES"
                  ),
                  labels = c(
                    "Species (A. insularis)", "Treatment (Hawk)", "Species × Treatment",
                    "Group Size", "Season (Winter)", "Hypotenuse*", "Playback*"
                  )
    )
  )


ggplot(all_data, aes(x = term, y = estimate, color = significant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, linewidth = 1) +
  

  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "#E15759")) +
  scale_y_continuous(breaks = c(-3, 0, 3), limits = c(-5, 5)) +
  scale_x_discrete(expand = expansion(add = c(0.5, 1.5))) +
  facet_wrap(~ model, scales = "fixed", ncol = 1) +
  theme_classic() +
  theme(
    legend.position  = "none",
    strip.text       = element_blank(),
    strip.background = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.spacing    = unit(0.5, "cm"),
    axis.text.x      = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y      = element_text(size = 14),
    axis.title       = element_blank()
  )

ggsave("fig3_flipped.png", last_plot(), width = 5, height = 10, dpi = 300)


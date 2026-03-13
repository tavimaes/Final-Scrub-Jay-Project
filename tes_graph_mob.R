pred_df_mob <- emmeans(m.mob, ~ SPECIES * TREATMENT) |>
  summary(type = "response") |>
  as.data.frame()

sjdf_clean <- read_csv("clean_data.csv")

raw_df <- sjdf_clean |>
  dplyr::select(SPECIES, TREATMENT, MOB)



ggplot() +
  stat_slab(data = sjdf_clean, show.legend = FALSE,
            aes(x = SPECIES,
                y = MOB,
                side = "both",
                fill = TREATMENT,
                group = TREATMENT),
            position = position_dodge(width = 0.5),
            alpha = 0.2,
            scale = 1) +
  
  geom_point(data = pred_df_mob,
             aes(x = SPECIES,
                 y = prob,
                 color = TREATMENT,
                 shape = TREATMENT),
             size = 5,
             position = position_dodge(width = 0.5)) +
  
  geom_errorbar(data = pred_df_mob, show.legend = FALSE,
                aes(x = SPECIES,
                    ymin = asymp.LCL,
                    ymax = asymp.UCL,
                    color = TREATMENT),
                width = 0.15,
                position = position_dodge(width = 0.5)) + 
  geom_line(data = pred_df_mob, show.legend = FALSE,
            aes(x = SPECIES,
                y = prob,
                color = TREATMENT,
                group = TREATMENT),
            position = position_dodge(width = 0.5),
            linewidth = 0.4) +
  scale_color_manual(values = c("CONTROL" = "gray60",
                                "HAWK" = "darkblue")) +
  scale_fill_manual(values = c("CONTROL" = "gray60",
                               "HAWK" = "darkblue")) +
  scale_shape_manual(values = c("CONTROL" = 15,
                                "HAWK" = 17)) +
  scale_x_discrete(labels = c(
    CASJ = expression(italic("A. californica")),
    ISSJ = expression(italic("A. insularis"))
  )) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  
  theme_classic(base_size = 14) +
  labs(x = "Species",
       y = "Probability of mobbing")


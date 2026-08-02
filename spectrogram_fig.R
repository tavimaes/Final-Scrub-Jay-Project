########################################################################
## Title: Spectrogram Visualization
## Input: One WAV file per species
## Output: Spectrogram figure (PDF)
## Date: 04/10/2026
########################################################################

# Load packages
library(warbleR)
library(tuneR)
library(here)
library(scico)
library(ggplot2)

# Define the four WAV files
wav_files <- c(
  here("WAV_files/casj.alarm.clipped.wav"),
  here("WAV_files/casj.mob.clipped.wav"),
  here("WAV_files/issj.alarm.clipped.wav"),
  here("WAV_files/issj.mob.clipped.wav")
)

# Function to generate spectrogram data frame
get_spectrogram_df <- function(path, wl = 512, ovlp = 90) {
  
  # Read audio file
  wav <- readWave(path)
  
  # Normalize amplitude
  wav <- normalize(wav, unit = "16")
  
  # Compute spectrogram
  sp <- spectro(
    wav,
    f = wav@samp.rate,
    wl = wl,
    ovlp = ovlp,
    wn = "hanning",
    plot = FALSE,
    scale = FALSE,
    osc = FALSE
  )
  
  # Convert spectrogram to data frame
  df <- expand.grid(
    time = sp$time,
    freq = sp$freq
  )
  
  df$amplitude <- as.vector(t(sp$amp))
  
  # Convert amplitude to decibels
  min_amp <- min(df$amplitude, na.rm = TRUE)
  df$amplitude_db <- 10 * log10(df$amplitude - min_amp + 1e-12)
  
  # Store file name for faceting
  df$file <- basename(path)
  
  return(df)
}

# Generate spectrogram data for all files
spectrogram_list <- lapply(wav_files, get_spectrogram_df)

# Combine into a single data frame
spectrogram_df <- do.call(rbind, spectrogram_list)

# Define color scale limits
color_limits <- quantile(
  spectrogram_df$amplitude_db,
  probs = c(0.05, 0.99),
  na.rm = TRUE
)

# Create spectrogram figure
spectrogram_plot <-
  ggplot(spectrogram_df,
         aes(x = time, y = freq, fill = amplitude_db)) +
  geom_raster() +
  facet_wrap(~file) +
  scale_fill_scico(
    palette = "lapaz",
    limits = color_limits,
    oob = scales::squish
  ) +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = function(x) {
      labs <- sprintf("%.2f", x)
      labs[x == 0] <- "0.0"
      labs[x == 1] <- "1.0"
      labs
    },
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 13),
    expand = c(0, 0),
    breaks = seq(0, 12, by = 3)
  ) +
  labs(
    x = "time (s)",
    y = "frequency (kHz)",
    fill = "dB"
  ) +
  theme_classic(base_size = 13) +
  theme_classic(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.spacing = unit(10, "mm"),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1.5
    )
  )

# Save figure
ggsave(
  "fig2.png",
  plot = spectrogram_plot,
  width = 200,
  height = 80,
  units = "mm",
  dpi = 300
)

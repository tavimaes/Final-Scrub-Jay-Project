# ============================================================
# PUBLICATION-STYLE TOPOGRAPHIC MAPS
# California + Monterey Bay + Santa Cruz Island
#
# Uses:
# - grayscale hillshade
# - contour lines
# - clean axes
# - high-res local topography appearance
# ============================================================

# INSTALL IF NEEDED:
# install.packages(c(
#   "sf",
#   "ggplot2",
#   "terra",
#   "geodata",
#   "rnaturalearth",
#   "rnaturalearthdata"
# ))

library(sf)
library(ggplot2)
library(terra)
library(geodata)
library(rnaturalearth)
library(rnaturalearthdata)

# ============================================================
# 1. CONVERT POINTS TO SF
# ============================================================

sf_mainland <- st_as_sf(
  mainland_all,
  coords = c("lon", "lat"),
  crs = 4326
)

sf_island <- st_as_sf(
  island_all,
  coords = c("lon", "lat"),
  crs = 4326
)

# ============================================================
# 2. CALIFORNIA OUTLINE
# ============================================================

states <- ne_states(
  country = "United States of America",
  returnclass = "sf"
)

california <- states[states$name == "California", ]

# ============================================================
# 3. DOWNLOAD ELEVATION DATA
# ============================================================

elev <- geodata::elevation_30s(
  country = "USA",
  path = tempdir()
)

# ============================================================
# 4. CROP TO STUDY REGIONS
# ============================================================

ca <- crop(
  elev,
  ext(-125, -114, 32, 42.5)
)

mb <- crop(
  elev,
  ext(-122.35, -121.65, 36.60, 37.20)
)

sci <- crop(
  elev,
  ext(-119.90, -119.55, 33.90, 34.10)
)

# ============================================================
# 5. CONVERT TO DATA FRAMES
# ============================================================

ca_df <- as.data.frame(ca, xy = TRUE, na.rm = TRUE)
mb_df <- as.data.frame(mb, xy = TRUE, na.rm = TRUE)
sci_df <- as.data.frame(sci, xy = TRUE, na.rm = TRUE)

# elevation column name
elev_col <- names(ca_df)[3]

# ============================================================
# 6. CONTOUR DATA
# ============================================================

ca_df  <- as.data.frame(ca,  xy = TRUE, na.rm = TRUE)
mb_df  <- as.data.frame(mb,  xy = TRUE, na.rm = TRUE)
sci_df <- as.data.frame(sci, xy = TRUE, na.rm = TRUE)

# elevation column name
elev_col <- names(ca_df)[3]

# ============================================================
# 7. COMMON THEME
# ============================================================

map_theme <- theme_bw() +
  
  theme(
    
    panel.grid.major = element_line(
      color = "grey85",
      linewidth = 0.25
    ),
    
    panel.grid.minor = element_blank(),
    
    axis.title = element_text(size = 10),
    
    axis.text = element_text(size = 8),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 12
    ),
    
    legend.position = "none"
  )

p_california <- ggplot() +
  
  geom_raster(
    data = ca_df,
    aes(
      x = x,
      y = y,
      fill = .data[[elev_col]]
    ),
    interpolate = TRUE
  ) +
  
  scale_fill_gradient(
    low = "white",
    high = "grey40"
  ) +
  
  geom_contour(
    data = ca_df,
    aes(
      x = x,
      y = y,
      z = .data[[elev_col]]
    ),
    color = "grey25",
    bins = 20,
    linewidth = 0.2,
    alpha = 0.5
  ) +
  
  geom_sf(
    data = california,
    fill = NA,
    color = "black",
    linewidth = 0.5
  ) +
  
  coord_sf(
    xlim = c(-125,-114),
    ylim = c(32,42.5),
    expand = FALSE
  ) +
  
  labs(
    title = "California",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  map_theme

p_monterey <- ggplot() +
  
  geom_raster(
    data = mb_df,
    aes(
      x = x,
      y = y,
      fill = .data[[elev_col]]
    ),
    interpolate = TRUE
  ) +
  
  scale_fill_gradient(
    low = "white",
    high = "grey40"
  ) +
  
  geom_contour(
    data = mb_df,
    aes(
      x = x,
      y = y,
      z = .data[[elev_col]]
    ),
    color = "grey25",
    bins = 35,
    linewidth = 0.2,
    alpha = 0.6
  ) +
  
  geom_sf(
    data = sf_mainland,
    color = "red3",
    size = 2
  ) +
  
  coord_sf(
    xlim = c(-122.35,-121.65),
    ylim = c(36.60,37.20),
    expand = FALSE
  ) +
  
  labs(
    title = "Monterey Bay",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  map_theme

p_sci <- ggplot() +
  
  geom_raster(
    data = sci_df,
    aes(
      x = x,
      y = y,
      fill = .data[[elev_col]]
    ),
    interpolate = TRUE
  ) +
  
  scale_fill_gradient(
    low = "white",
    high = "grey40"
  ) +
  
  geom_contour(
    data = sci_df,
    aes(
      x = x,
      y = y,
      z = .data[[elev_col]]
    ),
    color = "grey25",
    bins = 35,
    linewidth = 0.2,
    alpha = 0.6
  ) +
  
  geom_sf(
    data = sf_island,
    color = "red3",
    size = 2
  ) +
  
  coord_sf(
    xlim = c(-119.90,-119.55),
    ylim = c(33.90,34.10),
    expand = FALSE
  ) +
  
  labs(
    title = "Santa Cruz Island",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  map_theme

# ============================================================
# 11. VIEW MAPS
# ============================================================

print(p_california)
print(p_monterey)
print(p_sci)

# ============================================================
# 12. EXPORT
# ============================================================

ggsave(
  "california_map.png",
  p_california,
  width = 6,
  height = 6,
  dpi = 600
)

ggsave(
  "monterey_map.png",
  p_monterey,
  width = 6,
  height = 5,
  dpi = 600
)

ggsave(
  "santa_cruz_island_map.png",
  p_sci,
  width = 6,
  height = 5,
  dpi = 600
)
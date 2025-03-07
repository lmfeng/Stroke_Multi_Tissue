#### Figure 1 Radar Chart ####
# Data Preparation
library(dplyr)
library(scales)
library(tibble)
library(ggradar)

# Calculate tissue proportions
PBMC_ratio <- prop.table(table(gexpr$Cluster1, gexpr$condition), margin = 2) %>%
  as.data.frame() %>%
  mutate(percent = Freq*100) %>%
  rename(Cluster = Var1, Condition = Var2)

# Data loading (replace with actual file paths)
gexpr_percent <- read.table("gexpr_ratio.txt", sep="\t", header=TRUE)
gexpr_percent_all <- read.table("gexpr_ratio_minus.txt", sep="\t", header=TRUE)

#### Radar Plot Function ####
create_radar <- function(data, grid_params, colors, filename) {
  p <- ggradar(
    data,
    axis.labels = rep(NA, ncol(data)-1),
    grid.min = grid_params$min,
    grid.mid = grid_params$mid,
    grid.max = grid_params$max,
    values.radar = c(grid_params$min, grid_params$mid, grid_params$max),
    group.colours = colors,
    group.line.width = 1,
    group.point.size = 2,
    background.circle.colour = "white",
    gridline.mid.colour = grid_params$mid_color,
    legend.position = "none",
    plot.title = filename
  ) + theme(
    plot.background = element_blank(),
    panel.background = element_blank()
  )
  
  ggsave(paste0(filename, ".pdf"), plot = p, width = 5, height = 5)
  return(p)
}

#### Generate Radar Plots ####
# Brain Plot
brain_data <- gexpr_percent[c(2,3), ]
create_radar(
  brain_data,
  grid_params = list(min = 0, mid = 0.5, max = 1, mid_color = "#f9766c"),
  colors = c("#f9766e", "#00bfc4"),
  filename = "Brain_MCAO_vs_Sham"
)

# CBM Plot
cbm_data <- gexpr_percent[c(4,5), ]
create_radar(
  cbm_data,
  grid_params = list(min = 0, mid = 0.25, max = 0.5, mid_color = "#7cae00"),
  colors = c("#f9766e", "#00bfc4"),
  filename = "CBM_MCAO_vs_Sham"
)

# FBM Plot
fbm_data <- gexpr_percent[c(6,7), ]
create_radar(
  fbm_data,
  grid_params = list(min = 0, mid = 0.25, max = 0.5, mid_color = "#00c0c5"),
  colors = c("#f9766e", "#00bfc4"),
  filename = "FBM_MCAO_vs_Sham"
)

# PBMC Plot
pbmc_data <- gexpr_percent[c(8,9), ]
create_radar(
  pbmc_data,
  grid_params = list(min = 0, mid = 0.25, max = 0.5, mid_color = "#c77bff"),
  colors = c("#f9766e", "#00bfc4"),
  filename = "PBMC_MCAO_vs_Sham"
)
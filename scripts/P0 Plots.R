################################################################################
# Script: P0Plots.R
# Author: [Your Name]
#
# Purpose:
#   Generate:
#     1) Violin + beeswarm plots for lifetime reproductive success (LRS)
#     2) Violin + beeswarm plots for rate-sensitive fitness (Leslie λ)
#     3) Age-specific reproduction curves by treatment
#
# Expected project structure (example):
#   project/
#    ├ data/
#    │   └ P0repro.csv
#    └ scripts/
#        └ P0Plots.R
################################################################################

#############################
# 1. Load required packages #
#############################

library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(dplyr)
library(rlang)
library(scales)
library(ggpubr)
library(ggh4x)
library(cowplot)
library(legendry)
library(here)

################################################################################
# 2. Generic plotting function: violin + quasirandom + mean ± CI               #
################################################################################

plot_violin <- function(
    data,
    x,                     # grouping variable (factor/character)
    y,                     # numeric outcome (e.g., LRS, fitness)
    palette = NULL,        # vector of colors (one per group); NULL uses defaults
    y_lab = "Lifetime \n Reproductive Success",
    x_lab = NULL,
    title = NULL,
    subtitle = NULL,
    point_alpha = 0.4,
    point_size = 4,
    violin_alpha = 0.85,
    width_violin = 0.9
) {
  x <- enquo(x); y <- enquo(y)
  
  # Prepare grouping levels and default palette
  df <- data %>%
    filter(!is.na(!!x), !is.na(!!y)) %>%
    mutate(.group = factor(!!x))
  
  n_groups <- nlevels(df$.group)
  if (is.null(palette)) {
    palette <- c(
      "#D55E00",
      "#0072B2",
      "#009E73",
      "#CC79A7"
    )[seq_len(n_groups)]
  }
  
  p <- ggplot(df, aes(x = .group, y = !!y, fill = .group)) +
    # Violin for distribution
    geom_violin(
      width = width_violin,
      alpha = violin_alpha,
      color = NA,
      trim  = FALSE,
      scale = "width"
    ) +
    # Individual points (beeswarm)
    geom_quasirandom(
      color       = "black",
      dodge.width = 0.7,
      varwidth    = FALSE,
      alpha       = point_alpha,
      size        = point_size
    ) +
    # Mean point
    stat_summary(
      fun   = mean,
      geom  = "point",
      shape = 20,
      size  = 3,
      color = "black"
    ) +
    # Mean ± 95% CI (normal approximation)
    stat_summary(
      fun.data = mean_cl_normal,
      geom     = "errorbar",
      width    = 0.3,
      color    = "black",
      linewidth = 1.6
    ) +
    scale_fill_manual(values = palette) +
    scale_color_manual(values = palette) +
    labs(
      x        = x_lab %||% rlang::as_name(x),
      y        = y_lab,
      title    = title,
      subtitle = subtitle
    ) +
    guides(fill = "none", color = "none") +
    theme_classic(base_size = 11) +
    theme(
      plot.title.position = "plot",
      plot.title          = element_text(face = "bold"),
      axis.title.x        = element_blank(),
      axis.title.y        = element_text(margin = margin(r = 8), size = 30),
      axis.text           = element_text(color = "black", size = 25),
      axis.ticks          = element_line(color = "black"),
      strip.text          = element_text(size = 15),
      strip.placement     = "outside"
    ) +
    coord_cartesian(clip = "off")
  
  return(p)
}

################################################################################
# 3. Load and clean P0 data                                                    #
################################################################################

P0 <- read.csv(here("data", "P0repro.csv"))

# Reshape from wide to long (one row per Plate–Day)
data <- pivot_longer(
  P0,
  cols      = starts_with("D"),
  names_to  = "Day",
  values_to = "Reproduction"
)

# Convert Day to numeric, stripping "D" prefix
data$Day <- as.numeric(gsub("D", "", data$Day))

# Remove plates that have any NA in Reproduction
data <- data %>%
  group_by(Plate) %>%
  filter(!any(is.na(Reproduction))) %>%
  ungroup()

# Ensure Treatment is a factor with desired order
data <- data %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("Control", "Larval Starvation")
    )
  )

################################################################################
# 4. Lifetime Reproductive Success (LRS) plot                                  #
################################################################################

# Compute total reproduction per plate
LRS <- aggregate(
  Reproduction ~ Plate,
  data = data,
  FUN  = sum
)

# Merge in Treatment per plate
LRS <- merge(
  LRS,
  data %>% select(Plate, Treatment) %>% distinct(),
  by = "Plate"
)

# Ensure Treatment factor order is correct
LRS <- LRS %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("Control", "Larval Starvation")
    )
  )

# Plot LRS
LRS_plot <- plot_violin(
  data = LRS,
  x    = Treatment,
  y    = Reproduction,
  y_lab = "Lifetime \n Reproductive Success"
)

LRS_plot

################################################################################
# 5. Rate-sensitive fitness (Leslie λ) plot                                    #
################################################################################

# Compute dominant eigenvalue (λ) per plate
dominant_eigenvalues <- numeric(length(unique(data$Plate)))

for (i in unique(data$Plate)) {
  # Leslie matrix dimension = number of unique days
  Leslie <- matrix(
    0,
    nrow = length(unique(data$Day)),
    ncol = length(unique(data$Day))
  )
  
  # First row: reproduction for this plate across days
  Leslie[1, ] <- data$Reproduction[data$Plate == i]
  
  # Subdiagonal: survival (assumed 1)
  diag(Leslie[-1, ]) <- 1
  
  # Dominant eigenvalue
  eigen_res           <- eigen(Leslie)
  dominant_eigenvalue <- max(Re(eigen_res$values))
  dominant_eigenvalues[i] <- dominant_eigenvalue
}

# Make a data frame
fitness <- as.data.frame(dominant_eigenvalues)

# Unique Plate–Treatment combination
unique_combinations <- unique(data[c("Plate", "Treatment")])
udata <- data.frame(
  Plate     = unique_combinations$Plate,
  Treatment = unique_combinations$Treatment
)

# Remove any zero eigenvalues (if they occur)
fitness <- fitness %>%
  filter(dominant_eigenvalues != 0)

# Attach Plate and Treatment
fitness$Plate     <- udata$Plate
fitness$Treatment <- udata$Treatment

fitness <- fitness %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("Control", "Larval Starvation")
    )
  )

# Plot fitness
FitP0 <- plot_violin(
  data  = fitness,
  x     = Treatment,
  y     = dominant_eigenvalues,
  y_lab = "Fitness"
)

FitP0

################################################################################
# 6. Age-specific reproduction curves                                          #
################################################################################

# Optional: mean per Treatment–Day (not strictly needed for plotting)
# mean_values <- aggregate(
#   Reproduction ~ Treatment + Day,
#   data = data,
#   FUN  = mean
# )

Repro <- ggplot(data, aes(x = Day, y = Reproduction, colour = Treatment)) +
  # Raw observations (jittered points)
  geom_point(
    aes(fill = Treatment),
    position = position_jitter(width = 0.12, height = 0),
    pch      = 21,
    color    = "black",
    alpha    = 0.25,
    size     = 4,
    stroke   = 0
  ) +
  # Mean ± 95% CI per day and treatment
  stat_summary(
    aes(group = Treatment, colour = Treatment),
    fun.data = mean_cl_normal,
    geom     = "errorbar",
    width    = 0.25,
    linewidth = 1.5,
    alpha    = 0.5
  ) +
  # Mean line
  stat_summary(
    fun  = mean,
    geom = "line",
    linewidth = 2
  ) +
  # Mean points
  stat_summary(
    fun   = mean,
    geom  = "point",
    size  = 4,
    stroke = 0.2,
    fill  = "white",
    shape = 21,
    aes(color = Treatment)
  ) +
  # Scales
  scale_x_continuous(
    "Day",
    breaks = sort(unique(data$Day))
  ) +
  scale_y_continuous(
    "Reproduction",
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_colour_manual(
    values = c(
      "Control"          = "#D55E00",
      "Larval Starvation" = "#0072B2"
    ),
    name   = "Treatment",
    labels = c("Control", "Starvation")
  ) +
  scale_fill_manual(
    values = c(
      "Control"          = "#D55E00",
      "Larval Starvation" = "#0072B2"
    ),
    name   = "Treatment",
    labels = c("Control", "Starvation")
  ) +
  guides(
    fill   = "none",
    colour = guide_legend(override.aes = list(size = 3, alpha = 1))
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(margin = margin(r = 8), size = 30),
    axis.text    = element_text(color = "black", size = 25),
    axis.ticks   = element_line(color = "black"),
    legend.position = "right"
  )

Repro

################################################################################
# End of script                                                                #
################################################################################

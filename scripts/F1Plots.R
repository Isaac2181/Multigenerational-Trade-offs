################################################################################
# Script: F1Plots.R
# Author: Isaac Harriz
#
# Purpose:
#   Generate F1 figures:
#     1) Violin + quasirandom plot for lifetime reproductive success (LRS)
#     2) Violin + quasirandom plot for rate-sensitive fitness (Leslie λ)
#     3) Age-specific reproduction curves by combined treatment (P0–F1)
#
# Expected project structure (example):
#   project/
#    ├ data/
#    │   └ F1.csv
#    └ scripts/
#        └ F1Plots.R
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
# 2. Generic violin + quasirandom plotting function                            #
################################################################################

plot_violin_f1 <- function(
    data,
    x,                     # grouping variable (factor/character)
    y,                     # numeric outcome (e.g., LRS)
    palette = NULL,        # vector of colors (one per group); NULL uses defaults
    y_lab = "Lifetime \n Reproductive Success",
    x_lab = NULL,
    title = NULL,
    subtitle = NULL,
    point_alpha = 0.4,
    point_size = 4,
    violin_alpha = 0.85,
    box_alpha = 0.95,
    width_violin = 0.9,
    width_box = 0.2
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
      fill  = "black",
      color = "black"
    ) +
    # Mean ± 95% CI
    stat_summary(
      fun.data = mean_cl_normal,
      geom     = "errorbar",
      width    = 0.3,
      color    = "black",
      linewidth = 1.6
    ) +
    scale_fill_manual(values = palette) +
    scale_color_manual(values = palette) +
    # Basic x-labels (we’ll overlay brackets via legendry)
    scale_x_discrete(
      labels = c("Control",
                 "Starvation",
                 "Control",
                 "Starvation")
    ) +
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
# 3. Load and clean F1 data                                                    #
################################################################################

F1 <- read.csv(here("data", "F1.csv"))

# Long format (one row per Plate–Day)
data <- pivot_longer(
  F1,
  cols      = starts_with("D"),
  names_to  = "Day",
  values_to = "Reproduction"
)

# Convert Day to numeric
data$Day <- as.numeric(gsub("D", "", data$Day))

# Remove plates with any NA in Reproduction
data <- data %>%
  group_by(Plate) %>%
  filter(!any(is.na(Reproduction))) %>%
  ungroup()

################################################################################
# 4. Lifetime Reproductive Success (LRS)                                       #
################################################################################

# Total reproduction per plate
LRS <- aggregate(
  Reproduction ~ Plate,
  data = data,
  FUN  = sum
)

# Merge in Treatment code per plate
LRS <- merge(
  LRS,
  data %>% select(Plate, Treatment) %>% distinct(),
  by = "Plate"
)

# Treatment coding:
#   "C(C)" = Control (Control Parents)
#   "L(C)" = Starved (Control Parents)
#   "C(L)" = Control (Starved Parents)
#   "L(L)" = Starved (Starved Parents)
LRS <- LRS %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("C(C)", "L(C)", "C(L)", "L(L)")
    )
  )

levels(LRS$Treatment) <- c(
  "Control(Control Parents)",
  "Starved(Control Parents)",
  "Control(Starved Parents)",
  "Starved(Starved Parents)"
)

# LRS plot
LRSF1 <- plot_violin_f1(
  data = LRS,
  x    = Treatment,
  y    = Reproduction
)

# Add parental-condition brackets via {legendry}
LRSF1 <- LRSF1 +
  guides(
    x = compose_stack(
      "axis_base",
      primitive_bracket(
        key_range_manual(
          start = c(1, 3),
          end   = c(2, 4),
          name  = c("Control", "Starvation")
        )
      )
    )
  )

################################################################################
# 5. Rate-sensitive fitness (Leslie λ)                                         #
################################################################################

# Compute dominant eigenvalue (λ) per plate
dominant_eigenvalues <- numeric(length(unique(data$Plate)))

for (i in unique(data$Plate)) {
  Leslie <- matrix(
    0,
    nrow = length(unique(data$Day)),
    ncol = length(unique(data$Day))
  )
  
  Leslie[1, ] <- data$Reproduction[data$Plate == i]
  diag(Leslie[-1, ]) <- 1
  
  eigen_res           <- eigen(Leslie)
  dominant_eigenvalue <- max(Re(eigen_res$values))
  dominant_eigenvalues[i] <- dominant_eigenvalue
}

# Data frame of λ
fitness <- as.data.frame(dominant_eigenvalues)

# Unique Plate–Treatment combinations
unique_combinations <- unique(data[c("Plate", "Treatment", "Offspring.Treatment")])
udata <- data.frame(
  Plate     = unique_combinations$Plate,
  Treatment = unique_combinations$Treatment
)

# Filter out zero eigenvalues if any
fitness <- fitness %>%
  filter(dominant_eigenvalues != 0)

fitness$Plate     <- udata$Plate
fitness$Treatment <- udata$Treatment

fitness$Treatment <- factor(
  fitness$Treatment,
  levels = c("C(C)", "L(C)", "C(L)", "L(L)")
)

levels(fitness$Treatment) <- c(
  "Control(Control Parents)",
  "Starved(Control Parents)",
  "Control(Starved Parents)",
  "Starved(Starved Parents)"
)

# Fitness plot
FitF1 <- plot_violin_f1(
  data  = fitness,
  x     = Treatment,
  y     = dominant_eigenvalues,
  y_lab = "Fitness"
)

FitF1 <- FitF1 +
  guides(
    x = compose_stack(
      "axis_base",
      primitive_bracket(
        key_range_manual(
          start = c(1, 3),
          end   = c(2, 4),
          name  = c("Control", "Starvation")
        )
      )
    )
  )

################################################################################
# 6. Age-specific reproduction curves                                          #
################################################################################

# (Optional) mean values per Treatment–Day
# mean_values <- aggregate(
#   Reproduction ~ Treatment + Day,
#   data = data,
#   FUN  = mean
# )

ReproF1 <- ggplot(data, aes(x = Day, y = Reproduction, colour = Treatment)) +
  # Raw observations (light jitter)
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
      "C(C)" = "#D55E00",
      "L(C)" = "#0072B2",
      "C(L)" = "#009E73",
      "L(L)" = "#CC79A7"
    ),
    name   = "Treatment \n (P0 - F1/F3)",
    labels = c(
      "C(C)" = "Control - Control",
      "C(L)" = "Starvation - Control",
      "L(C)" = "Control - Starvation",
      "L(L)" = "Starvation - Starvation"
    )
  ) +
  scale_fill_manual(
    values = c(
      "C(C)" = "#D55E00",
      "L(C)" = "#0072B2",
      "C(L)" = "#009E73",
      "L(L)" = "#CC79A7"
    ),
    name   = "Treatment \n (P0 - F1/F3)",
    labels = c(
      "C(C)" = "Control - Control",
      "C(L)" = "Starvation - Control",
      "L(C)" = "Control - Starvation",
      "L(L)" = "Starvation - Starvation"
    )
  ) +
  theme_classic() +
  theme(
    axis.title       = element_text(size = 14),
    axis.title.x     = element_text(size = 30),
    axis.title.y     = element_text(margin = margin(r = 8), size = 30),
    axis.text        = element_text(color = "black", size = 25),
    axis.ticks       = element_line(color = "black"),
    legend.title     = element_text(size = 20, face = "bold"),
    legend.text      = element_text(size = 15)
  ) +
  guides(
    colour = guide_legend(override.aes = list(size = 3, alpha = 1)),
    fill   = "none"
  )

ReproF1

################################################################################
# End of script                                                                #
################################################################################


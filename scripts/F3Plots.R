################################################################################
# Script: F3Plots.R
# Author: Isaac Harris
#
# Purpose:
#   Generate F3 figures:
#     1) Violin + quasirandom plot for lifetime reproductive success (LRS)
#     2) Violin + quasirandom plot for rate-sensitive fitness (Leslie λ)
#     3) Age-specific reproduction curves by combined treatment (P0–F3)
#
# Expected project structure (example):
#   project/
#    ├ data/
#    │   └ F3.csv
#    └ scripts/
#        └ F3Plots.R
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
# 2. Define generic violin + quasirandom plotting function                     #
################################################################################

plot_violin_f3 <- function(
    data,
    x,                       
    y,                       
    palette = NULL,          
    y_lab = "Lifetime\nReproductive Success",
    x_lab = "F3 Treatment",
    title = NULL,
    subtitle = NULL,
    point_alpha = 0.3,
    point_size = 2.5,
    violin_alpha = 0.8,
    width_violin = 0.8
) {
  x <- enquo(x); y <- enquo(y)
  
  # 1. Prepare data
  df <- data %>%
    filter(!is.na(!!x), !is.na(!!y)) %>%
    mutate(.group = factor(!!x))
  
  # 2. Setup Palette
  n_groups <- nlevels(df$.group)
  if (is.null(palette)) {
    palette <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7")[seq_len(n_groups)]
  }
  
  # 3. Dynamic Bracket Positions
  max_y <- max(pull(df, !!y), na.rm = TRUE)
  
  # --- Calculate Heights ---
  # Increase these multipliers to push brackets higher
  bracket_height <- max_y * 1.45       # The horizontal line height
  text_height    <- max_y * 1.55       # The text height (above the line)
  whisker_bottom <- max_y * 1.20       # How far down the whiskers reach
  
  p <- ggplot(df, aes(x = .group, y = !!y, fill = .group)) +
    # --- Geometries ---
    geom_violin(
      width = width_violin, alpha = violin_alpha,
      color = NA, trim = FALSE, scale = "width"
    ) +
    geom_quasirandom(
      color = "black", dodge.width = 0.7, varwidth = FALSE,
      alpha = point_alpha, size = point_size, stroke = 0
    ) +
    stat_summary(
      fun = mean, geom = "point", shape = 21,
      size = 3, fill = "white", color = "black", stroke = 1
    ) +
    stat_summary(
      fun.data = mean_cl_normal, geom = "errorbar",
      width = 0.2, color = "black", linewidth = 1
    ) +
    
    # --- Parent Grouping Brackets ---
    
    # === Bracket 1: Control Ancestors ===
    # Horizontal line
    annotate("segment", x = 1, xend = 2, y = bracket_height, yend = bracket_height, linewidth = 2) +
    # Left Whisker
    annotate("segment", x = 1, xend = 1, y = bracket_height, yend = whisker_bottom, linewidth = 2) +
    # Right Whisker
    annotate("segment", x = 2, xend = 2, y = bracket_height, yend = whisker_bottom, linewidth = 2) +
    # Text Label
    annotate("text", x = 1.5, y = text_height, label = "Control Ancestors", size = 10, fontface = "bold") +
    
    # === Bracket 2: Starved Ancestors ===
    # Horizontal line
    annotate("segment", x = 3, xend = 4, y = bracket_height, yend = bracket_height, linewidth = 2) +
    # Left Whisker
    annotate("segment", x = 3, xend = 3, y = bracket_height, yend = whisker_bottom, linewidth = 2) +
    # Right Whisker
    annotate("segment", x = 4, xend = 4, y = bracket_height, yend = whisker_bottom, linewidth = 2) +
    # Text Label
    annotate("text", x = 3.5, y = text_height, label = "Starved Ancestors", size = 10, fontface = "bold") +
    
    # --- Scales & Theme ---
    scale_fill_manual(values = palette) +
    scale_color_manual(values = palette) +
    scale_x_discrete(
      labels = c("Control", "Starvation", "Control", "Starvation")
    ) +
    # Increased expansion at the top (0.3) to make room for the higher brackets
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    
    labs(
      x = x_lab,
      y = y_lab,
      title = title,
      subtitle = subtitle
    ) +
    ylim(0,text_height) +
    guides(fill = "none", color = "none") +
    theme_classic(base_size = 14) +
    theme(
      legend.key.size = unit(2,"cm"),
      plot.title.position = "plot",
      plot.title        = element_text(face = "bold"),
      axis.title.y      = element_text(margin = margin(r = 15), size = 20),
      axis.title.x      = element_text(margin = margin(t = 15), size = 26, face = "bold"),
      axis.text         = element_text(color = "black", size = 24),
      axis.ticks        = element_line(color = "black"),
      coord.clip        = "off"
    ) 
  
  return(p)
}
################################################################################
# 3. Load and clean F3 data                                                    #
################################################################################

F3 <- read.csv(here("data", "F3.csv"))

# Long format (one row per Plate–Day)
data <- pivot_longer(
  F3,
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

# Merge in Treatment per plate
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
    Treatment = fct_recode(Treatment,
                           "Control(Control Parents)" = "C(C)",
                           "Starved(Control Parents)" = "L(C)",
                           "Control(Starved Parents)" = "C(L)",
                           "Starved(Starved Parents)" = "L(L)"
    ),
    Treatment = fct_relevel(Treatment,
                            "Control(Control Parents)",  # P0 Control -> F1 Control
                            "Starved(Control Parents)",  # P0 Control -> F1 Starvation
                            "Control(Starved Parents)",  # P0 Starvation -> F1 Control
                            "Starved(Starved Parents)"   # P0 Starvation -> F1 Starvation
    )
  )
# LRS plot
LRSF3 <- plot_violin_f3(
  data = LRS,
  x    = Treatment,
  y    = Reproduction,
  x_lab = "F3 Treatment"
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

fitness <- fitness %>%
  mutate(
    Treatment = fct_recode(Treatment,
                           "Control(Control Parents)" = "C(C)",
                           "Starved(Control Parents)" = "L(C)",
                           "Control(Starved Parents)" = "C(L)",
                           "Starved(Starved Parents)" = "L(L)"
    ),
    Treatment = fct_relevel(Treatment,
                            "Control(Control Parents)",  # P0 Control -> F1 Control
                            "Starved(Control Parents)",  # P0 Control -> F1 Starvation
                            "Control(Starved Parents)",  # P0 Starvation -> F1 Control
                            "Starved(Starved Parents)"   # P0 Starvation -> F1 Starvation
    )
  )

# Fitness plot
FitF3 <- plot_violin_f3(
  data  = fitness,
  x     = Treatment,
  y     = dominant_eigenvalues,
  y_lab = "Fitness",
  x_lab = "F3 Treatment"
)


################################################################################
# 6. Age-specific reproduction curves                                          #
################################################################################


# Build age-specific reproduction plot
ReproF3 <- ggplot(data, aes(x = Day, y = Reproduction, colour = Treatment)) +
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
  # Axes
  scale_x_continuous(
    "Day",
    breaks = sort(unique(data$Day))
  ) +
  scale_y_continuous(
    "Reproduction",
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  # Manual colours and labels for combined treatments
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
    legend.key.size = unit(2,"cm"),
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

ReproF3

################################################################################
# End of script                                                                #
################################################################################

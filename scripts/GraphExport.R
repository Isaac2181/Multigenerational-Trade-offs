################################################################################
# Script: GraphExport.R
# Author: Isaac Harris
#
# Purpose:
#   Assemble and export all main and supplementary figures:
#     - Sup Fig 1: P0 survival + forest
#     - Sup Fig 2: F1 + F3 survival + forests
#     - Fig 2: P0 LRS, fitness, age-specific reproduction
#     - Fig 3: F1/F3 LRS, fitness, age-specific reproduction
#     - Fig 4: Simulation results
#
# Note:
#   This script assumes that the following objects already exist in the
#   environment (from previously sourced scripts):
#   P0, F1, F3         : ggsurvplot results (or plain ggplot objects)
#   forest, forestF1,
#   forestF3           : forest plot ggplot objects
#   LRS, LRSF1, LRSF3  : LRS plots
#   FitP0, FitF1, FitF3: fitness plots
#   Repro, ReproF1,
#   ReproF3            : age-specific reproduction plots
#   Plot, Plot2        : simulation plots
################################################################################

library(ggplot2)
library(patchwork)
library(here)

# Ensure export directory exists
export_dir <- here("plots", "Final", "Export")
if (!dir.exists(export_dir)) {
  dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
}

################################################################################
# Helper functions                                                             #
################################################################################

# General text-size bumping helper
bump_text <- function(p, base = 25, axis = 25, legend = 25,
                      add_grid = FALSE) {
  th <- theme(
    text         = element_text(size = base),
    axis.title   = element_text(size = base),
    axis.text    = element_text(size = axis),
    legend.title = element_text(size = legend),
    legend.text  = element_text(size = legend),
    strip.text   = element_text(size = base)
  )
  if (add_grid) {
    th <- th +
      theme(
        panel.grid.major = element_line(colour = "grey"),
        panel.grid.minor = element_line(colour = "grey")
      )
  }
  p + th
}

# Y-axis emphasis helper
y_enl <- function(p, title_size = 40, text_size = 40) {
  p + theme(
    axis.title.y = element_text(size = title_size),
    axis.text.y  = element_text(size = text_size)
  )
}

################################################################################
# Supplementary Figure 1: P0 survival + forest                                 #
################################################################################

# Extract ggplot from ggsurvplot object if needed
P0_plot <- if (inherits(P0splot, "ggsurvplot") || (is.list(P0splot) && !is.null(P0splot$plot))) {
  P0splot$plot
} else {
  P0splot
}

forest_sup1 <- bump_text(forest)
P0_plot_sup1 <- bump_text(P0_plot)

sup_fig1 <- P0_plot_sup1 / forest_sup1 +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 12))
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

sup_fig1

ggsave(
  filename = file.path(export_dir, "sup_fig1.svg"),
  plot     = sup_fig1,
  width    = 16,
  height   = 12,
  dpi      = 1200
)

################################################################################
# Supplementary Figure 2: F1 + F3 survival + forests                           #
################################################################################

F1_plot <- if (inherits(F1splot, "ggsurvplot") || (is.list(F1splot) && !is.null(F1splot$plot))) {
  F1splot$plot
} else {
  F1splot
}

F3_plot <- if (inherits(F3splot, "ggsurvplot") || (is.list(F3splot) && !is.null(F3splot$plot))) {
  F3splot$plot
} else {
  F3splot
}

# For these, we also add panel grids
F1_plot_sup2  <- bump_text(F1_plot,  add_grid = TRUE)
F3_plot_sup2  <- bump_text(F3_plot,  add_grid = TRUE)
forestF1_sup2 <- bump_text(forestF1, add_grid = TRUE)
forestF3_sup2 <- bump_text(forestF3, add_grid = TRUE)

sup_fig2 <- (forestF1_sup2 + forestF3_sup2) / F1_plot_sup2 / F3_plot_sup2 +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 12))
  ) &
  theme(legend.position = "right")

sup_fig2 <- sup_fig2 + plot_layout(guides = "collect")

sup_fig2
ggsave(
  filename = file.path(export_dir, "sup_fig2.svg"),
  plot     = sup_fig2,
  width    = 16,
  height   = 12,
  dpi      = 1200
)

################################################################################
# Figure 2: P0 LRS + fitness + age-specific reproduction                       #
################################################################################

Repro_big <- bump_text(Repro)
Repro_big <- Repro_big +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

FitP0_big <- bump_text(FitP0, axis = 24)
LRS_big   <- bump_text(LRS_plot)

fig2 <- (LRS_big + FitP0_big) / Repro_big +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

fig2

ggsave(
  filename = file.path(export_dir, "Fig2.svg"),
  plot     = fig2,
  width    = 16,
  height   = 12,
  dpi      = 1200
)

################################################################################
# Figure 3: F1/F3 LRS + fitness + age-specific reproduction                    #
################################################################################

# Layout design for patchwork
des <- ("aee
         bee
         cff
         dff")

# Bump sizes for F1/F3
ReproF1_big <- bump_text(ReproF1, base = 45, axis = 45, legend = 45)
ReproF3_big <- bump_text(ReproF3, base = 45, axis = 45, legend = 45)
FitF1_big   <- bump_text(FitF1, axis = 30)
FitF3_big   <- bump_text(FitF3, axis = 30)
LRSF1_big   <- bump_text(LRSF1, axis = 30)
LRSF3_big   <- bump_text(LRSF3, axis = 30)

ReproF1_big <- ReproF1_big +
  theme(
    axis.title.x = element_text(size = 45),
    axis.title.y = element_text(size = 45),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ReproF3_big <- ReproF3_big +
  theme(
    axis.title.x = element_text(size = 45),
    axis.title.y = element_text(size = 45),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

FitF1_big <- y_enl(FitF1_big)
FitF3_big <- y_enl(FitF3_big)
LRSF1_big <- y_enl(LRSF1_big)
LRSF3_big <- y_enl(LRSF3_big)

fig3 <- LRSF1_big + LRSF3_big + FitF1_big + FitF3_big + ReproF1_big + ReproF3_big +
  plot_annotation(tag_levels = "A") +
  plot_layout(design = des)

fig3 <- fig3 +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    plot.tag        = element_text(size = 40, face = "bold")
  )

fig3

ggsave(
  filename = file.path(export_dir, "Fig3.svg"),
  plot     = fig3,
  width    = 38,
  height   = 25,
  dpi      = 1200
)

################################################################################
# Figure 4: Simulation results                                                 #
################################################################################

# Simulation plots (Plot, Plot2) should come from your simulation script
Plot_big  <- bump_text(Plot,  add_grid = TRUE)
Plot2_big <- bump_text(Plot2, add_grid = TRUE)

fig4 <- (Plot_big / Plot2_big) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "inside")

fig4

ggsave(
  filename = file.path(export_dir, "Fig4.svg"),
  plot     = fig4,
  width    = 16,
  height   = 12,
  dpi      = 1200
)

################################################################################
# Supplementary Figure 3: Simulation parameters heatmap                        #
################################################################################
supfig3 <- ggplot(param_grid, aes(x = f1, y = f3)) +
  # The main heatmap with white outlines
  geom_tile(aes(fill = outcome), colour = "black", size = 0.2) +
  # The highlight square
  geom_rect(aes(xmin = orig_f1 - 0.01, xmax = orig_f1 + 0.01, 
                ymin = orig_f3 - 0.005, ymax = orig_f3 + 0.005, 
                color = "Simulation Parameters"), 
            fill = NA, size = 1.2) +
  scale_fill_gradient2(low = "#D55E00", mid = "white", high = "#0072B2", 
                       midpoint = 0.5, name = "Final Proportion of\nTransgenerational Strategy\n") +
  
  scale_color_manual(name = "", values = c("Simulation Parameters" = "red")) +
  # Visual tweaks
  labs(x = "Relative intergenerational fitness (F1)",
       y = "Relative transgenerational fitness (F3)") +
  theme_minimal() +
  theme(
    axis.title   = element_text(size = 30),
    axis.text    = element_text(size = 25),
    legend.key.size  = unit(1.5,"cm"),
    panel.grid.major = element_line(colour = "grey"),
    panel.grid.minor = element_line(colour = "grey"),
    legend.title = element_text(size = 25),
    legend.text  = element_text(size = 25)
  ) + guides(color = guide_legend(override.aes = list(fill = NA))) # Keeps the legend box hollow

ggsave(
  filename = file.path(r"(C:\Users\isaac\OneDrive\Desktop\OneDrive - University of East Anglia\Projects\Project- L1arrest\Github\Multigenerational-Trade-offs\plots\Final\Export)","supfig3.svg"),
  plot     = supfig3,
  width    = 18,
  height   = 12,
  dpi      = 1200
)
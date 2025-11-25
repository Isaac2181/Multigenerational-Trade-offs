################################################################################
# Script: P0Surv.R
# Author: Isaac Harris
#
# Purpose:
#   - Fit Cox mixed-effects model for P0 survival
#   - Check proportional hazards
#   - Create:
#       * Kaplan–Meier survival plot (P0)
#       * Forest-style plot of treatment effect (log hazard ratio)
#
# Expected project structure:
#   project/
#    ├ data/
#    │   └ P0.csv
#    └ scripts/
#        └ P0Surv.R
################################################################################

#############################
# Load required packages    #
#############################

library(here)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(coxme)

##########################################
# Read survival data (P0 generation)     #
##########################################

data <- read.csv(here("data", "P0GLMM.csv"))
# Columns required:
#   Age (numeric time)
#   Event (0/1 status)
#   Treatment (factor)
#   Plate (random effect for coxme)

##########################################
# Survival model                         #
##########################################

surv1 <- Surv(time = data$Age, event = data$Event)
fit1  <- survfit(surv1 ~ data$Treatment)

cox <- coxme(Surv(Age, Event) ~ Treatment + (1 | Plate), data = data)
summary(cox)

##########################################
# Helper: forest-style plot from coxme   #
##########################################

meforest <- function(cox, reference_label) {
  require(ggplot2)
  require(AICcmodavg)
  
  store <- matrix(nrow = length(cox$coefficients) + 1, ncol = 4)
  store[1, ] <- c(reference_label, 0, NA, NA)
  
  for (i in seq_along(cox$coefficients)) {
    m  <- cox$coefficients[i]
    se <- extractSE(cox)[i]
    L  <- m - 1.96 * se
    U  <- m + 1.96 * se
    store[i + 1, ] <- c(names(cox$coefficients[i]), m, L, U)
  }
  
  colnames(store) <- c("Treatment", "mean", "CIL", "CUL")
  store <- as.data.frame(store)
  store[, 2:4] <- lapply(store[, 2:4], \(x) as.numeric(as.character(x)))
  
  ggplot(store, aes(x = mean, xmin = CIL, xmax = CUL, y = Treatment)) +
    geom_point(size = 4, shape = 19, colour = c("#D55E00", "#0072B2")) +
    geom_errorbarh(height = 0.25, linewidth = 1.5, colour = c("#D55E00", "#0072B2")) +
    geom_vline(xintercept = 0, linetype = "dotted", size = 1.2) +
    xlab("Hazard Ratio (log scale)") +
    ylab("") +
    theme_classic() +
    theme(
      axis.text  = element_text(size = 18),
      axis.title = element_text(size = 20)
    ) +
    expand_limits(x = c(-1.5, 0.5))
}

##########################################
# Forest plot                            #
##########################################

forest <- meforest(cox, "F")   # "F" = reference Treatment
forest

##########################################
# Kaplan–Meier survival plot             #
##########################################

P0 <- ggsurvplot(
  fit1,
  data = data,
  palette = c("#D55E00", "#0072B2"),
  legend.labs = c("Control", "Starvation"),
  legend.title = "Treatment",
  size = 2,
  censor.size = 5,
  censor.colour = "black",
  xlab = "Day",
  ggtheme = theme_classic() +
    theme(
      plot.title.position = "plot",
      axis.title.x = element_text(margin = margin(r = 8), size = 30),
      axis.title.y = element_text(margin = margin(r = 8), size = 30),
      axis.text    = element_text(size = 25),
      axis.ticks   = element_line(color = "black"),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text  = element_text(size = 15),
      legend.position = "right"
    )
)

P0  # Plot is returned as ggsurvplot object

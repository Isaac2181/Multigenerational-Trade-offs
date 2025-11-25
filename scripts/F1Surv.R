################################################################################
# Script: F1Surv.R
# Purpose:
#   - Fit Cox mixed-effects model for F1 survival
#   - Generate:
#       * Kaplan–Meier survival curves by TreatmentCode
#       * Forest-style plot of log-hazard ratios
################################################################################

#############################
# 1. Load required packages #
#############################

library(here)
library(ggplot2)
library(survival)
library(survminer)
library(coxme)

##########################################
# 2. Read survival data (F1 generation)  #
##########################################

# F1.csv must live in the project's /data folder and contain:
#   Age            : numeric time
#   Event          : 0/1 event indicator
#   TreatmentCode  : factor with levels like F, cl, lc, ll
#   Plate          : grouping factor for random effect
data <- read.csv(here("data", "F1_survival.csv"))

##########################################
# 3. Survival object and Cox ME model    #
##########################################

surv1 <- Surv(time = data$Age, event = data$Event)
fit1  <- survfit(surv1 ~ data$TreatmentCode)

cox <- coxme(Surv(Age, Event) ~ TreatmentCode + (1 | Plate), data = data)
summary(cox)

##########################################
# 4. Forest-plot helper for coxme model  #
##########################################

# meforest():
#   - cox             : coxme model
#   - reference_label : label for the reference level (e.g. "F")
#
# Note:
#   Coefficients from coxme are on the log-hazard scale: log(HR).

meforest <- function(cox, reference_label) {
  require(AICcmodavg)
  require(ggplot2)
  
  # one row for reference + one per coefficient
  store <- matrix(nrow = length(cox$coefficients) + 1, ncol = 4)
  store[1, ] <- c(reference_label, 0, NA, NA)
  
  for (i in seq_along(cox$coefficients)) {
    mean <- cox$coefficients[i]
    se   <- extractSE(cox)[i]
    cil  <- mean - 1.96 * se
    cul  <- mean + 1.96 * se
    store[i + 1, ] <- c(names(cox$coefficients[i]), mean, cil, cul)
  }
  
  colnames(store) <- c("Treatment", "mean", "CIL", "CUL")
  store <- as.data.frame(store)
  
  store$mean <- as.numeric(as.character(store$mean))
  store$CIL  <- as.numeric(as.character(store$CIL))
  store$CUL  <- as.numeric(as.character(store$CUL))
  
  print(head(store))
  
  ggplot(
    store,
    aes(x = mean, xmin = CIL, xmax = CUL, y = Treatment, colour = Treatment)
  ) +
    geom_point(
      size   = 4,
      shape  = 19,
      colour = c("#D55E00", "#0072B2", "#009E73", "#CC79A7")
    ) +
    geom_errorbarh(
      height = 0.25,
      linewidth = 0.9,
      colour = c("#D55E00", "#0072B2", "#009E73", "#CC79A7")
    ) +
    theme_classic() +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.8) +
    xlab("Hazard Ratio") +
    ylab("") +
    theme(
      axis.text  = element_text(size = 18),
      axis.title = element_text(size = 20)
    ) +
    scale_y_discrete(
      limits = levels(store$Treatment),
      labels = c(
        "F"               = "",
        "TreatmentCodecl" = "",
        "TreatmentCodelc" = "",
        "TreatmentCodell" = ""
      )
    ) +
    expand_limits(x = c(-1.5, 0.5))
}

##########################################
# 5. Forest plot & PH diagnostics        #
##########################################

forestF1 <- meforest(cox, "F")

# PH diagnostics
inter <- ggcoxzph(cox.zph(cox))

# p-values for fixed effects 
coxp <- signif(summary(cox)$coefficients[, 5], digits = 3)

##########################################
# 6. Kaplan–Meier survival plot (F1)     #
##########################################

F1 <- ggsurvplot(
  fit1,
  data = data,
  legend.labs = c(
    "Control - Control",
    "Control - Larval Starvation",
    "Larval Starvation - Control",
    "Larval Starvation - Larval Starvation"
  ),
  legend.title = "P0 - F1/F3",
  size         = 2,
  censor.size  = 5,
  xlab         = "Day",
  palette      = c("#D55E00", "#0072B2", "#009E73", "#CC79A7"),
  ggtheme      = theme_classic() +
    theme(
      plot.title.position = "plot",
      plot.title          = element_text(face = "bold"),
      axis.title.x        = element_text(margin = margin(r = 8), size = 30),
      axis.title.y        = element_text(margin = margin(r = 8), size = 30),
      axis.text           = element_text(color = "black", size = 25),
      axis.ticks          = element_line(color = "black")
    )
)

F1  # prints the survival plot
################################################################################
# End of script                                                                #
################################################################################

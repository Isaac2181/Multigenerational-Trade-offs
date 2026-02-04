################################################################################
# Script: F3Surv.R
# Purpose:
#   - Fit Cox mixed-effects model for F3 survival
#   - Generate:
#       * Kaplan–Meier survival curve by TreatmentCode
#       * Forest-style plot of log-hazard ratios for TreatmentCode
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
# 2. Read survival data (F3 generation)  #
##########################################

# F3.csv should live in your project's /data directory and contain at least:
#   Age           : numeric time
#   Event         : 0/1 event indicator
#   TreatmentCode : factor for treatment combinations (F, cl, lc, ll, etc.)
#   Plate         : grouping factor for random effect
data <- read.csv(here("data", "F3_survival.csv"))

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
#   - cox : coxme model
#   - F   : label for the reference treatment
#
# Note:
#   Coefficients are on the log-hazard scale (log HR).

meforest <- function(cox, F) {  # Ed's forest-plot function
  require(AICcmodavg)
  require(ggplot2)
  
  store <- matrix(nrow = length(cox$coefficients) + 1, ncol = 4)
  ref   <- c(paste(F), 0, NA, NA)
  store[1, ] <- ref
  
  for (x in seq_along(cox$coefficients)) {
    y    <- x + 1
    mean <- cox$coefficients[x]
    CIL  <- cox$coefficients[x] - (1.96 * extractSE(cox)[x])
    CUL  <- cox$coefficients[x] + (1.96 * extractSE(cox)[x])
    store[y, 1:4] <- c(names(cox$coefficients[x]), mean, CIL, CUL)
  }
  
  colnames(store) <- c("Treatment", "mean", "CIL", "CUL")
  store <- as.data.frame(store)
  
  store$mean <- as.numeric(as.character(store$mean))
  store$CIL  <- as.numeric(as.character(store$CIL))
  store$CUL  <- as.numeric(as.character(store$CUL))
  print(head(store))
  
  forest <- ggplot(
    store,
    aes(x = mean, xmax = CUL, xmin = CIL, y = Treatment, colour = Treatment)
  ) +
    geom_point(
      size   = 4,
      shape  = 19,
      colour = c("#D55E00", "#0072B2", "#009E73", "#CC79A7")
    ) +
    geom_errorbarh(
      height = 0.25,
      size   = 0.9,
      colour = c("#D55E00", "#0072B2", "#009E73", "#CC79A7")
    ) +
    theme_classic() +
    geom_vline(xintercept = 0, linetype = "dotted", size = 0.8) +
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
  
  return(forest)
}

##########################################
# 5. Forest plot & PH diagnostics        #
##########################################

forestF3 <- meforest(cox, "F")
forestF3

# PH diagnostics (ggcoxzph object)
Trans <- ggcoxzph(cox.zph(cox))

# p-values for fixed effects (if needed)
coxp <- signif(summary(cox)$coefficients[, 5], digits = 3)

##########################################
# 6. Kaplan–Meier survival plot (F3)     #
##########################################

F3splot <- ggsurvplot(
  fit1,
  data = data,
  legend.labs = c(
    "Control - Control",
    "Starvation - Control",
    "Control - Starvation",
    "Starvation - Starvation"
  ),
  legend.title = "P0 - F1/F3",
  size         = 2,
  censor.size  = 5,
  xlab         = "Day",
  palette      = c("#D55E00", "#009E73", "#0072B2", "#CC79A7"),
  ggtheme      = theme_classic() +
    theme(
      plot.title.position = "plot",
      plot.title          = element_text(face = "bold"),
      axis.title.x        = element_text(margin = margin(r = 8), size = 30),
      axis.title.y        = element_text(margin = margin(r = 8), size = 30),
      axis.text           = element_text(color = "black", size = 25),
      axis.ticks          = element_line(color = "black"),
      legend.title        = element_blank(),
      legend.text         = element_text(size = 15)
    )
)

F3splot  # print KM plot

################################################################################
# End of script                                                                #
################################################################################

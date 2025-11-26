################################################################################
# Script: P0.R
# Author: Isaac Harris
#
# Purpose:
#   Analyse age-specific reproduction, rate-sensitive fitness (Leslie matrices),
#   lifetime reproductive success (LRS), and survival for P0 worms under
#   different treatments.
#
# Expected project structure (example):
#   project/
#    ├ data/
#    │   ├ P0repro.csv
#    │   └ P0GLMM.csv
#    └ scripts/
#        └ P0.R
#
# The script uses the {here} package for file paths, so it should be run
# from within the project directory
################################################################################

#############################
# 1. Load required packages #
#############################

library(tidyverse)   # data wrangling + pivot_longer
library(DHARMa)      # residual diagnostics for (G)LMMs
library(glmmTMB)     # generalized linear mixed models
library(car)         # Anova() for type-III tests
library(splines)     # ns() for spline smoothing
library(easystats)   # meta-package; used mainly via performance::
library(emmeans)     # estimated marginal means
library(here)        # reproducible file paths

#############################################
# 2. Read in data using here()              #
#############################################

# Age-specific reproduction data
P0 <- read.csv(here("data", "P0repro.csv"))

########################################
# 3. Tidy reproduction data (P0repro) #
########################################

# Reshape from wide (one column per day) to long format:
#   - 'Day' column: day of reproduction (D1, D2, ...)
#   - 'Reproduction' column: number of offspring
P0 <- pivot_longer(
  P0,
  cols      = starts_with("D"),
  names_to  = "Day",
  values_to = "Reproduction"
)

# Remove 'D' prefix from Day (e.g., D1 -> 1) and treat as factor
P0$Day <- as.factor(gsub("D", "", P0$Day))

# Coerce Treatment to factor
P0$Treatment <- as.factor(P0$Treatment)

# Remove any plates that have at least one NA in Reproduction across days
P0 <- P0 %>%
  group_by(Plate) %>%
  filter(!any(is.na(Reproduction)))

# Restrict to early days only (here days ≤ 5)
P0 <- P0 %>%
  group_by(Plate) %>%
  filter(as.numeric(Day) <= 5) %>%
  ungroup()


##########################################################
# 4. Age-specific reproduction: GLMM with spline on Day  #
##########################################################

# Model:
#   Reproduction ~ Treatment * f(Day) + (1 | Plate)
# where f(Day) is a natural cubic spline of Day with df = 4.
P0_age <- glmmTMB(
  Reproduction ~ Treatment * ns(as.numeric(Day), df = 4) + (1 | Plate),
  data   = P0,
  family = nbinom2
)

# Residual diagnostics: check dispersion and zero-inflation
simulationoutput <- simulateResiduals(fittedModel = P0_age, plot = TRUE)
testDispersion(simulationoutput)
testZeroInflation(simulationoutput)

# Allow dispersion to vary with Treatment and Day (spline)
P0_age_disp <- update(
  P0_age,
  dispformula = ~ Treatment * ns(as.numeric(Day), df = 4) + (1 | Plate),
  family      = nbinom2(),
  na.action   = na.fail
)

# Diagnostics for updated model
simulationoutput <- simulateResiduals(fittedModel = P0_age_disp, plot = TRUE)
testDispersion(simulationoutput)
testZeroInflation(simulationoutput)

# Model summary and type-III Anova
summary(P0_age_disp)
Anova(P0_age_disp, type = 3)

# Marginal means for Treatment
emmeans(P0_age_disp, ~ Treatment)

# Additional model checks (from easystats / performance)
performance::check_model(P0_age_disp)
performance::check_predictions(P0_age_disp)


##############################################
# 5. Rate-sensitive fitness (Leslie matrix)  #
##############################################

# Goal:
#   For each plate, construct a Leslie matrix using age-specific reproduction
#   and estimate rate-sensitive fitness as the dominant eigenvalue (λ).

# Unique plate IDs
plates <- unique(P0$Plate)

# Create a numeric vector to store plate-specific dominant eigenvalues
dominant_eigenvalues <- numeric(length(plates))

# Loop over plates to construct Leslie matrices and compute λ
for (p in seq_along(plates)) {
  
  plate_id <- plates[p]
  
  # Leslie matrix dimension = number of unique days
  n_age_classes <- length(unique(P0$Day))
  Leslie <- matrix(0, nrow = n_age_classes, ncol = n_age_classes)
  
  # First row: age-specific reproduction for this plate
  Leslie[1, ] <- P0$Reproduction[P0$Plate == plate_id]
  
  # Subdiagonal: survival from age x to x+1 (here assumed to be 1)
  diag(Leslie[-1, ]) <- 1
  
  # Eigen decomposition
  eigen_res <- eigen(Leslie)
  
  # Dominant eigenvalue (largest real part)
  dominant_eigenvalue <- max(Re(eigen_res$values))
  
  # Store plate-specific λ
  dominant_eigenvalues[p] <- dominant_eigenvalue
}

# Convert dominant eigenvalues to a data frame
fitness <- as.data.frame(dominant_eigenvalues)

# Create a data frame with unique Plate–Treatment combinations
unique_combinations <- unique(P0[c("Plate", "Treatment")])
udata <- data.frame(
  Plate     = unique_combinations$Plate,
  Treatment = unique_combinations$Treatment
)

# Filter out any zero eigenvalues (if they occur)
fitness <- fitness %>%
  filter(dominant_eigenvalues != 0)

# Attach Plate and Treatment to fitness dataframe
fitness$Plate     <- udata$Plate
fitness$Treatment <- udata$Treatment

########################################
# 6. Model for rate-sensitive fitness  #
########################################

# Model:
#   dominant_eigenvalues ~ Treatment + (1 | Plate)
#   with Gamma errors and log link.
P0_fit <- glmmTMB(
  dominant_eigenvalues ~ Treatment + (1 | Plate),
  data   = fitness,
  family = Gamma(link = "log")
)

# Diagnostics
simulationoutput <- simulateResiduals(fittedModel = P0_fit, plot = TRUE)
testDispersion(simulationoutput)
testZeroInflation(simulationoutput)

performance::check_model(P0_fit)
performance::check_predictions(P0_fit)

# Summary and treatment contrasts
summary(P0_fit)
contrast(emmeans(P0_fit, ~ Treatment))


###############################################
# 7. Lifetime reproductive success (LRS)      #
###############################################

# Aggregate reproduction across all days for each Plate–Treatment combination
LRSP0 <- aggregate(
  Reproduction ~ Plate + Treatment,
  data = P0,
  FUN  = sum
)

# Model:
#   LRS (total offspring) ~ Treatment + (1 | Plate)
#   with COM-Poisson family to allow for over- / under-dispersion.
P0_LRS <- glmmTMB(
  Reproduction ~ Treatment + (1 | Plate),
  data   = LRSP0,
  family = compois()
)

# Diagnostics
simulationoutput <- simulateResiduals(fittedModel = P0_LRS, plot = TRUE)
testDispersion(simulationoutput)
testZeroInflation(simulationoutput)

performance::check_model(P0_LRS)
performance::check_predictions(P0_LRS)

# Type-III Anova and treatment contrasts
Anova(P0_LRS, type = 3)
contrast(emmeans(P0_LRS, ~ Treatment))


########################
# 8. Survival analysis #
########################

# Read survival data
P0_surv <- read.csv(here("data", "P0GLMM.csv"))

# Convert Day to numeric, coercing non-numeric values to NA
P0_surv$Day <- suppressWarnings(as.numeric(as.character(P0_surv$Day)))

# Fill NA Day values by carrying forward and adding 1:
#   If Day[i] is NA but Day[i-1] is not, then Day[i] = Day[i-1] + 1.
for (i in 2:nrow(P0_surv)) {
  if (is.na(P0_surv$Day[i]) && !is.na(P0_surv$Day[i - 1])) {
    P0_surv$Day[i] <- P0_surv$Day[i - 1] + 1
  }
}

# Quick check of the first 10 rows
head(P0_surv, 10)

# Ensure Age is numeric
P0_surv$Age <- as.numeric(P0_surv$Age)

# Survival GLMM:
#   status (alive/dead) ~ Treatment + Age + (1 | Plate/Worm)
#   with binomial errors and zero-inflation modelled as a function of Age.
P0_sur <- glmmTMB(
  status ~ Treatment + Age + (1 | Plate/Worm),
  family    = binomial,
  ziformula = ~ Age,
  data      = P0_surv
)

# Diagnostics
simulationoutput <- simulateResiduals(fittedModel = P0_sur, plot = TRUE)
testDispersion(simulationoutput)
testZeroInflation(simulationoutput)

performance::check_predictions(P0_sur)


################################################################################
# End of script
################################################################################

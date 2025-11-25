################################################################################
# Script: F3.R
# Author: Isaac Harris
#
# Purpose:
#   Analyse age-specific reproduction, rate-sensitive fitness (Leslie matrices),
#   lifetime reproductive success (LRS), and survival for F3 worms under
#   different treatments and parental treatments.
#
# Expected project structure (example):
#   project/
#    ├ data/
#    │   ├ F3.csv             # reproduction data
#    │   └ F3_survival.csv    # survival data (was .../Mol/F3.csv)
#    └ scripts/
#        └ F3_reproduction_survival_analysis.R
################################################################################

#############################
# 1. Load required packages #
#############################

library(tidyverse)   # dplyr, tidyr, ggplot2, etc.
library(DHARMa)      # residual diagnostics for (G)LMMs
library(glmmTMB)     # GLMMs
library(car)         # Anova()
library(splines)     # ns()
library(emmeans)     # marginal means, if needed later
library(here)        # reproducible file paths

########################################
# 2. Read and tidy reproduction data   #
########################################

# Reproduction data (originally "F3.csv" in L1arrest/data)
F3 <- read.csv(here("data", "F3.csv"))

# Reshape from wide (one column per day) to long:
#   - 'Day' column: day of reproduction (D1, D2, ...)
#   - 'Reproduction' column: number of offspring
F3 <- pivot_longer(
  F3,
  cols      = starts_with("D"),
  names_to  = "Day",
  values_to = "Reproduction"
)

# Remove 'D' prefix from 'Day' and convert to factor
F3$Day <- as.factor(gsub("D", "", F3$Day))

# Ensure treatment variables are factors
F3$Parental.Treatment  <- as.factor(F3$Parental.Treatment)
F3$Offspring.Treatment <- as.factor(F3$Offspring.Treatment)
F3$Treatment           <- as.factor(F3$Treatment)

# Remove any plates that have at least one NA in Reproduction
F3 <- F3 %>%
  group_by(Plate) %>%
  filter(!any(is.na(Reproduction))) %>%
  ungroup()

# Optional observation ID (kept from original script)
F3$obsID <- seq_len(nrow(F3))


#########################################################
# 3. Age-specific reproduction model with splines on Day
#########################################################

# Model:
#   Reproduction ~ Treatment * f(Day) + (1 | Plate)
#   with negative binomial errors (nbinom2) and zero-inflation by Day.
F3_age <- glmmTMB(
  Reproduction ~ Treatment * ns(as.numeric(Day), df = 3) + (1 | Plate),
  family    = nbinom2(),
  ziformula = ~ Day,
  control   = glmmTMBControl(
    optimizer = optim,
    optArgs   = list(method = "BFGS")
  ),
  data      = F3
)

# Diagnostics
simulationoutput <- simulateResiduals(fittedModel = F3_age, plot = TRUE)
testDispersion(simulationoutput)
testZeroInflation(simulationoutput)
performance::check_predictions(F3_age)


##############################################
# 4. Rate-sensitive fitness (Leslie matrix)  #
##############################################

# Goal:
#   For each plate, construct a Leslie matrix using age-specific reproduction
#   and estimate rate-sensitive fitness as the dominant eigenvalue (λ).

plates <- unique(F3$Plate)

# Vector to store plate-specific dominant eigenvalues
dominant_eigenvalues <- numeric(length(plates))

# Loop over plates
for (p in seq_along(plates)) {
  
  plate_id <- plates[p]
  
  # Number of age classes = number of unique days
  n_age_classes <- length(unique(F3$Day))
  Leslie <- matrix(0, nrow = n_age_classes, ncol = n_age_classes)
  
  # First row: reproduction for this plate
  Leslie[1, ] <- F3$Reproduction[F3$Plate == plate_id]
  
  # Subdiagonal: survival from age x to x+1 (assumed 1)
  diag(Leslie[-1, ]) <- 1
  
  # Eigen decomposition and dominant eigenvalue
  eigen_res           <- eigen(Leslie)
  dominant_eigenvalue <- max(Re(eigen_res$values))
  dominant_eigenvalues[p] <- dominant_eigenvalue
}

# Convert eigenvalues to a data frame
fitness <- as.data.frame(dominant_eigenvalues)

# Unique Plate–Treatment combination data
unique_combinations <- unique(F3[c("Plate", "Parental.Treatment",
                                   "Offspring.Treatment", "Treatment")])

udata <- data.frame(
  Plate             = unique_combinations$Plate,
  ParentalTreatment = unique_combinations$Parental.Treatment,
  OffspringTreatment= unique_combinations$Offspring.Treatment,
  Treatment         = unique_combinations$Treatment
)

# Remove any zero eigenvalues (if they occur)
fitness <- fitness %>%
  filter(dominant_eigenvalues != 0)

# Attach plate and treatment info to fitness dataframe
fitness$Plate              <- udata$Plate
fitness$ParentalTreatment  <- udata$ParentalTreatment
fitness$OffspringTreatment <- udata$OffspringTreatment
fitness$Treatment          <- udata$Treatment


########################################
# 5. Model for rate-sensitive fitness  #
########################################

# Model:
#   dominant_eigenvalues ~ Treatment + (1 | Plate)
#   with Gamma errors and log link.
F3_fit <- glmmTMB(
  dominant_eigenvalues ~ Treatment + (1 | Plate),
  data        = fitness,
  family      = Gamma(link = "log"),
  control     = glmmTMBControl(
    optimizer = optim,
    optArgs   = list(method = "BFGS")
  ),
  dispformula = ~ (1 | Plate)
)

# Diagnostics
simulationoutput <- simulateResiduals(fittedModel = F3_fit, plot = TRUE)
testDispersion(simulationoutput)
testZeroInflation(simulationoutput)

performance::check_predictions(F3_fit)


###############################################
# 6. Lifetime reproductive success (LRS)      #
###############################################

# Aggregate reproduction across all days for each Plate–Treatment combination
F3LRS <- aggregate(
  Reproduction ~ Plate + Offspring.Treatment + Parental.Treatment + Treatment,
  data = F3,
  FUN  = sum
)

# Model:
#   LRS (total offspring) ~ Treatment + (1 | Plate)
#   with COM-Poisson family (to allow over-/under-dispersion).
F3_LRS <- glmmTMB(
  Reproduction ~ Treatment + (1 | Plate),
  data      = F3LRS,
  family    = compois(),    # added () for clarity/consistency
  na.action = na.fail
)

# Diagnostics
simulationoutput <- simulateResiduals(fittedModel = F3_LRS, plot = TRUE)
testDispersion(simulationoutput)
testZeroInflation(simulationoutput)

performance::check_predictions(F3_LRS)


########################
# 7. Survival analysis #
########################

# Survival data (originally "C:/.../Mol/F3.csv")
# Suggestion: place this file in /data/ as "F3_survival.csv"
LSF3 <- read.csv(here("data", "F3_survival.csv"))

## EH analysis without matricides censored

# Order by age
LSF3 <- LSF3 %>%
  arrange(Age)

# Wide form: one row per worm, one column per age
EH_LSF3 <- LSF3 %>%
  pivot_wider(
    names_from  = Age,
    values_from = Event
  )

range(LSF3$Age)

# Back to long format: one row per worm-age
EH_LSF3 <- EH_LSF3 %>%
  pivot_longer(
    cols      = "1":"32",
    names_to  = "Age",
    values_to = "Event"
  )

# Keep worms that have at least one event (1) or NA
EH_LSF3 <- EH_LSF3 %>%
  group_by(Worm) %>%
  filter(any(Event == 1 | is.na(Event)))

# Helper: replace NA with 0 before the first 1 for each worm
replace_na_before_first_1 <- function(x) {
  # Find first index of 1
  first_1_index <- which(x == 1)[1]
  
  if (is.na(first_1_index)) {
    # If no 1 is found, return x as-is
    return(x)
  }
  
  # Replace NAs before first 1 with 0
  x[seq_len(first_1_index - 1)][is.na(x[seq_len(first_1_index - 1)])] <- 0
  return(x)
}

# Apply per worm
EH_LSF3 <- EH_LSF3 %>%
  group_by(Worm) %>%
  mutate(status = replace_na_before_first_1(Event))

# Drop original Event column
EH_LSF3$Event <- NULL

# Remove rows with remaining NAs
EH_LS3F3 <- na.omit(EH_LSF3)

# Ensure Treatment and (if present) Parental.Treatment are factors
EH_LS3F3 <- EH_LS3F3 %>%
  mutate(
    Treatment          = as.factor(Treatment),
    Parental.Treatment = as.factor(Parental.Treatment)
  )

# Ensure Age is numeric
EH_LS3F3 <- EH_LS3F3 %>%
  mutate(Age = as.integer(Age))

# Fill in missing days (Age) for each worm:
#   - complete sequence of ages
#   - status defaults to 0 for missing rows
#   - fill Treatment, Plate, Cause down/up
EH_LS3F3 <- EH_LS3F3 %>%
  group_by(Worm) %>%
  complete(Age = full_seq(Age, 1), fill = list(status = 0)) %>%
  fill(Treatment, Plate, Cause, Parental.Treatment, .direction = "downup") %>%
  ungroup()

# Survival GLMM:
#   status (alive/dead) ~ Parental.Treatment * Treatment + Age + (1 | Plate / Worm)
#   with binomial family and zero-inflation by Age.
F3_surv <- glmmTMB(
  status ~ Parental.Treatment * Treatment + Age + (1 | Plate/Worm),
  family    = binomial,
  ziformula = ~ Age,
  data      = EH_LS3F3
)

# Diagnostics
simulationoutput <- simulateResiduals(fittedModel = F3_surv, plot = TRUE)
testDispersion(simulationoutput)
testZeroInflation(simulationoutput)

performance::check_predictions(F3_surv)

################################################################################
# End of script
################################################################################


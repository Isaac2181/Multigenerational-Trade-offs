################################################################################
# Script: F1.R
# Author: Isaac Harris
#
# Purpose:
#   Analyse age-specific reproduction, rate-sensitive fitness (Leslie matrices),
#   lifetime reproductive success (LRS), and survival for F1 worms under
#   different treatments and parental treatments.
#
# Expected project structure (example):
#   project/
#    ├ data/
#    │   ├ F1.csv             # reproduction data
#    │   └ F1_survival.csv    # survival data (was in .../Mol/F1.csv)
#    └ scripts/
#        └ F1_reproduction_survival_analysis.R
################################################################################

#############################
# 1. Load required packages #
#############################

library(tidyverse)   # dplyr, tidyr, ggplot2, etc.
library(DHARMa)      # residual diagnostics
library(glmmTMB)     # GLMMs
library(car)         # Anova()
library(splines)     # ns()
library(emmeans)     # marginal means
library(here)        # reproducible file paths

########################################
# 2. Read and tidy reproduction data   #
########################################

# Reproduction data (originally "F1.csv" in L1arrest/data)
F1 <- read.csv(here("data", "F1.csv"))

# Reshape from wide (one column per day) to long:
#   - 'Day' column: day of reproduction (D1, D2, ...)
#   - 'Reproduction' column: number of offspring
F1 <- pivot_longer(
  F1,
  cols      = starts_with("D"),
  names_to  = "Day",
  values_to = "Reproduction"
)

# Remove 'D' prefix from 'Day' and convert to factor
F1$Day <- as.factor(gsub("D", "", F1$Day))

# Ensure treatment variables are factors
F1$Parental.Treatment  <- as.factor(F1$Parental.Treatment)
F1$Offspring.Treatment <- as.factor(F1$Offspring.Treatment)
F1$Treatment           <- as.factor(F1$Treatment)

# Remove any plates that have at least one NA in Reproduction
F1 <- F1 %>%
  group_by(Plate) %>%
  filter(!any(is.na(Reproduction))) %>%
  ungroup()

# Optional observation ID (not used below but kept from original script)
F1$obsID <- seq_len(nrow(F1))

summary(F1)


#########################################################
# 3. Age-specific reproduction model with splines on Day #
#########################################################

# Model:
#   Reproduction ~ Treatment * f(Day) + (1 | Plate)
#   with generalized Poisson family and zero-inflation by Day.
F1_age <- glmmTMB(
  Reproduction ~ Treatment * ns(as.numeric(Day), df = 3) + (1 | Plate),
  family    = genpois(),
  ziformula = ~ Day,
  control   = glmmTMBControl(optimizer = optim,
                             optArgs   = list(method = "BFGS")),
  data      = F1
)

AIC(F1_age)

# Diagnostics
simulationoutput <- simulateResiduals(fittedModel = F1_age, plot = TRUE)
plot(simulationoutput)
testDispersion(simulationoutput)
testZeroInflation(simulationoutput)
performance::check_predictions(F1_age_dis)


##############################################
# 4. Rate-sensitive fitness (Leslie matrix)  #
##############################################

# Goal:
#   For each plate, construct a Leslie matrix using age-specific reproduction
#   and estimate rate-sensitive fitness as the dominant eigenvalue (λ).

plates <- unique(F1$Plate)

# Vector to store plate-specific dominant eigenvalues
dominant_eigenvalues <- numeric(length(plates))

# Loop over plates
for (p in seq_along(plates)) {
  
  plate_id <- plates[p]
  
  # Number of age classes = number of unique days
  n_age_classes <- length(unique(F1$Day))
  Leslie <- matrix(0, nrow = n_age_classes, ncol = n_age_classes)
  
  # First row: reproduction for this plate
  Leslie[1, ] <- F1$Reproduction[F1$Plate == plate_id]
  
  # Subdiagonal: survival from age x to x+1 (assumed 1)
  diag(Leslie[-1, ]) <- 1
  
  # Eigen decomposition and dominant eigenvalue
  eigen_res            <- eigen(Leslie)
  dominant_eigenvalue  <- max(Re(eigen_res$values))
  dominant_eigenvalues[p] <- dominant_eigenvalue
}

# Convert eigenvalues to a data frame
fitness <- as.data.frame(dominant_eigenvalues)

# Unique Plate–Treatment combination data
unique_combinations <- unique(F1[c("Plate", "Parental.Treatment",
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
fitness$Plate             <- udata$Plate
fitness$ParentalTreatment <- udata$ParentalTreatment
fitness$OffspringTreatment<- udata$OffspringTreatment
fitness$Treatment         <- udata$Treatment


########################################
# 5. Model for rate-sensitive fitness  #
########################################

# Model:
#   dominant_eigenvalues ~ Treatment + (1 | Plate)
#   with Gamma errors and log link.
F1_fit <- glmmTMB(
  dominant_eigenvalues ~ Treatment + (1 | Plate),
  data   = fitness,
  family = Gamma(link = "log")
)

# Diagnostics
simulationoutput <- simulateResiduals(fittedModel = F1_fit, plot = TRUE)
testDispersion(simulationoutput)
testZeroInflation(simulationoutput)

performance::check_predictions(F1_fit)

# ANOVA and marginal means
Anova(F1_fit)
emmeans(F1_fit, pairwise ~ Treatment)


###############################################
# 6. Lifetime reproductive success (LRS)      #
###############################################

# Aggregate reproduction across all days for each Plate–Treatment combination
F1LRS <- aggregate(
  Reproduction ~ Plate + Offspring.Treatment + Parental.Treatment + Treatment,
  data = F1,
  FUN  = sum
)

# Model:
#   LRS (total offspring) ~ Treatment + (1 | Plate)
#   COM-Poisson family allows for over-/under-dispersion.
F1_LRS <- glmmTMB(
  Reproduction ~ Treatment + (1 | Plate),
  data        = F1LRS,
  family      = compois(),                 # <- added () for clarity
  dispformula = ~ Treatment + (1 | Plate),
  na.action   = na.fail
)

# Diagnostics
simulationoutput <- simulateResiduals(fittedModel = F1_LRS, plot = TRUE)
testDispersion(simulationoutput)
testZeroInflation(simulationoutput)

performance::check_predictions(F1_LRS)

summary(F1_LRS)
emmeans(F1_LRS, pairwise ~ Treatment)


########################
# 7. Survival analysis #
########################

# Survival data (originally "C:/.../Mol/F1.csv")
# Suggestion: put this file in /data/ as, e.g., "F1_survival.csv"
LSF1 <- read.csv(here("data", "F1_survival.csv"))

## EH analysis without matricides censored

# Order by age
LSF1 <- LSF1 %>%
  arrange(Age)

# Put data in wide form: one row per worm, one column per age
EH_LSF1 <- LSF1 %>%
  pivot_wider(
    names_from  = Age,
    values_from = Event
  )

# (Optional) remove unnecessary columns here if needed
# EH_LSF1 <- EH_LSF1 %>%
#   select(-Death_date, -D1)

range(LSF1$Age)

# Back to long format: one row per worm-age combination
EH_LSF1 <- EH_LSF1 %>%
  pivot_longer(
    cols      = "2":"31",
    names_to  = "Age",
    values_to = "Event"
  )

# Keep worms that have at least one event (1) or NA
EH_LSF1 <- EH_LSF1 %>%
  group_by(Worm) %>%
  filter(any(Event == 1 | is.na(Event)))

# Helper: replace NA with 0 before the first 1 for each worm
replace_na_before_first_1 <- function(x) {
  first_event_row <- which(x == 1)[1]  # get first index of 1
  
  if (is.na(first_event_row)) {
    return(x)  # no 1 found, return as-is
  }
  
  x[seq_len(first_event_row - 1)][is.na(x[seq_len(first_event_row - 1)])] <- 0
  return(x)
}


# Apply per worm
EH_LSF1 <- EH_LSF1 %>%
  group_by(Worm) %>%
  mutate(status = replace_na_before_first_1(Event))

# Drop original Event column
EH_LSF1$Event <- NULL

# Remove rows with any remaining NAs
EH_LS3F1 <- na.omit(EH_LSF1)

# Ensure Treatment and Parental.Treatment are factors
EH_LS3F1 <- EH_LS3F1 %>%
  mutate(
    Treatment          = as.factor(Treatment),
    Parental.Treatment = as.factor(Parental.Treatment)
  )

# Ensure Age is numeric
EH_LS3F1 <- EH_LS3F1 %>%
  mutate(Age = as.integer(Age))

# Fill in missing days (Age) for each worm:
#   - complete sequence of ages
#   - status defaults to 0 for missing rows
#   - fill Treatment, Parental.Treatment, Plate, Cause down/up
EH_LS3F1 <- EH_LS3F1 %>%
  group_by(Worm) %>%
  complete(Age = full_seq(Age, 1), fill = list(status = 0)) %>%
  fill(Treatment, Parental.Treatment, Plate, Cause, .direction = "downup") %>%
  ungroup()

# Survival GLMM:
#   status (alive/dead) ~ Treatment * Parental.Treatment + Age + (1 | Plate)
#   with binomial family and zero-inflation by Age.
F1_Sur <- glmmTMB(
  status ~ Treatment * Parental.Treatment + Age + (1 | Plate),
  family    = binomial,
  control   = glmmTMBControl(optimizer = optim,
                             optArgs   = list(method = "BFGS")),
  ziformula = ~ Age,
  data      = EH_LS3F1
)

# Diagnostics
simulationoutput <- simulateResiduals(fittedModel = F1_Sur, plot = TRUE)
testDispersion(simulationoutput)
testZeroInflation(simulationoutput)

performance::check_predictions(F1_Sur)

################################################################################
# End of script
################################################################################



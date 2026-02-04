
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

# Reproduction data
transgen <- read.csv(here("data", "transgen.csv"))

# Reshape from wide (one column per day) to long:
#   - 'Day' column: day of reproduction (D1, D2, ...)
#   - 'Reproduction' column: number of offspring
transgen <- pivot_longer(
  transgen,
  cols      = starts_with("D"),
  names_to  = "Day",
  values_to = "Reproduction"
)

# Remove 'D' prefix from 'Day' and convert to factor
transgen$Day <- as.factor(gsub("D", "", transgen$Day))

# Ensure treatment variables are factors
transgen$Parental.Treatment  <- as.factor(transgen$Parental.Treatment)
transgen$Offspring.Treatment <- as.factor(transgen$Offspring.Treatment)
transgen$Treatment           <- as.factor(transgen$Treatment)

# Remove any plates that have at least one NA in Reproduction
transgen <- transgen %>%
  group_by(Plate) %>%
  filter(!any(is.na(Reproduction))) %>%
  ungroup()

# Optional observation ID
#transgen$obsID <- seq_len(nrow(transgen))


###############################################
# 6. Lifetime reproductive success (LRS)      #
###############################################

# Aggregate reproduction across all days for each Plate–Treatment combination
transgenLRS <- aggregate(
  Reproduction ~ Plate + Offspring.Treatment + Parental.Treatment + Treatment + Generation,
  data = transgen,
  FUN  = sum
)

transgen_LRS <- glmmTMB(
  Reproduction ~ Parental.Treatment*Offspring.Treatment*Generation + (1 | Plate),
  dispformula = ~ (1 | Plate),
  data      = transgenLRS,
  family    = nbinom2(),
  na.action = na.fail
)

# Diagnostics
simulationoutput <- simulateResiduals(fittedModel = transgen_LRS, plot = TRUE)
testDispersion(simulationoutput)
testZeroInflation(simulationoutput)

performance::check_predictions(transgen_LRS)

summary(transgen_LRS)

emm_LRS <- emmeans(transgen_LRS, ~ Parental.Treatment * Offspring.Treatment * Generation)

contrast(emm_LRS, interaction = c("pairwise", "pairwise", "pairwise"))

emmip(emm_LRS, Generation ~ Offspring.Treatment | Parental.Treatment, 
      type = "response",
      CIs = TRUE) + 
  labs(y = "Lifetime Reproductive Success", 
       title = "Intergenerational vs. Transgenerational Effects") +
  theme_bw()

##############################################
# 4. Rate-sensitive fitness (Leslie matrix)  #
##############################################

# Goal:
#   For each plate, construct a Leslie matrix using age-specific reproduction
#   and estimate rate-sensitive fitness as the dominant eigenvalue (λ).

plates <- unique(transgen$Plate)

# Vector to store plate-specific dominant eigenvalues
dominant_eigenvalues <- numeric(length(plates))

# Loop over plates
for (p in seq_along(plates)) {
  
  plate_id <- plates[p]
  
  # Number of age classes = number of unique days
  n_age_classes <- length(unique(transgen$Day))
  Leslie <- matrix(0, nrow = n_age_classes, ncol = n_age_classes)
  
  # First row: reproduction for this plate
  Leslie[1, ] <- transgen$Reproduction[transgen$Plate == plate_id]
  
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
unique_combinations <- unique(transgen[c("Plate", "Parental.Treatment",
                                   "Offspring.Treatment", "Treatment", "Generation")])

udata <- data.frame(
  Plate             = unique_combinations$Plate,
  ParentalTreatment = unique_combinations$Parental.Treatment,
  OffspringTreatment= unique_combinations$Offspring.Treatment,
  Generation        = unique_combinations$Generation,
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
fitness$Generation         <- udata$Generation


########################################
# 5. Model for rate-sensitive fitness  #
########################################

# Model:
#   dominant_eigenvalues ~ Treatment + (1 | Plate)
#   with Gamma errors and log link.
transgen_fit <- glmmTMB(
  dominant_eigenvalues ~ ParentalTreatment*OffspringTreatment*Generation + (1 | Plate),
  data        = fitness,
  family      = Gamma(link = "log"),
  control     = glmmTMBControl(
    optimizer = optim,
    optArgs   = list(method = "BFGS")
  ),
  dispformula = ~ (1 | Plate)
)

# Diagnostics
simulationoutput <- simulateResiduals(fittedModel = transgen_fit, plot = TRUE)
testDispersion(simulationoutput)
testZeroInflation(simulationoutput)

performance::check_predictions(transgen_fit)

summary(transgen_fit)

emm_fit <- emmeans(transgen_fit, ~ ParentalTreatment * OffspringTreatment * Generation)

contrast(emm_fit, interaction = c("pairwise", "pairwise", "pairwise"))

emmip(emm_fit, Generation ~ OffspringTreatment | ParentalTreatment, 
      type = "response",
      CIs = TRUE) + 
  labs(y = "Fitness (Dominant Eigenvalues)", 
       title = "Intergenerational vs. Transgenerational Effects") +
  theme_bw()

################################################################################
# run_all.R
#
# Master pipeline runner for generating model outputs and plots used in "Evolutionary trade-offs between intergenerational and transgenerational fitness effects"
#
# This script:
#   1) Sets up the project root for here::here()
#   2) Runs P0, F1, F3 analyses
#   3) Runs the simulation
#   4) Fits survival models and builds survival plots
#   5) Builds and exports all figures (main + supplementary)
#   6) Exports model output tables (gt -> .docx)
#   7) Saves sessionInfo() for reproducibility
#
# Usage (from project root):
#   source("run_all.R")
################################################################################

#############################
# 0. Basic setup            #
#############################

# Make sure the key packages exist; individual scripts will load their own.
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

# OPTIONAL: pin the project root for 'here' (only needed once per machine)
# This will create a small .here file in the root directory.
getwd()  # should now be your repo root
here::i_am("run_all.R")
here::here()

# Create standard output folders if they don't exist
dirs_to_create <- c(
  here("plots", "Final", "Export"),
  here("model_output")
)

for (d in dirs_to_create) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

#############################
# 1. Core analyses          #
#############################

cat(">>> Running P0 analysis...\n")
source(here("scripts", "P0.R"))

cat(">>> Running F1 analysis...\n")
source(here("scripts", "F1.R"))

cat(">>> Running F3 analysis...\n")
source(here("scripts", "F3.R"))

#############################
# 2. Simulation             #
#############################

cat(">>> Running simulation script...\n")
source(here("scripts", "WormSimulation.R"))

#############################
# 3. plots         #
#############################

cat(">>> Running P0 survival analysis and plots...\n")
source(here("scripts", "P0Surv.R"))

cat(">>> Running F1 survival analysis and plots...\n")
source(here("scripts", "F1Surv.R"))

cat(">>> Running F3 survival analysis and plots...\n")
source(here("scripts", "F3Surv.R"))

cat(">>> Running P0 plots...\n")
source(here("scripts", "P0 Plots.R"))

cat(">>> Running F1 plots...\n")
source(here("scripts", "F1Plots.R"))

cat(">>> Running F3 plots...\n")
source(here("scripts", "F3Plots.R"))
#############################
# 4. Figures (main & supp)  #
#############################

# This script should:
#   - construct Fig2, Fig3, Fig4, etc. from objects defined in the analysis
#   - export them as SVG to plots/Final/Export/
cat(">>> Building and exporting figures...\n")
source(here("scripts", "GraphExport.R"))

#############################
# 5. Model output tables    #
#############################

# This script should:
#   - use make_mixed_model_table() and make_survival_table()
#   - export .docx tables to model_output/
cat(">>> Exporting model output tables...\n")
source(here("scripts", "ModelOutputs.R"))

#############################
# 6. Session info           #
#############################

cat(">>> Saving sessionInfo() to model_output/sessionInfo.txt ...\n")
sink(here("model_output", "sessionInfo.txt"))
sessionInfo()
sink()

cat(">>> All done! ✅\n")
################################################################################
# End of run_all.R                                                             #
################################################################################

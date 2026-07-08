################################################################################
### Decision for dialysis versus conservative care
### Run all
################################################################################

# remove history
rm(list=ls(all.names=TRUE))

# prepare data
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/00 - Data preparation.R")

# apply eligibility criteria
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/01 - Apply eligibility criteria.R")

# create covariates and outcomes
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/02 - Covariate and outcome derivation.R")

# fit transportability and propensity score model, and create weights
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/03 - Compute weights.R")

# create descriptive statistics
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/04 - Descriptives.R")

# calculate overall average treatment effects
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/05 - ATE.R")

# Effect modification analysis
# Warning: The `size` argument of `element_line()` is deprecated comes from ggthemes
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/06 - Continuous HTE.R")

# Example patient in methods
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/07 - Example Supplemental Methods.R")

# Sensitivity analysis for positivity
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/08 - SA positivity.R")

# Time-to-dialysis
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/09 - Competing risk time-to-dialysis.R")

# Sensitivity analysis for unmeasured confounding
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/10 - SA unmeasured confounding.R")

# Illustation relative and absolute HTE
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/11 - Supplemental Figure HTE.R")

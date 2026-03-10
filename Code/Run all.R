################################################################################
### Decision for dialysis versus conservative care
### Run all
################################################################################

# remove history
rm(list=ls(all.names=TRUE))

# prepare data
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/0 - Data preparation.R")

# apply eligibility criteria
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/1 - Apply eligibility criteria.R")

# create covariates and outcomes
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/2 - Covariate and outcome derivation.R")

# fit transportability and propensity score model, and create weights
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/3 - Compute weights.R")

# create descriptive statistics
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/4 - Descriptives.R")

# calculate overall average treatment effects
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/5 - ATE.R")

# Effect modification analysis
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/6 - Continuous HTE.R")

# Sensitivity analysis for positivity
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/7 - SA positivity.R")

# Sensitivity analysis for positivity
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/8 - Competing risk time-to-dialysis.R")

# Example patient in methods
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/9 - Example Supplemental Methods.R")

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
# TODO: add CRP to weighting model?
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/3 - Compute weights.R")

# create descriptive statistics
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/4 - Descriptives.R")

# calculate overall average treatment effects
# TODO: check censoring tick marks
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/5 - ATE.R")

# Effect modification analysis
# TODO: The `size` argument of `element_line()` is deprecated 
# TODO:  In coxph.fit(X, Y, istrat, offset, init, control, weights = weights,  : Loglik converged before variable  3 ; coefficient may be infinite. 
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/6 - Continuous HTE.R")

# Example patient in methods
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/7 - Example Supplemental Methods.R")

# Sensitivity analysis for positivity
# TODO: `geom_errorbarh()` was deprecated in ggplot2 4.0.0. Please use the `orientation` argument of `geom_errorbar()` instead.
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/8 - SA positivity.R")

# Time-to-dialysis
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/9 - Competing risk time-to-dialysis.R")

# Sensitivity analysis for unmeasured confounding
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/10 - SA unmeasured confounding.R")

# Illustate relative and absolute HTE
source("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Code/11 - Supplemental Figure HTE.R")

################################################################################
### External 2-year binary risk model
################################################################################
# extract predictors
predictors <- baseline[, .(A = trt, age, egfr2021, cancer, dm, ihd, vhd, pvd, female, albumin)]

# median CRP in development is 76 for dialysis and 52 for CC
predictors$crp <- ifelse(predictors$A == 1, 76, 52)

# use external model to make predictions
baseline$prob_ext_D <- prob_D_and_CC(X = predictors, arm = "D")
baseline$prob_ext_CC <- prob_D_and_CC(X = predictors, arm = "CC")

# define 5-year time-to-event outcome
Surv_5y <- survival::Surv(baseline$time2event_death_5y, baseline$event_death_5y)

# define 2-year time-to-event outcome
Surv_2y <- Surv_5y
Surv_2y[Surv_5y[, 1] > 2 * 365, 1] <- 2 * 365 # time horizon set to 2 years
Surv_2y[Surv_5y[, 1] > 2 * 365, 2] <- 0     # censor if event after 2 years
baseline$time2event_death_2y <- Surv_2y[, 1]
baseline$death_2y <- Surv_2y[, 2]

# validate external model
# rms::val.prob(p=baseline$prob_ext_D, y=baseline$death_2y)
# rms::val.prob(p=baseline$prob_ext_CC, y=baseline$death_2y)
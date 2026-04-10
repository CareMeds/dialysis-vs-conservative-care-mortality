################################################################################
### Decision for dialysis versus conservative care
### PART 5 - Average treatment effect
################################################################################

# remove history
rm(list = ls(all.names = TRUE))
knitr::opts_knit$set(root.dir = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/")

# set directory
setwd(
  "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/"
)
results_path <- "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Results/"

# load libraries
library(data.table)
library(patchwork) # combine figures
library(foreach)   # parallel computation
library(doRNG)     # handle parallel seeds
set.seed(1)        # set seed for parallel backend

# load functions
source("Code/utils/plots.R")
source("Code/utils/tables.R")
source("Code/utils/weighting.R")
source("Code/utils/data_manipulation.R")
source("Code/utils/compute_absolute_relative_risks.R")

################################################################################
### Load data ##################################################################
################################################################################
load("Data/cohort_with_weights.Rdata")

################################################################################
### Median follow-up time ######################################################
################################################################################
# Create reversed KM model
print(quantile(
  prodlim::prodlim(
    prodlim::Hist(time2event_death_inf / 30.5, 1 - event_death_inf) ~ 1,
    data = baseline,
    reverse = TRUE
  )
))
# Median follow-up time is 28.20 (11.41;55.31)
print(quantile(
  prodlim::prodlim(
    prodlim::Hist(time2event_death_inf / 30.5, 1 - event_death_inf) ~ 1,
    data = baseline[trt == 1],
    reverse = TRUE
  )
))
# Dialysis: 37.28 (16.03;66.00)
print(quantile(
  prodlim::prodlim(
    prodlim::Hist(time2event_death_inf / 30.5, 1 - event_death_inf) ~ 1,
    data = baseline[trt == 0],
    reverse = TRUE
  )
))
# Conservative care: 14.56 (6.30;30.00)

################################################################################
### Combine absolute and relative results in one Table and make KM plot
################################################################################

# Define labels
outcome_label <- "2-year all-cause mortality"
outcome_var <- "event_death_2y"
time2outcome_var <- "time2event_death_2y"
competing_events_var <- "none"

# set unit for dRMST
unit <- "months"

# Loop through each outcome variable
n_bootstraps <- 1000
horizon <- 2 * 365
# Loop through each weighting method
for (w_meth in w_meths) {
  # Set weights: use 1 for unweighted, otherwise use specified weights
  if (w_meth == "unweighted") {
    weights_meth <- rep(1, nrow(baseline))
  } else {
    weights_meth <- baseline[[paste0("sw_", w_meth)]]
  }
  
  # compute estimates
  out_est <- compute_estimates_with_CI(
    data = baseline,
    unit = unit,
    horizon = horizon,
    elig_cohort = elig_cohort,
    model_PS = model_PS,
    model_S = model_S,
    event_var = outcome_var,
    competing_event_var = competing_events_var,
    time2event_var = time2outcome_var,
    trt_var = trt_var,
    w_meth = w_meth,
    trim_meth = "no_trimming",
    catvar = catvar,
    contvar = contvar,
    n_bootstraps = n_bootstraps,
    bootstrap_seed = 1
  )
  assign(paste0("out_est_", outcome_var, "_", w_meth), out_est)
  
  # Kaplan-Meier curves
  out_KM <- create_KM_plot(
    data = baseline,
    trt_var = trt_var,
    event_var = outcome_var,
    time2event_var = time2outcome_var,
    w_meth = w_meth,
    out_est = out_est,
    horizon = horizon,
    unit = unit,
    manual_colors = manual_colors,
    trt_labels = c(control_label, treatment_label)
  )
  assign(paste0("out_KM_", outcome_var, "_", w_meth), out_KM)
}

# for sensitivity analysis for unmeasured confounding
surv_out <- out_KM_event_death_2y_IPTW$KM_plot$data
cat(
  " Risk of all-cause mortality among patients who chose dialysis at 6 months         :",
  (1 - as.numeric(surv_out[surv_out$time == 182 &
                             surv_out$strata == 1, "surv"])) * 100,
  "\n",
  "Risk of all-cause mortality among patients who chose conservative care at 6 months:",
  (1 - as.numeric(surv_out[surv_out$time == 182 &
                             surv_out$strata == 0, "surv"])) * 100,
  "\n"
)

# Create Figures
ggplot2::ggsave(
  plot = ggpubr::ggarrange(
    out_KM_event_death_2y_IPTW$KM_plot,
      # ggplot2::ggtitle(outcome_label) +
      # ggplot2::theme(axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = -10))),
    out_KM_event_death_2y_unweighted$KM_table +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA)),
    ncol = 1,
    heights = c(3, 1),
    align = "v"
  ),
  filename = paste0(results_path, "Main/Figure_1.png"),
  width = 6,
  height = 6,
  dpi = 600
)

# Supplemental Table for all outcomes and horizons
for (w_meth in w_meths) {
  # Initialize the results table with one row for the outcome name
  results_df <- data.frame(Conservative_care = outcome_label, Dialysis =
                             "")
  rownames(results_df) <- "Outcome"
  
  # Add a header row for the weighting method
  results_df[ifelse(w_meth == "unweighted",
                    "Unweighted",
                    paste("Weighting", w_meth)), ] <- rep("", 2)
  
  # obtain sample size
  results_df["Sample size", ] <- c(sum(baseline$trt == 0), sum(baseline$trt == 1))
  
  # Add raw number of events to the table
  results_df["Number of events", ] <- c(sum(baseline[[outcome_var]] == 1 &
                                                baseline$trt == 0),
                                          sum(baseline[[outcome_var]] == 1 &
                                                baseline$trt == 1))
  
  # obtain correct estimates
  out_est <- eval(parse(text = paste0("out_est_", outcome_var, "_", w_meth)))
  
  # fill absolute risks
  results_df[paste("Risk, % (95% CI)", w_meth), ] <-
    c(
      fmt_ci(out_est$R0 * 100, out_est$R0_lower * 100, out_est$R0_upper * 100),
      fmt_ci(out_est$R1 * 100, out_est$R1_lower * 100, out_est$R1_upper * 100)
    )
  
  # fill risk difference
  results_df[paste("Risk difference, % (95% CI)", w_meth), ] <-
    c("Reference",
      fmt_ci(out_est$RD * 100, out_est$RD_lower * 100, out_est$RD_upper * 100))
  
  # fill risk ratio
  results_df[paste("Risk ratio (95% CI)", w_meth), ] <-
    c("Reference",
      fmt_ci(out_est$RR, out_est$RR_lower, out_est$RR_upper, 2))
  
  # fill RMST
  results_df[paste("RMST, ", unit, " (95% CI)", w_meth), ] <-
    c(
      fmt_ci(out_est$RMST0, out_est$RMST0_lower, out_est$RMST0_upper),
      fmt_ci(out_est$RMST1, out_est$RMST1_lower, out_est$RMST1_upper)
    )
  
  # fill RMST difference
  results_df[paste("\u0394RMST, ", unit, " (95% CI)", w_meth), ] <-
    c("Reference",
      fmt_ci(out_est$dRMST, out_est$dRMST_lower, out_est$dRMST_upper))
  
  # fill hazard ratio
  results_df[paste("HR (95% CI)", w_meth), ] <-
    c("Reference",
      fmt_ci(out_est$HR, out_est$HR_lower, out_est$HR_upper, 2))
  
  assign(paste0("results_df_", w_meth), results_df)
}

# Save the combined results table for all outcomes to Excel
openxlsx::write.xlsx(
  results_df_IPTW,
  rowNames = TRUE,
  file = paste0(results_path, "Main/Table_2.xlsx")
)
openxlsx::write.xlsx(
  results_df_IPSW_IPTW,
  rowNames = TRUE,
  file = paste0(results_path, "Supplemental/Table_S9_IPSW.xlsx")
)
openxlsx::write.xlsx(
  rbind(results_df_IPTW, results_df_SMR_ATT[-c(1, 2), ], results_df_SMR_ATU[-c(1, 2), ]),
  rowNames = TRUE,
  file = paste0(results_path, "Supplemental/Table_S7_ATT_ATU.xlsx")
)

# save cohort
save(
  id_name,
  listvar,
  listvar_main,
  catvar,
  contvar,
  non_normal_vars,
  outcome_var,
  time2outcome_var,
  competing_events_var,
  treatment_label,
  control_label,
  baseline,
  elig_cohort,
  model_PS,
  coef_PS_overall,
  elig_cohort,
  w_meths,
  trt_var,
  horizon,
  unit,
  manual_colors,
  n_bootstraps,
  file = file.path("Data/cohort_with_models.Rdata")
)

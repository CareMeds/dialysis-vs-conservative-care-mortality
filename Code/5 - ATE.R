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
library(patchwork) # combine figures
library(data.table)
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
outcome_labels <- c("1-year all-cause mortality", "2-year all-cause mortality")
outcome_vars <- paste0("event_death_", 1:2, "y")
time2outcome_vars <- paste0("time2event_death_", 1:2, "y")
competing_events_vars <- rep("none", 2)

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
# Median follow-up time is 28.20 (11.48;55.31)
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
# # Conservative care: 14.62 (6.30;30.03)

################################################################################
### Combine absolute and relative results in one Table and make KM plot
################################################################################

# Define treatment labels
treatment <- c(control_label, treatment_label)

# set unit for dRMST
unit <- "months"

# Loop through each outcome variable
n_bootstraps <- 1000
# for (outcome_var in 2){
for (outcome_var in 1:length(outcome_vars)) {
  # only for all-cause mortality predict for 1 and 2 years
  if (outcome_var == 1) {
    horizon <- 365
    w_meths_ATE <- w_meths[c(1, 2, 6, 7)]
  } else{
    horizon <- 2 * 365
    w_meths_ATE <- w_meths
  }
  
  # Loop through each weighting method
  for (w_meth in w_meths_ATE) {
    # for (w_meth in "IPTW") {
    cat("Outcome", outcome_var, "years, weighting:", w_meth, "\n")
    # Set weights: use 1 for unweighted, otherwise use specified weights
    if (w_meth == "") {
      weights_meth <- rep(1, nrow(baseline))
    } else {
      weights_meth <- baseline[[paste0("sw_", w_meth)]]
    }
    
    # Kaplan-Meier curves
    out_KM <- create_KM_plot(
      data = baseline,
      elig_cohort = elig_cohort,
      horizon = horizon,
      unit = unit,
      model_PS = model_PS,
      model_S = model_S,
      w_meth = w_meth,
      weights_meth = weights_meth,
      catvar = catvar,
      contvar = contvar,
      event_var = outcome_vars[outcome_var],
      time2event_var = time2outcome_vars[outcome_var],
      trt_var = "trt",
      competing_event_var = competing_events_vars[outcome_var],
      n_bootstraps = n_bootstraps,
      bootstrap_seed = 1,
      plotTitle = outcome_labels[outcome_var],
      plotColors = manual_colors[1:2],
      annotate_figure = FALSE
    )
    assign(paste0("out_KM_", outcome_vars[outcome_var], "_", w_meth), out_KM)
  }
}

# for sensitivity analysis for unmeasured confounding
surv_out <- out_KM_event_death_2y_IPTW$KM_plot$data
cat(
  "Risk of all-cause mortality among patients who chose dialysis at 6 months:",
  (1 - surv_out[surv_out$time == 183 &
                  surv_out$trt == 1, "surv"]) * 100,
  "\n",
  "Risk of all-cause mortality among patients who chose conservative care at 6 months:",
  (1 - surv_out[surv_out$time == 183 &
                  surv_out$trt == 0, "surv"]) * 100,
  "\n"
)

# Create Figures
ggplot2::ggsave(
  plot = ggpubr::ggarrange(
    out_KM_event_death_2y_IPTW$KM_plot +
      ggplot2::ggtitle(outcome_labels[2]) +    # add title
      ggplot2::theme(axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(r = -10)    # make y-axis label closer
      )),
    out_KM_event_death_2y_IPTW$KM_table,
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
for (outcome_var in 1:length(outcome_vars)) {
  # Initialize the results table with one row for the outcome name
  rows <- c("Outcome")
  results_df <- data.frame(matrix(nrow = length(rows), ncol = length(treatment)))
  rownames(results_df) <- rows
  colnames(results_df) <- treatment
  
  # Add the outcome name to the first row
  results_df["Outcome", ] <- c(outcome_labels[outcome_var], "")
  
  # Add sample size and percentage for each treatment group
  results_df["Sample size", ] <- c(sum(baseline$trt == 0), sum(baseline$trt == 1))
  
  # Loop through each weighting method
  if (outcome_var == 1) {
    w_meths_ATE <- w_meths[c(1, 2, 6, 7)]
  } else{
    w_meths_ATE <- w_meths
  }
  
  for (w_meth in w_meths_ATE) {
    # obtain correct IR and KM
    out_KM <- eval(parse(text = paste0("out_KM_", outcome_vars[outcome_var], "_", w_meth)))
    
    # Add a header row for the weighting method
    results_df[ifelse(w_meth == "", "Unweighted", paste("Weighting", w_meth)), ] <- rep("", 2)
    
    # Add number of events and incidence rates to the table
    if (w_meth == "") {
      results_df["Number of events", ] <- c(sum(baseline[[outcome_vars[outcome_var]]] ==
                                                  1 & baseline$trt == 0),
                                            sum(baseline[[outcome_vars[outcome_var]]] ==
                                                  1 & baseline$trt == 1))
    }
    
    # fill absolute risks
    results_df[paste("Risk, % (95% CI)", w_meth), ] <-
      c(
        fmt_ci(out_KM$R0 * 100, out_KM$R0_lower * 100, out_KM$R0_upper * 100),
        fmt_ci(out_KM$R1 * 100, out_KM$R1_lower * 100, out_KM$R1_upper * 100)
      )
    
    # fill risk difference
    results_df[paste("Risk difference, % (95% CI)", w_meth), ] <-
      c("Reference",
        fmt_ci(out_KM$RD * 100, out_KM$RD_lower * 100, out_KM$RD_upper * 100))
    
    # fill risk ratio
    results_df[paste("Risk ratio (95% CI)", w_meth), ] <-
      c("Reference",
        fmt_ci(out_KM$RR, out_KM$RR_lower, out_KM$RR_upper, 2))
    
    # fill RMST
    results_df[paste("RMST, ", unit, " (95% CI)", w_meth), ] <-
      c(
        fmt_ci(out_KM$RMST0, out_KM$RMST0_lower, out_KM$RMST0_upper),
        fmt_ci(out_KM$RMST1, out_KM$RMST1_lower, out_KM$RMST1_upper)
      )
    
    # fill RMST difference
    results_df[paste("\u0394RMST, ", unit, " (95% CI)", w_meth), ] <-
      c("Reference",
        fmt_ci(out_KM$dRMST, out_KM$dRMST_lower, out_KM$dRMST_upper))
    
    # fill hazard ratio
    results_df[paste("HR (95% CI)", w_meth), ] <-
      c("Reference",
        fmt_ci(out_KM$HR, out_KM$HR_lower, out_KM$HR_upper, 2))
  }
  
  # Store results for current outcome
  assign(paste0("results_df_", outcome_var), results_df)
}

# Save the combined results table for all outcomes to Excel
openxlsx::write.xlsx(
  results_df_2[c(2, 4, 12:17), ],
  # 2-year IPTW
  rowNames = TRUE,
  file = paste0(results_path, "Main/Table_2.xlsx")
)
openxlsx::write.xlsx(
  results_df_2[c(4, 33:38), ],
  # 2-year IPSW + IPTW
  rowNames = TRUE,
  file = paste0(results_path, "Supplemental/Table_S8_IPSW.xlsx")
)
openxlsx::write.xlsx(
  results_df_2[c(
    2,
    4,
    11:17,
    # 2-year IPTW
    39:46,
    # 2-year ATT
    47:52
  ), ],
  # 2-year ATU
  rowNames = TRUE,
  file = paste0(results_path, "Supplemental/Table_S9_ATT_ATU.xlsx")
)

# save cohort
save(
  id_name,
  listvar,
  catvar,
  contvar,
  non_normal_vars,
  outcome_vars,
  time2outcome_vars,
  competing_events_vars,
  treatment_label,
  control_label,
  baseline,
  elig_cohort,
  model_PS,
  coef_PS_overall,
  elig_cohort,
  w_meths,
  horizon,
  unit,
  manual_colors,
  n_bootstraps,
  file = file.path("Data/cohort_with_models.Rdata")
)
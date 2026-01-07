################################################################################
### Decision for dialysis versus conservative care
### PART 5 - Average treatment effect
################################################################################

# remove history
rm(list = ls(all.names = TRUE))
knitr::opts_knit$set(root.dir = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/")
set.seed(1)

# set directory
setwd("P:/SCREAM2/SCREAM2_Research/Carolien Maas/")
results_path <- "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Results/"

# load libraries
library(patchwork) # combine figures
library(data.table)
library(foreach)

# load functions
source("Code/utils/plots.R")
source("Code/utils/tables.R")
source("Code/utils/weighting.R")
source("Code/utils/data_manipulation.R")
source("Code/utils/outcomes_absolute_risks.R")
source("Code/utils/outcomes_incidence_rates.R")
source("Code/utils/outcomes_hazard_ratios.R")

################################################################################
### Load data ##################################################################
################################################################################
load("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/cohort_with_weights.Rdata")
outcome_labels <- c(
  "1-year all-cause mortality",
  "2-year all-cause mortality"
)
outcome_vars <- paste0("event_death_", 1:2, "y")
time2outcome_vars <- paste0("time2event_death_", 1:2, "y")
competing_events_vars <- rep("none", 2)

################################################################################
### Median follow-up time ######################################################
################################################################################
# Create reversed KM model
print(quantile(prodlim::prodlim(
  prodlim::Hist(time2event_death_inf / 30.5, 1 - event_death_inf) ~ 1,
  data = baseline,
  reverse = TRUE
)))
# Median follow-up time is 29.18 (11.61;56.66)
print(quantile(prodlim::prodlim(
  prodlim::Hist(time2event_death_inf / 30.5, 1 - event_death_inf) ~ 1,
  data = baseline[trt == 1],
  reverse = TRUE
)))
# Dialysis: 37.87 (16.20;67.51)
print(quantile(prodlim::prodlim(
  prodlim::Hist(time2event_death_inf / 30.5, 1 - event_death_inf) ~ 1,
  data = baseline[trt == 0],
  reverse = TRUE
)))
# Conservative care: 14.75 (6.36;30.98)

################################################################################
### Combine absolute and relative results in one Table and make KM plot
################################################################################

# Define treatment labels
treatment <- c(control_label, treatment_label)

# set unit for dRMST
unit <- "months"

# Loop through each outcome variable
n_bootstraps <- 1000
# for (outcome_var in 1:length(outcome_vars)) {
for (outcome_var in 2){
  # only for all-cause mortality predict for 1 and 2 years
  if (outcome_var == 1) {
    horizon <- 365
    w_meths_ATE <- w_meths[c(1, 2, 6, 7)]
  } else{
    horizon <- 2 * 365
    w_meths_ATE <- w_meths # [c(1:3, 6:7)]
  }
  
  # Loop through each weighting method
  # for (w_meth in w_meths_ATE) {
  for (w_meth in "IPTW") {
    cat("Outcome", outcome_var, "years, weighting:", w_meth, "\n")
    # Set weights: use 1 for unweighted, otherwise use specified weights
    if (w_meth == "") {
      weights_meth <- rep(1, nrow(baseline))
    } else {
      weights_meth <- baseline[[paste0("sw_", w_meth)]]
    }
    
    # incidence rate
    IR <- get_incidence_rate_stratified(data = baseline,
                                        id_name = id_name,
                                        trt = baseline$trt,
                                        outcome = baseline[[outcome_vars[outcome_var]]],
                                        time2outcome = baseline[[time2outcome_vars[outcome_var]]],
                                        weights = weights_meth)
    assign(paste0("IR_", outcome_vars[outcome_var], "_", w_meth), IR)
    
    # Kaplan-Meier curves
    out_KM <- create_KM_plot(
      data = baseline,
      elig_cohort = elig_cohort,
      horizon = horizon,
      unit = unit,
      model_PS = model_PS,
      w_meth = w_meth,
      weights_meth = weights_meth,
      catvar = catvar,
      contvar = contvar,
      event_var = outcome_vars[outcome_var],
      time2event_var = time2outcome_vars[outcome_var],
      trt_var = "trt",
      competing_event_var = competing_events_vars[outcome_var],
      n_bootstraps = n_bootstraps, # ifelse(w_meth=="IPSW", 10, n_bootstraps),
      bootstrap_seed = 1,
      plotTitle = paste0(outcome_labels[outcome_var], "\n", plot_title(w_meth))
    )
    assign(paste0("out_KM_", outcome_vars[outcome_var], "_", w_meth), out_KM)
  }
}

# Create Figures
ggplot2::ggsave(
  plot = ggpubr::ggarrange(
    out_KM_event_death_2y_IPTW$KM_plot +
      ggplot2::ggtitle(outcome_labels[2]) +    # add title
      ggplot2::theme(
        axis.title.y = ggplot2::element_text(
          margin = ggplot2::margin(r = -10)    # make y-axis label closer
        )
      ),
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
  results_df["Sample size", ] <- c(sum(baseline$trt == 0), 
                                   sum(baseline$trt == 1))
  
  # Loop through each weighting method
  if (outcome_var == 1) {
    w_meths_ATE <- w_meths[c(1, 2, 6, 7)]
  } else{
    w_meths_ATE <- w_meths # [c(1:3, 6:7)]
  }
  
  for (w_meth in w_meths_ATE) {
    # obtain correct IR and KM
    out_KM <- eval(parse(text = paste0("out_KM_", outcome_vars[outcome_var], "_", w_meth)))
    IR <- eval(parse(text = paste0("IR_", outcome_vars[outcome_var], "_", w_meth)))
    
    # Add a header row for the weighting method
    results_df[ifelse(w_meth == "", "Unweighted",
                      paste("Weighting", w_meth)), ] <- rep("", 2)
    
    # Add number of events and incidence rates to the table
    results_df[paste("Number of events", w_meth), ] <- c(IR$event[1], IR$event[2])
      # c(IR$raw_table["Conservative care", "No. of events"], IR$raw_table["Dialysis", "No. of events"])
    # results_df[paste("IR per 100PY", w_meth), ] <-
    #   c(IR$raw_table["Conservative care", "IR per 100PY (95% CI)"], IR$raw_table["Dialysis", "IR per 100PY (95% CI)"])
    
    # fill absolute risks
    results_df[paste("Risk, % (95% CI)", w_meth), ] <- 
      c(fmt_ci(out_KM$R0 * 100, out_KM$R0_lower * 100, out_KM$R0_upper * 100),
        fmt_ci(out_KM$R1 * 100, out_KM$R1_lower * 100, out_KM$R1_upper * 100))
    
    # fill risk difference
    results_df[paste("Risk difference, % (95% CI)", w_meth), ] <- 
      c("Reference", fmt_ci(out_KM$RD * 100, out_KM$RD_lower * 100, out_KM$RD_upper * 100))
    
    # fill risk ratio
    results_df[paste("Risk ratio (95% CI)", w_meth), ] <- 
      c("Reference", fmt_ci(out_KM$RR, out_KM$RR_lower, out_KM$RR_upper, 2))
    
    # fill RMST
    results_df[paste("RMST, ", unit, " (95% CI)", w_meth), ] <- 
      c(fmt_ci(out_KM$RMST0, out_KM$RMST0_lower, out_KM$RMST0_upper),
        fmt_ci(out_KM$RMST1, out_KM$RMST1_lower, out_KM$RMST1_upper))
    
    # fill RMST difference
    results_df[paste("\u0394RMST, ", unit, " (95% CI)", w_meth), ] <- 
      c("Reference", fmt_ci(out_KM$dRMST, out_KM$dRMST_lower, out_KM$dRMST_upper))
    
    # fill hazard ratio
    results_df[paste("HR (95% CI)", w_meth), ] <- 
      c("Reference", fmt_ci(out_KM$HR, out_KM$HR_lower, out_KM$HR_upper, 2))
  }
  
  # Store results for current outcome
  assign(paste0("results_df_", outcome_var), results_df)
}

# Save the combined results table for all outcomes to Excel
openxlsx::write.xlsx(
  results_df_2[c(2, 4, 13:18), ], # 2-year unweighted and IPTW
  rowNames = TRUE,
  file = paste0(results_path, "Main/Table_2.xlsx")
)
openxlsx::write.xlsx(
  results_df_2[c(2, 4, 3, 5:10,   # 2-year unweighted
                 27, 29:34,       # 2-year IPSW
                 11, 13:18,       # 2-year IPTW
                 35, 37:42), ],   # 2-year IPSW + IPTW
  rowNames = TRUE,
  file = paste0(results_path, "Supplemental/Table_S6.xlsx")
)
openxlsx::write.xlsx(
  results_df_2[c(2, 4, 11, 13:18, # 2-year IPTW
                 43, 45:50,       # 2-year ATT
                 51, 53:58), ],   # 2-year ATU
  rowNames = TRUE,
  file = paste0(results_path, "Supplemental/Table_S8.xlsx")
)

# save cohort
save(
  id_name, 
  listvar,
  varnames,
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
  n_bootstraps,
  out_KM_event_death_2y_IPTW,
  file = file.path(
    "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/cohort_with_models.Rdata"
  )
)

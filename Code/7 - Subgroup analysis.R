################################################################################
### Decision for dialysis versus conservative care
### PART 7 - Subgroup analysis
################################################################################

# Clear the R environment to avoid conflicts with previous variables
rm(list = ls(all.names = TRUE))

# Set the knitting root directory (used for relative file paths in reports)
knitr::opts_knit$set(root.dir = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/")

# Set seed for reproducibility of any random processes
set.seed(1)

# Set working directory
setwd("P:/SCREAM2/SCREAM2_Research/Carolien Maas/")
# Define path to save results
results_path <- "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Results/"

# load libraries
library(data.table)
library(survival)
library(foreach)   # parallel computation
library(rms)

# Load custom plotting and helper functions
source("Code/utils/plots.R")
source("Code/utils/tables.R")
source("Code/utils/weighting.R")
source("Code/utils/data_manipulation.R")
source("Code/utils/outcomes_absolute_risks.R")

# Load cohort data with pre-calculated models
load("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/cohort_with_models.Rdata")

################################################################################
### Define subgroups
################################################################################
# define risk bounds
lower_risk_bounds <- c(0, 55, 65)
upper_risk_bounds <- c(55, 65, 100)

# Create separate data subsets for each subgroup
subgroup_1 <- baseline                                  # All patients
subgroup_2 <- baseline[age < 80]                        # 65 <= Age < 80
subgroup_3 <- baseline[age >= 80]                       # Age >= 80
subgroup_4 <- baseline[egfr2021 >= 15]
subgroup_5 <- baseline[egfr2021 < 15]
subgroup_6 <- baseline[aki == 0]
subgroup_7 <- baseline[aki == 1]
# subgroup_6 <- baseline[cancer == 0]
# subgroup_7 <- baseline[cancer == 1]
# subgroup_8 <- baseline[dm == 0]
# subgroup_9 <- baseline[dm == 1]
# subgroup_10 <- baseline[ihd == 0]
# subgroup_11 <- baseline[ihd == 1]
# subgroup_12 <- baseline[vhd == 0]
# subgroup_13 <- baseline[vhd == 1]
# subgroup_14 <- baseline[pvd == 0]
# subgroup_15 <- baseline[pvd == 1]
# subgroup_16 <- baseline[female == 1]
# subgroup_17 <- baseline[female == 0]
# subgroup_18 <- baseline[albumin >= 35]
# subgroup_19 <- baseline[albumin < 35]
# subgroup_8 <- baseline[pred_risk > lower_risk_bounds[1]/100 & 
#                          pred_risk < upper_risk_bounds[1]/100]
# subgroup_9 <- baseline[pred_risk > lower_risk_bounds[2]/100 & 
#                          pred_risk < upper_risk_bounds[2]/100]
# subgroup_10 <- baseline[pred_risk > lower_risk_bounds[3]/100 & 
#                          pred_risk < upper_risk_bounds[3]/100]

# Initialize dataframe to store results for all subgroups
nr.subgroups <- 7
subgroup_results <- data.frame(
  analysis_name = c(
    "Total population",
    "Age < 80",
    "Age >= 80",
    "eGFR 15 - 20",
    "eGFR < 15",
    "No history of acute kidney injury",
    "History of acute kidney injury"
    # paste("No malignancy"),
    # paste("Malignancy"),
    # paste("No diabetes mellitus"),
    # paste("Diabetes mellitus"),
    # paste("No ischemic heart disease"),
    # paste("Ischemic heart disease"),
    # paste("No valvular heart disease"),
    # paste("Valvular heart disease"),
    # paste("No peripheral heart disease"),
    # paste("Peripheral heart disease"),
    # paste("Females"),
    # paste("Males"),
    # paste("Albumin >= 35"),
    # paste("Albumin < 35"),
    # paste0("Predicted mortality risk ", lower_risk_bounds[1], 
    #        "-", upper_risk_bounds[1], "%"),
    # paste0("Predicted mortality risk ", lower_risk_bounds[2], 
    #        "-", upper_risk_bounds[2], "%"),
    # paste0("Predicted mortality risk ", lower_risk_bounds[3], 
    #        "-", upper_risk_bounds[3], "%")
  ),
  N = rep(NA, nr.subgroups),
  Nonoverlap = rep(NA, nr.subgroups),
  Imbalance = rep(NA, nr.subgroups),
  RD = rep(NA, nr.subgroups),
  RD_lower = rep(NA, nr.subgroups),
  RD_upper = rep(NA, nr.subgroups),
  dRMST = rep(NA, nr.subgroups),
  dRMST_lower = rep(NA, nr.subgroups),
  dRMST_upper = rep(NA, nr.subgroups),
  HR = rep(NA, nr.subgroups),
  HR_lower = rep(NA, nr.subgroups),
  HR_upper = rep(NA, nr.subgroups)
)

# Loop over each subgroup to calculate survival statistics
for (subgroup_nr in 1:nr.subgroups) {
  cat("Currently subgroup", subgroup_nr, "\n")
  # Dynamically select subgroup dataset
  subgroup_data <- eval(parse(text = paste0("subgroup_", subgroup_nr)))
  
  # nonoverlap
  nonoverlap <- nrow(
    trim_propensity_scores(
      data = subgroup_data,
      PS_varname = "ps",
      trt_varname = "trt",
      trim_meth = "common_range"
    )$nonoverlap
  )
  
  # list the variables to be excluded because of sample size issues in categorical variables
  samp_threshold <- 0
  omit_var <- c()
  for (var in catvar){
    table_var <- table(subgroup_data[[var]], subgroup_data[, trt])
    tab_present <- table_var[table_var > 0]
    if (length(tab_present) <= 2 | any(tab_present<samp_threshold)) {
      omit_var <- c(omit_var, var)
    }
  }
  
  # update propensity score model
  if (is.null(omit_var)){
    # keep full model
    model_PS_sub <- model_PS
  } else{
    # remove categorical values as their own name or as.factor() name
    remove_terms <- paste(
      c(paste0("- ", omit_var), paste0("- as.factor(", omit_var, ")")),
      collapse = " "
    )
    
    model_PS_sub <- update(
      model_PS,
      as.formula(paste("~ .", remove_terms))
    )
  }
  
  # check characteristics before weighting
  table_one <- create_baseline_table(
    data = subgroup_data,
    id_name = "LOPNR",
    vars = listvar,
    categoricalVars = catvar,
    IQRVars = non_normal_vars,
    treatmentColumn = "trt",
    treatmentLabel = treatment_label,
    controlLabel = control_label,
    tableCaption = paste("Subgroup", subgroup_nr)
  )
  if (subgroup_nr != 1) {
    openxlsx::write.xlsx(
      table_one$raw_table,
      rowNames = TRUE,
      file = paste0(
        results_path,
        "Other/Descriptives_unweighted_subgroup_analysis_",
        subgroup_nr,
        ".xlsx"
      )
    )
  }
  
  # re-estimation of weights
  subgroup_data_reestimated <- create_weights(
    data = subgroup_data,
    model_PS = model_PS_sub,
    w_meth = "IPTW",
    verbose = FALSE
  )
  
  # Create PS distribution plot
  PS_title_hist <- paste0(2*subgroup_nr - 1, "A. ", 
                          subgroup_results[subgroup_nr, "analysis_name"], "\n")
  PS_title_scaled_hist <- paste0(2*subgroup_nr - 1, "B. ", 
                                 subgroup_results[subgroup_nr, 
                                                  "analysis_name"], "\n")
  unweighted_fig_ps <- create_ps_distribution_plot(
    data = subgroup_data,
    PS_varname = "ps", 
    trt_varname = "trt",
    PS_title_hist = PS_title_hist,
    PS_title_scaled_hist = PS_title_scaled_hist,
    TextSize = 22,
    xlab = ifelse(subgroup_nr == nr.subgroups, TRUE, FALSE),
    ylab = TRUE
  )
  assign(paste0("unweighted_fig_ps_", subgroup_nr),
         unweighted_fig_ps)
  
  # check weights
  subgroup_weights <- subgroup_data_reestimated$data$w
  
  # trim weights
  subgroup_weights_trimmed <- subgroup_weights
  lower_trunc <- as.numeric(quantile(subgroup_weights, 0.001))
  upper_trunc <- as.numeric(quantile(subgroup_weights, 0.999))
  subgroup_weights_trimmed[subgroup_weights_trimmed<lower_trunc] <- lower_trunc
  subgroup_weights_trimmed[subgroup_weights_trimmed>upper_trunc] <- upper_trunc
  
  # Create PS distribution plot
  PS_title_hist <- paste0(2*subgroup_nr, "A. ", 
                          subgroup_results[subgroup_nr, "analysis_name"], "\n")
  PS_title_scaled_hist <- paste0(2*subgroup_nr, "B. ",
                                 subgroup_results[subgroup_nr, "analysis_name"], "\n")
  weighted_fig_ps <- create_ps_distribution_plot(
    data = subgroup_data,
    PS_varname = "ps", 
    trt_varname = "trt",
    weights = subgroup_weights_trimmed,
    PS_title_hist = PS_title_hist,
    PS_title_scaled_hist = PS_title_scaled_hist,
    TextSize = 22,
    xlab = ifelse(subgroup_nr == nr.subgroups, TRUE, FALSE),
    ylab = FALSE
  )
  assign(paste0("weighted_fig_ps_", subgroup_nr),
         weighted_fig_ps)
  
  # check SMD
  table_one_weighted <- create_baseline_table(
    data = subgroup_data_reestimated$data,
    id_name = "LOPNR", 
    weights = subgroup_weights_trimmed,
    vars = listvar,
    categoricalVars = catvar,
    IQRVars = non_normal_vars,
    treatmentColumn = "trt",
    treatmentLabel = treatment_label,
    controlLabel = control_label,
    tableCaption = paste("Subgroup", subgroup_nr)
  )
  if (subgroup_nr != 1) {
    openxlsx::write.xlsx(
      table_one_weighted$raw_table,
      rowNames = TRUE,
      file = paste0(
        results_path,
        "Other/Descriptives_weighted_subgroup_analysis_",
        subgroup_nr,
        ".xlsx"
      )
    )
  }
  
  # add weighted SMD
  SMDs <- table_one_weighted$smd_table
  nr_large_SMDs <- sum(SMDs > 0.1)
  
  # Generate Kaplan-Meier plot and calculate metrics using custom function
  out_KM <- create_KM_plot(
    data = subgroup_data_reestimated$data,
    horizon = horizon,
    unit = unit,
    model_PS = model_PS,
    w_meth = "IPTW",
    weights_meth = subgroup_weights_trimmed,
    event_var = outcome_vars[2],
    time2event_var = time2outcome_vars[2],
    trt_var = "trt",
    competing_event_var = competing_events_vars[2],
    n_bootstraps = n_bootstraps,
    bootstrap_seed = 1,
    plotTitle = subgroup_results[subgroup_nr, "analysis_name"]
  )
  assign(paste0("KM_", subgroup_nr), out_KM$KM_plot)
  
  # Save risk difference (RD) and confidence intervals
  subgroup_results[subgroup_nr, "N"] <- nrow(subgroup_data)
  subgroup_results[subgroup_nr, "Nonoverlap"] <- nonoverlap
  subgroup_results[subgroup_nr, "Imbalance"] <- sum(SMDs > 0.1)
  subgroup_results[subgroup_nr, "RD"] <- as.numeric(out_KM$RD * 100)
  subgroup_results[subgroup_nr, "RD_lower"] <- as.numeric(out_KM$RD_lower * 100)
  subgroup_results[subgroup_nr, "RD_upper"] <- as.numeric(out_KM$RD_upper * 100)
  
  # Save difference in restricted mean survival time (dRMST)
  subgroup_results[subgroup_nr, "dRMST"] <- as.numeric(out_KM$dRMST)
  subgroup_results[subgroup_nr, "dRMST_lower"] <- as.numeric(out_KM$dRMST_lower)
  subgroup_results[subgroup_nr, "dRMST_upper"] <- as.numeric(out_KM$dRMST_upper)
  
  # Save hazard ratio (HR) and confidence intervals
  subgroup_results[subgroup_nr, "HR"] <- as.numeric(out_KM$HR)
  subgroup_results[subgroup_nr, "HR_lower"] <- as.numeric(out_KM$HR_lower)
  subgroup_results[subgroup_nr, "HR_upper"] <- as.numeric(out_KM$HR_upper)
}

# PS distributions
plot_names <- c(paste0("unweighted_fig_ps_", 1:nr.subgroups),
                paste0("weighted_fig_ps_", 1:nr.subgroups))
ggplot2::ggsave(
  plot = patchwork::wrap_plots(combine_PS_plots(mget(plot_names)),
                               ncol = 2, nrow = nr.subgroups, byrow = FALSE),
  filename = paste0(results_path, "Supplemental/Figure_S7.png"),
  width = 30,
  height = 30,
  dpi = 300
)

# Kaplan-Meier plots
km_names <- paste0("KM_", 2:nr.subgroups)
ggplot2::ggsave(
  plot = patchwork::wrap_plots(mget(km_names), 
                               ncol = 2, 
                               nrow = 3),
  filename = paste0(results_path, "Supplemental/Figure_S8.png"),
  width = 12,
  height = 15,
  dpi = 300
)

# create forest plots of subgroup analysis
subgroup_forest <- create_forest_plot_all_measures(dt = setDT(subgroup_results),
                                                   print_metrics = c("N", 
                                                                     "Nonoverlap",
                                                                     "Imbalance"))
ggplot2::ggsave(
  plot = subgroup_forest$combined_plot,
  filename = paste0(results_path, "Main/Figure_2.png"),
  width = 11,
  height = 10,
  dpi = 300
)

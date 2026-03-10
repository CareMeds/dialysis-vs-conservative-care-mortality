################################################################################
### Decision for dialysis versus conservative care
### PART 6 - Sensitivity analysis for positivity
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
source("Code/utils/weighting.R")
source("Code/utils/plots.R")
source("Code/utils/tables.R")
source("Code/utils/data_manipulation.R")
source("Code/utils/compute_absolute_relative_risks.R")

################################################################################
### Load data ##################################################################
################################################################################
load("Data/cohort_with_models.Rdata")

# Define main survival model formula
survival_model <- survival::Surv(time2event_death_2y, event_death_2y) ~ trt

################################################################################
### Save untrimmed data
################################################################################
baseline_untrimmed <- baseline
rm(baseline)

################################################################################
### Sensitivity anaysis for non-positivity, investigate non-overlapping region of PS
################################################################################
trim_meths <- c("", "IPTW", "overlap", "Crump", "Stürmer", "Walker")

# Initialize summary table
summary_table <- data.frame(
  analysis_name = character(),
  N = numeric(),
  Excluded = numeric(),
  Nonoverlap = numeric(),
  Imbalance = numeric(),
  RD = numeric(),
  RD_lower = numeric(),
  RD_upper = numeric(),
  dRMST = numeric(),
  dRMST_lower = numeric(),
  dRMST_upper = numeric(),
  HR = numeric(),
  HR_lower = numeric(),
  HR_upper = numeric(),
  stringsAsFactors = FALSE
)

for (trim_meth in trim_meths) {
  cat("Trimming method:", trim_meth, "\n")
  # Use untrimmed data for common_range, IPTW, and overlap methods
  if (trim_meth %in% c("Crump", "Stürmer", "Walker")) {
    # Apply trimming method to remove extreme PS values
    baseline_trimmed <- trim_propensity_scores(
      data = baseline_untrimmed,
      PS_varname = "ps",
      trt_varname = "trt",
      trim_meth = trim_meth
    )$overlap
  } else{
    baseline_trimmed <- baseline_untrimmed
  }
  
  # Identify non-overlap patients
  nonoverlap <- trim_propensity_scores(
    data = baseline_trimmed,
    PS_varname = "ps",
    trt_varname = "trt",
    trim_meth = "common_range"
  )$nonoverlap
  
  # extract weights
  if (trim_meth == "") {
    # extract weights
    weights_meth <- NULL
    
    # use untrimmed sample
    baseline <- copy(baseline_untrimmed)
    
  } else if (trim_meth == "IPTW" | trim_meth == "overlap") {
    # extract weights
    weights_meth <- baseline_untrimmed[[paste0("sw_", trim_meth)]]
    
    # use untrimmed sample
    baseline <- copy(baseline_untrimmed)
    
  } else if (trim_meth %in% c("Crump", "Stürmer", "Walker")) {
    # re-estimate PS in the trimmed population
    baseline_reestimated <- create_weights(
      data = baseline_trimmed,
      model_PS = model_PS,
      w_meth = "IPTW",
      verbose = FALSE
    )$data
    
    # extract weights
    weights_meth <- baseline_reestimated$w
    
    # use trimmed and re-estimated PS
    baseline <- copy(baseline_reestimated)
  }
  
  # create baseline table of non overlap
  table_one <- create_baseline_table(
    data = baseline,
    id_name = "LOPNR",
    weights = weights_meth,
    vars = listvar,
    categoricalVars = catvar,
    IQRVars = non_normal_vars,
    treatmentColumn = "trt",
    treatmentLabel = treatment_label,
    controlLabel = control_label,
    tableCaption = ""
  )
  
  # Create PS distribution plot
  unweighted_fig_ps <- create_ps_distribution_plot(
    data = baseline,
    PS_varname = "ps",
    trt_varname = "trt",
    PS_title_hist = trim_meth,
    PS_title_scaled_hist = trim_meth,
    xlab = ifelse(trim_meth == tail(trim_meths, 1), TRUE, FALSE),
    ylab = TRUE,
    palette = manual_colors[1:2]
  )
  assign(paste0("unweighted_fig_ps_", trim_meth), unweighted_fig_ps)
  
  weighted_fig_ps <- create_ps_distribution_plot(
    data = baseline,
    PS_varname = "ps",
    trt_varname = "trt",
    weights = weights_meth,
    PS_title_hist = trim_meth,
    PS_title_scaled_hist = trim_meth,
    xlab = ifelse(trim_meth == tail(trim_meths, 1), TRUE, FALSE),
    ylab = TRUE,
    palette = manual_colors[1:2]
  )
  assign(paste0("weighted_fig_ps_", trim_meth), weighted_fig_ps)
  
  # obtain KM, no weighting
  out_KM <- create_KM_plot(
    data = baseline,
    horizon = horizon,
    unit = unit,
    model_PS = model_PS,
    w_meth = ifelse(
      trim_meth == "",
      "",
      ifelse(trim_meth == "overlap", trim_meth, "IPTW")
    ),
    weights_meth = weights_meth,
    event_var = outcome_vars[2],
    time2event_var = time2outcome_vars[2],
    trt_var = "trt",
    competing_event_var = competing_events_vars[2],
    n_bootstraps = n_bootstraps,
    bootstrap_seed = 1,
    plotColors = manual_colors[1:2]
  )
  
  # Append summary statistics for current trimming method
  summary_table <- rbind(
    summary_table,
    data.frame(
      analysis_name = trim_meth,
      N = nrow(baseline),
      Excluded = nrow(baseline_untrimmed) - nrow(baseline),
      Nonoverlap = nrow(nonoverlap),
      Imbalance = sum(table_one$smd_table > 0.1),
      RD = out_KM$RD * 100,
      RD_lower = out_KM$RD_lower * 100,
      RD_upper = out_KM$RD_upper * 100,
      dRMST = out_KM$dRMST,
      dRMST_lower = out_KM$dRMST_lower,
      dRMST_upper = out_KM$dRMST_upper,
      HR = out_KM$HR,
      HR_lower = out_KM$HR_lower,
      HR_upper = out_KM$HR_upper
    )
  )
}

################################################################################
### Combine histograms and scaled histograms of propensity score distributions
### across different trimming methods
################################################################################
# Define the layout: 
# Each letter represents a plot area. 
# # represents an empty space.
# We define a 4-column grid.
layout_design <- "
  ABCD
  ##EF
  GHIJ
  KLMN
  OPQR
"

# Combine the plots into a list to pass to wrap_plots
plot_list <- list(
  # Row 1
  unweighted_fig_ps_$hist + ggplot2::labs(title = "1A. Total population"),
  unweighted_fig_ps_$scaled_hist + ggplot2::labs(title = "1B. Total population"),
  weighted_fig_ps_IPTW$hist + ggplot2::labs(title = "2A. IPTW"),
  weighted_fig_ps_IPTW$scaled_hist + ggplot2::labs(title = "2B. IPTW"),
  
  # Row 2 (Starts at column 3 based on design ##EF)
  weighted_fig_ps_overlap$hist + ggplot2::labs(title = "3A. Overlap weighting"),
  weighted_fig_ps_overlap$scaled_hist + ggplot2::labs(title = "3B. Overlap weighting"),
  
  # Row 3
  unweighted_fig_ps_Crump$hist + ggplot2::labs(title = "4A. Crump trimming"),
  unweighted_fig_ps_Crump$scaled_hist + ggplot2::labs(title = "4B. Crump trimming"),
  weighted_fig_ps_Crump$hist + ggplot2::labs(title = "5A. Crump trimming"),
  weighted_fig_ps_Crump$scaled_hist + ggplot2::labs(title = "5B. Crump trimming"),
  
  # Row 4
  unweighted_fig_ps_Stürmer$hist + ggplot2::labs(title = "6A. Stürmer trimming"),
  unweighted_fig_ps_Stürmer$scaled_hist + ggplot2::labs(title = "6B. Stürmer trimming"),
  weighted_fig_ps_Stürmer$hist + ggplot2::labs(title = "7A. Stürmer trimming"),
  weighted_fig_ps_Stürmer$scaled_hist + ggplot2::labs(title = "7B. Stürmer trimming"),
  
  # Row 5
  unweighted_fig_ps_Walker$hist + ggplot2::labs(title = "8A. Walker trimming"),
  unweighted_fig_ps_Walker$scaled_hist + ggplot2::labs(title = "8B. Walker trimming"),
  weighted_fig_ps_Walker$hist + ggplot2::labs(title = "9A. Walker trimming"),
  weighted_fig_ps_Walker$scaled_hist + ggplot2::labs(title = "9B. Walker trimming")
)

# Generate final plot
final_plot <- patchwork::wrap_plots(plot_list, design = layout_design) + 
  patchwork::plot_layout(guides = "collect")

# Save
ggplot2::ggsave(
  plot = final_plot,
  filename = paste0(results_path, "Supplemental/Figure_S7.png"),
  width = 20,
  height = 20,
  dpi = 300
)

# create forest plots of subgroup analysis
summary_table <- summary_table[-1, ]
summary_table$analysis_name <- c(
  "Main analysis",
  "Overlap weighting",
  "Crump trimming",
  "Stürmer trimming",
  "Walker trimming"
)
nonoverlap_forest <- create_forest_plot_all_measures(dt = setDT(summary_table),
                                                     print_metrics = c("N", "Nonoverlap"))
ggplot2::ggsave(
  plot = nonoverlap_forest$combined_plot,
  filename = paste0(results_path, "Supplemental/Figure_S8.png"),
  width = 11,
  height = 8,
  dpi = 300
)

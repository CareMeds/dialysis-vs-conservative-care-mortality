################################################################################
### Decision for dialysis versus conservative care
### PART 6 - Sensitivity analysis for positivity
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
library(foreach)   # parallel computation

# load functions
source("Code/utils/weighting.R")
source("Code/utils/plots.R")
source("Code/utils/tables.R")
source("Code/utils/data_manipulation.R")
source("Code/utils/outcomes_absolute_risks.R")

################################################################################
### Load data ##################################################################
################################################################################
load("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/cohort_with_models.Rdata")

# Define main survival model formula
survival_model <- survival::Surv(time2event_death_2y, event_death_2y) ~ trt

################################################################################
### Save untrimmed data
################################################################################
baseline_untrimmed <- baseline
rm(baseline)

# hist(baseline_untrimmed$prob_int_cox)
quantile(baseline_untrimmed$prob_int_cox)

################################################################################
### Sensitivity anaysis for non-positivity, investigate non-overlapping region of PS
################################################################################
no_trimming_methods <- c("common_range",
                         "overlap")
trimming_methods <- c("Crump",
                      "Stürmer",
                      "Walker")
trim_meths <- c(
  no_trimming_methods,
  trimming_methods
)

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

n_bootstraps <- 2
for (trim_meth in trim_meths) {
  # Use untrimmed data for common_range, IPTW, and overlap methods
  if (trim_meth %in% trimming_methods) {
    # Apply trimming method to remove extreme PS values
    baseline_trimmed <- trim_propensity_scores(
      data = baseline_untrimmed,
      PS_varname = "ps",
      trt_varname = "trt",
      trim_meth = trim_meth
    )$overlap
  } else{
    # unweighted, IPTW, overlap does not require trimming
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
  if (trim_meth == "common_range") {
    # extract weights
    weights_meth <- NULL
    
    # use untrimmed sample
    baseline <- copy(baseline_untrimmed)
    
  } else if (trim_meth == "IPTW" | trim_meth == "overlap") {
    # extract weights
    weights_meth <- baseline_untrimmed[[paste0("sw_", trim_meth)]]
    
    # use untrimmed sample
    baseline <- copy(baseline_untrimmed)
    
  } else if (trim_meth %in% trimming_methods){
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
    tableCaption = "")
  
  # Create PS distribution plot
  nr.meth <- which(trim_meth == trim_meths)
  PS_title_hist <- paste0(ifelse(trim_meth=="common_range", 1, 2*nr.meth-3),
                          "A. ", plot_title(trim_meth, weighting = FALSE), "\n")
  PS_title_scaled_hist <- paste0(ifelse(trim_meth=="common_range", 1, 2*nr.meth-3), 
                                 "B. ", plot_title(trim_meth, weighting = FALSE), "\n")
  unweighted_fig_ps <- create_ps_distribution_plot(
    data = baseline,
    PS_varname = "ps", 
    trt_varname = "trt",
    PS_title_hist = PS_title_hist,
    PS_title_scaled_hist = PS_title_scaled_hist,
    xlab = ifelse(trim_meth == tail(trim_meths, 1), TRUE, FALSE),
    ylab = TRUE
  )
  assign(paste0("unweighted_fig_ps_", trim_meth), unweighted_fig_ps)
  
  PS_title_hist <- paste0(2*nr.meth-2, "A. ", plot_title(trim_meth), "\n")
  PS_title_scaled_hist <- paste0(2*nr.meth-2, "B. ", plot_title(trim_meth), "\n")
  weighted_fig_ps <- create_ps_distribution_plot(
    data = baseline,
    PS_varname = "ps", 
    trt_varname = "trt",
    weights = weights_meth,
    PS_title_hist = PS_title_hist,
    PS_title_scaled_hist = PS_title_scaled_hist,
    xlab = ifelse(trim_meth == tail(trim_meths, 1), TRUE, FALSE),
    ylab = FALSE
  )
  assign(paste0("weighted_fig_ps_", trim_meth), weighted_fig_ps)
  
  # obtain KM, no weighting
  out_KM <- create_KM_plot(
    data = baseline,
    horizon = horizon,
    unit = unit,
    model_PS = model_PS,
    w_meth = ifelse(trim_meth=="common_range", "", 
                    ifelse(trim_meth == "overlap", trim_meth, "IPTW")),
    weights_meth = weights_meth,
    event_var = outcome_vars[2],
    time2event_var = time2outcome_vars[2],
    trt_var = "trt",
    competing_event_var = competing_events_vars[2],
    n_bootstraps = n_bootstraps,
    bootstrap_seed = 1
  )
  
  # remove x-axis for all plots but the last one
  KM_plot <- out_KM$KM_plot +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank()
    ) +
    ggplot2::ggtitle(paste(nr.meth, ".", plot_title(trim_meths[nr.meth])))
  
  if (trim_meth!=length(trim_meths)) {
    KM_table <- out_KM$KM_table +
      ggplot2::xlab("Time in years")
  } else{
    KM_table <- out_KM$KM_table +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank()
      )
  }
  
  KM_plot_table <- ggpubr::ggarrange(KM_plot,
                                     KM_table, 
                                     heights=c(3, 1),
                                     nrow=2,
                                     ncol=1, 
                                     align="v")
  assign(paste0("KM_", trim_meth), KM_plot_table)
  
  # Append summary statistics for current trimming method
  summary_table <- rbind(
    summary_table,
    data.frame(
      analysis_name = trim_meth,
      N = nrow(baseline),
      Excluded = nrow(baseline_untrimmed) - nrow(baseline),
      Nonoverlap = nrow(nonoverlap),
      Imbalance = sum(table_one$smd_table > 0.1),
      RD = out_KM$RD*100,
      RD_lower = out_KM$RD_lower*100,
      RD_upper = out_KM$RD_upper*100,
      dRMST = out_KM$dRMST,
      dRMST_lower = out_KM$dRMST_lower,
      dRMST_upper = out_KM$dRMST_upper,
      HR = out_KM$HR,
      HR_lower = out_KM$HR_lower,
      HR_upper = out_KM$HR_upper
    )
  )
}

# Combine histograms and scaled histograms of propensity score distributions
# across different trimming methods
# empty_plot <- ggplot2::ggplot() + ggplot2::theme_void()
# empty_hists <- list(hist = empty_plot, scaled_hist = empty_plot)
plot_names <- c("unweighted_fig_ps_common_range", 
                "unweighted_fig_ps_Crump",
                "unweighted_fig_ps_Stürmer",
                "unweighted_fig_ps_Walker",
                "weighted_fig_ps_overlap",
                "weighted_fig_ps_Crump",
                "weighted_fig_ps_Stürmer",
                "weighted_fig_ps_Walker")
plot_names
ggplot2::ggsave(
  plot = patchwork::wrap_plots(combine_PS_plots(mget(plot_names)),
                               ncol = 2, nrow = 4, byrow = FALSE),
  filename = paste0(results_path, "Supplemental/Figure_S10.png"),
  width = 20,
  height = 20,
  dpi = 300
)

# Combine Kaplan-Meier survival curves and corresponding summary tables
# for each trimming methods
km_names <- paste0("KM_", trim_meths)
ggplot2::ggsave(
  plot = patchwork::wrap_plots(mget(km_names), ncol = 2, nrow = 3),
  filename = paste0(results_path, "Other/KM_positivity.png"),
  width = 15,
  height = 15,
  dpi = 300
)

# create forest plots of subgroup analysis
summary_table$analysis_name <- c(
  "Unweighted",
  "Main analysis",
  "Overlap weighting",
  "Crump trimming",
  "Stürmer trimming",
  "Walker trimming"
)
nonoverlap_forest <- create_forest_plot_all_measures(dt = setDT(summary_table[-1, ]),
                                                     print_metrics = c("N", "Nonoverlap"))
ggplot2::ggsave(
  plot = nonoverlap_forest$combined_plot,
  filename = paste0(results_path, "Supplemental/Figure_S11.png"),
  width = 11,
  height = 8,
  dpi = 300
)
print(nonoverlap_forest$dt)

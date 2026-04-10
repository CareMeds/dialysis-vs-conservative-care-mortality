################################################################################
### Decision for dialysis versus conservative care
### PART 8 - Continuous heterogeneous treatment effect estimation
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
library(survival)
library(rms)
library(foreach)   # parallel computation
library(doRNG)     # handle parallel seeds
set.seed(1)        # set seed for parallel backend

# load functions
source("Code/utils/weighting.R")
source("Code/utils/plots.R")
source("Code/utils/tables.R")
source("Code/utils/data_manipulation.R")
source("Code/utils/compute_absolute_relative_risks.R")

# perform internal validation 
validate <- TRUE

################################################################################
### Load data ##################################################################
################################################################################
load("Data/cohort_with_models.Rdata")

################################################################################
### Internal 2-year time-to-event risk model
################################################################################
# make outcome for elig cohort
elig_Surv <- survival::Surv(elig_cohort$time2event_death_2y, elig_cohort$event_death_2y)

# use external predictors identified by Chava
elig_cohort[, log_crp := log(crp + 1)]
baseline[, log_crp := log(crp + 1)]
risk_model_cox <- survival::coxph(
  elig_Surv ~ age + egfr2021 + cancer + dm + ihd +
    vhd + pvd + female + albumin + log_crp,
  data = elig_cohort,
  method = "breslow"
)

# save table
table_risk <- risk_model_table(
  model_cox  = risk_model_cox,
  predictor_labels = c(
    "   Age (per 1 year)",
    "   eGFR (per 1 ml/min/1.73 m²)",
    "   Malignancies",
    "   Diabetes mellitus",
    "   Ischemic heart disease",
    "   Valvular heart disease",
    "   Peripheral vascular disease",
    "   Female vs. male",
    "   Albumin (per 1 g/L)",
    "   CRP (log-transformed, per unit increase)"
  ),
  horizon      = horizon
)

# create plot on hazard scale
dd <- rms::datadist(elig_cohort)
options(datadist = "dd")
risk_model <- rms::cph(
  elig_Surv ~ age + egfr2021 + cancer + dm + ihd +
    vhd + pvd + female + albumin + log(crp + 1),
  data = elig_cohort,
  method = c("breslow"),
  y = TRUE,
  x = TRUE
)
png(
  filename = "Results/Supplemental/Figure_M1_Risk_model.png",
  width = 2000,
  height = 2000,
  res = 300
)
plot(Predict(risk_model, fun = exp), ylab = "Hazard Ratio")
dev.off()

# ################################################################################
# ### Internal validation
# ################################################################################
# --- Initialise storage --------------------------------------------------
if (!validate){
  n_bootstraps_elig <- 1
} else {
  n_bootstraps_elig <- 1000
}
n_iter <- n_bootstraps_elig + 1L          # iteration 1 = original sample

metrics <- list(
  orig = data.frame(
    elig_Intercept = numeric(n_iter),
    elig_Slope     = numeric(n_iter),
    elig_AUC       = numeric(n_iter),
    trt_Intercept  = numeric(n_iter),
    trt_Slope      = numeric(n_iter),
    trt_AUC        = numeric(n_iter)
  ),
  boot = data.frame(
    elig_Intercept = numeric(n_iter),
    elig_Slope     = numeric(n_iter),
    elig_AUC       = numeric(n_iter),
    trt_Intercept  = numeric(n_iter),
    trt_Slope      = numeric(n_iter),
    trt_AUC        = numeric(n_iter)
  )
)

cal_plots <- list()

# --- Main loop -----------------------------------------------------------
for (B in seq_len(n_iter)) {
  
  is_original <- (B == 1L)
  
  # 1. Draw samples -------------------------------------------------------
  if (is_original) {
    boot_elig    <- elig_cohort
    boot_trt_dec <- baseline
  } else {
    boot_elig    <- elig_cohort[sample(nrow(elig_cohort), replace = TRUE), ]
    boot_trt_dec <- baseline[sample(nrow(baseline),      replace = TRUE), ]
  }
  
  # 2. Fit risk model on bootstrap sample ---------------------------------
  boot_risk_model <- rms::cph(
    survival::Surv(time2event_death_2y, event_death_2y) ~
      age + egfr2021 + cancer + dm + ihd + vhd + pvd + female + albumin,
    data   = boot_elig,
    method = "breslow",
    y      = TRUE,
    x      = TRUE
  )
  
  # 3. Evaluate on original sample (apparent for B=1, optimism for B>1) --
  elig_orig <- compute_measures(boot_risk_model, data = elig_cohort,
                                plot = is_original)
  trt_orig  <- compute_measures(boot_risk_model, data = baseline,
                                plot = is_original)
  
  if (is_original) {
    cal_plots$elig <- calibration_plot(elig_orig)
    cal_plots$trt  <- calibration_plot(trt_orig)
  }
  
  # 4. Evaluate on bootstrap sample (needed only for optimism correction) -
  if (!is_original) {
    elig_boot <- compute_measures(boot_risk_model, data = boot_elig)
    trt_boot  <- compute_measures(boot_risk_model, data = boot_trt_dec)
  }
  
  # 5. Store results ------------------------------------------------------
  metrics$orig[B, ] <- c(extract_metrics(elig_orig),
                         extract_metrics(trt_orig))
  
  if (!is_original) {
    metrics$boot[B, ] <- c(extract_metrics(elig_boot),
                           extract_metrics(trt_boot))
  }
}

# --- Compute optimism (rows 2:n = bootstrap iterations) -----------------
# metrics$boot[1, ] is never filled (all zeros) — optimism correctly uses [-1, ]
optimism   <- colMeans(metrics$orig[-1, ] - metrics$boot[-1, ])
apparent   <- metrics$orig[1, ]
orig_boots <- metrics$orig[-1, ]

# --- Build and save plots -----------------------------------------------
for (cohort in list(list(prefix = "elig", label = "elig"),
                    list(prefix = "trt",  label = "trt"))) {
  
  annotated <- annotate_cal_plot(
    cal_plot_obj  = cal_plots[[cohort$prefix]],
    apparent_vals = cohort_apparent(cohort$prefix),
    boot_vals     = cohort_boots(cohort$prefix),
    optimism_vals = cohort_optimism(cohort$prefix)
  )
  
  save_cal_plot(
    cal_plot_obj   = cal_plots[[cohort$prefix]],
    annotated_plot = annotated,
    filename       = paste0("Figure_M1_calibration_", cohort$label, ".png")
  )
}

################################################################################
### Explore linearity of age
################################################################################
fit_age_linear <- survival::coxph(
  survival::Surv(time2event_death_2y, event_death_2y) ~
    trt * age,
  data = baseline,
  method = c("breslow"),
  weights = baseline$sw_IPTW,
  x = TRUE,
  y = TRUE
)

# save table
ITE_model_age <- risk_model_table(
  model_cox  = fit_age_linear,
  predictor_labels = c(
    "   Dialysis",
    "   Age",
    "   Dialysis * age"
  ),
  horizon      = horizon
)

################################################################################
### Explore linearity of risk
################################################################################
baseline$lp_risk <- predict(risk_model, newdata = baseline)
fit_lp_linear <- survival::coxph(
  survival::Surv(time2event_death_2y, event_death_2y) ~
    trt * lp_risk,
  data = baseline,
  method = c("breslow"),
  weights = baseline$sw_IPTW,
  x = TRUE,
  y = TRUE
)

# save table
ITE_model_lp <- risk_model_table(
  model_cox  = fit_lp_linear,
  predictor_labels = c(
    "   Dialysis",
    "   Linear predictor",
    "   Dialysis * Linear predictor"
  ),
  horizon      = horizon
)

# save to xlsx
summmarized_table <- rbind(
  c("Risk model", "", "", ""),
  table_risk$risk_model_table,
  c("ITE model for age", "", "", ""),
  ITE_model_age$risk_model_table,
  c("ITE model for linear predictor", "", "", ""),
  ITE_model_lp$risk_model_table
)
openxlsx::write.xlsx(summmarized_table, 
                     file = paste0(results_path, "Supplemental/Table_S4_risk_ITE.xlsx"), 
                     rowNames = FALSE)

################################################################################
### Continuous HTE
################################################################################
data_sets <- c("baseline", "baseline[Davies_score >= 2]")
# probabilities between 40% and 90% mortality risk
p_range <- c(0.2, 0.4, 0.9)
for (nr_analysis in 1:2) {
  cat("Perform analysis on",
      ifelse(nr_analysis == 1, "full data\n", "Davies >= 2 data\n"))
  # determine LP and survival probability using risk model
  analysis_data <- eval(parse(text = data_sets[nr_analysis]))
  analysis_data$lp_risk <- predict(risk_model, newdata = analysis_data)
  analysis_data$pred_risk <- PredictionTools::fun.event(h0 = table_risk$h0, 
                                                        lp = analysis_data$lp_risk)
  
  # compute estimates across 100 risk points
  for (B in 1:(n_bootstraps + 1)) {
    if (B == 1) {
      # Original sample
      bootstrap <- analysis_data
    } else{
      # create bootstrap sample
      bootstrap <- analysis_data[sample(1:nrow(analysis_data), replace = TRUE), ]
    }
    
    # re-estimate weights
    bootstrap_reestimated <- create_weights(
      data = bootstrap,
      model_PS = model_PS,
      w_meth = "IPTW",
      catvar = catvar,
      contvar = contvar,
      verbose = FALSE
    )
    bootstrap$sw_IPTW <- bootstrap_reestimated$data$w
    
    # check SMDs
    if (B == 1) {
      table_one_weighted <- create_baseline_table(
        data = bootstrap_reestimated$data,
        id_name = "LOPNR",
        weights = bootstrap_reestimated$data$w,
        vars = listvar,
        categoricalVars = catvar,
        IQRVars = non_normal_vars,
        treatmentColumn = trt_var,
        treatmentLabel = treatment_label,
        controlLabel = control_label,
        tableCaption = paste("Subgroup", 2)
      )
      cat("Number of SMDs > 0.1",
          sum(table_one_weighted$smd_table > 0.1),
          "\n")
    }
    
    # compute HTE across age
    show_test <- ifelse(B == 1, TRUE, FALSE)
    age_df <- compute_HTE(
      data = bootstrap,
      unit = unit,
      horizon = horizon,
      event_var = outcome_var,
      time2event_var = time2outcome_var,
      effect_modifier = "age",
      effect_modifier_range = seq(65, 95, 1),
      add_interaction = TRUE,
      test_relative_HTE = show_test
    )
    
    # compute HTE across predicted risk
    # get corresponding linear predictor
    lp_range <- log(-log(1 - p_range) / table_risk$h0)
    pred_risk_df <- compute_HTE(
      data = bootstrap,
      unit = unit,
      horizon = horizon,
      event_var = outcome_var,
      time2event_var = time2outcome_var,
      effect_modifier = "lp_risk",
      effect_modifier_range = sort(c(
        seq(lp_range[1], lp_range[3], length.out = 99), lp_range[2]
      )),
      add_interaction = TRUE,
      test_relative_HTE = show_test
    )
    
    # save results
    if (B == 1) {
      # report effect modifier using probabilities and not linear predictor
      pred_risk_df$effect_modifier_range <- PredictionTools::fun.event(h0 = table_risk$h0,
                                                                       lp = pred_risk_df$effect_modifier_range)
      summary(pred_risk_df$effect_modifier_range)
      
      # original sample
      age_RD <- age_df[, c("effect_modifier_range", "RD")]
      age_RR <- age_df[, c("effect_modifier_range", "RR")]
      age_dRMST <- age_df[, c("effect_modifier_range", "dRMST")]
      age_HR <- age_df[, c("effect_modifier_range", "HR")]
      p_value_age <- unique(age_df$p_for_HTE)
      p_value_age <- ifelse(p_value_age < 0.001, "<0.001", sprintf("%.3f", p_value_age))
      
      pred_risk_RD <- pred_risk_df[, c("effect_modifier_range", "RD")]
      pred_risk_RR <- pred_risk_df[, c("effect_modifier_range", "RR")]
      pred_risk_dRMST <- pred_risk_df[, c("effect_modifier_range", "dRMST")]
      pred_risk_HR <- pred_risk_df[, c("effect_modifier_range", "HR")]
      p_value_pred_risk <- unique(pred_risk_df$p_for_HTE)
      p_value_pred_risk <- ifelse(p_value_pred_risk < 0.001, "<0.001", sprintf("%.3f", p_value_pred_risk))
      
    } else{
      # bootstrapped sample
      age_RD[, paste0("boot_", B - 1)] <- age_df$RD
      age_RR[, paste0("boot_", B - 1)] <- age_df$RR
      age_dRMST[, paste0("boot_", B - 1)] <- age_df$dRMST
      age_HR[, paste0("boot_", B - 1)] <- age_df$HR
      
      pred_risk_RD[, paste0("boot_", B - 1)] <- pred_risk_df$RD
      pred_risk_RR[, paste0("boot_", B - 1)] <- pred_risk_df$RR
      pred_risk_dRMST[, paste0("boot_", B - 1)] <- pred_risk_df$dRMST
      pred_risk_HR[, paste0("boot_", B - 1)] <- pred_risk_df$HR
    }
  }
  
  assign(paste0("age_RD_", nr_analysis), age_RD)
  assign(paste0("age_RR_", nr_analysis), age_RR)
  assign(paste0("age_dRMST_", nr_analysis), age_dRMST)
  assign(paste0("age_HR_", nr_analysis), age_HR)
  assign(paste0("p_value_age_", nr_analysis), p_value_age)
  assign(paste0("pred_risk_RD_", nr_analysis), pred_risk_RD)
  assign(paste0("pred_risk_RR_", nr_analysis), pred_risk_RR)
  assign(paste0("pred_risk_dRMST_", nr_analysis), pred_risk_dRMST)
  assign(paste0("pred_risk_HR_", nr_analysis), pred_risk_HR)
  assign(paste0("p_value_pred_risk_", nr_analysis), p_value_pred_risk)
}

# compute figures and tables
HTE_table <- data.frame()
for (nr_analysis in 1:2) {
  for (effect_modifier in c("age", "pred_risk")) {
    # create histogram of variable stratified by treatment
    analysis_data <- eval(parse(text = data_sets[nr_analysis]))
    analysis_data$lp_risk <- predict(risk_model, newdata = analysis_data)
    analysis_data$pred_risk <- PredictionTools::fun.event(h0 = table_risk$h0, 
                                                          lp = analysis_data$lp_risk)
    
    hist_stratified <- create_histogram_stratified(
      dt = analysis_data,
      var_name = effect_modifier,
      trt_name = trt_var,
      manual_colors = manual_colors
    )
    
    for (measure in c("RD", "dRMST", "RR", "HR")) {
      # legend of histogram
      if (effect_modifier == "pred_risk" & measure == "dRMST") {
        # add legend title to dRMST plot
        hist_stratified <- hist_stratified +
          ggplot2::theme(legend.text = ggplot2::element_text(size = 18))
      } else {
        # remove legends from other plots
        hist_stratified <- hist_stratified +
          ggplot2::theme(legend.position = "none")
      }
      
      # Retrieve the object from the environment
      df_name <- paste0(effect_modifier, "_", measure, "_", nr_analysis)
      estimates_df <- get(df_name)
      
      # Calculate quantiles
      # Assuming the first column(s) are the bootstrap replicates to be excluded
      estimates_df[, paste0(measure, "_lower")] <- apply(estimates_df[, -c(1:2)],
                                                         1,
                                                         quantile,
                                                         probs = 0.025,
                                                         na.rm = TRUE)
      estimates_df[, paste0(measure, "_upper")] <- apply(estimates_df[, -c(1:2)],
                                                         1,
                                                         quantile,
                                                         probs = 0.975,
                                                         na.rm = TRUE)
      
      # make effect plot
      effect_plot_out <- effect_plot(
        estimates_df = estimates_df,
        y_middle = ifelse(measure == "HR" | measure == "RR", 1, 0),
        measure = measure,
        y_min_RD = -70,
        y_max_RD = ifelse(nr_analysis == 1, 21, 56),
        y_min_RR = 0,
        y_max_RR = ifelse(nr_analysis == 1, 1.3, 1.8),
        y_min_dRMST = ifelse(nr_analysis == 1, -3, -8),
        y_max_dRMST = 10,
        y_min_HR = 0,
        y_max_HR = ifelse(nr_analysis == 1, 1.3, 1.8)
      )
      
      # control x-axis of effect plot
      if (effect_modifier == "age") {
        effect_plot_out <- effect_plot_out +
          ggplot2::scale_x_continuous(limits = c(60, 100),
                                      breaks = seq(60, 100, 5))
      } else{
        effect_plot_out <- effect_plot_out +
          ggplot2::scale_x_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, 0.1),
            labels = scales::percent_format(accuracy = 1)
          )
      }
      
      # combine the effect and histogram
      combined_effect_hist <- (effect_plot_out / hist_stratified) +
        plot_layout(heights = c(1, 0.2))
      
      assign(paste0(effect_modifier, "_", measure, "_plot"),
             combined_effect_hist)
      
      # save table
      if (effect_modifier == "age") {
        sel_rows <- c(
          which(estimates_df$effect_modifier_range == 65),
          which(estimates_df$effect_modifier_range == 90)
        )
        range <- estimates_df[sel_rows, "effect_modifier_range"]
      } else{
        sel_rows <- c(
          which(estimates_df$effect_modifier_range == p_range[2]),
          which(estimates_df$effect_modifier_range == p_range[3])
        )
        range <- sprintf("%.2f", estimates_df[sel_rows, "effect_modifier_range"])
      }
      est <- sprintf(ifelse(measure == "RR" |
                              measure == "HR", "%.2f", "%.1f"),
                     estimates_df[sel_rows, measure])
      CI <- paste0(
        "(",
        sprintf(
          ifelse(measure == "RR" | measure == "HR", "%.2f", "%.1f"),
          estimates_df[sel_rows, paste0(measure, "_lower")]
        ),
        ", ",
        sprintf(
          ifelse(measure == "RR" | measure == "HR", "%.2f", "%.1f"),
          estimates_df[sel_rows, paste0(measure, "_upper")]
        ),
        ")"
      )
      
      if (measure == "RD") {
        HTE_tab <- cbind(range, est, CI)
      } else{
        HTE_tab <- cbind(HTE_tab, est, CI)
      }
    }
    HTE_table <- rbind(HTE_table, HTE_tab)
  }
  
  # save figure
  ggplot2::ggsave(
    plot = (age_RD_plot | age_RR_plot | age_dRMST_plot | age_HR_plot) /
      (
        pred_risk_RD_plot | 
          pred_risk_RR_plot |
          pred_risk_dRMST_plot | 
          pred_risk_HR_plot
      ),
    filename = paste0(
      results_path,
      ifelse(
        nr_analysis == 1,
        "Main/Figure_3.png",
        "Supplemental/Figure_S7_DCS.png"
      )
    ),
    width = 20,
    height = 15,
    dpi = 300
  )
}

# save table
HTE_table_final <- data.frame(
  HTE_table,
  c(p_value_age_1, "", p_value_pred_risk_1, "",
    p_value_age_2, "", p_value_pred_risk_2, "")
)
colnames(HTE_table_final) <- c("Effect modifier",
                               "RD",
                               "95% CI",
                               "dRMST",
                               "95% CI",
                               "RR",
                               "95% CI",
                               "HR",
                               "95% CI",
                               "p-value")
openxlsx::write.xlsx(
  HTE_table_final,
  rowNames = FALSE,
  file = paste0(results_path, "Supplemental/Table_S6_HTE.xlsx")
)

# save variables
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
  model_PS,
  coef_PS_overall,
  elig_cohort,
  w_meths,
  trt_var,
  horizon,
  unit,
  n_bootstraps,
  manual_colors,
  estimates_df,
  metrics,
  table_risk,
  ITE_model_lp,
  ITE_model_age,
  file = file.path("Data/cohort_with_prob.Rdata")
)

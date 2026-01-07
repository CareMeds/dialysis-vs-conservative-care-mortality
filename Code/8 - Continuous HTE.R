################################################################################
### Decision for dialysis versus conservative care
### PART 8 - Continuous heterogeneous treatment effect estimation
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
library(survival)
library(foreach)
library(rms)

# load functions
source("Code/utils/weighting.R")
source("Code/utils/plots.R")
source("Code/utils/tables.R")
source("Code/utils/data_manipulation.R")
source("Code/utils/outcomes_absolute_risks.R")
source("Code/utils/risk_model.R")

################################################################################
### Load data ##################################################################
################################################################################
load("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/cohort_with_models.Rdata")

################################################################################
### Internal 2-year time-to-event risk model
################################################################################
# make outcome for elig cohort
elig_Surv <- survival::Surv(elig_cohort$time2event_death_2y, elig_cohort$event_death_2y)

# use external predictors identified by Chava
dd <- rms::datadist(elig_cohort)
options(datadist = "dd")
int_cox_model <- rms::cph(
  elig_Surv ~ age + egfr2021 + cancer + dm + ihd +
    vhd + pvd + female + albumin,
  data = elig_cohort,
  y = TRUE,
  x = TRUE
)
# Wald test
anova(int_cox_model)
# log relative hazard scale
plot(Predict(int_cox_model))
# hazard scale
plot(Predict(int_cox_model, fun = exp))

# Hazard ratios from model
sum_df <- as.data.frame(summary(int_cox_model))
wald_test <- anova(int_cox_model)[-nrow(anova(int_cox_model)), 1]
pred_df <- data.frame(
  short = rownames(sum_df[sum_df$Type == 1, ]),
  long = c(
    paste("Age", quantile(baseline$age)[2], "vs", quantile(baseline$age)[4], "years"),
    paste("eGFR", round(quantile(baseline$egfr2021)[2]), "vs", round(quantile(baseline$egfr2021)[4])),
    "Maligancies",
    "Diabetes mellitus",
    "Ischemic heart disease",
    "Heart valve disease",
    "Primary vascular disease",
    "Female versus male",
    paste("Albumin", quantile(baseline$albumin)[2], "vs", quantile(baseline$albumin)[4])
  )
)
hr_table <- data.frame(
  Predictor = pred_df[rownames(sum_df[sum_df$Type == 1, ]) %in% pred_df$short, "long"],
  HR_CI = fmt_ci(sum_df[sum_df$Type == 2, "Effect"],
                 sum_df[sum_df$Type == 2, "Lower 0.95"],
                 sum_df[sum_df$Type == 2, "Upper 0.95"], digits = 2),
  Wald = fmt(wald_test)
)
openxlsx::write.xlsx(
  hr_table,
  rowNames = FALSE,
  file = paste0(results_path, "Supplemental/Table_S7.xlsx")
)

# backward selection from full model
# Davies_comorbidities is cancer, ihd, pvd, hf, dm, scvd, copd, cirr, psycho, hiv
int_cox_model_full <- rms::cph(
  elig_Surv ~ age + female + log(Davies_score+1) + 
    region + clinic_level + calendar_year_cat +
    rcs(egfr2021, 3) + pol(sbp, 2) + dbp +
    rcs(calcium_total, 3) + phosphate + albumin + pol(hb, 2) + prd_cat +
    cancer + ihd + pvd + hf + dm + # scvd +   # exclude for multi-colinearity
    copd + cirr + acs + psycho + 
    hyperten + vhd + cevd + af + arrh + lung + thrombo + liver + fracture + aki +
    bblock + calblock + diuretic + rasi + lipid + phosbinder + esa + vitamind +
    digoxin + vasodilator + antiplatelet + anticoag + iron_cat +
    log(n_hospital+1) + log(n_cvd_hospital+1) + edu,
  data = elig_cohort,
  y = TRUE,
  x = TRUE
)
anova(int_cox_model_full)
# plot(Predict(int_cox_model_full, Davies_score), main = "Davies_score")
rms::fastbw(int_cox_model_full, k.aic = log(nrow(elig_cohort))) # BIC criterion
# age, Davies_score, phosphate, albumin, n_hospital  
# age, phosphate, albumin, hf, n_hospital, edu

bw_model <- rms::cph(
  elig_Surv ~ age + log(Davies_score+1) + phosphate + albumin + log(n_hospital+1),
  # elig_Surv ~ age + phosphate + albumin + hf + log(n_hospital+1) + edu,
  data = elig_cohort,
  y = TRUE,
  x = TRUE
)
BIC(bw_model)
# log relative hazard scale
plot(Predict(bw_model))
anova(bw_model)
# int_cox_model <- bw_model

# obtain baseline hazard from model
basehaz <- survival::basehaz(int_cox_model)
h0_cox_model <- basehaz$hazard[basehaz$time == max(basehaz$time[basehaz$time <= horizon])]
cat("Baseline hazard at", horizon, " ", unit, ":", h0_cox_model, "\n")

# CHECK average prediction is equal to 1 - KM at 5 years
lp_cox_elig <- predict(int_cox_model, type = "lp")
pred_risk_elig <- PredictionTools::fun.event(lp = lp_cox_elig, h0 = h0_cox_model)
mean(pred_risk_elig)
KM <- survival::survfit(survival::Surv(time2event_death_2y / 365, event_death_2y) ~ 1,
                        data = elig_cohort)
1 - KM$surv[KM$time == max(KM$time)]

# obtain linear predictor from model
lp_cox_model <- predict(
  int_cox_model,
  newdata = baseline, # predict for analysis cohort
  type = "lp"
)

# make probabilities
baseline$lp <- lp_cox_model
baseline$pred_risk <- PredictionTools::fun.event(lp = lp_cox_model, 
                                                 h0 = h0_cox_model)

# validate internal model
survival::concordance(int_cox_model)$concordance
# TODO: survival load is required
baseline_Surv <- survival::Surv(baseline$time2event_death_2y, baseline$event_death_2y)
png(filename = paste0(results_path, "Other/Performance_Model.png"))
PredictionTools::val.surv.mi(p = as.matrix(baseline$pred_risk), y = baseline_Surv)
dev.off()

cat(
  " C-index in full eligible cohort: ",
  timeROC::timeROC(
    T = elig_cohort$time2event_death_2y,
    delta = elig_cohort$event_death_2y,
    marker = predict(int_cox_model, type = "lp"),
    times = 2 * 365,
    cause = 1
  )$AUC[2],
  "\n",
  "C-index in treatment decision cohort: ",
  timeROC::timeROC(
    T = baseline$time2event_death_2y,
    delta = baseline$event_death_2y,
    marker = predict(int_cox_model, newdata = baseline, type = "lp"),
    times = 2 * 365,
    cause = 1
  )$AUC[2],
  "\n"
)

################################################################################
### Continuous HTE
################################################################################
dd <- rms::datadist(baseline)
options(datadist = "dd")

# explore linearity of risk
fit_linear <- rms::cph(
  survival::Surv(time2event_death_2y, event_death_2y) ~ 
    trt * pred_risk,
  data = baseline,
  weights = baseline$sw_IPTW,
  x = TRUE, 
  y = TRUE
)
fit_pol_2 <- rms::cph(
  survival::Surv(time2event_death_2y, event_death_2y) ~ 
    trt * rms::pol(pred_risk, 2),
  data = baseline,
  weights = baseline$sw_IPTW,
  x = TRUE, 
  y = TRUE
)
fit_pol_3 <- rms::cph(
  survival::Surv(time2event_death_2y, event_death_2y) ~ 
    trt * rms::pol(pred_risk, 3),
  data = baseline,
  weights = baseline$sw_IPTW,
  x = TRUE, 
  y = TRUE
)
fit_rcs_3 <- rms::cph(
  survival::Surv(time2event_death_2y, event_death_2y) ~ 
    trt * rms::rcs(pred_risk, 3),
  data = baseline,
  weights = baseline$sw_IPTW,
  x = TRUE, 
  y = TRUE
)
fit_rcs_4 <- rms::cph(
  survival::Surv(time2event_death_2y, event_death_2y) ~ 
    trt * rms::rcs(pred_risk, 4),
  data = baseline,
  weights = baseline$sw_IPTW,
  x = TRUE, 
  y = TRUE
)
BIC(fit_linear)
BIC(fit_pol_2)
BIC(fit_pol_3)
BIC(fit_rcs_3)
BIC(fit_rcs_4)

# compute estimates across 100 risk points
nr_pred <- 100
estimates_df <- compute_HTE(data = baseline,
                            horizon = horizon,
                            nr_pred = nr_pred, 
                            print_test_HTE = TRUE)

# perform bootstrapping to obtain confidence intervals
# n_bootstraps <- 100
RD_boot <- data.frame(pred_risk = estimates_df$pred_risk)
dRMST_boot <- data.frame(pred_risk = estimates_df$pred_risk)
HR_boot <- data.frame(pred_risk = estimates_df$pred_risk)
for (B in 1:n_bootstraps) {
  # create bootstrap sample
  bootstrap <- baseline[sample(1:nrow(baseline), replace = TRUE), ]
  
  # estimate RD, dRMST, and HR on each bootstrap
  estimates_df_boot <- compute_HTE(bootstrap, 
                                         horizon = horizon,
                                         nr_pred = nr_pred)
  
  # save estimates
  RD_boot[, paste0("boot_", B)] <- estimates_df_boot$RD
  dRMST_boot[, paste0("boot_", B)] <- estimates_df_boot$dRMST
  HR_boot[, paste0("boot_", B)] <- estimates_df_boot$HR
}

# compute confidence intervals
estimates_df[, "RD_lower"] <- apply(RD_boot[, -1], 1, quantile, probs = 0.025)
estimates_df[, "RD_upper"] <- apply(RD_boot[, -1], 1, quantile, probs = 0.975)
estimates_df[, "dRMST_lower"] <- apply(dRMST_boot[, -1], 1, quantile, probs = 0.025)
estimates_df[, "dRMST_upper"] <- apply(dRMST_boot[, -1], 1, quantile, probs = 0.975)
estimates_df[, "HR_lower"] <- apply(HR_boot[, -1], 1, quantile, probs = 0.025)
estimates_df[, "HR_upper"] <- apply(HR_boot[, -1], 1, quantile, probs = 0.975)

# Plot effect of risk on hazard ratio
manual_colors <- c("#9161BD", "#2CA02C", "#EC7F12")
RD_plot <- effect_plot(dt = estimates_df,
                       y_middle = 0,
                       measure = "RD")
RD_plot
dRMST_plot <- effect_plot(dt = estimates_df,
                          y_middle = 0,
                          measure = "dRMST")
dRMST_plot
HR_plot <- effect_plot(dt = estimates_df,
                       y_middle = 1,
                       measure = "HR")
HR_plot

# Histogram of predicted risk
hist_predicted_risk <- ggplot2::ggplot(baseline, ggplot2::aes(x = pred_risk)) +
  ggplot2::geom_histogram(binwidth = 0.01,
                          fill = manual_colors[2],
                          color = "white") +
  ggplot2::labs(x = "Predicted 2-year mortailty risk (%)", y = "Count") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    text = ggplot2::element_text(size = 18),
    axis.ticks.y = ggplot2::element_line(color = "black", size = 0.5),
    axis.line.y = ggplot2::element_line(color = "black", size = 0.5)
  )

# Treatment-stratified histogram of predicted risk
hist_predicted_risk_stratified <- ggplot2::ggplot(baseline, ggplot2::aes(
  x = pred_risk,
  fill = as.factor(trt)        # color by trt
)) +
  ggplot2::geom_histogram(
    binwidth = 0.01,
    color = NA,
    position = "dodge"           # side-by-side
  ) +
  ggplot2::labs(x = "Predicted mortality risk (%)", y = "Count") +
  ggplot2::scale_fill_manual(
    values = c("0" = manual_colors[3], "1" = manual_colors[1]),
    labels = c("Conservative", "Dialysis")
  ) +
  ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.1),
                              labels = scales::percent_format(accuracy = 1)) +
  ggplot2::coord_cartesian(xlim = c(0, 1)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.grid = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 18),
    axis.ticks.x = ggplot2::element_line(color = "black", size = 0.5),
    axis.line.x = ggplot2::element_line(color = "black", size = 0.5),
    axis.ticks.y = ggplot2::element_line(color = "black", size = 0.5),
    axis.line.y = ggplot2::element_line(color = "black", size = 0.5),
    legend.position = "bottom",          # move legend below plot
    legend.title = ggplot2::element_blank(),
    legend.direction = "horizontal"      # horizontal layout
  )
hist_predicted_risk_stratified

# save figures
combined_plot <- (RD_plot / hist_predicted_risk_stratified + 
                    ggplot2::theme(legend.position = "none")) +
  plot_layout(heights = c(1, 0.2))  |
  (dRMST_plot / hist_predicted_risk_stratified + 
     ggplot2::theme(legend.text = ggplot2::element_text(size = 18))) +
  plot_layout(heights = c(1, 0.2))  |
  (HR_plot / hist_predicted_risk_stratified + 
     ggplot2::theme(legend.position = "none")) +
  plot_layout(heights = c(1, 0.2))
combined_plot
ggplot2::ggsave(
  plot = combined_plot,
  filename = paste0(results_path, "Supplemental/Figure_S9.png"),
  width = 20,
  height = 10,
  dpi = 300
)

save(
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
  model_PS,
  coef_PS_overall,
  elig_cohort,
  w_meths,
  horizon,
  unit,
  n_bootstraps,
  file = file.path(
    "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/cohort_with_prob.Rdata"
  )
)

################################################################################
### abstract ERA figure
################################################################################
common_tag_theme <- ggplot2::theme(
  plot.tag = ggplot2::element_text(size = 20,
                                   face = "bold")
)
KM_plot <- ggpubr::ggarrange(
  out_KM_event_death_2y_IPTW$KM_plot + 
    ggplot2::labs(title = "", tag = "A") + 
    common_tag_theme + 
    ggplot2::theme(
      legend.position = "none",
      axis.text.y = ggplot2::element_text(size = 16),
      axis.title.y = ggplot2::element_text(
        size = 18,
        margin = ggplot2::margin(r = -10)    # make y-axis label closer
      )
    ),
  out_KM_event_death_2y_IPTW$KM_table + 
    ggplot2::labs(tag = "B") + 
    common_tag_theme,
  nrow = 2, 
  ncol = 1,
  heights = c(7, 3),
  align = "v"
)
HTE_plot <- ggpubr::ggarrange(
  HR_plot + 
    ggplot2::labs(tag = "C") +
    common_tag_theme,
  hist_predicted_risk_stratified + 
    ggplot2::labs(tag = "D") + 
    common_tag_theme + 
    ggplot2::theme(legend.position = "none"),
  nrow = 2, 
  ncol = 1,
  heights = c(7, 3),
  align = "v"
)
ggplot2::ggsave(
  plot = KM_plot | HTE_plot,
  filename = paste0(results_path, "Abstract_ERA.png"),
  width = 12,
  height = 8,
  dpi = 600
)

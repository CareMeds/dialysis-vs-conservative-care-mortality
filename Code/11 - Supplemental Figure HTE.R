################################################################################
### Decision for dialysis versus conservative care
### PART 11 - Illustration relationship between RD, RR, RMST, and HR
################################################################################

# remove history
rm(list = ls(all.names = TRUE))

# Define baseline risk
baseline_risk <- seq(0.1, 0.9, by = 0.01)

tau <- 1  # time horizon (set to match the period your baseline risk covers)

# ---- RMST difference function (exponential survival assumption) ----
dRMST <- function(baseline_risk, hr, tau = 1) {
  lambda      <- -log(1 - baseline_risk) / tau
  rmst_ctrl   <- (1 - exp(-lambda * tau)) / lambda
  lambda_t    <- lambda * hr                          # treated hazard
  rmst_treat  <- (1 - exp(-lambda_t * tau)) / lambda_t
  return(rmst_treat - rmst_ctrl)
}

# ---- CONSTANT HAZARD RATIO (HR = 0.66) ----
hr_const <- 0.66
data_constant_HR <- data.frame(
  baseline_risk = baseline_risk,
  HR            = hr_const,
  RR            = (1 - (1 - baseline_risk)^hr_const) / baseline_risk,
  RD            = (1 - (1 - baseline_risk)^hr_const) - baseline_risk,
  Scenario      = "Constant HR (0.66)"
)
data_constant_HR$dRMST <- dRMST(data_constant_HR$baseline_risk, data_constant_HR$HR, tau)

# ---- INCREASING HAZARD RATIO (0.4 to 0.9) ----
hr_start <- 0.4
hr_end   <- 0.9
data_inc_HR <- data.frame(baseline_risk = baseline_risk)
data_inc_HR$HR       <- hr_start + (hr_end - hr_start) * (data_inc_HR$baseline_risk - min(baseline_risk)) / (max(baseline_risk) - min(baseline_risk))
data_inc_HR$RR       <- (1 - (1 - data_inc_HR$baseline_risk)^data_inc_HR$HR) / data_inc_HR$baseline_risk
data_inc_HR$RD       <- (1 - (1 - data_inc_HR$baseline_risk)^data_inc_HR$HR) - data_inc_HR$baseline_risk
data_inc_HR$Scenario <- "Increasing HR"
data_inc_HR$dRMST <- dRMST(data_inc_HR$baseline_risk, data_inc_HR$HR, tau)

# ---- CONSTANT RISK DIFFERENCE (RD = 0.10) ----
rd_const   <- -0.10
br_filtered <- baseline_risk[baseline_risk > rd_const]
data_constant_RD <- data.frame(
  baseline_risk = br_filtered,
  RD            = rd_const,
  RR            = (br_filtered + rd_const) / br_filtered,
  HR            = log(1 - (br_filtered + rd_const)) / log(1 - br_filtered),
  Scenario      = "Constant RD (-0.10)"
)
data_constant_RD$dRMST <- dRMST(data_constant_RD$baseline_risk, data_constant_RD$HR, tau)

# ---- Plot function ----
plot_metric <- function(df,
                        y_var,
                        y_lab,
                        x_lab = "",
                        y_lim = c(0, 1)) {
  ggplot2::ggplot(df, ggplot2::aes_string(x = "baseline_risk", y = y_var)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    ggplot2::ylim(y_lim) +
    ggplot2::labs(y = y_lab, x = x_lab) +
    ggplot2::theme_classic()
}

# ---- Y-axis limits ----
y_lim_rd   <- c(-0.2, 0)
y_lim_rmst <- c(0.15, 0)

# ---- Row 1: Constant HR ----
p_RD_constant_HR  <- plot_metric(data_constant_HR, "RD", "RD", y_lim = y_lim_rd)
p_RR_constant_HR  <- plot_metric(data_constant_HR, "RR", "RR")
p_dRMST_constant_HR  <- plot_metric(data_constant_HR, "dRMST", "dRMST", y_lim = y_lim_rmst)
p_HR_constant_HR  <- plot_metric(data_constant_HR, "HR", "HR")

# ---- Row 2: Increasing HR ----
p_RD_inc_HR  <- plot_metric(data_inc_HR, "RD", "RD", y_lim = y_lim_rd)
p_RR_inc_HR  <- plot_metric(data_inc_HR, "RR", "RR")
p_dRMST_inc_HR  <- plot_metric(data_inc_HR, "dRMST", "dRMST", y_lim = y_lim_rmst)
p_HR_inc_HR  <- plot_metric(data_inc_HR, "HR", "HR")

# ---- Row 3: Constant RD ----
p_RD_constant_RD  <- plot_metric(data_constant_RD, "RD", "RD", "Baseline Risk", y_lim = y_lim_rd)
p_RR_constant_RD <- plot_metric(data_constant_RD, "RR", "RR", "Baseline Risk")
p_dRMST_constant_RD <- plot_metric(data_constant_RD, "dRMST", "dRMST", "Baseline Risk", y_lim = y_lim_rmst)
p_HR_constant_RD <- plot_metric(data_constant_RD, "HR", "HR", "Baseline Risk")

# ---- Combine into 3x4 grid ----
final_plot <- ggpubr::ggarrange(
  p_RD_constant_HR,
  p_RR_constant_HR,
  p_dRMST_constant_HR,
  p_HR_constant_HR,
  p_RD_inc_HR,
  p_RR_inc_HR,
  p_dRMST_inc_HR,
  p_HR_inc_HR,
  p_RD_constant_RD,
  p_RR_constant_RD,
  p_dRMST_constant_RD,
  p_HR_constant_RD,
  nrow = 3,
  ncol = 4
)

# Warning: first value is NA due to log(0)
data_constant_RD[1, ]

# ---- Save ----
ggplot2::ggsave(
  plot     = final_plot,
  filename = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Results/Supplemental/Figure_HTE_simulation.png",
  width    = 10,
  height   = 5,
  dpi      = 600
)

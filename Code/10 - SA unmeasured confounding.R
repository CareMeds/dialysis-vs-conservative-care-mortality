################################################################################
### Decision for dialysis versus conservative care
### PART 10 - Sensitivity analysis for unmeasured confounding
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
library(patchwork)
set.seed(1)        # set seed for parallel backend

################################################################################
### Load data ##################################################################
################################################################################
load("Data/cohort_with_models.Rdata")

################################################################################
### Adjusted Kaplan-Meier curve
################################################################################
# 1. Run the weighted KM on the baseline dataset in years (time / 365)
KM_fit <- survival::survfit(
  survival::Surv(time2event_death_2y, event_death_2y) ~ trt,
  data = baseline,
  weights = baseline$sw_IPTW
)

# 2. Extract ALL KM time points (do not restrict to eval_times_years)
KM_curve <- summary(KM_fit, times = seq(0, 730, 1))

results_primary <- data.table(
  time_years = KM_curve$time / 365,
  time_days = KM_curve$time,
  trt = as.numeric(KM_curve$strata) - 1,
  surv = KM_curve$surv, 
  lower = KM_curve$lower,  # TODO: bootstrap?
  upper = KM_curve$upper,
  n.censor = KM_curve$n.censor
)

# 3. Constant bias factor
bf_fixed <- 0.55

# 4. Create the 'Adjusted Dialysis' entry using the interpolated bias
adj_dialysis <- results_primary[trt == 1, .(
  time = time_days,
  surv = surv^(1 / bf_fixed),
  lower = lower^(1 / bf_fixed),
  upper = upper^(1 / bf_fixed),
  strata = "Dialysis (Bias-Adjusted)",
  n.censor = n.censor
)]

# 5. Prepare the Observed data
plot_data_long <- results_primary[, .(
  time = time_days,
  surv = surv,
  lower = lower,
  upper = upper,
  strata = factor(
    trt,
    levels = c(0, 1),
    labels = c("Conservative Care", "Dialysis (Observed)")
  ),
  n.censor = n.censor
)]

# 6. Combine and Plot
KM_plot_data_final <- rbind(plot_data_long, adj_dialysis)
KM_plot_data_final[, strata := factor(
  strata,
  levels = c(
    "Conservative Care",
    "Dialysis (Observed)",
    "Dialysis (Bias-Adjusted)"
  )
)]

# censor data
censor_data <- KM_plot_data_final[n.censor>0]

# 7. Plot adjusted KM curve
adj_KM <- ggplot2::ggplot(
  KM_plot_data_final,
  ggplot2::aes(
    x = time,
    y = surv,
    color = strata,
    fill = strata,
    linetype = strata
  )
) +
  ggplot2::geom_step(linewidth = 1) +
  pammtools::geom_stepribbon(
    ggplot2::aes(ymin = lower, ymax = upper),
    alpha = 0.2,
    # Transparency for the shading
    color = NA   # Remove the outline from the ribbon itself
  ) +
  ggplot2::geom_point(data = censor_data,  # add censoring
                      ggplot2::aes(x = time, y = surv),
                      shape = 3, 
                      size = 2,
                      stroke = 0.8,
                      show.legend = FALSE) +
  ggplot2::labs(x = "Time (months)", y = "Survival probability (%)") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(color = "black"),
    axis.ticks = ggplot2::element_line(color = "black"),
    text = ggplot2::element_text(size = 14),
    plot.background = ggplot2::element_rect(fill = "white", color = NA)
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),
    labels = seq(0, 100, by = 10),
    expand = c(0, 0)
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, horizon, by = 365 / 2),
    # 2-year horizon
    labels = function(x)
      round(x / 30) # years -> months
  ) +
  ggplot2::scale_color_manual(values = manual_colors) +
  ggplot2::scale_fill_manual(values = manual_colors) +
  ggplot2::scale_linetype_manual(values = c("solid", "solid", "solid"))
show(adj_KM)

################################################################################
### Risk table
################################################################################
# Define the time points in days
time_points <- seq(0, horizon, by = 365 / 2)

# Extract counts from the raw baseline data
counts_list_C <- lapply(time_points, function(t) {
  # CC group at risk: trt 0 and follow-up >= t
  cc_count <- baseline[trt == 0 & time2event_death_2y >= t, .N]
  
  # Dialysis group at risk: trt 1 and follow-up >= t
  # Note: This is the 'Observed' group, so we only look at death follow-up
  dialysis_count <- baseline[trt == 1 &
                               time2event_death_2y >= t, .N]
  
  data.table(time = t,
             CC = cc_count,
             Dialysis = dialysis_count)
})

# Combine and reshape for plotting
table_data_C <- rbindlist(counts_list_C)
table_data_C_long <- melt(
  table_data_C,
  id.vars = "time",
  variable.name = "strata",
  value.name = "n_at_risk"
)

# Set factor levels to match the KM plot legend
table_data_C_long[, strata := factor(
  strata,
  levels = c("CC", "Dialysis"),
  labels = c("Conservative Care", "Dialysis (Observed)")
)]

# Create the Table Plot for Panel C
# Define your labels with HTML color tags
# This matches the colors to your specific strata manually
labels_with_color <- c(
  glue::glue(
    "<span style='color:{manual_colors[1]};'>Conservative Care</span>"
  ),
  glue::glue(
    "<span style='color:{manual_colors[2]};'>Dialysis (Observed)</span>"
  )
)

# 2. Create the Table Plot
risk_table_C <- ggplot2::ggplot(table_data_C_long,
                                ggplot2::aes(x = time, y = strata, label = n_at_risk)) +
  ggplot2::geom_text(size = 4) +
  ggplot2::scale_x_continuous(limits = c(0, horizon), breaks = time_points) +
  # Map the custom HTML labels to the y-axis
  ggplot2::scale_y_discrete(labels = labels_with_color) +
  ggplot2::theme_void() +
  ggplot2::theme(
    # Use element_markdown to render the HTML/CSS colors
    axis.text.y = ggtext::element_markdown(size = 10, hjust = 1),
    plot.margin = ggplot2::margin(5, 5, 5, 5)
  )
show(risk_table_C)

# save figure
ggplot2::ggsave(
  plot = adj_KM / risk_table_C + patchwork::plot_layout(heights = c(4, 1)),
  filename = paste0(results_path, "Supplemental/Figure_M2_KM.png"),
  width = 7,
  height = 7,
  dpi = 600
)

################################################################################
### 1D plot
################################################################################

# 1. Parameters
bf_fixed <- 0.55
p0_fixed <- 0.5
n <- 100
rr_seq <- seq(2, 10, length.out = n)

# 2. Calculate the required P_C1 to maintain BF = 0.58
# Formula rearranged: P_C1 = (BF * (P0*(RR-1) + 1) - 1) / (RR-1)
p1_required <- (bf_fixed * (p0_fixed * (rr_seq - 1) + 1) - 1) / (rr_seq - 1)

# 3. Open PNG Device
png(
  paste0(
    results_path,
    "Supplemental/Figure_M2_Fixed_Bias_Isoline.png"
  ),
  width = 1200,
  height = 900,
  res = 150
)

# 4. Create 2D Plot
plot(
  rr_seq,
  p1_required,
  type = "l",
  lwd = 3,
  col = "darkred",
  ylim = c(0, 1),
  xlab = "Confounder Strength (RR_CD)",
  ylab = "Prevalence in Dialysis Group (P_C1)",
  main = paste0("Prevalence gap needed to nullify the effect (HR = ", bf_fixed, ")"),
  sub = "Assumes Prevalence in Conservative Care (P_C0) = 0.5"
)

grid(lty = "dotted", col = "gray")
abline(h = 0.5, col = "blue", lty = 2) # Reference line for balance
text(8, 0.52, "P_C0 = 0.5", col = "blue", cex = 0.8)

# 5. Save
dev.off()

################################################################################
### 3D plot
################################################################################

# 1. Parameters
p1_seq <- seq(0, 0.5, length.out = n)

# 2. Calculate Required P_C0
calc_p0 <- function(rr, p1) {
  num <- p1 * (rr - 1) + 1 - bf_fixed
  den <- bf_fixed * (rr - 1)
  p0 <- num / den
  
  # Ensure P_C0 stays within probability bounds [0, 1]
  p0[p0 < 0 | p0 > 1] <- NA
  return(p0)
}

z_p0_matrix <- outer(rr_seq, p1_seq, calc_p0)

# 3. Color Mapping (Green to Yellow to White)
# We map colors to the Z-axis (P_C0) values
nbcol <- n
color_palette <- colorRampPalette(c("blue", "purple", "green4"))(nbcol)

# Calculate color levels for each facet
z_facet <- (z_p0_matrix[-1, -1] + z_p0_matrix[-1, -ncol(z_p0_matrix)] +
              z_p0_matrix[-nrow(z_p0_matrix), -1] + z_p0_matrix[-nrow(z_p0_matrix), -ncol(z_p0_matrix)]) / 4
facet_col <- color_palette[cut(z_facet, nbcol)]

# 4. Save and Plot
png(
  paste0(
    results_path,
    "Supplemental/Figure_M2_Schneeweiss_3D_Plot.png"
  ),
  width = 1200,
  height = 1000,
  res = 150
)

# Adjust margins to fit the legend later if needed
par(mar = c(2, 2, 4, 2))

res <- persp(
  rr_seq,
  p1_seq,
  z_p0_matrix,
  theta = 310,
  phi = 20,
  expand = 0.8,
  col = facet_col,
  lwd = 0.2,
  ticktype = "detailed",
  border = "black",
  xlab = "Confounder-outcome strength (RR_CD)",
  ylab = "Prevalence dialysis (P_C1)",
  zlab = "Prevalence conservative care (P_C0)",
  main = paste0("Sensitivity Analysis: BF = ", bf_fixed),
  cex.main = 1.2,
  cex.axis = 0.8,
  cex.lab = 1
)

dev.off()

################################################################################
### Table
################################################################################
rr_seq <- c(2, 5, 10)
p1_required <- (bf_fixed * (p0_fixed * (rr_seq - 1) + 1) - 1) / (rr_seq - 1)
data.frame(
  RR_CD = rr_seq,
  P_C0 = p0_fixed * 100,
  P_C1 = round(p1_required * 100, 1),
  Gap = round((p0_fixed - p1_required) * 100, 1)
)

KM_plot_data_final[KM_plot_data_final$time == 730, c("strata", "surv")]
cat("Survival probability observed dialysis:",
    (1-as.numeric(KM_plot_data_final[KM_plot_data_final$time == 730 & KM_plot_data_final$strata=="Dialysis (Observed)", "surv"]))*100,
    "\nSurvival probability adjusted dialysis:",
    (1-as.numeric(KM_plot_data_final[KM_plot_data_final$time == 730 & KM_plot_data_final$strata=="Dialysis (Bias-Adjusted)", "surv"]))*100,
    "\nSurvival probability conservative care:",
    (1-as.numeric(KM_plot_data_final[KM_plot_data_final$time == 730 & KM_plot_data_final$strata=="Conservative Care", "surv"]))*100,
    "\nObserved risk difference              :",
    (1-as.numeric(KM_plot_data_final[KM_plot_data_final$time == 730 & KM_plot_data_final$strata=="Dialysis (Observed)", "surv"]))*100-
      (1-as.numeric(KM_plot_data_final[KM_plot_data_final$time == 730 & KM_plot_data_final$strata=="Conservative Care", "surv"]))*100,
    "\nAdjusted risk difference              :",
    (1-as.numeric(KM_plot_data_final[KM_plot_data_final$time == 730 & KM_plot_data_final$strata=="Dialysis (Bias-Adjusted)", "surv"]))*100-
      (1-as.numeric(KM_plot_data_final[KM_plot_data_final$time == 730 & KM_plot_data_final$strata=="Conservative Care", "surv"]))*100)

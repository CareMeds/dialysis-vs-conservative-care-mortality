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
library(tidyverse)
set.seed(1)        # set seed for parallel backend

################################################################################
### 2D plot, fix p0
################################################################################
# Parameters
p0_values       <- c(0.4, 0.6, 0.8, 0.9)
n               <- 100
rr_cd           <- seq(2, 10, length.out = n)
rr_confounded   <- 0.43
rr_unconfounded <- 1

# Single source of truth for factor levels
p0_labels <- paste0("P_C0 = ", p0_values)

# Build long data frame
df <- map_dfr(p0_values, function(p0) {
  p1 <- (rr_confounded / rr_unconfounded * (p0 * (rr_cd - 1) + 1) - 1) / (rr_cd - 1)
  tibble(
    rr_cd       = rr_cd,
    p1_required = p1,
    p0_fixed    = factor(paste0("P_C0 = ", p0), levels = p0_labels)
  )
}) %>%
  filter(p1_required >= 0, p1_required <= 1)

# Reference lines
ref_df <- tibble(
  p0_fixed   = factor(p0_labels, levels = p0_labels),
  yintercept = p0_values
)

# Colours keyed to the same p0_labels vector
my_colours <- setNames(
  c("#1D9E75", "#378ADD", "#D85A30", "#8B0000"),
  p0_labels
)

# Plot
p <- ggplot(df, aes(x = rr_cd, y = p1_required, colour = p0_fixed)) +
  geom_line(linewidth = 1) +
  geom_hline(
    data      = ref_df,
    aes(yintercept = yintercept, colour = p0_fixed),
    linetype  = "dashed",
    linewidth = 0.6
  ) +
  geom_text(
    data        = ref_df,
    aes(x = 9.5, y = yintercept + 0.025, label = as.character(p0_fixed), colour = p0_fixed),
    size        = 3,
    hjust       = 1,
    show.legend = FALSE
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = seq(2, 10, 2)) +
  scale_colour_manual(values = my_colours) +
  labs(
    x      = "Confounder strength (RR_CD)",
    y      = "Prevalence in dialysis group (P_C1)",
    title  = "Prevalence gap needed to nullify the effect",
    colour = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linetype = "dotted", colour = "gray80"),
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 14)
  )
p

ggsave(
  file.path(results_path, "Supplemental/Figure_M2_fix_p0.png"),
  plot   = p,
  width  = 8,
  height = 6,
  dpi    = 150
)

################################################################################
### 2D plot, fix RR_CD
################################################################################
# Parameters
rr_cd_values    <- c(2, 3, 4, 5)          # fixed RR_CD scenarios
n               <- 100
p0              <- seq(0.01, 0.99, length.out = n)
rr_confounded   <- 0.43
rr_unconfounded <- 1

# Single source of truth for factor levels
rr_cd_labels <- paste0("RR_CD = ", rr_cd_values)

# Build long data frame
df2 <- map_dfr(rr_cd_values, function(rcd) {
  p1 <- (rr_confounded / rr_unconfounded * (p0 * (rcd - 1) + 1) - 1) / (rcd - 1)
  tibble(
    p0_fixed    = p0,
    p1_required = p1,
    rr_cd_fixed = factor(paste0("RR_CD = ", rcd), levels = rr_cd_labels)
  )
}) %>%
  filter(p1_required >= 0, p1_required <= 1)

# Reference line: p1 = p0 (perfect balance, no confounding)
ref_line <- tibble(x = c(0, 1), y = c(0, 1))

# Colours keyed to rr_cd_labels
my_colours2 <- setNames(
  c("#1D9E75", "#378ADD", "#D85A30", "#8B0000"),
  rr_cd_labels
)

# Plot
p2 <- ggplot(df2, aes(x = p0_fixed, y = p1_required, colour = rr_cd_fixed)) +
  geom_line(linewidth = 1) +
  geom_abline(
    slope     = 1,
    intercept = 0,
    linetype  = "dashed",
    colour    = "gray50",
    linewidth = 0.5
  ) +
  annotate(
    "text",
    x     = 0.85,
    y     = 0.82,
    label = "P_C1 = P_C0 (no gap)",
    size  = 3,
    colour = "gray50"
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = my_colours2) +
  labs(
    x      = "Prevalence in reference group (P_C0)",
    y      = "Prevalence in dialysis group (P_C1)",
    title  = "Prevalence gap needed to nullify the effect",
    subtitle = paste0("RR_confounded = ", rr_confounded, ", RR_unconfounded = ", rr_unconfounded),
    colour = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linetype = "dotted", colour = "gray80"),
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 14)
  )
p2

ggsave(
  file.path(results_path, "Supplemental/Figure_M3_fix_rr_cd.png"),
  plot   = p2,
  width  = 8,
  height = 6,
  dpi    = 150
)

################################################################################
### 3D plot
################################################################################

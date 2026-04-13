################################################################################
### Decision for dialysis versus conservative care
### PART 10 - Sensitivity analysis for unmeasured confounding
################################################################################

rm(list = ls(all.names = TRUE))
knitr::opts_knit$set(root.dir = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/")

setwd("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/")
results_path <- "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Results/"

library(tidyverse)
set.seed(1)

################################################################################
### Shared parameters
################################################################################
rr_confounded   <- 0.43
rr_unconfounded <- 1
n               <- 100

################################################################################
### 2D plot: fixed p_c0 — confounder strength (RR_CD) vs required p_c1
################################################################################
p_c0            <- c(0.4, 0.6, 0.8, 0.9)
rr_cd           <- seq(2, 10, length.out = n)

p_c0_labels       <- paste0("p_c0 = ", p_c0)
my_colours      <- setNames(c("#1D9E75", "#378ADD", "#D85A30", "#8B0000"), p_c0_labels)

df_fixed_p_c0 <- map_dfr(p_c0, function(p_c0) {
  p_c1 <- (rr_confounded / rr_unconfounded * (p_c0 * (rr_cd - 1) + 1) - 1) / (rr_cd - 1)
  tibble(
    rr_cd       = rr_cd,
    p_c1        = p_c1,
    p_c0_fixed    = factor(paste0("p_c0 = ", p_c0), levels = p_c0_labels)
  )
}) %>%
  filter(p_c1 >= 0, p_c1 <= 1)

ref_df_fixed_p_c0 <- tibble(
  p_c0_fixed   = factor(p_c0_labels, levels = p_c0_labels),
  yintercept = p_c0
)

plot_fixed_p_c0 <- ggplot(df_fixed_p_c0, aes(x = rr_cd, y = p_c1, colour = p_c0_fixed)) +
  geom_line(linewidth = 1) +
  geom_hline(
    data      = ref_df_fixed_p_c0,
    aes(yintercept = yintercept, colour = p_c0_fixed),
    linetype  = "dashed",
    linewidth = 0.6
  ) +
  geom_text(
    data        = ref_df_fixed_p_c0,
    aes(x = 9.5, y = yintercept + 0.025, label = as.character(p_c0_fixed), colour = p_c0_fixed),
    size        = 3,
    hjust       = 1,
    show.legend = FALSE
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = seq(2, 10, 2)) +
  scale_colour_manual(values = my_colours) +
  labs(
    x      = "Confounder strength (RR_CD)",
    y      = "Prevalence in dialysis group (p_c1)",
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
plot_fixed_p_c0

ggsave(
  file.path(results_path, "Supplemental/Figure_M2_fix_p_c0.png"),
  plot   = plot_fixed_p_c0,
  width  = 8,
  height = 6,
  dpi    = 150
)

################################################################################
### 2D plot: fixed RR_CD — reference group prevalence (p_c0) vs required p_c1
################################################################################
rr_cd         <- c(2, 3, 4, 5)
p_c0          <- seq(0.01, 0.99, length.out = n)

rr_cd_labels  <- paste0("RR_CD = ", rr_cd)
my_colours2   <- setNames(c("#1D9E75", "#378ADD", "#D85A30", "#8B0000"), rr_cd_labels)

df_fixed_rr_cd <- map_dfr(rr_cd, function(rcd) {
  p1 <- (rr_confounded / rr_unconfounded * (p_c0 * (rcd - 1) + 1) - 1) / (rcd - 1)
  tibble(
    p_c0_fixed  = p_c0,
    p_c1 = p1,
    rr_cd = factor(paste0("RR_CD = ", rcd), levels = rr_cd_labels)
  )
}) %>%
  filter(p_c1 >= 0, p_c1 <= 1)

plot_fixed_rr_cd <- ggplot(df_fixed_rr_cd, aes(x = p_c0_fixed, y = p_c1, colour = rr_cd)) +
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
    x      = 0.85,
    y      = 0.82,
    label  = "p_c1 = p_c0 (no gap)",
    size   = 3,
    colour = "gray50"
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = my_colours2) +
  labs(
    x        = "Prevalence in reference group (p_c0)",
    y        = "Prevalence in dialysis group (p_c1)",
    title    = "Prevalence gap needed to nullify the effect",
    subtitle = paste0("RR_confounded = ", rr_confounded, ", RR_unconfounded = ", rr_unconfounded),
    colour   = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linetype = "dotted", colour = "gray80"),
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 14)
  )
plot_fixed_rr_cd

ggsave(
  file.path(results_path, "Supplemental/Figure_M2_fix_rr_cd.png"),
  plot   = plot_fixed_rr_cd,
  width  = 8,
  height = 6,
  dpi    = 150
)

################################################################################
### 3D plot
################################################################################
rr_cd       <- seq(2, 10, length.out = n)
p_c1        <- seq(0, 0.5, length.out = n)
p_c0_matrix <- outer(
  rr_cd, p_c1,
  function(rr_cd, p_c1) {
    p_c0 <- (p_c1 * (rr_cd - 1) + 1 - rr_confounded) / (rr_confounded * (rr_cd - 1))
    p_c0[p_c0 < 0 | p_c0 > 1] <- NA
    return(p_c0)
  }
)

# Color mapping: blue -> purple -> green4 along the p_c0 (z) axis
nbcol         <- n
color_palette <- colorRampPalette(c("blue", "purple", "green4"))(nbcol)

z_facet <- (p_c0_matrix[-1, -1] +
              p_c0_matrix[-1, -ncol(p_c0_matrix)] +
              p_c0_matrix[-nrow(p_c0_matrix), -1] +
              p_c0_matrix[-nrow(p_c0_matrix), -ncol(p_c0_matrix)]) / 4
facet_col <- color_palette[cut(z_facet, nbcol)]

png(
  file.path(results_path, paste0("Supplemental/Figure_M2_3D_Plot.png")),
  width  = 1200,
  height = 1000,
  res    = 150
)
par(mar = c(2, 2, 4, 2))
persp(
  rr_cd, p_c0, p_c0_matrix,
  theta     = 310,
  phi       = 20,
  expand    = 0.8,
  col       = facet_col,
  lwd       = 0.2,
  ticktype  = "detailed",
  border    = "black",
  xlab      = "Confounder-outcome strength (RR_CD)",
  ylab      = "Prevalence dialysis (p_c1)",
  zlab      = "Prevalence conservative care (p_c0)",
  main      = "Sensitivity Analysis",
  cex.main  = 1.2,
  cex.axis  = 0.8,
  cex.lab   = 1
)
dev.off()

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
library(patchwork) # combine figures
library(foreach)   # parallel computation
library(doRNG)     # handle parallel seeds
set.seed(1)        # set seed for parallel backend

# load functions
source("Code/utils/plots.R")
source("Code/utils/tables.R")
source("Code/utils/competing_risk.R")
source("Code/utils/data_manipulation.R")

################################################################################
### Load data ##################################################################
################################################################################
load("Data/cohort_with_prob.Rdata")

################################################################################
### 2D plot, fix p0
################################################################################
# Parameters
p0_fixed        <- 0.9 # TODO: 0.4, 0.6, 0.8
n               <- 100
rr_cd           <- seq(2, 10, length.out = n)
rr_confounded   <- 0.43
rr_unconfounded <- 1
p1_required     <- ((rr_confounded / rr_unconfounded) * (p0_fixed * (rr_cd - 1 ) + 1) - 1) / (rr_cd - 1)
  
# Create 2D Plot
png(
  paste0(
    results_path,
    "Supplemental/Figure_M2_fix_p0.png"
  ),
  width = 1200,
  height = 900,
  res = 150
)
plot(
  rr_cd,
  p1_required,
  type = "l",
  lwd = 3,
  col = "darkred",
  ylim = c(0, 1),
  xlab = "Confounder Strength (RR_CD)",
  ylab = "Prevalence in Dialysis Group (P_C1)",
  main = paste0(
    "Prevalence gap needed to nullify the effect" 
  )
)
grid(lty = "dotted", col = "gray")
abline(h = p0_fixed, col = "blue", lty = 2) # Reference line for balance
text(8,
     p0_fixed + 0.02,
     paste("P_C0 =", p0_fixed),
     col = "blue",
     cex = 0.8)
dev.off()

################################################################################
### 3D plot
################################################################################
# Parameters
p1_seq      <- seq(0, 0.5, length.out = n)
p0_seq      <- (p1 * (rr - 1) + 1) / (rr - 1)
z_p0_matrix <- outer(rr_seq, p1_seq, calc_p0)

# Color Mapping (Green to Yellow to White)
# We map colors to the Z-axis (P_C0) values
nbcol         <- n
color_palette <- colorRampPalette(c("blue", "purple", "green4"))(nbcol)

# Calculate color levels for each facet
z_facet   <- (z_p0_matrix[-1, -1] +
                z_p0_matrix[-1, -ncol(z_p0_matrix)] +
                z_p0_matrix[-nrow(z_p0_matrix), -1] +
                z_p0_matrix[-nrow(z_p0_matrix), -ncol(z_p0_matrix)]) / 4
facet_col <- color_palette[cut(z_facet, nbcol)]

# Save and Plot
png(
  paste0(
    results_path,
    "Supplemental/Figure_M2_3D_Plot_",
    landmark,
    ".png"
  ),
  width = 1200,
  height = 1000,
  res = 150
)
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
  main = paste0("Sensitivity Analysis"),
  cex.main = 1.2,
  cex.axis = 0.8,
  cex.lab = 1
)
dev.off()
################################################################################
### Decision for dialysis versus conservative care
### PART - Competing risk analysis for time-to-dialysis
################################################################################

# remove history
rm(list = ls(all.names = TRUE))
knitr::opts_knit$set(root.dir = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/")

# set directory
setwd("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/")
results_path <- "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Results/"

# load libraries
library(mstate)
library(patchwork)
library(foreach)   # parallel computation
library(doRNG)     # handle parallel seeds
set.seed(1)        # set seed for parallel backend

# load functions
source("Code/utils/data_manipulation.R")
source("Code/utils/weighting.R")
source("Code/utils/competing_risk.R")
source("Code/utils/plots.R")
source("Code/utils/compute_absolute_relative_risks.R")

# load data
load("Data/merged_ckd.Rdata")
load("Data/cohort_with_models.Rdata")

# Two patients have their KRT exactly at the time horizon, we can ignore warnings
baseline[LOPNR==287508019, .(LOPNR, 
                             event_KRT_2y, time2event_KRT_2y,
                             event_KRT_inf, time2event_KRT_inf,
                             event_death_2y, time2event_death_2y,
                             event_death_inf, time2event_death_inf)]
baseline[LOPNR==366058177, .(LOPNR, 
                             event_KRT_2y, time2event_KRT_2y,
                             event_KRT_inf, time2event_KRT_inf,
                             event_death_2y, time2event_death_2y,
                             event_death_inf, time2event_death_inf)]

# prepare data with patients who chose dialysis in long format
cmp_dt <- baseline[trt == 1, c(
  id_name,
  "event_KRT_2y",
  "time2event_KRT_2y",
  "event_death_2y",
  "time2event_death_2y",
  "sw_IPTW"
), with = FALSE]

# estimate state probabilities
cmp_result <- state_probabilities(dt = cmp_dt, horizon = horizon)
cmp_result

# create data frame for CI plot
state_prob_df <- data.frame(
  time = rep(cmp_result$times, 3)/365,
  state_prob = c(cmp_result$state_prob$P11, 
                 cmp_result$state_prob$P12, 
                 cmp_result$state_prob$P13),
  Event = rep(c("Event_free",
                "Dialysis", 
                "Death"),
              each = nrow(cmp_result$state_prob))
)

### Use bootstrapping for 95% CI
# 1. Make a bootstrap sample of the full baseline dt (dia and CCC)
# 2. Fit a IPTW model on each bootstrap sample
# 3. Fit Cox model using bootstrapped model and bootstrapped IPTW
n_bootstraps <- 2  # TODO!!

state_prob_boot <- data.frame(time = state_prob_df$time)
RMST_boot <- data.frame(RMST = c("Event_free", "Dialysis", "Death"))
for (B in 1:n_bootstraps){
  # create bootstrap sample
  bootstrap <- baseline[sample(1:nrow(baseline), replace = TRUE), ]
  
  # compute IPTW on bootstrap
  out_weights <- create_weights(
    data = bootstrap,
    model_PS = model_PS,
    w_meth = "IPTW",
    verbose = FALSE
  )
  bootstrap[["sw_IPTW"]] <- as.numeric(out_weights$data$w)
  
  # prepare data with patients who chose dialysis in long format
  cmp_boot_dt <- bootstrap[trt == 1, c(
    id_name,
    "event_KRT_2y",
    "time2event_KRT_2y",
    "event_death_2y",
    "time2event_death_2y",
    "sw_IPTW"
  ), with = FALSE]
  
  # estimate state probabilities
  cmp_boot_result <- state_probabilities(cmp_boot_dt, horizon)
  
  # save state probabilities for each bootstrap
  state_prob_boot[, paste0("boot_", B)] <- 
    c(cmp_boot_result$state_prob$P11, 
      cmp_boot_result$state_prob$P12, 
      cmp_boot_result$state_prob$P13)
  RMST_boot[, paste0("boot_", B)] <-
    unlist(cmp_boot_result$RMST)
}
state_prob_df$state_prob_lower <- apply(state_prob_boot[, -1], 1, quantile, probs = 0.025)
state_prob_df$state_prob_upper <- apply(state_prob_boot[, -1], 1, quantile, probs = 0.975)
RMST_lower <- apply(RMST_boot[, -1], 1, quantile, probs = 0.025)
RMST_upper <- apply(RMST_boot[, -1], 1, quantile, probs = 0.975)

################################################################################
### SANITY CHECK
################################################################################
# RMST_death has to be the same as the KM RMST
fit <- survival::survfit(survival::Surv(time2event_death_2y, 
                                        event_death_2y) ~ 1,
                         data = baseline[trt==1, ])
RMST <- summary(fit, rmean = horizon)$table["rmean"]
# P13 RMST "decision -> death" should be 19.7
RMST/30.5
# 19.71609 
cmp_result$RMST$RMST_13
# 19.27663 

################################################################################
### Create plot
################################################################################
# --- Left: State probability curves ---
state_colors <- c(manual_colors[3], manual_colors[2], manual_colors[1])
p1 <- ggplot2::ggplot(state_prob_df, ggplot2::aes(x = time, 
                                                  y = state_prob, 
                                                  fill = Event,
                                                  color = Event)) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = state_prob_lower,
                                    ymax = state_prob_upper,
                                    group = Event,
                                    fill = Event),
                       alpha = 0.2,
                       colour = NA) +
  ggplot2::scale_fill_manual(
    values = c("Event_free" = state_colors[1],
               "Dialysis" = state_colors[2],
               "Death" = state_colors[3])
  ) +
  ggplot2::scale_color_manual(
    labels = c("Alive dialysis-free",
               "Transition to dialysis",
               "Transition to death"),
    values = c("Event_free" = state_colors[1],
               "Dialysis" = state_colors[2],
               "Death" = state_colors[3]),
    breaks = c("Event_free",
               "Dialysis",
               "Death")
  ) +
  ggplot2::guides( 
    color = ggplot2::guide_legend(override.aes = list(fill = NA)),
    fill = "none"
    ) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 1, by = 0.1)) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    x = "Time (years)",
    y = "Probability (%)",
    color = "Transition"
  ) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 legend.position = "bottom",        
                 legend.direction = "horizontal",   
                 legend.title = ggplot2::element_blank()
                 ) +
  ggplot2::annotate(
    "text",
    x = 0.5,
    y = 1,
    label = paste0(
      "At 2 years\nAlive dialysis-free\nTransition to dialysis\nTransition to death"),
    hjust = 0,
    vjust = 1,
    size = 3
  ) +
  ggplot2::annotate(
    "text",
    x = 1.1,
    y = 1,
    label = paste0(
      "%, 95% CI\n",
      fmt(as.numeric(state_prob_df |> 
                         dplyr::filter(Event == "Event_free" & 
                                         time == horizon/365) |>
                         dplyr::select(state_prob)) * 100, 1),
      " (",
      fmt(as.numeric(state_prob_df |> 
                         dplyr::filter(Event == "Event_free" & 
                                         time == horizon/365) |>
                         dplyr::select(state_prob_lower)) * 100, 1),
      ", ",
      fmt(as.numeric(state_prob_df |> 
                         dplyr::filter(Event == "Event_free" & 
                                         time == horizon/365) |>
                         dplyr::select(state_prob_upper)) * 100, 1),
      ")\n",
      fmt(as.numeric(state_prob_df |> 
                         dplyr::filter(Event == "Dialysis" & 
                                         time == horizon/365) |>
                         dplyr::select(state_prob)) * 100, 1),
      " (",
      fmt(as.numeric(state_prob_df |> 
                         dplyr::filter(Event == "Dialysis" & 
                                         time == horizon/365) |>
                         dplyr::select(state_prob_lower)) * 100, 1),
      ", ",
      fmt(as.numeric(state_prob_df |> 
                         dplyr::filter(Event == "Dialysis" & 
                                         time == horizon/365) |>
                         dplyr::select(state_prob_upper)) * 100, 1),
      ")\n",
      fmt(as.numeric(state_prob_df |> 
                         dplyr::filter(Event == "Death" & 
                                         time == horizon/365) |>
                         dplyr::select(state_prob)) * 100, 1),
      " (",
      fmt(as.numeric(state_prob_df |> 
                         dplyr::filter(Event == "Death" & 
                                         time == horizon/365) |>
                         dplyr::select(state_prob_lower)) * 100, 1),
      ", ",
      fmt(as.numeric(state_prob_df |> 
                         dplyr::filter(Event == "Death" & 
                                         time == horizon/365) |>
                         dplyr::select(state_prob_upper)) * 100, 1),
      ")"
    ),
    hjust = 0,
    vjust = 1,
    size = 3
  ) + 
  ggplot2::annotate(
    "text",
    x = 1.6,
    y = 1,
    label = paste0(
      "months, 95% CI\n",
      fmt(cmp_result$RMST$RMST_11, 1),
      " (",
      fmt(RMST_lower[1], 1),
      ", ",
      fmt(RMST_upper[1], 1),
      ")\n",
      fmt(cmp_result$RMST$RMST_12, 1),
      " (",
      fmt(RMST_lower[2], 1),
      ", ",
      fmt(RMST_upper[2], 1),
      ")\n",
      fmt(cmp_result$RMST$RMST_13, 1),
      " (",
      fmt(RMST_lower[3], 1),
      ", ",
      fmt(RMST_upper[3], 1),
      ")"
    ),
    hjust = 0,
    vjust = 1,
    size = 3
  )

# --- Right: Stacked state probabilities ---
# Long format for stacked plot
df_long <- state_prob_df[, 1:3] |>
  tidyr::pivot_longer(
    cols = state_prob,
    names_to = "State",
    values_to = "Probability"
  ) |>
  dplyr::mutate(Event = factor(Event,
                               levels = c("Dialysis", "Event_free", "Death")))
p2 <- ggplot2::ggplot(df_long, ggplot2::aes(x = time, y = Probability, fill = Event)) +
  ggplot2::geom_area(color = "black",
                     linewidth = 0.2,
                     alpha = 0.8) +
  ggplot2::scale_fill_manual(values = c(
    "Event_free" = state_colors[1],
    "Dialysis" = state_colors[2],
    "Death" = state_colors[3]
  ),
  labels = c("Transition to dialysis",
             "Alive dialysis-free",
             "Transition to death")
  ) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 1, by = 0.1)) +
  ggplot2::labs(x = "Time (years)") +
  ggplot2::theme_minimal() +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 legend.position = "bottom",          # move legend below plot
                 legend.direction = "horizontal",      # horizontal layout
                 legend.title = ggplot2::element_blank()
                 )

# --- Arrange plots side by side ---
ggplot2::ggsave(
  plot = gridExtra::grid.arrange(p1, p2, ncol = 2),
  filename = paste0(results_path, "Main/Figure_2.png"),
  width = 10,
  height = 5,
  dpi = 300
)

### Median follow-up time
# Extract at which time CI == 0.5
# Median time from decision to dialysis is 1.3 years, while you do not die
# Probability of having started dialysis in the next 1.3 year is 0.5%
# Probability of having died in the next 1.3 year is 0.15%

################################################################################
### Sensitivity analysis for confounding
################################################################################
KM_fit <- survival::survfit(
  survival::Surv(time2event_death_2y / 365, event_death_2y) ~ trt,
  data = baseline,
  weights = baseline$sw_IPTW
)
KM_curve <- summary(KM_fit, times = unique(state_prob_df$time))
S_dia <- data.table::data.table(time = KM_curve$time[KM_curve$strata=="trt=1"],
                                surv = KM_curve$surv[KM_curve$strata=="trt=1"],
                                lower = KM_curve$lower[KM_curve$strata=="trt=1"],
                                upper = KM_curve$upper[KM_curve$strata=="trt=1"],
                                RR_confounded = (1-KM_curve$surv[KM_curve$strata=="trt=1"]) / 
                                  (1-KM_curve$surv[KM_curve$strata=="trt=0"]))
p_dia <- state_prob_df[state_prob_df$Event=="Dialysis", 
                       c("time", "state_prob")]
S_dia <- merge(S_dia,
               p_dia,
               by = "time",
               all.x = TRUE)
S_dia[, RR_unconfounded := 1 - state_prob]
S_dia[, BF := RR_confounded / RR_unconfounded]
head(S_dia)

# check
print("Transition probability at 6 months")
print(S_dia[time > 0.498 & time < 0.501, ])

# update survival of dialysis to non-starters
KM_curve$surv[KM_curve$strata=="trt=1"] <- S_dia$surv ^ (1/S_dia$BF)
KM_curve$lower[KM_curve$strata=="trt=1"] <- S_dia$lower ^ (1/S_dia$BF)
KM_curve$upper[KM_curve$strata=="trt=1"] <- S_dia$upper ^ (1/S_dia$BF)

# make KM plot
KM_plot_data <- data.frame(time = KM_curve$time,
                           strata = KM_curve$strata,
                           surv = KM_curve$surv,
                           lower = KM_curve$lower,
                           upper = KM_curve$upper,
                           n.censor = KM_curve$n.censor)
censor_data <- KM_plot_data[KM_plot_data$n.censor > 0,]
adj_KM <- ggplot2::ggplot(KM_plot_data, 
                ggplot2::aes(x = time,
                             y = surv,
                             fill = strata,
                             color = strata)) +
  ggplot2::geom_step(linewidth = 1) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                       alpha = 0.2,
                       color = NA) +
  ggplot2::geom_point(data = censor_data,
                      ggplot2::aes(x = time, y = surv),
                      shape = 3, 
                      size = 2,
                      stroke = 0.8) +
  ggplot2::labs(x = "Time (months)",
                y = "Survival probability (%)") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "none",
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(color = "black"),
    axis.ticks = ggplot2::element_line(color = "black"),
    text = ggplot2::element_text(size = 14),
    plot.background = ggplot2::element_rect(fill = "white")
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),
    labels = seq(0, 100, by = 10),
    expand = c(0, 0)     # remove padding below 0
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, horizon, by = 0.5),
    labels = function(x)
      round(x * 12)            # convert years → months
  ) +
  ggplot2::scale_fill_manual(values = manual_colors[1:2])
show(adj_KM)

ggplot2::ggsave(
  plot = adj_KM,
  filename = paste0(results_path, "Supplemental/Figure_M2.png"),
  width = 5,
  height = 5,
  dpi = 600
)

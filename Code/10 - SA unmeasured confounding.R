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
source("Code/utils/tables.R")
source("Code/utils/competing_risk.R")
source("Code/utils/data_manipulation.R")

################################################################################
### Load data ##################################################################
################################################################################
load("Data/cohort_with_prob.Rdata")

################################################################################
### Perform sensitivity analysis at three landmark points
################################################################################
n_bootstraps <- 2
# extract risks at one, three, and six months
landmarks <- c(30, 90, 183)
bf_dt <- data.frame(
  "Bootstrap" = numeric(),
  "Landmark" = numeric(),
  "Risk - conservative care (landmark, %)" = numeric(),
  "Risk - dialysis (landmark, observed, %)" = numeric(),
  "Confounded RR" = numeric(),
  "Probability of initiating dialysis (%)" = numeric(),
  "No. patients initiating dialysis" = numeric(),
  "Risk - initiators (%)" = numeric(),
  "Unconfounded risk - dialysis (%)" = numeric(),
  "Unconfounded RR" = numeric(),
  "Bias factor" = numeric(),
  "Risk - conservative care (%)" = numeric(),
  "Risk - dialysis (observed, %)" = numeric(),
  "Risk - adjusted (adjusted, %)" = numeric(),
  "Risk difference (observed, %)" = numeric(),
  "Risk difference (adjusted, %)" = numeric(),
  "Risk ratio (observed)" = numeric(),
  "Risk ratio (adjusted)" = numeric(),
  "RMST - conservative care (months)" = numeric(),
  "RMST - dialysis (observed, months)" = numeric(),
  "RMST - dialysis (adjusted, months)" = numeric(),
  "dRMST (observed, months)" = numeric(),
  "dRMST (adjusted, months)" = numeric(),
  check.names = FALSE
)
QBA_examples <- c()
for (B in 1:(n_bootstraps + 1)) {
  if (B == 1) {
    # Original sample
    bootstrap <- baseline
  } else{
    # create bootstrap sample
    bootstrap <- baseline[sample(1:nrow(baseline), replace = TRUE), ]
  }
  
  # run for each landmark separately
  for (landmark in landmarks) {
    # 1. Run the weighted KM on the baseline dataset in years (time / 365) ------#
    KM_fit <- survival::survfit(
      survival::Surv(time2event_death_2y, event_death_2y) ~ trt,
      data = bootstrap,
      weights = bootstrap$sw_IPTW
    )
    
    # 2. Extract risks at landmark ----------------------------------------------#
    KM_landmark   <- summary(KM_fit, times = landmark)
    R_CC          <- 1 - KM_landmark$surv[KM_landmark$strata == "trt=0"]
    R_dialysis    <- 1 - KM_landmark$surv[KM_landmark$strata == "trt=1"]
    RR_confounded <- R_dialysis / R_CC
    
    # 3. Estimate Probability of having initiated dialysis at landmark ----------#
    vars <- c(
      id_name,
      "event_KRT_inf",
      "time2event_KRT_inf",
      "event_death_inf",
      "time2event_death_inf",
      "sw_IPTW"
    )
    cmp_dt         <- bootstrap[trt == 1, ..vars, with = FALSE]
    state_prob_out <- state_probabilities(dt = cmp_dt, horizon = landmark)
    p_initiated    <- tail(state_prob_out$state_prob$P12, 1)
    
    # 4. Estimate risk of those who initiated dialysis --------------------------#
    # 4.1 select those that initially chose dialysis
    # 4.2 select those who started KRT within landmark period
    bootstrap[, initiated := trt == 1 &
                event_KRT_inf == 1 &
                time2event_KRT_inf <= landmark]
    initiated_patients <- bootstrap[(initiated)]
    N_initiated <- nrow(initiated_patients)
    
    # obtain weights P(intiatiated before or at landmark)
    dialysis_df <- bootstrap[trt == 1, ]
    model_initiated <- glm(
      update(model_PS, "initiated ~ ."),
      family = stats::binomial(),
      data = dialysis_df
    )
    dialysis_df$w_initiated <- 1 / predict(model_initiated, newdata = dialysis_df, type = "response")
    initiated_patients$w_initiated <- 1 / predict(model_initiated, newdata = initiated_patients, type = "response")
    
    # describe patients who initiated and those who did not
    if (B == 1) {
      table_one <- create_baseline_table(
        data = dialysis_df,
        id_name = "LOPNR",
        vars = listvar,
        categoricalVars = catvar,
        IQRVars = non_normal_vars,
        treatmentColumn = "initiated",
        treatmentLabel = "Initiated",
        controlLabel = "Not initiated",
        tableCaption = paste0("Stratified by initiation at landmark", landmark)
      )
      openxlsx::write.xlsx(
        table_one$raw_table,
        rowNames = TRUE,
        file = paste0(
          results_path,
          "Other/Descriptives_initiated_at_landmark_",
          landmark,
          ".xlsx"
        )
      )
      table_one_weighted <- create_baseline_table(
        data = dialysis_df,
        weights = dialysis_df$sw_IPTW * dialysis_df$w_initiated,
        id_name = "LOPNR",
        vars = listvar,
        categoricalVars = catvar,
        IQRVars = non_normal_vars,
        treatmentColumn = "initiated",
        treatmentLabel = "Initiated",
        controlLabel = "Not initiated",
        tableCaption = paste0("Stratified by initiation at landmark", landmark)
      )
      openxlsx::write.xlsx(
        table_one_weighted$raw_table,
        rowNames = TRUE,
        file = paste0(
          results_path,
          "Other/Descriptives_initiated_at_landmark_",
          landmark,
          "_weighted.xlsx"
        )
      )
    }
    
    # 4.3 obtain survival probability at landmark
    KM_initiated <- survival::survfit(
      survival::Surv(time2event_death_2y, event_death_2y) ~ 1,
      data = initiated_patients,
      weights = sw_IPTW * w_initiated
    )
    R_initiated <- 1 - summary(KM_initiated, times = landmark)$surv
    
    # 5. Estimate unconfounded risk of dialysis, risk ratio and bias factor
    R_dialysis_unconfounded <- p_initiated * R_initiated + (1 - p_initiated) * R_CC
    RR_unconfounded         <- R_dialysis_unconfounded / R_CC
    bias_factor             <- RR_confounded / RR_unconfounded
    
    if (B == 1) {
      # 6. Unadjusted Kaplan-Meier curve ------------------------------------------#
      # Extract ALL KM time points (do not restrict to eval_times_years)
      KM_curve <- summary(KM_fit, times = seq(0, horizon, 1))
      results_primary <- data.table(
        time_years = KM_curve$time / 365,
        time_days = KM_curve$time,
        trt = as.numeric(KM_curve$strata) - 1,
        surv = KM_curve$surv,
        lower = KM_curve$lower,
        upper = KM_curve$upper,
        n.censor = KM_curve$n.censor
      )
      obs_KM <- results_primary[, .(
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
      
      # 7. Adjusted KM curve ------------------------------------------------------#
      adj_dialysis <- results_primary[trt == 1, .(
        time = time_days,
        surv = surv^(1 / bias_factor),
        lower = lower^(1 / bias_factor),
        upper = upper^(1 / bias_factor),
        strata = "Dialysis (Bias-Adjusted)",
        n.censor = n.censor
      )]
      
      # Combine plot data
      KM_plot_data_final <- rbind(obs_KM, adj_dialysis)
      KM_plot_data_final[, strata := factor(
        strata,
        levels = c(
          "Conservative Care",
          "Dialysis (Observed)",
          "Dialysis (Bias-Adjusted)"
        )
      )]
      
      # Censor data
      censor_data <- KM_plot_data_final[n.censor > 0]
      
      # Plot adjusted KM curve
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
        ggplot2::geom_point(
          data = censor_data,
          # add censoring
          ggplot2::aes(x = time, y = surv),
          shape = 3,
          size = 2,
          stroke = 0.8,
          show.legend = FALSE
        ) +
        ggplot2::labs(
          x = "Time (months)",
          y = "Survival probability (%)",
          title = paste0(
            "Bias factor of ",
            round(bias_factor, 2),
            " determined at ",
            round(landmark / 30),
            " months"
          )
        ) +
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
      assign(paste0("adj_KM_", landmark), adj_KM)
      
      # 8. Risk table --------------------------------------------------------------#
      # Define the time points in days
      time_points <- seq(0, horizon, by = 365 / 2)
      
      # Extract counts from the raw bootstrap data
      count_list <- lapply(time_points, function(t) {
        # CC group at risk: trt 0 and follow-up >= t
        cc_count <- bootstrap[trt == 0 &
                                time2event_death_2y >= t, .N]
        
        # Dialysis group at risk: trt 1 and follow-up >= t
        # Note: This is the 'Observed' group, so we only look at death follow-up
        dialysis_count <- bootstrap[trt == 1 &
                                      time2event_death_2y >= t, .N]
        
        data.table(time = t,
                   CC = cc_count,
                   Dialysis = dialysis_count)
      })
      
      # Combine and reshape for plotting
      risk_table      <- rbindlist(count_list)
      risk_table_long <- melt(
        risk_table,
        id.vars = "time",
        variable.name = "strata",
        value.name = "n_at_risk"
      )
      
      # Set factor levels to match the KM plot legend
      risk_table_long[, strata := factor(
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
      
      # Create the Table Plot
      risk_table_plot <- ggplot2::ggplot(risk_table_long,
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
      
      # 9. 1D plot ----------------------------------------------------------------#
      # Parameters
      p0_fixed <- 0.9
      n        <- 100
      rr_seq   <- seq(2, 10, length.out = n)
      
      # Calculate the required P_C1 to maintain BF = 0.58
      # Formula rearranged: P_C1 = (BF * (P0*(RR-1) + 1) - 1) / (RR-1)
      p1_required <- (bias_factor * (p0_fixed * (rr_seq - 1) + 1) - 1) / (rr_seq - 1)
      
      # Create 2D Plot
      png(
        paste0(
          results_path,
          "Supplemental/Figure_M2_1D_Plot_",
          landmark,
          ".png"
        ),
        width = 1200,
        height = 900,
        res = 150
      )
      plot(
        rr_seq,
        p1_required,
        type = "l",
        lwd = 3,
        col = "darkred",
        ylim = c(0, 1),
        xlab = "Confounder Strength (RR_CD)",
        ylab = "Prevalence in Dialysis Group (P_C1)",
        main = paste0(
          "Prevalence gap needed to nullify the effect (HR = ",
          round(bias_factor, 2),
          ")"
        ),
        sub = "Assumes Prevalence in Conservative Care (P_C0) = 0.5"
      )
      grid(lty = "dotted", col = "gray")
      abline(h = p0_fixed,
             col = "blue",
             lty = 2) # Reference line for balance
      text(
        8,
        p0_fixed + 0.02,
        paste("P_C0 =", p0_fixed),
        col = "blue",
        cex = 0.8
      )
      dev.off()
      
      # 10. 3D plot ---------------------------------------------------------------#
      # Parameters
      p1_seq      <- seq(0, 0.5, length.out = n)
      z_p0_matrix <- outer(rr_seq, p1_seq, calc_p0)
      # TODO: direct formula for p0
      
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
        main = paste0("Sensitivity Analysis: BF = ", round(bias_factor, 2)),
        cex.main = 1.2,
        cex.axis = 0.8,
        cex.lab = 1
      )
      dev.off()
      
      # 11. Examples of combinations of parameters --------------------------------#
      rr_seq <- rep(c(2, 3), each = 11)
      p0_seq <- rep(seq(0, 1, 0.1), 2)
      p1_required <- (bias_factor * (p0_seq * (rr_seq - 1) + 1) - 1) / (rr_seq - 1)
      prevalence_gap <- (p0_seq - p1_required) * 100
      QBA_example <- data.frame(
        landmark = round(landmark / 30),
        RR_confounded = sprintf("%.2f", RR_confounded),
        bias_factor = sprintf("%.2f", bias_factor),
        RR_unconfounded = sprintf("%.2f", RR_unconfounded),
        RR_CD = rr_seq,
        P_C0 = paste0(p0_seq * 100, "%"),
        P_C1 = paste0(
          ifelse(p1_required < 0, "Impossible (", ""),
          paste0(sprintf("%.1f", p1_required * 100), "%"),
          ifelse(p1_required < 0, ")", "")
        ),
        Prevalence_gap = ifelse(p1_required < 0, "", paste0(
          sprintf("%.1f", prevalence_gap), "%"
        ))
      )
      QBA_examples <- rbind(QBA_examples, QBA_example)
    }
    
    # Adjusted RMST -----------------------------------------------------------#
    # Extract full KM curve at event times only (not seq)
    KM_curve <- summary(KM_fit, times = seq(0, horizon, 1))
    
    # Step-function RMST = sum of surv[i] * (time[i+1] - time[i])
    rmst_step <- function(times, surv, horizon) {
      # Prepend time 0 with surv = 1
      times <- c(0, times)
      surv  <- c(1, surv)
      
      # Step function: surv[i] held constant until times[i+1]
      dt <- diff(times)
      sum(surv[-length(surv)] * dt)
    }
    
    # set unit
    if (unit == "years") {
      div.fact <- 365.25
    } else if (unit == "months") {
      div.fact <- 30.5
    } else if (unit == "days") {
      div.fact <- 1
    }
    
    # Adjusted
    surv_dial_obs <- KM_curve$surv[KM_curve$strata == "trt=1"]
    surv_dial_adj <- surv_dial_obs^(1 / bias_factor)
    times_dial    <- KM_curve$time[KM_curve$strata == "trt=1"]
    RMST_dial_adj <- rmst_step(times_dial, surv_dial_adj, horizon) / div.fact
    
    # Extract KM results at horizon -------------------------------------------#
    KM_horizon <- summary(KM_fit, times = horizon)
    risk_CC <- 1 - KM_horizon$surv[KM_horizon$strata == "trt=0"]
    risk_dia <- 1 - KM_horizon$surv[KM_horizon$strata == "trt=1"]
    risk_dia_adj <- 1 - KM_horizon$surv[KM_horizon$strata == "trt=1"]^(1 / bias_factor)
    RMST_horizon <- summary(KM_fit, rmean = horizon)$table[, "rmean"] / div.fact
    
    # save results ------------------------------------------------------------#
    summarized_results <- data.frame(
      Bootstrap = B - 1,
      Landmark = landmark,
      `Risk - conservative care (landmark, %)` = R_CC * 100,
      `Risk - dialysis (landmark, observed, %)` = R_dialysis * 100,
      `Confounded RR` = RR_confounded,
      `Probability of initiating dialysis (%)` = p_initiated * 100,
      `No. patients initiating dialysis` = N_initiated,
      `Risk - initiators (%)` = R_initiated * 100,
      `Unconfounded risk - dialysis (%)` = R_dialysis_unconfounded * 100,
      `Unconfounded RR` = RR_unconfounded,
      `Bias factor` = bias_factor,
      `Risk - conservative care (%)` = risk_CC * 100,
      `Risk - dialysis (observed, %)` = risk_dia * 100,
      `Risk - adjusted (adjusted, %)` = risk_dia_adj * 100,
      `Risk difference (observed, %)` = (risk_dia - risk_CC) * 100,
      `Risk difference (adjusted, %)` = (risk_dia_adj - risk_CC) * 100,
      `Risk ratio (observed)` = risk_dia / risk_CC,
      `Risk ratio (adjusted)` = risk_dia_adj / risk_CC,
      `RMST - conservative care (months)` = RMST_horizon["trt=0"],
      `RMST - dialysis (observed, months)` = RMST_horizon["trt=1"],
      `RMST - dialysis (adjusted, months)` = RMST_dial_adj,
      `dRMST (observed, months)` = RMST_horizon["trt=1"] - RMST_horizon["trt=0"],
      `dRMST (adjusted, months)` = RMST_dial_adj - RMST_horizon["trt=0"],
      check.names = FALSE
    )
    bf_dt <- rbind(bf_dt, summarized_results)
    combined_plot <- adj_KM / risk_table_plot + patchwork::plot_layout(heights = c(4, 1))
    assign(
      paste0("combined_plot_", landmark),
      adj_KM / risk_table_plot + patchwork::plot_layout(heights = c(4, 1))
    )
  }
  
  if (B == 1) {
    openxlsx::write.xlsx(QBA_examples,
                         file = paste0(results_path, "Supplemental/Table_M2_QBA_examples.xlsx"))
    ggplot2::ggsave(
      plot = combined_plot_30 | combined_plot_90 | combined_plot_183,
      filename = paste0(results_path, "Supplemental/Figure_M2_KM.png"),
      width = 21,
      height = 7,
      dpi = 600
    )
  }
}

# Point estimates (original sample)
point_estimates <- bf_dt |> dplyr::filter(Bootstrap == 0)

# Bootstrap CIs (percentile method)
boot_ci <- bf_dt |>
  dplyr::filter(Bootstrap > 0) |>
  dplyr::group_by(Landmark) |>
  dplyr::summarise(dplyr::across(
    where(is.numeric) & !c(Bootstrap),
    list(
      lower = \(x) quantile(x, 0.025, na.rm = TRUE),
      upper = \(x) quantile(x, 0.975, na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  ))

# Format as "estimate (lower, upper)" using fmt_ci
measures <- bf_dt |> dplyr::select(where(is.numeric) &
                                     !c(Bootstrap, Landmark)) |> names()
combined <- point_estimates |> dplyr::left_join(boot_ci, by = "Landmark")
digits_map <- c(
  "Confounded RR" = 2,
  "Unconfounded RR" = 2,
  "Bias factor" = 2,
  "No. patients initiating dialysis" = 0
)
combined_results <- combined |>
  dplyr::mutate(dplyr::across(
    dplyr::all_of(measures),
    \(x) fmt_ci(x, combined[[paste0(dplyr::cur_column(), "_lower")]], combined[[paste0(dplyr::cur_column(), "_upper")]], digits = dplyr::coalesce(digits_map[dplyr::cur_column()], 1L))
  )) |>
  dplyr::select(-dplyr::ends_with("_lower"), -dplyr::ends_with("_upper"))
final_table <- t(combined_results)
final_table_with_headers <- rbind(
  final_table[2, ],
  rep("", 3),
  # short-term mortality risk
  final_table[3:5, ],
  rep("", 3),
  # dialysis initiators at landmark
  final_table[6:8, ],
  rep("", 3),
  # bias adjustment
  final_table[9:11, ],
  rep("", 3),
  # 2-year outcomes
  final_table[12:23, ]
)

# save to table
openxlsx::write.xlsx(
  final_table_with_headers,
  rowNames = TRUE,
  colNames = FALSE,
  file = paste0(results_path, "Supplemental/Table_M2_QBA_with_CI.xlsx")
)

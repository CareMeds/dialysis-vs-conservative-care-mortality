create_ps_distribution_plot <- function(data,
                                        PS_varname,
                                        trt_varname,
                                        weights = NULL,
                                        PS_title_hist,
                                        PS_title_scaled_hist,
                                        titleSize = 22,
                                        TextSize = 26,
                                        xlab = TRUE,
                                        ylab = TRUE,
                                        palette = c("blue", "darkgreen"),
                                        x_axis_text = "Propensity score") {
  # extract ps and trt
  data <- copy(data)
  data$ps <- data[, PS_varname, with = FALSE][[1]]
  data$trt <- data[, trt_varname, with = FALSE][[1]]
  
  # calculate unadjusted AUC
  AUC <- pROC::auc(pROC::roc(
    predictor = data$ps,
    response = data$trt,
    quiet = TRUE
  ))
  # AUC <- WeightedROC::WeightedAUC(WeightedROC::WeightedROC(guess=data$ps,
  #                                                          label=data$trt,
  #                                                          weight=weights))
  
  # create histogram for propensity scores
  # !! added to make sure to use weights from data
  hist <- ggplot2::ggplot(data, ggplot2::aes(
    x = ps,
    fill = as.factor(trt),
    weight = !!weights
  )) +
    ggplot2::geom_histogram(
      alpha = 0.6,
      position = "dodge",
      binwidth = 0.01,
      boundary = 0
    ) +
    ggplot2::scale_x_continuous(x_axis_text, limits = c(0, 1)) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::theme(
      text = ggplot2::element_text(size = TextSize),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = titleSize
      )
    ) +
    ggplot2::ggtitle(PS_title_hist)
  
  # create scaled density for propensity scores
  scaled_hist <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = ps,
      fill = as.factor(trt),
      weight = !!weights,
      y = ggplot2::after_stat(scaled)
    )
  ) +
    ggplot2::geom_density(alpha = 0.6, bw = 0.05) +
    ggplot2::scale_x_continuous(x_axis_text, limits = c(0, 1)) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::theme(
      text = ggplot2::element_text(size = TextSize),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = titleSize
      )
    ) +
    ggplot2::ggtitle(PS_title_scaled_hist)
  
  if (!xlab) {
    hist <- hist +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
    scaled_hist <- scaled_hist +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }
  
  if (!ylab) {
    hist <- hist +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
    scaled_hist <- scaled_hist +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  }
  
  return(list(
    hist = hist,
    scaled_hist = scaled_hist,
    AUC = AUC
  ))
}

create_age_distribution_plot <- function(data,
                                         # containing ps and trt
                                         title,
                                         treatmentLabel,
                                         controlLabel,
                                         titleSize = 26,
                                         TextSize = 26,
                                         palette = c("blue", "darkgreen"),
                                         weights = NULL) {
  # create histogram for age
  hist <- ggplot2::ggplot(data, ggplot2::aes(
    x = age,
    fill = factor(trt),
    weight = !!weights
  )) + # !! added to make sure to use weights from data
    ggplot2::geom_histogram(alpha = 0.6,
                            position = "dodge",
                            binwidth = 1) +
    ggplot2::scale_fill_manual(
      name = "Treatment",
      values = palette,
      breaks = c("0", "1"),
      labels = c(controlLabel, treatmentLabel)
    ) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(
      text = ggplot2::element_text(size = TextSize),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      legend.justification = c(1, 1),
      legend.position = "top",
      legend.key = ggplot2::element_rect(colour = NA),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = titleSize
      )
    )
  
  return(hist)
}

create_histogram <- function(data,
                             title,
                             x_label = "Months",
                             y_label = "Frequency",
                             TextSize = 26) {
  # Remove NA from data
  data <- data[!is.na(data)]
  
  # --- Define bin width as 1 month (≈30.44 days) ---
  max_days <- max(data)
  month_bins <- seq(0, max_days + 30.44, by = 30.44) # extend by one month to include max
  
  # Create histogram manually using monthly bins
  h <- hist(data, breaks = month_bins, plot = FALSE)
  
  # Put into df format to work with ggplot
  h_df <- data.frame(
    breaks = h$breaks[-length(h$breaks)],
    # left edges of bins
    counts = h$counts
  )
  
  # --- Create ggplot ---
  h_plot <- ggplot2::ggplot(h_df, ggplot2::aes(x = breaks, y = counts)) +
    ggplot2::geom_bar(stat = "identity",
                      fill = "skyblue",
                      color = "black") +
    ggplot2::labs(title = title, x = x_label, y = y_label) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, max(h_df$breaks), by = 12 * 30.44),
      # tick every 6 months
      labels = function(x)
        round(x / 365.25)                                  # label in months
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text = ggplot2::element_text(size = TextSize),
      axis.text = ggplot2::element_text(size = TextSize)
    )
  return(h_plot)
}

# create love plot to compare SMD before and after weighting
love_plot <- function(SMDs_dt,
                      SMD_names,
                      plot_title,
                      plotColors,
                      xlab_title = "Standardized mean difference",
                      xmax,
                      titleSize = 22,
                      TextSize = 26) {
  # make one love plot
  love.plot <- ggplot2::ggplot(data = SMDs_dt, ggplot2::aes(y = factor(rownames(SMDs_dt), levels = rev(
    rownames(SMDs_dt)
  )))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid") +
    ggplot2::geom_vline(xintercept = 0.1, linetype = "dashed") +
    ggplot2::geom_point(x = SMDs_dt[, SMD_names[1]],
                        colour = plotColors[1],
                        size = 4) +
    ggplot2::geom_point(x = SMDs_dt[, SMD_names[2]],
                        colour = plotColors[2],
                        size = 4) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::xlab(xlab_title) +
    ggplot2::ylab("") +
    ggplot2::scale_x_continuous(limits = c(0, max(max(SMDs_dt), xmax)), breaks =
                                  seq(0, max(max(SMDs_dt), xmax), 0.1)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = titleSize
      ),
      text = ggplot2::element_text(size = TextSize),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank()
    )
  
  if (length(SMD_names) == 4) {
    love.plot <- love.plot +
      ggplot2::geom_point(x = SMDs_dt[, SMD_names[3]],
                          colour = plotColors[1],
                          size = 4) +
      ggplot2::geom_point(x = SMDs_dt[, SMD_names[4]],
                          colour = plotColors[2],
                          size = 4)
  }
  return(love.plot)
}

# create Kaplan-Meier plot with effect measures
create_KM_plot <- function(data,
                           elig_cohort = NULL,
                           horizon,
                           unit = "years",
                           model_PS,
                           model_S = NULL, 
                           w_meth,
                           weights_meth = NULL,
                           catvar = NULL,
                           contvar = NULL,
                           event_var,
                           time2event_var,
                           trt_var,
                           competing_event_var,
                           n_bootstraps = 100,
                           bootstrap_seed = 1,
                           plotTitle = "",
                           plotColors,
                           TextSize = 28,
                           MetricsText = 4,
                           annotate_figure = TRUE) {
  # unadjusted analysis
  model_formula <- as.formula(paste0(
    "survival::Surv(",
    time2event_var,
    ", ",
    event_var,
    ") ~ ",
    trt_var
  ))
  fit_args_unadjusted <- list(formula = model_formula,
                              data = data,
                              robust = TRUE)

  # Create a list of arguments for the survfit call
  if (w_meth != "") {
    keep <- weights_meth > 0
    fit_args <- list(
      formula = model_formula,
      data = data[keep],
      robust = TRUE,
      weights = weights_meth[keep]
    )
  } else{
    fit_args <- fit_args_unadjusted
  }

  # Use do.call to ensure the function call is constructed correctly
  fit_unadjusted <- do.call(survival::survfit, fit_args_unadjusted)
  fit <- do.call(survival::survfit, fit_args)
  
  # Create KM plot
  KM_plot <- survminer::ggsurvplot(
    fit = fit,
    data = data,
    risk.table = FALSE,
    conf.int = TRUE,
    legend.labs = c("Conservative", "Dialysis"),
    legend.title = "",
    xlab = "Time (months)",
    ylab = "Survival probability (%)",
    xscale = 365 / 2,
    break.time.by = 365 / 2,
    palette = c(plotColors[1], plotColors[2]),
    legend = "none"
  )$plot
  
  # Create KM table
  KM_table <- survminer::ggsurvplot(
    fit = fit_unadjusted,
    data = data,
    risk.table = TRUE,
    conf.int = TRUE,
    legend.labs = c("Conservative", "Dialysis"),
    xlab = "",
    xscale = 365 / 2,
    break.time.by = 365 / 2,
    palette = c(plotColors[1], plotColors[2])
  )$table
  
  # set unit for RMST
  if (unit == "years") {
    div.fact <- 365.25
  } else if (unit == "months") {
    div.fact <- 30.5
  } else if (unit == "days") {
    div.fact <- 1
  }
  
  # only needed for bootstraps to compute new model_PS on each bootstrap
  est <- compute_absolute_relative_risks(
    data = data,
    horizon = horizon,
    event_var = event_var,
    competing_event_var = competing_event_var,
    time2event_var = time2event_var,
    trt = "trt",
    w_meth = w_meth,
    weights_meth = weights_meth,
    catvar = catvar,
    contvar = contvar
  ) |>
    dplyr::filter(time == round(horizon))
  
  est_CI <- risks_boots(
    data = data,
    horizon = horizon,
    elig_cohort = elig_cohort,
    model_PS = model_PS,
    model_S = model_S, 
    event_var = event_var,
    competing_event_var = competing_event_var,
    time2event_var = time2event_var,
    trt_var = trt_var,
    w_meth = w_meth,
    weights_meth = weights_meth,
    catvar = catvar,
    contvar = contvar,
    n_bootstraps = n_bootstraps,
    bootstrap_seed = bootstrap_seed
  ) |>
    dplyr::filter(time == round(horizon))
  
  low_CI <- est_CI |>
    dplyr::filter(name == "conf.low")
  high_CI <- est_CI |>
    dplyr::filter(name == "conf.high")
  
  R0 <- as.numeric(est$R0)
  R0_lower <- as.numeric(low_CI$R0)
  R0_upper <- as.numeric(high_CI$R0)
  
  R1 <- as.numeric(est$R1)
  R1_lower <- as.numeric(low_CI$R1)
  R1_upper <- as.numeric(high_CI$R1)
  
  RD <- as.numeric(est$RD)
  RD_lower <- as.numeric(low_CI$RD)
  RD_upper <- as.numeric(high_CI$RD)
  
  RR <- as.numeric(est$RR)
  RR_lower <- as.numeric(low_CI$RR)
  RR_upper <- as.numeric(high_CI$RR)
  
  RMST0 <- as.numeric(est$RMST0) / div.fact
  RMST0_lower <- as.numeric(low_CI$RMST0) / div.fact
  RMST0_upper <- as.numeric(high_CI$RMST0) / div.fact
  
  RMST1 <- as.numeric(est$RMST1) / div.fact
  RMST1_lower <- as.numeric(low_CI$RMST1) / div.fact
  RMST1_upper <- as.numeric(high_CI$RMST1) / div.fact
  
  dRMST <- as.numeric(est$dRMST) / div.fact
  dRMST_lower <- as.numeric(low_CI$dRMST) / div.fact
  dRMST_upper <- as.numeric(high_CI$dRMST) / div.fact
  
  HR <- as.numeric(est$HR)
  HR_lower <- as.numeric(low_CI$HR)
  HR_upper <- as.numeric(high_CI$HR)
  
  # remove only x and y scales; keep color scale intact
  KM_plot$scales$scales <- Filter(
    f = function(s)
      ! inherits(s, "ScaleContinuousPosition"),
    x = KM_plot$scales$scales
  )
  
  # modify plot
  KM_plot <- KM_plot +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      text = ggplot2::element_text(size = TextSize)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, by = 0.1),
      labels = seq(0, 100, by = 10),
      expand = c(0, 0)     # remove padding below 0
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 365 * horizon, by = 365 / 2),
      labels = function(x)
        round(x / 30)            # convert days → months
    ) +
    ggplot2::ggtitle(plotTitle)
  
  if (annotate_figure) {
    KM_plot <- KM_plot +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0.15,
        hjust = 0,
        size = MetricsText,
        label = paste0(
          "N = ",
          nrow(data),
          "\n",
          "Risk difference, % = ",
          fmt_ci(RD * 100, RD_lower * 100, RD_upper * 100),
          "\n",
          "\u0394RMST, ",
          unit,
          " = ",
          fmt_ci(dRMST, dRMST_lower, dRMST_upper),
          "\n",
          "Hazard ratio = ",
          fmt_ci(HR, HR_lower, HR_upper, 2)
        )
      )
  }
  # store tables
  KM_table <- KM_table +
    ggplot2::theme(
      plot.title   = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks   = ggplot2::element_blank(),
      axis.line    = ggplot2::element_blank()
    )
  
  return(
    list(
      KM_plot = KM_plot,
      KM_table = KM_table,
      R0 = R0,
      R0_lower = R0_lower,
      R0_upper = R0_upper,
      R1 = R1,
      R1_lower = R1_lower,
      R1_upper = R1_upper,
      RR = RR,
      RR_lower = RR_lower,
      RR_upper = RR_upper,
      RD = RD,
      RD_lower = RD_lower,
      RD_upper = RD_upper,
      RMST0 = RMST0,
      RMST0_lower = RMST0_lower,
      RMST0_upper = RMST0_upper,
      RMST1 = RMST1,
      RMST1_lower = RMST1_lower,
      RMST1_upper = RMST1_upper,
      dRMST = dRMST,
      dRMST_lower = dRMST_lower,
      dRMST_upper = dRMST_upper,
      HR = HR,
      HR_lower = HR_lower,
      HR_upper = HR_upper
    )
  )
  
}

# combine histograms on propensity score distributions
combine_PS_plots <- function(plot_list) {
  lapply(seq_along(plot_list), function(i) {
    plot_list[[i]]$hist + plot_list[[i]]$scaled_hist
  })
}

get_metric_positions <- function(n) {
  if (n == 1) {
    return(c(0.50))
  } else if (n == 2) {
    return(c(0.40, 0.60))
  } else if (n == 3) {
    return(c(0.37, 0.49, 0.65))
  } else {
    stop("Only 1, 2, or 3 metrics supported.")
  }
}

create_forest_plot_all_measures <- function(dt,
                                            print_metrics = c("N", "Nonoverlap", "Imbalance")) {
  # Loop over each measure to generate forest plots
  plot_list <- list()
  for (measure in c("RD", "dRMST", "HR")) {
    # create label
    dt[, paste0(measure, "_label") :=
         fmt_ci(get(measure),
                get(paste0(measure, "_lower")),
                get(paste0(measure, "_upper")),
                digits = ifelse(measure == "RD", 1, 2))]
    
    # fix the order
    dt$analysis_name <- factor(dt$analysis_name, levels = rev(dt$analysis_name))
    
    #------------------------------------------------------------
    # Base header
    #------------------------------------------------------------
    header_table <- ggplot2::ggplot(data.frame(y = 0), ggplot2::aes(y = y)) +
      ggplot2::geom_text(
        x = 0,
        label = "Subgroup",
        hjust = 0,
        vjust = 0,
        fontface = "bold"
      ) +
      ggplot2::geom_text(
        x = 1,
        label = ifelse(
          measure == "RD",
          "2-year RD,\n% (95% CI)",
          ifelse(
            measure == "dRMST",
            "2-year \u0394RMST,\nmonths (95% CI)",
            "2-year HR\n(95% CI)"
          )
        ),
        hjust = 1,
        vjust = 0,
        fontface = "bold"
      ) +
      ggplot2::theme_void() +
      ggplot2::scale_y_continuous(limits = c(0, 1))
    
    
    #------------------------------------------------------------
    # Base table
    #------------------------------------------------------------
    table <- ggplot2::ggplot(dt, ggplot2::aes(y = analysis_name)) +
      ggplot2::geom_text(ggplot2::aes(x = 0, label = analysis_name), hjust = 0) +
      ggplot2::geom_text(ggplot2::aes(x = 1, label = .data[[paste0(measure, "_label")]]), hjust = 1) +
      ggplot2::theme_void() +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = 0,
          xend = 1,
          y = as.numeric(analysis_name) + 0.45,
          yend = as.numeric(analysis_name) + 0.45
        ),
        color = "gray90"
      )
    
    
    #------------------------------------------------------------
    # Dynamic metric columns
    #------------------------------------------------------------
    metric_positions <- get_metric_positions(length(print_metrics))
    names(metric_positions) <- print_metrics
    # metric_positions <- c(N = 0.36, Nonoverlap = 0.48, Imbalance = 0.65)
    metric_labels    <- c(N = "Patients,\nNo.",
                          Nonoverlap = "Non overlap,\nNo.",
                          Imbalance = "SMD > 0.1,\nNo.")
    
    for (m in print_metrics) {
      # Header label
      header_table <- header_table +
        ggplot2::geom_text(
          x = metric_positions[m],
          label = metric_labels[m],
          hjust = 0,
          vjust = 0,
          fontface = "bold"
        )
      
      # Table column
      table <- table +
        ggplot2::geom_text(
          data = dt,
          mapping = ggplot2::aes(y = analysis_name, label = .data[[m]]),
          x = metric_positions[m],
          hjust = 0,
          inherit.aes = FALSE
        )
    }
    
    # header plot
    header_forest <- ggplot2::ggplot(data.frame(y = 0), ggplot2::aes(y = y)) +
      ggplot2::geom_text(
        x = 0.4,
        label = ifelse(measure == "dRMST", "Favor conservative", "Favor dialysis"),
        hjust = 1,
        vjust = 0,
        fontface = "bold"
      ) +
      ggplot2::geom_text(
        x = 0.6,
        label = ifelse(measure == "dRMST", "Favor dialysis", "Favor conservative"),
        hjust = 0,
        vjust = 0,
        fontface = "bold"
      ) +
      ggplot2::theme_void() +
      ggplot2::scale_y_continuous(limits = c(0, 1))
    
    # forest plot
    if (measure == "HR") {
      center <- 1
    } else {
      center <- 0  # default fallback
    }
    
    x_dev <- max(abs(dt[[paste0(measure, "_lower")]] - center), 
                 abs(dt[[paste0(measure, "_upper")]] - center), 
                 na.rm = TRUE)
    x_min <- ifelse(measure=="RD", -75, 
                    ifelse(measure=="dRMST", -14, 0))
    x_max <- ifelse(measure=="RD", 75, 
                    ifelse(measure=="dRMST", 14, 2))
    x_break <- ifelse(measure=="RD", 10,
                      ifelse(measure=="dRMST", 2, 0.2))
    forest <- ggplot2::ggplot(dt, ggplot2::aes(y = analysis_name)) +
      ggplot2::geom_point(ggplot2::aes(x = .data[[measure]]),
                          shape = 15,
                          size = 2) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data[[paste0(measure, "_lower")]], 
                                           xmax = .data[[paste0(measure, "_upper")]]),
                              height = 0) +
      ggplot2::geom_vline(
        xintercept = center,
        linetype = "dashed",
        color = "gray50"
      ) +
      # x-axis line in padding area
      ggplot2::geom_hline(
        yintercept = 0.5,
        # below the first row (padding area)
        color = "black"
      ) +
      ggplot2::coord_cartesian(clip = "off") +  # allow drawing in the padding area
      ggplot2::scale_x_continuous(limits = c(x_min, x_max),
                                  breaks = seq(x_min, x_max, x_break),
                                  labels = seq(x_min, x_max, x_break)) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     plot.margin = ggplot2::margin(0, 0, -10, 0))
    
    plot_list[[paste0("header_table_", measure)]] <- header_table
    plot_list[[paste0("table_", measure)]] <- table
    plot_list[[paste0("header_forest_", measure)]] <- header_forest
    plot_list[[paste0("forest_", measure)]] <- forest
  }
  
  # Save plot to results folder
  blank_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  combined_plot <- cowplot::plot_grid(
    plot_list$header_table_RD,
    plot_list$header_forest_RD,
    plot_list$table_RD,
    plot_list$forest_RD,
    plot_list$header_table_dRMST,
    plot_list$header_forest_dRMST,
    plot_list$table_dRMST,
    plot_list$forest_dRMST,
    plot_list$header_table_HR,
    plot_list$header_forest_HR,
    plot_list$table_HR,
    plot_list$forest_HR,
    blank_plot,
    blank_plot,
    ncol = 2,
    nrow = 7,
    rel_heights = c(0.4, 1, 0.4, 1, 0.4, 1, 0.1),
    rel_widths = c(1, 0.5)
  ) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA))
  
  return(list(combined_plot = combined_plot, dt = dt))
}

# create effect plot of benefit versus risk
effect_plot <- function(estimates_df, 
                        y_middle, 
                        measure,
                        y_max_RD,
                        y_min_dRMST,
                        y_max_HR) {
  # add padding on y-axis
  padding_y <- ifelse(measure == "HR", 0.1, ifelse(measure == "RD", 5, 1))
  
  # Define y-axis breaks based on measure
  if (measure == "HR") {
    y_min <- 0 
    y_max <- y_max_HR
    y_breaks <- sort(c(seq(y_min, y_max, by = 0.2), y_middle))
  } else if (measure == "RD") {
    y_min <- -70
    y_max <- y_max_RD
    y_breaks <- sort(c(seq(y_min, y_max, by = 10), y_middle))
  } else {
    y_min <- y_min_dRMST
    y_max <- 10
    y_breaks <- sort(c(seq(y_min, y_max, by = 1), y_middle))
  }
  
  # x-axis minimum
  if (min(estimates_df$effect_modifier_range) < 1) {
    x_min <- 0
    x_min_text <- 0.05
  } else{
    x_min <- 60
    x_min_text <- 61.5
  }
  
  # effect plot
  plot <- ggplot2::ggplot(estimates_df,
                          ggplot2::aes(x = effect_modifier_range, y = get(measure))) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = get(paste0(measure, "_lower")), 
                                      ymax = get(paste0(measure, "_upper"))), 
                         alpha = 0.1) +
    ggplot2::geom_hline(
      yintercept = y_middle,
      linetype = "dashed",
      color = "darkgrey"
    ) +
    ggplot2::scale_y_continuous(limits = c(y_min, y_max),
                                breaks = y_breaks) +
    ggplot2::labs(x = NULL,
                  y = ifelse(
                    measure == "RD",
                    "Risk difference in %",
                    ifelse(measure == "dRMST", "\u0394RMST in months", 
                           "Hazard ratio")
                  )) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      text = ggplot2::element_text(size = 18),
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.5)
    ) +
    ggplot2::annotate(
      "text",
      x = x_min_text,
      y = ifelse(measure=="dRMST", y_middle + padding_y, y_middle - padding_y),
      label = "Favor dialysis",
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6
    ) +
    ggplot2::annotate(
      "segment",
      x = x_min,
      xend = x_min,
      y = ifelse(measure=="dRMST", y_middle + padding_y, y_middle - padding_y),
      yend = ifelse(measure=="dRMST", y_max, y_min),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
      color = "black"
    ) 
  
  # reverse y-axis for dRMST
  if (measure=="dRMST"){
    plot <- plot +
      ggplot2::scale_y_reverse(limits = c(y_max, y_min),
                               breaks = rev(y_breaks))
  }
  
  return(plot)
}

# histogram of effect modifier
create_histogram_stratified <- function(dt, var_name, trt_name, manual_colors) {
  # histogram plot
  hist_stratified <- ggplot2::ggplot(dt, ggplot2::aes(x = get(var_name), fill = as.factor(get(trt_name)))) +
    ggplot2::geom_histogram(
      binwidth = ifelse(var_name == "pred_risk", 0.01, 1),
      color = NA,
      position = "dodge"
    ) +
    ggplot2::labs(
      x = ifelse(
        var_name == "pred_risk",
        "Predicted 2-year mortality risk (%)",
        "Age in years"
      ),
      y = "Count"
    ) +
    ggplot2::scale_fill_manual(
      values = c("0" = manual_colors[1], "1" = manual_colors[2]),
      labels = c("Conservative", "Dialysis")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid = ggplot2::element_blank(),
      text = ggplot2::element_text(size = 18),
      axis.ticks.x = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.direction = "horizontal"
    )
  
  # histogram x-axis limits
  if (var_name == "age") {
    hist_stratified <- hist_stratified +
      ggplot2::scale_x_continuous(breaks = seq(60, 100, 5),
                                  labels = seq(60, 100, 5)) +
      ggplot2::coord_cartesian(xlim = c(60, 100))
  } else if (var_name == "pred_risk") {
    hist_stratified <- hist_stratified +
      ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.1),
                                  labels = scales::percent_format(accuracy = 1)) +
      ggplot2::coord_cartesian(xlim = c(0, 1))
  }
  
  return(hist_stratified)
}

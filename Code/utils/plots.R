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
    ggplot2::scale_fill_manual(
      values = palette
    ) +
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
  scaled_hist <- ggplot2::ggplot(data,
                                 ggplot2::aes(
                                   x = ps,
                                   fill = as.factor(trt),
                                   weight = !!weights,
                                   y = ggplot2::after_stat(scaled)
                                 )) +
    ggplot2::geom_density(alpha = 0.6, bw = 0.05) +
    ggplot2::scale_x_continuous(x_axis_text, limits = c(0, 1)) +
    ggplot2::scale_fill_manual(
      values = palette
    ) +
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
  
  if (!xlab){
    hist <- hist +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
    scaled_hist <- scaled_hist +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
  }
  
  if (!ylab){
    hist <- hist +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())
    scaled_hist <- scaled_hist +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())
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
    breaks = h$breaks[-length(h$breaks)],            # left edges of bins
    counts = h$counts
  )
  
  # --- Create ggplot ---
  h_plot <- ggplot2::ggplot(h_df, ggplot2::aes(x = breaks, y = counts)) +
    ggplot2::geom_bar(stat = "identity",
                      fill = "skyblue",
                      color = "black") +
    ggplot2::labs(title = title, x = x_label, y = y_label) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, max(h_df$breaks), by = 12 * 30.44),  # tick every 6 months
      labels = function(x)
        round(x / 365.25)                                  # label in months
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(size = TextSize),
                   axis.text = ggplot2::element_text(size = TextSize))
  return(h_plot)
}

# create love plot to compare SMD before and after weighting
love_plot <- function(SMDs_dt,
                      SMD_names,
                      plot_title,
                      xlab_title = "Standardized mean difference",
                      xmax,
                      titleSize = 22,
                      TextSize = 26) {
  # make one love plot
  love.plot <- ggplot2::ggplot(data = SMDs_dt, ggplot2::aes(y = factor(
    rownames(SMDs_dt), levels = rev(rownames(SMDs_dt))
  ))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid") +
    ggplot2::geom_vline(xintercept = 0.1, linetype = "dashed") +
    ggplot2::geom_point(x = SMDs_dt[, SMD_names[1]],
                        colour = "#EC7F12",
                        size = 4) +
    ggplot2::geom_point(x = SMDs_dt[, SMD_names[2]],
                        colour = "#9161BD",
                        size = 4) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::xlab(xlab_title) +
    ggplot2::ylab("") +
    ggplot2::scale_x_continuous(limits = c(0, max(max(SMDs_dt), xmax)), breaks =
                                  seq(0, max(max(SMDs_dt), xmax), 0.1)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, 
                                         face = "bold",
                                         size = titleSize),
      text = ggplot2::element_text(size = TextSize),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank()
    )
  
  if (length(SMD_names)==4){
    love.plot <- love.plot + 
      ggplot2::geom_point(x = SMDs_dt[, SMD_names[3]],
                          colour = "#2CA02C",
                          size = 4) +
      ggplot2::geom_point(x = SMDs_dt[, SMD_names[4]],
                          colour = "#77B1D4",
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
                           TextSize = 28) {
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
  if (w_meth!="") {
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
    ylab= "Survival probability (%)",
    xscale = 365/2,
    break.time.by = 365/2,
    palette = c("#EC7F12", "#9161BD"),
    legend = c(0.15, 0.4)
  )$plot
  
  # Create KM table
  KM_table <- survminer::ggsurvplot(
    fit = fit_unadjusted,
    data = data,
    risk.table = TRUE,
    conf.int = TRUE,
    legend.labs = c("Conservative", "Dialysis"),
    xlab = "",
    xscale = 365/2,
    break.time.by = 365/2,
    palette = c("#EC7F12", "#9161BD")
  )$table
  
  # Calculate relative risk
  Cox_fit <- do.call(survival::coxph, fit_args)
  HR <- exp(coef(Cox_fit))
  HR_CI <- exp(confint(Cox_fit))
  HR_lower <- HR_CI[[1]]
  HR_upper <- HR_CI[[2]]
  
  # set unit for RMST
  if (unit == "years") {
    div.fact <- 365.25
  } else if (unit == "months") {
    div.fact <- 30.5
  } else if (unit == "days") {
    div.fact <- 1
  }
  
  if (w_meth!=""){
    # only needed for bootstraps to compute new model_PS on each bootstrap
    KM_est <- compute_absolute_risks(
      data = data,
      event = event_var,
      competing_event = competing_event_var,
      time2event = time2event_var,
      trt = "trt",
      w_meth = w_meth,
      weights_meth = weights_meth,
      catvar = catvar,
      contvar = contvar
    ) |>
      dplyr::filter(time == round(horizon))
    
    KM_CI <- risks_boots(
      data = data,
      elig_cohort = elig_cohort,
      model_PS = model_PS,
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
    
    low_CI <- KM_CI |>
      dplyr::filter(name == "conf.low")
    high_CI <- KM_CI |>
      dplyr::filter(name == "conf.high")
    
    R0 <- as.numeric(KM_est$R0)
    R0_lower <- as.numeric(low_CI$R0)
    R0_upper <- as.numeric(high_CI$R0)
    
    R1 <- as.numeric(KM_est$R1)
    R1_lower <- as.numeric(low_CI$R1)
    R1_upper <- as.numeric(high_CI$R1)
    
    RD <- as.numeric(KM_est$rd)
    RD_lower <- as.numeric(low_CI$rd)
    RD_upper <- as.numeric(high_CI$rd)
    
    RR <- as.numeric(KM_est$rr)
    RR_lower <- as.numeric(low_CI$rr)
    RR_upper <- as.numeric(high_CI$rr)
    
    RMST0 <- as.numeric(KM_est$rmst0) / div.fact
    RMST0_lower <- as.numeric(low_CI$rmst0) / div.fact
    RMST0_upper <- as.numeric(high_CI$rmst0) / div.fact
    
    RMST1 <- as.numeric(KM_est$rmst1) / div.fact
    RMST1_lower <- as.numeric(low_CI$rmst1) / div.fact
    RMST1_upper <- as.numeric(high_CI$rmst1) / div.fact
    
    dRMST <- as.numeric(KM_est$rmst_diff) / div.fact
    dRMST_lower <- as.numeric(low_CI$rmst_diff) / div.fact
    dRMST_upper <- as.numeric(high_CI$rmst_diff) / div.fact
  } else{
    # calculate absolute risk
    KM_summary <- summary(fit, times = horizon, extend = TRUE)
    
    # calculate RMST
    RMST_summary <- summary(fit, rmean = horizon, extend = TRUE)$table
    
    # Delta-method
    p0 <- KM_summary$surv[1]
    R0 <- 1 - p0
    R0_se <- KM_summary$std.err[1]
    R0_lower <- R0 - 1.96 * sqrt(sum(R0_se^2))
    R0_upper <- R0 + 1.96 * sqrt(sum(R0_se^2))
    
    p1 <- KM_summary$surv[2]
    R1 <- 1 - p1
    R1_se <- KM_summary$std.err[2]
    R1_lower <- R1 - 1.96 * sqrt(sum(R1_se^2))
    R1_upper <- R1 + 1.96 * sqrt(sum(R1_se^2))
    
    RD <- R1 - R0
    RD_se <- sqrt(sum((KM_summary$std.err)^2))
    RD_lower <- RD - 1.96 * sqrt(sum(RD_se^2))
    RD_upper <- RD + 1.96 * sqrt(sum(RD_se^2))
    
    RR <- R1 / R0
    RR_se_log <- sqrt((R1_se^2)/(R1^2) + (R0_se^2)/(R0^2))
    RR_lower <- exp(log(RR) - 1.96 * RR_se_log)
    RR_upper <- exp(log(RR) + 1.96 * RR_se_log)
    
    RMST0 <- RMST_summary[1, "rmean"] / div.fact
    RMST0_se <- RMST_summary[1, "se(rmean)"] / div.fact
    RMST0_lower <- RMST0 - 1.96 * RMST0_se 
    RMST0_upper <- RMST0 + 1.96 * RMST0_se 
    
    RMST1 <- RMST_summary[2, "rmean"] / div.fact
    RMST1_se <- RMST_summary[2, "se(rmean)"] / div.fact
    RMST1_lower <- RMST1 - 1.96 * RMST1_se 
    RMST1_upper <- RMST1 + 1.96 * RMST1_se 
    
    dRMST <- as.numeric(diff(RMST_summary[, "rmean"])) / div.fact
    dRMST_se <- sqrt(sum((RMST_summary[, "se(rmean)"])^2)) / div.fact
    dRMST_lower <- dRMST - 1.96 * sqrt(sum(dRMST_se^2))
    dRMST_upper <- dRMST + 1.96 * sqrt(sum(dRMST_se^2))
  }
  
  # remove only x and y scales; keep color scale intact
  KM_plot$scales$scales <- Filter(
    f = function(s) !inherits(s, "ScaleContinuousPosition"),
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
      breaks = seq(0, 365*horizon, by = 365/2),   # same as break.time.by
      labels = function(x) round(x / 30)            # convert days → months
    ) +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0.15,
      hjust = 0,
      label = paste0(
        "N = ",
        nrow(data),
        "\n",
        "Risk difference, % = ",
        fmt_ci(RD * 100, RD_lower * 100, RD_upper * 100),
        "\n",
        "\u0394RMST, ", unit, " = ",
        fmt_ci(dRMST, dRMST_lower, dRMST_upper),
        "\n",
        "Hazard ratio = ",
        fmt_ci(HR, HR_lower, HR_upper, 2)
      )
    ) +
    ggplot2::ggtitle(plotTitle)
  
  # store tables
  KM_table <- KM_table +
    ggplot2::theme(plot.title   = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.x  = ggplot2::element_blank(),
                   axis.ticks   = ggplot2::element_blank(),
                   axis.line    = ggplot2::element_blank())
  
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

# set title of plot accordign to weighting method
plot_title <- function(method, weighting = TRUE) {
  # store plots
  if (method == "" | method == "common_range") {
    paste0("Total population")
  } else if (method == "IPTW" | method == "SMR") {
    method
  } else if (method == "ATE" | method == "ATT") {
    paste0("Fine Stratification ", method)
  } else if (method == "transportability") {
    "Transportability weighting"
  } else if (method == "Crump" |
             method == "Stürmer" | method == "Walker") {
    paste0(method, " trimming")
  } else if (method == "overlap") {
    paste0("Overlap weighting")
  } 
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
                                            print_metrics = c("N", "Nonoverlap", "Imbalance")
                                            ){
  # Loop over each measure to generate forest plots
  plot_list <- list()
  for (measure in c("RD", "dRMST", "HR")) {
    # create label
    dt[, paste0(measure, "_label") :=
         fmt_ci(
           get(measure),
           get(paste0(measure, "_lower")),
           get(paste0(measure, "_upper")),
           digits = ifelse(measure=="RD", 1, 2)
         )
    ]
    
    # fix the order
    dt$analysis_name <- factor(dt$analysis_name, levels = rev(dt$analysis_name))
    
    #------------------------------------------------------------
    # Base header
    #------------------------------------------------------------
    header_table <- ggplot2::ggplot(data.frame(y = 0), ggplot2::aes(y = y)) +
      ggplot2::geom_text(x = 0, label = "Subgroup", hjust = 0, vjust = 0, fontface = "bold") +
      ggplot2::geom_text(
        x = 1,
        label = ifelse(
          measure == "RD", "2-year RD,\n% (95% CI)",
          ifelse(measure == "dRMST", "2-year \u0394RMST,\nmonths (95% CI)", "2-year HR\n(95% CI)")
        ),
        hjust = 1, vjust = 0, fontface = "bold"
      ) +
      ggplot2::theme_void() +
      ggplot2::scale_y_continuous(limits = c(0, 1))
    
    
    #------------------------------------------------------------
    # Base table
    #------------------------------------------------------------
    table <- ggplot2::ggplot(dt, ggplot2::aes(y = analysis_name)) +
      ggplot2::geom_text(ggplot2::aes(x = 0, label = analysis_name), hjust = 0) +
      ggplot2::geom_text(ggplot2::aes(x = 1, label = .data[[ paste0(measure, "_label") ]]), hjust = 1) +
      ggplot2::theme_void() +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = 0, xend = 1,
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
          mapping = ggplot2::aes(
            y = analysis_name,
            label = .data[[m]]
          ),
          x = metric_positions[m],
          hjust = 0,
          inherit.aes = FALSE
        )
    }
    
    # header plot
    header_forest <- ggplot2::ggplot(data.frame(y = 0), ggplot2::aes(y = y)) +
      ggplot2::geom_text(
        x = 0.4,
        label = ifelse(measure=="dRMST", "Favor conservative", "Favor dialysis"),
        hjust = 1,
        vjust = 0, 
        fontface = "bold"
      ) +
      ggplot2::geom_text(
        x = 0.55,
        label = ifelse(measure=="dRMST", "Favor dialysis", "Favor conservative"),
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
    x_min <- center - x_dev
    x_max <- center + x_dev
    forest <- ggplot2::ggplot(dt, ggplot2::aes(y = analysis_name)) +
      ggplot2::geom_point(ggplot2::aes(x = .data[[measure]]), 
                          shape = 15, size = 2) +
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
        yintercept = 0.5,        # below the first row (padding area)
        color = "black"
      ) +
      # x-axis tick and label
      ggplot2::annotate(
        "text",
        x = center,
        y = 0.3,                  # slightly below the axis line
        label = as.character(center),
        vjust = 1
      ) +
      ggplot2::theme_void() +
      ggplot2::coord_cartesian(clip = "off") +  # allow drawing in the padding area
      ggplot2::scale_x_continuous(limits = c(x_min, x_max))
    
    plot_list[[paste0("header_table_", measure)]] <- header_table
    plot_list[[paste0("table_", measure)]] <- table
    plot_list[[paste0("header_forest_", measure)]] <- header_forest
    plot_list[[paste0("forest_", measure)]] <- forest
  }
  
  # Save plot to results folder
  blank_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  combined_plot <- cowplot::plot_grid(
    plot_list$header_table_RD, plot_list$header_forest_RD,
    plot_list$table_RD, plot_list$forest_RD,
    plot_list$header_table_dRMST, plot_list$header_forest_dRMST,
    plot_list$table_dRMST, plot_list$forest_dRMST,
    plot_list$header_table_HR, plot_list$header_forest_HR,
    plot_list$table_HR, plot_list$forest_HR,
    blank_plot, blank_plot,
    ncol = 2, nrow = 7, 
    rel_heights = c(0.4, 1,
                    0.4, 1, 
                    0.4, 1, 0.1),
    rel_widths = c(1, 0.5)
  ) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA))
  
  return(list(combined_plot = combined_plot, dt = dt))
}

# create effect plot of benefit versus risk
effect_plot <- function(dt, y_middle, measure) {
  padding_y <- ifelse(measure=="HR", 0.1, ifelse(measure=="RD", 5, 1))
  
  # Define y-axis breaks based on measure
  y_breaks <- if(measure == "HR") {
    seq(0, 2, by = 0.2)        # more intermediate steps for HR
  } else if(measure == "RD") {
    seq(-60, 60, by = 10)      # more intermediate steps for RD
  } else {  # dRMST
    seq(-10, 10, by = 1)       # more intermediate steps for dRMST
  }
  
  plot <- ggplot2::ggplot(dt, ggplot2::aes(x = pred_risk, y = get(measure))) +
    ggplot2::geom_line(color = manual_colors[1]) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = get(paste0(measure, "_lower")),
                                      ymax = get(paste0(measure, "_upper"))),
                         alpha = 0.1) +
    ggplot2::geom_hline(
      yintercept = y_middle,
      linetype = "dashed",
      color = "darkgrey"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks
    ) +
    ggplot2::labs(x = "Predicted Risk",
                  y = ifelse(measure=="RD", "Risk difference in %", 
                             ifelse(measure=="dRMST", "\u0394RMST in months", 
                                    "Hazard ratio"))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      text = ggplot2::element_text(size = 18),
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_line(color = "black", size = 0.5),
      axis.line.y = ggplot2::element_line(color = "black", size = 0.5)
    ) +
    ggplot2::annotate(
      "text",
      x = 0,
      y = y_middle - padding_y,
      label = paste("Favor", ifelse(measure=="dRMST", "conservative", "dialysis")),
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6
    ) +
    ggplot2::annotate(
      "text",
      x = 0,
      y = y_middle + padding_y,
      label = paste("Favor", ifelse(measure=="dRMST", "dialysis", "conservative")),
      angle = 90,
      hjust = 0,
      vjust = 0.5,
      size = 6
    ) +
    ggplot2::annotate(
      "segment",
      x = 0.05,
      xend = 0.05,
      y = y_middle + padding_y,
      yend = ifelse(measure=="HR", 2, ifelse(measure=="RD", 60, 10)),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
      color = "black"
    ) +
    ggplot2::annotate(
      "segment",
      x = 0.05,
      xend = 0.05,
      y = y_middle - padding_y,
      yend = ifelse(measure=="HR", 0, ifelse(measure=="RD", -60, -10)),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
      color = "black"
    )
  
  return(plot)
}

compute_HTE <- function(data, horizon, nr_pred = 100, print_test_HTE = FALSE) {
  
  # 1. Fit model using rms (ensures rcs is handled correctly)
  # Using 'x=TRUE, y=TRUE' is good practice for rms objects
  fit <- rms::cph(
    survival::Surv(time2event_death_2y, event_death_2y) ~ 
      trt * pred_risk,
    data = data,
    weights = data$sw_IPTW,
    x = TRUE, 
    y = TRUE
  )
  if (print_test_HTE){
    cat(" BIC:", BIC(fit), "\n",
        "Test for HTE p-value = ", anova(fit)["trt * pred_risk  (Factor+Higher Order Factors)", 3], "\n")
  }
  
  # 2. Define prediction grid
  pred_risk_range <- seq(min(data$pred_risk, na.rm = TRUE)+0.05, 
                         max(data$pred_risk, na.rm = TRUE)-0.05, 
                         length.out = nr_pred)
  
  data_0 <- data.frame(trt = 0, pred_risk = pred_risk_range)
  data_1 <- data.frame(trt = 1, pred_risk = pred_risk_range)
  
  # 3. Generate Survival Objects
  sf_0 <- survival::survfit(fit, newdata = data_0)
  sf_1 <- survival::survfit(fit, newdata = data_1)
  
  # 4. Absolute Risk (at specific horizon)
  # 1 - S(t) = Probability of event occurring
  risk_0 <- as.numeric(1 - summary(sf_0, times = horizon)$surv)
  risk_1 <- as.numeric(1 - summary(sf_1, times = horizon)$surv)
  
  # 5. Restricted Mean Survival Time (RMST)
  # Extracting the restricted mean (area under S(t))
  s0_tab <- summary(sf_0, rmean = horizon)$table
  s1_tab <- summary(sf_1, rmean = horizon)$table
  
  # Handle potential column naming variance (*rmean vs rmean)
  RMST_0 <- as.numeric(s0_tab[, "rmean"]) / 30.5
  RMST_1 <- as.numeric(s1_tab[, "rmean"]) / 30.5
  
  # 6. Hazard Ratio (HR)
  # For interaction models with splines, it is safest to predict the 
  # log(HR) as the difference in Linear Predictors (LP)
  lp_0 <- predict(fit, newdata = data_0, type = "lp")
  lp_1 <- predict(fit, newdata = data_1, type = "lp")
  HR   <- exp(lp_1 - lp_0) 
  
  # 7. Consolidate results
  return(
    data.frame(
      pred_risk = pred_risk_range,
      risk_0    = risk_0,
      risk_1    = risk_1,
      RD        = (risk_1 - risk_0)*100,  # Risk Difference (Absolute Risk Reduction)
      dRMST     = RMST_1 - RMST_0,  # Difference in RMST (Time Gain)
      HR        = HR                # Point-wise Hazard Ratio
    )
  )
}

compute_estimates_with_CI <- function(data,
                                      unit = "months",
                                      horizon,
                                      elig_cohort = NULL,
                                      model_PS,
                                      model_S = NULL,
                                      event_var,
                                      competing_event_var,
                                      time2event_var,
                                      trt_var,
                                      w_meth,
                                      trim_meth = NULL,
                                      catvar,
                                      contvar,
                                      n_bootstraps = 100,
                                      bootstrap_seed = 1) {
  # compute point estimates at multiple times
  est_full <- compute_risk_at_multiple_times(
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
    trim_meth = trim_meth,
    catvar = catvar,
    contvar = contvar
  )
  
  # compute point estimat
  est <- est_full |>
    dplyr::filter(time == round(horizon))
  
  # compute 95% CI at multiple time points
  est_CI_full <- compute_CI_at_multiple_times(
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
    trim_meth = trim_meth,
    catvar = catvar,
    contvar = contvar,
    n_bootstraps = n_bootstraps,
    bootstrap_seed = bootstrap_seed
  ) 
  
  # compute 95% CI at horizon
  est_CI <- est_CI_full |>
    dplyr::filter(time == round(horizon))
  low_CI <- est_CI |>
    dplyr::filter(name == "conf.low")
  high_CI <- est_CI |>
    dplyr::filter(name == "conf.high")
  
  # set unit for RMST
  if (unit == "years") {
    div.fact <- 365.25
  } else if (unit == "months") {
    div.fact <- 30.5
  } else if (unit == "days") {
    div.fact <- 1
  }
  
  return(
    list(
      R0 = as.numeric(est$R0),
      R0_lower = as.numeric(low_CI$R0),
      R0_upper = as.numeric(high_CI$R0),
      R1 = as.numeric(est$R1),
      R1_lower = as.numeric(low_CI$R1),
      R1_upper = as.numeric(high_CI$R1),
      RD = as.numeric(est$RD),
      RD_lower = as.numeric(low_CI$RD),
      RD_upper = as.numeric(high_CI$RD),
      RR = as.numeric(est$RR),
      RR_lower = as.numeric(low_CI$RR),
      RR_upper = as.numeric(high_CI$RR),
      RMST0 = as.numeric(est$RMST0) / div.fact,
      RMST0_lower = as.numeric(low_CI$RMST0) / div.fact,
      RMST0_upper = as.numeric(high_CI$RMST0) / div.fact,
      RMST1 = as.numeric(est$RMST1) / div.fact,
      RMST1_lower = as.numeric(low_CI$RMST1) / div.fact,
      RMST1_upper = as.numeric(high_CI$RMST1) / div.fact,
      dRMST = as.numeric(est$dRMST) / div.fact,
      dRMST_lower = as.numeric(low_CI$dRMST) / div.fact,
      dRMST_upper = as.numeric(high_CI$dRMST) / div.fact,
      HR = as.numeric(est$HR),
      HR_lower = as.numeric(low_CI$HR),
      HR_upper = as.numeric(high_CI$HR),
      est_full = est_full,
      est_CI_full = est_CI_full
    )
  )
}

# compute absolute and relative risks at multiple time points
compute_risk_at_multiple_times <- function(data,
                                           horizon,
                                           elig_cohort,
                                           model_PS,
                                           model_S,
                                           event_var,
                                           competing_event_var = "none",
                                           time2event_var,
                                           trt_var,
                                           w_meth,
                                           trim_meth = NULL,
                                           catvar,
                                           contvar) {
  ### Trimming #################################################################
  if (trim_meth!="no_trimming") {
    data <- trim_propensity_scores(
      data = data,
      trt_var = trt_var,
      trim_meth = trim_meth,
      w_meth = "IPTW",
      model_PS = model_PS,
      catvar = catvar,
      contvar = contvar
    )$overlap$data
  }
  
  ### Compute absolute risks ###################################################
  dat.table <- compute_absolute_risk_at_multiple_times(
    data = data,
    elig_cohort = elig_cohort,
    event_var = event_var,
    competing_event_var = competing_event_var,
    time2event_var = time2event_var,
    trt_var = trt_var,
    model_PS = model_PS,
    model_S = model_S,
    w_meth = w_meth,
    catvar = catvar,
    contvar = contvar,
    estimate = "absolute_risk"
  ) |>
    # Pivot data to wide format, only keeping strata and estimate
    dplyr::select(strata, time, estimate) |>
    tidyr::pivot_wider(names_from = strata, values_from = estimate)
  
  # set column names
  colnames(dat.table) <- c("time", "R0", "R1")
  
  ### Compute relative risks ###################################################
  # define model variable
  model_formula <- as.formula(paste0(
    "survival::Surv(",
    time2event_var,
    ", ",
    event_var,
    ") ~ ",
    trt_var
  ))
  
  if (w_meth == "unweighted") {
    # unadjusted analysis
    fit_args <- list(formula = model_formula,
                     data = data,
                     robust = TRUE)
  } else{
    # adjusted analysis
    weights_meth <- create_weights(
      data = data,
      elig_cohort = elig_cohort,
      model_PS = model_PS,
      model_S = model_S,
      w_meth = w_meth,
      catvar = catvar,
      contvar = contvar,
      verbose = FALSE
    )$data$w
    
    # omit those with weights equal to zero
    keep <- weights_meth > 0
    
    # fit weighted model
    fit_args <- list(
      formula = model_formula,
      data = data[keep],
      robust = TRUE,
      weights = weights_meth[keep]
    )
  }
  
  # Fit cox model
  Cox_fit <- do.call(survival::coxph, fit_args)
  HR <- exp(coef(Cox_fit))
  
  # Compute risk differences, risk ratios, RMST, and HR
  dat.table <- dat.table |>
    dplyr::mutate(
      R1 = R1,
      R0 = R0,
      RD = R1 - R0,
      RR = R1 / R0,
      RMST0 = cumsum(1 - R0),
      RMST1 = cumsum(1 - R1),
      dRMST = RMST1 - RMST0,
      HR = HR
    ) 
    # dplyr::filter(time %in% horizon)
  
  # If RR is undefined, set to zero
  dat.table[is.na(dat.table$RR), "RR"] <- 0
  
  return(dat.table)
}

# compute absolute stimates at multiple time points
compute_absolute_risk_at_multiple_times <- function(data,
                                                    elig_cohort,
                                                    event_var,
                                                    competing_event_var,
                                                    time2event_var,
                                                    trt_var,
                                                    model_PS,
                                                    model_S,
                                                    w_meth,
                                                    catvar,
                                                    contvar,
                                                    estimate) {
  # Set some variables depending on whether to estimate absolute risk
  # or numbers at risk
  if (estimate == "absolute_risk") {
    state_selection <- "1"
    estimate_colname <- "estimate"
    baseline_estimate <- 0
  } else if (estimate == "numbers_at_risk") {
    state_selection <- "(s0)"
    estimate_colname <- "n.risk"
    baseline_estimate <- NA
  } else {
    print("Type of estimate is not supported!")
  }
  
  # Change column names for competing events (ce)
  dat.ce <- data
  dat.ce$event_main <- data[[event_var]]
  dat.ce$time2event_ce <- data[[time2event_var]]
  dat.ce$trt_ce <- data[[trt_var]]
  
  if (w_meth != "unweighted") {
    out_weights <- create_weights(
      data = data,
      elig_cohort = elig_cohort,
      model_PS = model_PS,
      model_S = model_S,
      w_meth = w_meth,
      catvar = catvar,
      contvar = contvar,
      verbose = FALSE
    )
    dat.ce$weights_ce <- out_weights$data$w
  }
  
  # Change event variable to represent censoring. 0 is censoring (lowest is seen as censoring), 1 the event, and 2 competing event.
  # If competing event is none, we don't need to do this to compute a normal cumulative incidence curve.
  if (competing_event_var != "none") {
    dat.ce$event_competing <- data[[competing_event_var]]
    dat.ce <- dat.ce |>
      dplyr::mutate(event_ce = ifelse(event_main == 1, 1, ifelse(event_competing == 1, 2, 0)))
  } else {
    dat.ce <- dat.ce |> dplyr::mutate(event_ce = event_main)
  }
  
  # Fit the survival model using the mstate package for multi-state survival models.
  # Individuals can transition between different states over time. The type argument in the Surv() function is set
  # to "mstate" to indicate that the analysis should be performed using multi-state survival models.
  if (w_meth == "unweighted") {
    model <- survival::survfit(
      survival::Surv(time2event_ce, as.factor(event_ce), type = "mstate") ~
        as.factor(trt_ce),
      data = dat.ce
    )
  } else {
    model <- survival::survfit(
      survival::Surv(time2event_ce, as.factor(event_ce), type = "mstate") ~
        as.factor(trt_ce),
      weights = dat.ce$weights_ce,
      data = dat.ce
    )
  }
  
  # Create preliminary dataset with estimates
  dat.risk <- model |>
    # Change fit2 list to dataframe
    broom::tidy() |>
    # Keep only events
    dplyr::filter(state == state_selection) |>
    # Rename column of estimate
    dplyr::mutate(estimate = !!dplyr::sym(estimate_colname)) |>
    # Select only relevant variables
    dplyr::select(time, strata, estimate)
  
  # Get time 0 numbers
  time.0 <- data.frame(
    strata = unique(dat.risk$strata),
    time = 0,
    estimate = baseline_estimate
  )
  
  # For loop to get closets number at risk to each timepoint on the x axis
  dat.table <- dat.risk |>
    # Add time 0 to the dataset
    dplyr::bind_rows(time.0) |>
    # Group by strata
    dplyr::group_by(strata) |>
    # Sort on ascending time and keep only latest timepoint for each chunk
    dplyr::arrange(time) |>
    # Fill in missing chunks (sometimes no events are observed in a chunk)
    dplyr::right_join(tidyr::expand(dat.risk, strata, time = seq(0, max(time))),
                      by = c("strata", "time")) |>
    # For added empty intervals, use estimates of previous timepoint
    dplyr::arrange(time) |>
    tidyr::fill(estimate, .direction = "down") |>
    dplyr::ungroup()
  
  return(dat.table)
}

compute_CI_at_multiple_times <- function(data,
                                         horizon,
                                         elig_cohort,
                                         model_PS,
                                         model_S = NULL,
                                         event_var,
                                         competing_event_var,
                                         time2event_var,
                                         trt_var,
                                         w_meth,
                                         weights_meth,
                                         trim_meth = "no_trimming",
                                         catvar,
                                         contvar,
                                         n_bootstraps,
                                         bootstrap_seed = 123) {
  # Define columns to check
  cols_to_check <- c("R0", "R1", "RD", "RR", "RMST0", "RMST1", "dRMST", "HR")
  
  # 1. Register a parallel backend
  num_cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  
  # 2. SET THE SEED (Crucial for doRNG)
  # This ensures every bootstrap sample is reproducible
  registerDoRNG(seed = bootstrap_seed)
  
  # 3. Execute the bootstrap loop
  # Creates bootsamples and runs the model to each sample
  # perform in parallel
  bootsamps <- foreach::foreach(
    boot = 1:n_bootstraps,
    .packages = c("dplyr", "data.table"),
    .export = c(
      "filter_terms_from_formula",
      "create_weights",
      "trim_propensity_scores",
      "compute_risk_at_multiple_times",
      "compute_absolute_risk_at_multiple_times"
    )
  ) %dorng% {
    tryCatch({
      # Sample from dataset
      d <- sample(1:nrow(data),
                  size = nrow(data),
                  replace = T)
      
      # select patients
      ds_b <- setDT(data[d, ])
      
      # Compute absolute risks in weighted subpopulation
      output <- compute_risk_at_multiple_times(
        data = ds_b,
        horizon = horizon,
        elig_cohort = elig_cohort,
        model_PS = model_PS,
        model_S = model_S,
        event_var = event_var,
        competing_event_var = competing_event_var,
        time2event_var = time2event_var,
        trt_var = trt_var,
        w_meth = w_meth,
        trim_meth = trim_meth,
        catvar = catvar,
        contvar = contvar
      )
      
      # Check for non-finite values (NA, NaN, Inf) across all existing key columns
      # We use 'any(!is.finite(.))' inside 'summarise' to get a logical check
      is_invalid <- output |>
        dplyr::summarise(dplyr::across(dplyr::all_of(cols_to_check), ~ any(!is.finite(.x)))) |>
        any()
      if (is_invalid) {
        return(list(error = "Calculation produced non-finite values (NA, NaN, or Inf)."))
      }
      
      # Return the clean dplyr table of estimates at multiple time points
      output
      
    }, error = function(e) {
      # This captures the actual R error message
      return(list(error = paste0("Error in boot ", boot, ": ", e$message)))
    })
  }
  
  # 4. Stop the parallel backend and clean up
  parallel::stopCluster(cl)
  registerDoSEQ() # Returns R to sequential processing mode
  
  ##############################################################################
  ### Error reporting
  ##############################################################################
  # Separate the successes from the errors
  errors <- Filter(function(x)
    "error" %in% names(x), bootsamps)
  successes <- Filter(function(x)
    ! "error" %in% names(x), bootsamps)
  
  # Print Detailed Report
  n_failed <- length(errors)
  if (n_failed > 0) {
    message(sprintf(
      "Warning: %d of %d bootstraps failed.",
      n_failed,
      n_bootstraps
    ))
    
    # Group and count unique error messages
    error_summary <- table(sapply(errors, `[[`, "error"))
    print(error_summary)
  }
  
  # Get dataframe of bootstrap results at multiple time points
  totalboot <- dplyr::bind_rows(bootsamps) |>
    dplyr::arrange(time)
  
  # remove failed bootstraps
  totalboot <- totalboot[!is.na(totalboot$time),]
  
  # Calculate confidence intervals at multiple time points
  totalboot <- totalboot |>
    dplyr::group_by(time) |>
    dplyr::reframe(dplyr::across(.cols = dplyr::all_of(cols_to_check), ~ quantile(
      .x, probs = c(0.025, 0.975), na.rm = TRUE
    )),
    name = c("conf.low", "conf.high")) |>
    dplyr::ungroup()
  
  return(totalboot)
}

compute_measures <- function(risk_model, 
                             data,
                             plot = FALSE) {
  # psuedovalues
  Score <- riskRegression::Score(
    object = list("model" = risk_model),
    formula = Surv(time2event_death_2y, event_death_2y) ~ 1,
    cens.method = "pseudo",
    data = data,
    times = horizon,
    conf.int = TRUE,
    plots = "calibration"
  )
  
  # discriminiation
  AUC <- Score$AUC$score[["AUC"]]
  
  # calibration intercept
  pseudos <- data.frame(Score$Calibration$plotframe)
  pseudos$cll_pred <- log(-log(1 - pseudos$risk))
  fit_cal_int <- geepack::geese(
    pseudovalue ~ offset(cll_pred),
    data = pseudos,
    id = riskRegression_ID,
    scale.fix = TRUE,
    family = gaussian,
    mean.link = "cloglog",
    corstr = "independence",
    jack = TRUE
  )
  Intercept <- summary(fit_cal_int)$mean$estimate
  
  # calibration slope
  fit_cal_slope <- geepack::geese(
    pseudovalue ~ offset(cll_pred) + cll_pred,
    data = pseudos,
    id = riskRegression_ID,
    scale.fix = TRUE,
    family = gaussian,
    mean.link = "cloglog",
    corstr = "independence",
    jack = TRUE
  )
  Slope <- 1 + summary(fit_cal_slope)$mean["cll_pred", ]$estimate
  
  # smooth pseudo values for calibration plot
  if (plot){
    smooth_pseudos <- predict(stats::loess(
      pseudovalue ~ risk,
      data = pseudos,
      degree = 1,
      span = 0.5
    ),
    se = TRUE)
  } else {
    smooth_pseudos <- NA
  }
  
  return(
    list(
      Intercept = Intercept,
      Slope = Slope,
      AUC = AUC,
      smooth_pseudos = smooth_pseudos,
      pseudos = pseudos
    )
  )
}

compute_HTE <- function(data,
                        horizon,
                        effect_modifier = "lp_risk",
                        effect_modifier_range,
                        print_test_HTE = FALSE) {
  # 0. Set formula
  data <- copy(data)
  data[, effect_modifier := get(effect_modifier)]
  data[, trt := as.numeric(trt) - 1]
  
  # 1. Fit model using rms (ensures rcs is handled correctly)
  fit <- survival::coxph(
    survival::Surv(time2event_death_2y, event_death_2y) ~
      trt * effect_modifier,
    data = data,
    weights = data$sw_IPTW,
    robust = TRUE,
    x = TRUE,
    y = TRUE
  )
  
  if (print_test_HTE) {
    cat(" BIC:",
        BIC(fit),
        "\n",
        "Test for HTE p-value = ",
        summary(fit)$coefficients[3, "Pr(>|z|)"],
        "\n")
  }
  
  # 2. Define prediction grid
  pred_grid_0 <- data.frame(
    trt = 0,
    effect_modifier = effect_modifier_range,
    `trt*effect_modifier` = 0
  )
  pred_grid_1 <- data.frame(
    trt = 1,
    effect_modifier = effect_modifier_range,
    `trt*effect_modifier` = effect_modifier_range
  )
  
  # 3. Generate Survival Objects
  sf_0 <- survival::survfit(fit, newdata = pred_grid_0)
  sf_1 <- survival::survfit(fit, newdata = pred_grid_1)
  
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
  lp_0 <- predict(fit, newdata = pred_grid_0, type = "lp")
  lp_1 <- predict(fit, newdata = pred_grid_1, type = "lp")
  HR   <- exp(lp_1 - lp_0)
  
  # 7. Consolidate results
  return(
    data.frame(
      effect_modifier_range = effect_modifier_range,
      risk_0 = risk_0,
      risk_1 = risk_1,
      RD     = (risk_1 - risk_0) * 100,
      dRMST  = RMST_1 - RMST_0,
      HR     = HR
    )
  )
}

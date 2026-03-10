get_absolute_risks <- function(data,
                               elig_cohort,
                               model_PS = NULL,
                               model_S = NULL,
                               event_vars,
                               competing_event_vars = NULL,
                               time2event_vars,
                               trt_var,
                               # extracts weights already obtained for cohort
                               w_meth,
                               weights_meth,
                               catvar,
                               contvar,
                               bootstrap = FALSE,
                               n_bootstraps = 1000,
                               # only needed for bootstraps to compute new model_PS on each bootstrap
                               bootstrap_seed = 123) {
  # Loop over each outcome
  results = list()
  for (i in seq_along(event_vars)) {
    event_var = event_vars[i]
    time2event = time2event_vars[i]
    competing_event_var = ifelse(is.null(competing_event_vars),
                             "none",
                             competing_event_vars[i])
    
    # Calculate absolute risks
    df_risks <- compute_absolute_relative_risks(
      data = data,
      horizon = horizon,
      elig_cohort = elig_cohort,
      event_var = event_var,
      competing_event_var = competing_event_var,
      time2event_var = time2event_var,
      trt_var = trt_var,
      w_meth = w_meth, 
      weights_meth = weights_meth,
      catvar = catvar,
      contvar = contvar
    ) |>
      dplyr::mutate(name = "estimate")
    
    # Compute confidence intervals if bootstrap is requested
    if (bootstrap) {
      df_risks_CI <- risks_boots(
        data = data,
        model_PS = model_PS,
        model_S = model_S, 
        event_var = event_var,
        competing_event_var = competing_event_var,
        time2event_var = time2event_var,
        trt_var = trt_var,
        w_meth = w_meth,
        weights_meth = weights_meth,
        n_bootstraps = n_bootstraps,
        bootstrap_seed = bootstrap_seed
      )
      
      df_risks <- df_risks |>
        dplyr::bind_rows(df_risks_CI)
    }
    
    # Some reformatting
    df_risks <- df_risks |>
      dplyr::rename(A0 = "as.factor(trt_ce)=0") |>
      dplyr::rename(A1 = "as.factor(trt_ce)=1") |>
      tidyr::pivot_wider(
        names_from = name,
        values_from = c(A0, A1, RD, RR, RMST0, RMST1, dRMST)
      )
    
    # Compute numbers at risk
    df_numbers_at_risk <- compute_numbers_at_risk(
      data = data,
      event_var = event_var,
      competing_event_var = competing_event_var,
      time2event_var = time2event_var,
      trt_var = trt_var,
      weights_meth = weights_meth
    ) |>
      # Some reformatting
      dplyr::mutate(strata = ifelse(strata == "as.factor(trt_ce)=0", "N0", strata)) |>
      dplyr::mutate(strata = ifelse(strata == "as.factor(trt_ce)=1", "N1", strata)) |>
      tidyr::pivot_wider(names_from = strata, values_from = estimate) |>
      # Set the initial numbers at risk
      dplyr::mutate(N0 = ifelse(is.na(N0), sum(1 - data[[trt_var]]), N0)) |>
      dplyr::mutate(N1 = ifelse(is.na(N1), sum(data[[trt_var]]), N1))
    
    df <- df_risks |>
      dplyr::left_join(df_numbers_at_risk, by = "time")
    
    results[[i]] = df
  }
  
  return(results)
}

compute_numbers_at_risk <- function(data,
                                    event_var,
                                    competing_event_var = "none",
                                    time2event_var,
                                    trt_var,
                                    w_meth, 
                                    weights_meth) {
  tab <- compute_estimates(
    data = data,
    event_var = event_var,
    competing_event_var = competing_event_var,
    time2event_var = time2event_var,
    trt_var = trt_var,
    w_meth = w_meth, 
    weights_meth = weights_meth,
    estimate = "numbers_at_risk"
  )
  
  return(tab)
}

compute_absolute_relative_risks <- function(data,
                                   horizon, 
                                   elig_cohort,
                                   event_var,
                                   competing_event_var = "none",
                                   time2event_var,
                                   trt_var,
                                   w_meth, 
                                   weights_meth,
                                   catvar,
                                   contvar,
                                   fit_args) {
  ### Compute absolute risks ###################################################
  dat.table <- compute_estimates(
    data = data,
    elig_cohort = elig_cohort,
    event_var = event_var,
    competing_event_var = competing_event_var,
    time2event_var = time2event_var,
    trt_var = trt_var,
    w_meth = w_meth,
    weights_meth = weights_meth,
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
  
  if (w_meth==""){
    # unadjusted analysis
    fit_args <- list(formula = model_formula,
                     data = data,
                     robust = TRUE)
  } else{
    # adjusted analysis
    keep <- weights_meth > 0
    fit_args <- list(
      formula = model_formula,
      data = data[keep],
      robust = TRUE,
      weights = weights_meth[keep]
    )
  }
  
  # fit cox model
  Cox_fit <- do.call(survival::coxph, fit_args)
  HR <- exp(coef(Cox_fit))
  
  # Compute risk differences, risk ratios, RMST, and HR
  dat.table <- dat.table |>
    dplyr::mutate(
      RD = R1 - R0,
      RR = R1 / R0,
      RMST0 = cumsum(1 - R0),
      RMST1 = cumsum(1 - R1),
      dRMST = RMST1 - RMST0,
      HR = HR
    ) |>
    dplyr::filter(time %in% horizon)
  
  # If RR is undefined, set to zero
  dat.table[is.na(dat.table$RR), "RR"] <- 0
  
  return(dat.table)
}

compute_estimates <- function(data,
                              elig_cohort,
                              event_var,
                              competing_event_var,
                              time2event_var,
                              trt_var,
                              w_meth, 
                              weights_meth,
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
  dat.ce$weights_ce <- weights_meth
  
  # Change event variable to represent censoring. 0 is censoring (lowest is seen as censoring), 1 the event, and 2 competing event.
  # If competing event is none, we don't need to do this to compute a normal cumulative incidence curve.
  if (competing_event_var != "none") {
    dat.ce$event_competing <- data[[competing_event_var]]
    dat.ce <- dat.ce |>
      dplyr::mutate(event_ce = ifelse(event_main == 1, 1, 
                                      ifelse(event_competing == 1, 2, 0)))
  } else {
    dat.ce <- dat.ce |> dplyr::mutate(event_ce = event_main)
  }
  
  # Fit the survival model using the mstate package for multi-state survival models.
  # Individuals can transition between different states over time. The type argument in the Surv() function is set
  # to "mstate" to indicate that the analysis should be performed using multi-state survival models.
  if (w_meth=="") {
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

risks_boots <- function(data,
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
                        catvar,
                        contvar,
                        n_bootstraps,
                        bootstrap_seed = 123) {
  
  # Get all categorical variables from propensity mode
  categorical_vars <- extract_categorical_vars_from_formula(model_PS,
                                                            data)
  categorical_vars <- c(names(categorical_vars[categorical_vars==TRUE]))
  
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
      "categorical_vars",
      "extract_categorical_vars_from_formula",
      "filter_terms_from_formula",
      "create_weights",
      "compute_absolute_relative_risks",
      "compute_estimates",
      "encode_factors",
      "safe_sbw"
    )
  ) %dorng% {
    tryCatch({
      # Sample from dataset
      d <- sample(1:nrow(data),
                  size = nrow(data),
                  replace = T)
      
      # select patients
      ds_b <- setDT(data[d, ]) # id := .I]
      
      # check if there are any categorical variables with only one level
      check_single_level <- apply(ds_b[, ..categorical_vars], 2, 
                                  function(col) length(unique(col)) == 1)
      var_single_level <- c(names(check_single_level[check_single_level == TRUE]))
      
      # Remove categorical variables with a single level from the ps formula
      model_PS <- filter_terms_from_formula(model_PS, var_single_level)
      
      # add weights to data
      cw <- create_weights(
        data = ds_b,
        elig_cohort = elig_cohort,
        model_PS = model_PS,
        model_S = model_S,
        w_meth = w_meth,
        catvar = catvar,
        contvar = contvar,
        verbose = FALSE
      )
      
      # if no weights can be obtained, failed bootstrap
      if (is.null(cw)) {
        cat("Failed bootstrap sample \n")
        return(NULL)  # skip this bootstrap draw
      }
      
      ds_b <- cw$data
      
      # Compute absolute risks in weighted subpopulation
      output <- compute_absolute_relative_risks(
        data = ds_b,
        horizon = horizon,
        elig_cohort = elig_cohort,
        event_var = event_var,
        competing_event_var = competing_event_var,
        time2event_var = time2event_var,
        trt_var = trt_var,
        w_meth = w_meth, 
        weights_meth = ds_b$w,
        catvar = catvar,
        contvar = contvar
      )
      
      # return output
      output
    })
  }
  
  # 4. Stop the parallel backend and clean up
  parallel::stopCluster(cl)
  registerDoSEQ() # Returns R to sequential processing mode
  
  # Remove NULL from bootsamps
  print(paste0(
    "Number of failed bootstraps: ",
    sum(sapply(bootsamps, is.null)),
    " out of ",
    n_bootstraps
  ))
  bootsamps <- bootsamps[!sapply(bootsamps, is.null)]
  
  # Get dataframe of bootstrap results
  totalboot <- dplyr::bind_rows(bootsamps) |>
    dplyr::arrange(time)
  
  # Calculate confidence intervals
  totalboot <- totalboot |>
    dplyr::group_by(time) |>
    dplyr::reframe(
      dplyr::across(.cols = everything(),
                    ~ quantile(.x, probs = c(0.025, 0.975))),
      name = c("conf.low", "conf.high")
    ) |>
    dplyr::ungroup()
  
  return(totalboot)
}

compute_HTE <- function(data, 
                        horizon, 
                        effect_modifier = "pred_risk", 
                        effect_modifier_range,
                        print_test_HTE = FALSE) {
  # 0. Set formula
  data <- copy(data)
  data[, effect_modifier := get(effect_modifier)]
  
  # 1. Fit model using rms (ensures rcs is handled correctly)
  # Using 'x=TRUE, y=TRUE' is good practice for rms objects
  fit <- rms::cph(
    survival::Surv(time2event_death_2y, event_death_2y) ~ 
      trt * effect_modifier,
    data = data,
    weights = data$sw_IPTW,
    x = TRUE, 
    y = TRUE
  )
  
  # set datadist
  dd <- rms::datadist(baseline)
  options(datadist = "dd")
  
  if (print_test_HTE){
    print(fit)
    cat(" BIC:", BIC(fit), "\n",
        "Test for HTE p-value = ",
        anova(fit)["trt * effect_modifier  (Factor+Higher Order Factors)", 3], "\n")
  }
  
  # 2. Define prediction grid
  data_0 <- data.frame(trt = 0, effect_modifier = effect_modifier_range)
  data_1 <- data.frame(trt = 1, effect_modifier = effect_modifier_range)
  
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
      effect_modifier_range = effect_modifier_range,
      risk_0 = risk_0,
      risk_1 = risk_1,
      RD     = (risk_1 - risk_0)*100,
      dRMST  = RMST_1 - RMST_0,
      HR     = HR          
    )
  )
}

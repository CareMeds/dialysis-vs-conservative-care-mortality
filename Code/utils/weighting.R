# Extract variables from formula that are categorical
extract_categorical_vars_from_formula <- function(formula, data) {
  # Get all variables in the formula
  all_labels <- attributes(terms(formula))$term.labels
  
  # Create model frame using the data to evaluate variables
  mf <- model.frame(formula, data)
  
  # Check if the variables are factors
  factor_vars <- sapply(all_labels, function(var) {
    # Check if the variable exists in the data and is a factor
    if (var %in% colnames(mf)) {
      return(is.factor(mf[[var]]))
    } else {
      return(FALSE) # Variable not found in data
    }
  })
  
  # Return a named logical vector indicating which variables are factors
  names(factor_vars) <- all_labels
  
  return(factor_vars)
}

# Function to filter terms
filter_terms_from_formula <- function(formula, exclude_vars) {
  # Extract the LHS and RHS of the formula
  lhs <- formula[[2]]
  rhs <- xfun::attr2(terms(formula), "term.labels")
  
  # Filter RHS terms
  filtered_rhs <- rhs[!sapply(rhs, function(term) {
    any(sapply(exclude_vars, function(var)
      grepl(var, term)))
  })]
  
  # Recreate the formula with the original LHS and filtered RHS
  stats::reformulate(filtered_rhs, response = lhs)
}

# whenever no weights can be found, i.e., in bootstrap, continue
safe_sbw <- function(...) {
  tryCatch(
    sbw::sbw(...),
    error = function(e) NULL
  )
}

# obtain weights
create_weights <- function(data,
                           elig_cohort = NULL,
                           model_PS,
                           model_S = NULL,
                           w_meth = "IPTW",
                           catvar = NULL,
                           contvar = NULL,
                           verbose = TRUE) {
  # Remove categorical variables with a single level from the ps formula
  categorical_vars <- extract_categorical_vars_from_formula(model_PS, data)
  
  vars_with_single_level <- colnames(data[, ..categorical_vars][, which(sapply(.SD, function(col)
    length(unique(col)) == 1))])
  
  if (length(vars_with_single_level) != 0) {
    cat("Excluded variables from propensity score model: ",
        vars_with_single_level,
        "\n")
  }
  model_PS <- filter_terms_from_formula(model_PS, vars_with_single_level)
  
  # fit logistic regression for PS
  denom.fit <- suppressWarnings(stats::glm(model_PS, family = stats::binomial(), data = data))
  coef_ps <- coef(denom.fit)          # return coefficients of propensity score model
  ci_ps <- confint.default(denom.fit) # Wald intervals beta_hat +/- 1.96*SE(beta_hat) sufficient for large samples
  
  # variance inflation factor
  # print(car::vif(denom.fit, type="terms"))
  
  # attach ps
  data[, ps := as.numeric(stats::predict(denom.fit, type = "response"))]
  
  ### INVERSE PROPENSITY TREATMENT WEIGHTING
  if (w_meth == "") {
    data$w <- rep(1, nrow(data))
  } else if (w_meth == "IPTW" | w_meth == "IPSW" | w_meth == "IPSW_IPTW") {
    # numerator of weights
    numer.fit <- stats::glm(trt ~ 1, family = binomial(), data = data)
    pred_num <- stats::predict(numer.fit, type = "response")
    
    # calculation of Stablized weights (mean of all weights == 1)
    # weighted pseudopopulation will be as large as original population
    IPTW_weights <- ifelse(data$trt == 0,
                           ((1 - pred_num) / (1 - data$ps)), 
                           (pred_num / data$ps))
  
    # predict probability of being in eligible cohort for each individual in the analysis data
    if (w_meth == "IPSW") {
      # fit selection model (as you already did)
      fit_S <- stats::glm(model_S, family = binomial(), data = elig_cohort)
      
      # get P(S = 1 | X) for everyone
      pred_S <- stats::predict(fit_S, type = "response")
      
      # use weights only to those with the decision (S == 1)
      PS_sel <- pred_S[elig_cohort$S == 1]
      
      # create Dahabreh's stabilized generalizability weights 
      IPSW_weights <- 1 / PS_sel
    }
    
    if (w_meth == "IPTW"){
      data$w <- IPTW_weights
    } else if (w_meth == "IPSW"){
      data$w <- IPSW_weights
    } else if (w_meth == "IPSW_IPTW"){
      data$w <- data$sw_IPTW * data$sw_IPSW
    }
  } else if (w_meth == "SMR_ATT") {
    # STANDARDIZED MORTALITY RATIO FOR ATT
    data$w <- ifelse(data$trt == 0, (data$ps / (1 - data$ps)), 1)
  } else if (w_meth == "SMR_ATU") {
    # STANDARDIZED MORTALITY RATIO FOR ATU
    data$w <- ifelse(data$trt == 0, 1, ((1 - data$ps) / data$ps))
  } else if (w_meth == "Fine_ATE") {
    # FINE STRATIFICATION ATE
    psfinestrat <- WeightIt::weightit(
      stats::as.formula(trt ~ 1),
      data = data,
      method = NULL,
      estimand = "ATE",
      stabilize = TRUE,
      focal = NULL,
      by = NULL,
      s.weights = NULL,
      ps = data$ps,
      moments = NULL,
      int = FALSE,
      subclass = 50,
      missing = NULL,
      verbose = FALSE,
      include.obj = FALSE
    )
    
    data$w <- psfinestrat[["weights"]]
  } else if (w_meth == "Fine_ATT") {
    # FINE STRATIFICATION ATT
    psfinestrat <- WeightIt::weightit(
      stats::as.formula(trt ~ 1),
      data = data,
      method = NULL,
      estimand = "ATT",
      stabilize = TRUE,
      focal = NULL,
      by = NULL,
      s.weights = NULL,
      ps = data$ps,
      moments = NULL,
      int = FALSE,
      subclass = 50,
      missing = NULL,
      verbose = FALSE,
      include.obj = FALSE
    )
    
    data$w <- psfinestrat[["weights"]]
  } else if (w_meth == "overlap") {
    # overlap weights
    data$w <- ifelse(
      data$trt == 0,
      # control
      data$ps,
      # treated
      1 - data$ps 
    )
  }
  
  if (w_meth != "" & verbose) {
    # distribution of weights
    h_overall = create_histogram(data$w, title = paste("Overall distribution of", w_meth))
    
    data0 <- data[trt == 0, "w"]
    h0 = create_histogram(data0$w, title = paste("Distribution of", w_meth, "within A0"))
    
    data1 <- data[trt == 1, "w"]
    h1 = create_histogram(data1$w, title = paste("Distribution of", w_meth, "within A1"))
    
    # Combine plots
    dist <- h_overall | h0 | h1
    
    # create table of summary statistics
    print("Summary statistics of weights")
    
    # Overall summary
    results_overall <- data[, .(
      Min = min(w, na.rm = TRUE),
      Q1 = quantile(w, 0.25, na.rm = TRUE),
      Median = median(w, na.rm = TRUE),
      Mean = mean(w, na.rm = TRUE),
      Q3 = quantile(w, 0.75, na.rm = TRUE),
      Max = max(w, na.rm = TRUE)
    )]
    results_overall[, trt := "Overall"]
    
    # Stratified summary
    results_strat <- data[, .(
      Min = min(w, na.rm = TRUE),
      Q1 = quantile(w, 0.25, na.rm = TRUE),
      Median = median(w, na.rm = TRUE),
      Mean = mean(w, na.rm = TRUE),
      Q3 = quantile(w, 0.75, na.rm = TRUE),
      Max = max(w, na.rm = TRUE)
    ), by = trt]
    
    # Combine and reorder columns so 'trt' is first
    results <- data.table::rbindlist(list(results_overall, results_strat),
                                     use.names = TRUE,
                                     fill = TRUE)
    data.table::setcolorder(results, c("trt", setdiff(names(results), "trt")))
    print(results)
  } else{
    dist <- NA
  }
  
  return(list(
    data = data,
    dist = dist,
    coef_ps = coef_ps,
    lower_ps = ci_ps[, 1],
    upper_ps = ci_ps[, 2]
  ))
}

trim_propensity_scores <- function(data,
                                   PS_varname, 
                                   trt_varname, 
                                   trim_meth = "common_range") {
  # extract data
  data <- copy(data)
  data$ps <- data[, PS_varname, with = FALSE]
  data$trt <- data[, trt_varname, with = FALSE]
  
  if (trim_meth == "common_range") {
    # Get min and max PS for each group
    ps_control <- data[trt == 0]
    # w_control <- eval(parse(text = paste0("ps_control$", PS_varname)))
    control_min <- min(ps_control$ps)
    control_max <- max(ps_control$ps)
    
    ps_treatment <- data[trt == 1]
    # w_treatment <- eval(parse(text = paste0("ps_treatment$", PS_varname)))
    treatment_min <- min(ps_treatment$ps)
    treatment_max <- max(ps_treatment$ps)
    
    # Find non-overlapping region
    low_region_control <- ps_control[ps < max(treatment_min, control_min)]
    low_region_treatment <- ps_treatment[ps < max(treatment_min, control_min)]
    high_region_control <- ps_control[ps > min(treatment_max, control_max)]
    high_region_treatment <- ps_treatment[ps > min(treatment_max, control_max)]
    nonoverlap <- dplyr::bind_rows(
      low_region_control,
      low_region_treatment,
      high_region_control,
      high_region_treatment
    )
    
    # Trimming non-overlapping region
    overlap <- data[(ps >= max(control_min, treatment_min)) &
                      (ps <= min(control_max, treatment_max))]
  } else if (trim_meth == "Crump") {
    # PS between 0.1 and 0.9
    nonoverlap <- data[ps <= 0.1 | ps >= 0.9]
    overlap <- data[ps > 0.1 & ps < 0.9]
  } else if (trim_meth == "Stürmer") {
    # fifth PS percentile in the treated
    PS_5_treated <- quantile(data$ps[data$trt == 1], prob = 0.05)
    # 95th PS percentile in the untreated
    PS_95_untreated <- quantile(data$ps[data$trt == 0], prob = 0.95)
    
    # define nonoverlap and overlap region
    nonoverlap <- data[ps <= PS_5_treated | ps >= PS_95_untreated]
    overlap <- data[ps > PS_5_treated & ps < PS_95_untreated]
  } else if (trim_meth == "Walker") {
    # preference score is defined as the logit of the PS minus the logit of the treatment prevalence
    trt_prev <- sum(data$trt==1) / nrow(data)
    data$pref_score <- plogis(qlogis(data$ps) - qlogis(trt_prev))
    
    # define nonoverlap and overlap region
    nonoverlap <- data[pref_score <= 0.3 | pref_score >= 0.7]
    overlap <- data[pref_score > 0.3 & pref_score < 0.7]
  }
  
  return(list(overlap = overlap, nonoverlap = nonoverlap))
}

# compute statistics on weights
summarize.weights <- function(data) {
  DT <- data.table::as.data.table(data)
  cols <- setdiff(names(DT), "trt")
  
  # helper function to compute summary stats
  summarize_col <- function(x) {
    c(
      Min    = min(x),
      Q1     = quantile(x, 0.25),
      Median = median(x),
      Mean   = mean(x),
      Q3     = quantile(x, 0.75),
      Max    = max(x)
    )
  }
  
  # overall stats: list of numeric vectors
  overall <- DT[, lapply(.SD, summarize_col), .SDcols = cols]
  
  # stratified stats
  strat <- DT[, lapply(.SD, summarize_col), by = trt, .SDcols = cols]
  
  # flatten each column properly
  Overall <- sprintf("%.2f", as.numeric(overall[[1]]))
  Control <- sprintf("%.2f", as.numeric(strat[trt == 0, .SD, .SDcols = cols][[1]]))
  Treated <- sprintf("%.2f", as.numeric(strat[trt == 1, .SD, .SDcols = cols][[1]]))
  
  # combine into a clean data.table
  result <- data.table::data.table(
    Statistic = c("Min", "Q1", "Median", "Mean", "Q3", "Max"),
    Overall = Overall,
    Control = Control,
    Treated = Treated
  )
  
  return(result)
}

# quantify extreme PS
extreme_PS <- function(nonoverlap, ps) {
  ps_df <- data.frame(stats = c(
    nrow(nonoverlap),
    sum(ps < 0.01),
    sum(ps < 0.05),
    sum(ps > 0.95),
    sum(ps > 0.99)
  ))
  rownames(ps_df) <- c("nonoverlap", "# PS < 0.01", "# PS < 0.05", "# PS > 0.95", "# PS > 0.99")
  return(ps_df)
}
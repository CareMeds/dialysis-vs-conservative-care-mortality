get_incidence_rate <- function(data,
                               id_name, 
                               outcome,
                               time2outcome,
                               weights = NA,
                               .decimals = 1,
                               .per_person_years = 100) {
  # Extract IDs
  data <- copy(data)
  setnames(data, id_name, "ID")
  
  # Convert time from days to years
  time2outcome_in_years <- time2outcome / 365.25
  
  # Estimate weighted incidence rate per 100 person-years
  if (any(is.na(weights))) {
    model <- geepack::geeglm(
      outcome ~ offset(log(time2outcome_in_years)),
      family = poisson,
      data = data,
      id = ID
    )
  } else {
    model <- geepack::geeglm(
      outcome ~ offset(log(time2outcome_in_years)),
      family = poisson,
      data = data,
      id = ID,
      weights = weights
    )
  }
  
  # Extract model statistics
  stats <- generics::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  incid <- stats$estimate[1] * .per_person_years
  lower_95_ci <- stats$conf.low[1] * .per_person_years
  upper_95_ci <- stats$conf.high[1] * .per_person_years
  
  # Create output
  if (any(is.na(weights))) {
    weights <- 1
  }
  incid_ovr <- data.frame(incid) |>
    dplyr::mutate(
      event = sum(outcome * weights),
      tstop = sum(time2outcome_in_years * weights),
      rate = paste0(
        sprintf(paste0("%.", .decimals, "f"), incid),
        " (",
        sprintf(paste0("%.", .decimals, "f"), lower_95_ci),
        ", ",
        sprintf(paste0("%.", .decimals, "f"), upper_95_ci),
        ")"
      )
    ) |>
    dplyr::select(event, tstop, rate)
  return(incid_ovr)
}

get_incidence_rate_stratified <- function(data,
                                          id_name, 
                                          trt,
                                          outcome,
                                          time2outcome,
                                          weights = NA,
                                          .decimals = 1,
                                          .per_person_years = 100) {
  incid_str <- data.frame(
    trt = integer(),
    event = numeric(),
    tstop = numeric(),
    rate = character()
  )
  
  for (trt_val in c(0, 1)) {
    subset_indices <- trt == trt_val
    subset_data <- data[subset_indices, ]
    subset_outcome <- outcome[subset_indices]
    subset_time2outcome <- time2outcome[subset_indices]
    if (length(weights) > 1) {
      subset_weight <- weights[subset_indices]
    } else {
      subset_weight <- NA
    }
    
    result <- get_incidence_rate(
      data = subset_data,
      id_name = id_name,
      outcome = subset_outcome,
      time2outcome = subset_time2outcome,
      weights = subset_weight,
      .decimals = .decimals,
      .per_person_years = .per_person_years
    )
    incid_str <- rbind(incid_str, data.frame(trt = trt_val, result))
  }
  
  rownames(incid_str) <- paste("trt", incid_str$trt, sep = "=")
  incid_str <- incid_str |> dplyr::select(-trt)
  return(incid_str)
}

create_incidence_table <- function(data,
                                   id_name, 
                                   outcome_vars,
                                   time2outcome_vars,
                                   outcome_labels,
                                   treatment_label,
                                   control_label,
                                   weights = NA,
                                   table_caption = paste0("Incidence rates per ", .per_person_years, " person-years"),
                                   .decimals = 1,
                                   .include_overall = FALSE,
                                   .per_person_years = 100) {
  # Create empty table
  incid_outcomes <- matrix("", nrow = length(outcome_vars) * 4, ncol = 3)
  colnames(incid_outcomes) <- c("No. of events",
                                "Person Years",
                                paste0("IR per ", .per_person_years, "PY (95% CI)"))
  row_labels <- c(rbind(outcome_labels, matrix(
    c(
      paste0("Overall"),
      paste0(control_label),
      paste0(treatment_label)
    ),
    nrow = 3,
    ncol = length(outcome_labels)
  )))
  row.names(incid_outcomes) <- row_labels
  
  # Fill table with incidence rates, both overall and stratified by treatment
  for (i in 1:length(outcome_vars)) {
    incid_ovr <- get_incidence_rate(
      data = data,
      id_name = id_name, 
      outcome = data[[outcome_vars[[i]]]],
      time2outcome = data[[time2outcome_vars[[i]]]],
      weights = weights,
      .decimals = .decimals,
      .per_person_years = .per_person_years
    )
    incid_str <- get_incidence_rate_stratified(
      data,
      id_name = id_name, 
      data$trt,
      data[[outcome_vars[[i]]]],
      data[[time2outcome_vars[[i]]]],
      weights = weights,
      .decimals = .decimals,
      .per_person_years = .per_person_years
    )
    incid_outcomes[((i - 1) * 4 + 2), ] <- as.matrix(incid_ovr)
    incid_outcomes[((i - 1) * 4 + 3):((i - 1) * 4 + 4), ] <- as.matrix(incid_str)
  }
  
  # Round No. of events
  incid_outcomes[, 1] <- ifelse(incid_outcomes[, 1] == "", "", round(as.numeric(incid_outcomes[, 1])))
  
  # Round Person Years to requested number of .decimals
  incid_outcomes[, 2] <- ifelse(incid_outcomes[, 2] == "", "", round(as.numeric(incid_outcomes[, 2]), .decimals))
  
  # Remove rows with Overall values if overall is not requested
  stride <- 4
  if (!.include_overall) {
    rows_to_remove <- grep("Overall", rownames(incid_outcomes))
    incid_outcomes <- incid_outcomes[-rows_to_remove, ]
    stride <- 3
  }
  
  # Create table
  formatted_table <- knitr::kable(
    incid_outcomes,
    align = "c",
    caption =  as.character(table_caption),
    escape = FALSE
  )
  
  return(list(formatted_table = formatted_table, raw_table = incid_outcomes))
}
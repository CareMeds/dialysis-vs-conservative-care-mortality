get_HR_from_cox <- function(cox_model,
                            robust_se = TRUE,
                            .decimals = 2) {
  coefficients_summary <- summary(cox_model)$coefficients
  
  hr <- round(exp(coefficients_summary[, "coef"]), .decimals)
  if (robust_se) {
    lower_95_ci <- round(exp(coefficients_summary[, "coef"] -
                               1.96 * coefficients_summary[, "robust se"]),
                         .decimals)
    upper_95_ci <- round(exp(coefficients_summary[, "coef"] +
                               1.96 * coefficients_summary[, "robust se"]),
                         .decimals)
  } else {
    lower_95_ci <- round(exp(coefficients_summary[, "coef"] -
                               1.96 * coefficients_summary[, "se(coef)"]),
                         .decimals)
    upper_95_ci <- round(exp(coefficients_summary[, "coef"] +
                               1.96 * coefficients_summary[, "se(coef)"]),
                         .decimals)
  }
  
  hr_formatted <- paste0(hr, " (", lower_95_ci, ", ", upper_95_ci, ")")
  return(hr_formatted)
}

create_hr_table <- function(data,
                            survival_models,
                            outcome_labels,
                            treatment_label,
                            control_label,
                            table_caption,
                            weights,
                            .decimals = 2) {
  # Add weights to dataframe in temp variable
  # Note: This is necessary because the coxph function does not accept weights as a vector
  if (!any(is.na(weights))) {
    data <- data |> mutate(weights_temp = weights)
  }
  
  # Create a function to fit the Cox model, with or without weights
  with_weights = !any(is.na(weights))
  fit_cox_model <- ifelse(with_weights, function(model) {
    survival::coxph(model,
                    data = data,
                    weights = weights_temp,
                    robust = TRUE)
  }, function(model) {
    survival::coxph(model, data = data)
  })
  
  # Fit the Cox model for each outcome
  cox_models <- lapply(survival_models, fit_cox_model)
  
  # Extract the HRs and 95% CIs
  hr <- lapply(cox_models,
               get_HR_from_cox,
               robust_se = with_weights,
               .decimals = .decimals)
  
  # Format the table by creating row labels and adding slots for reference values
  ref_values <- matrix(c("", "-"),
                       nrow = 2,
                       ncol = length(outcome_labels))
  table <- as.matrix(c(rbind(ref_values, hr)))
  
  # Add row and column labels
  row.names(table) <- c(rbind(outcome_labels, matrix(
    c(paste0(control_label), paste0(treatment_label)),
    nrow = 2,
    ncol = length(outcome_labels)
  )))
  colnames(table) <- ifelse(with_weights, "Weighted HR (95% CI)", "Unweighted HR (95% CI)")
  
  formatted_table <- knitr::kable(
    table,
    align = "c",
    caption =  as.character(table_caption),
    escape = FALSE
  )
  
  return(list(formatted_table = formatted_table, raw_table = table))
}
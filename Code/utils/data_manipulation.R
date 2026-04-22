# make variables that are binary numeric and multi-level factors dummies
encode_factors <- function(dt, catvar, contvar, expand = TRUE){
  if (expand){
    # count number of levels for each variable
    n_levels <- sapply(dt[, ..catvar], nlevels)
    binary_fac <- catvar[n_levels == 2]
    
    # create dummies for variables with more than two levels
    multi_fac  <- catvar[n_levels > 2]
    dummy_list <- lapply(multi_fac, function(v) {
      mm <- model.matrix(~ get(v) - 1, data = dt)
      colnames(mm) <- paste0(v, "_", sub("^get\\(v\\)", "", colnames(mm)))
      as.data.table(mm)
    })
  } else{
    binary_fac <- catvar
    dummy_list <- list()
  }
  
  # make binary variables numeric
  bin_dt <- dt[, lapply(.SD, function(x) as.numeric(x) - 1),
               .SDcols = binary_fac]
  
  # combine continuous, binary, and multiple levels variables
  dt_num <- cbind(
    dt[, ..contvar],
    bin_dt,
    do.call(cbind, dummy_list)
  )
  
  return(dt_num)
}

# Function to calculate age based on birthdate and visit_date
calculate_age <- function(dt, birthdate_col, visitdate_col) {
  dt[, age := as.integer(floor(lubridate::time_length(lubridate::interval(get(birthdate_col), get(visitdate_col)), "years")))]
  return(dt)
}

# remove decimals from counts
fix_counts <- function(x) {
  sapply(x, function(cell) {
    # Fix "370.9 (11.3)" -> "371 (11.3)"
    m <- regmatches(cell, regexpr("^\\d+\\.\\d+(?=\\s*\\()", cell, perl = TRUE))
    if (length(m) == 1) {
      cell <- sub("^\\d+\\.\\d+", round(as.numeric(m)), cell)
    }
    cell
  }, USE.NAMES = FALSE)
}

# For each visit, retrieve information from the past up to 'max_roll_days' days
retrieve_past_info <- function(dt,
                               dictionary,
                               id_name,
                               date_name,
                               var_name,
                               lookback_months = 12,
                               fill_with_zero = FALSE) {
  
  # Subset dictionary to valid entries for the specific variable
  dict_sub <- dictionary[!is.na(get(var_name)), .(
    id = get(id_name),
    dict_date = get(date_name),
    value = get(var_name)
  )]
  
  # Sort dictionary to ensure the most recent record is selected via mult = 'last'
  data.table::setkeyv(dict_sub, c("id", "dict_date"))
  
  # Identify records with missing values
  missing_sub <- dt[is.na(get(var_name)), .(
    id = get(id_name),
    visit_date = get(date_name)
  )]
  
  # Calculate the calendar-aware start date
  if (is.infinite(lookback_months)) {
    missing_sub[, start_date := as.Date("1900-01-01")]
  } else {
    # Using period() avoids the namespace export issue with months()
    # and handles the rollback logic correctly
    missing_sub[, start_date := lubridate::add_with_rollback(
      visit_date,
      lubridate::period(as.integer(-lookback_months), units = "months"),
      roll_to_first = FALSE
    )]
  }
  
  # Non-equi join to find the last known value within the window
  filled <- dict_sub[missing_sub,
                     on = .(id, dict_date >= start_date,
                            dict_date <= visit_date),
                     mult = "last"]
  
  # Update the original data.table by reference
  dt[is.na(get(var_name)), (var_name) := filled$value]
  
  # Set remaining NA values to 0 if requested
  if (fill_with_zero) {
    dt[is.na(get(var_name)), (var_name) := 0]
  }
  
  return(dt)
}

# functions for calculation of eGFR
ckd_epi_2021_cr <- function(creatinine, age, female) {
  k <- ifelse(female == 1, 62, 80)
  alpha <- ifelse(female == 1, -0.241, -0.302)
  return(ifelse(
    female == 1,
    142 *
      (pmin(creatinine / k, 1)^alpha) *
      (pmax(creatinine / k, 1)^(-1.200)) *
      (0.9938^age) *
      1.012,
    142 *
      (pmin(creatinine / k, 1)^alpha) *
      (pmax(creatinine / k, 1)^(-1.200)) *
      (0.9938^age)
  ))
}

# add visit_date at xth birthday date if patient has visits around xth bday
add_bday_rows <- function(dt,
                          bday_year,
                          birthdate_name,
                          id_name,
                          date_name,
                          death_name, 
                          study_end) {
  # ensure unified names
  dt <- copy(dt)
  setnames(dt, id_name, "id")
  setnames(dt, date_name, "visit_date")
  setnames(dt, birthdate_name, "birthdate")
  setnames(dt, death_name, "date_of_death")
  
  # define xth bday
  dt[, bday := data.table::as.IDate(lubridate::add_with_rollback(birthdate,
                                                                 lubridate::years(bday_year)))]
  
  # identify patients with visits before xth birthday
  eligible_bday <- dt[, .(
    bday_val     = bday[1],
    visit_before = any(visit_date < bday[1]),
    visit_after  = any(visit_date > bday[1]),
    died_after   = any(is.na(date_of_death) | date_of_death > bday[1])
  ), by = id][
    bday_val <= study_end & ((visit_before & visit_after) | (visit_before & died_after)), 
    id
  ]
  
  # add xth birthday row
  bday_rows <- dt[id %in% eligible_bday, 
                  if (!any(visit_date == bday[1]))
                    .(
                      visit_date = bday[1],
                      bday = bday[1],
                      birthdate = birthdate[1],
                      date_of_death = date_of_death[1]
                      ), by = id]
  
  # add rows to original data
  dt_added_bday <- rbind(dt, bday_rows, fill = TRUE)
  data.table::setkey(dt_added_bday, id, visit_date)
  
  # update age
  dt_added_bday <- calculate_age(dt = dt_added_bday,
                                 birthdate_col = "birthdate",
                                 visitdate_col = "visit_date")
  
  # revert names
  setnames(dt_added_bday, "id", id_name)
  setnames(dt_added_bday, "bday", paste0("bday", bday_year))
  setnames(dt_added_bday, "visit_date", date_name)
  setnames(dt_added_bday, "birthdate", birthdate_name)
  setnames(dt_added_bday, "date_of_death", death_name)
  
  
  return(dt_added_bday)
}

# determine end of eligibility
eligibility_end <- function(visit_date, next_visit, date_of_death, study_end) {
  
  # Calculate 1 year forward using robust calendar logic
  # Using period() to avoid potential namespace issues with years()
  plus_1_yr <- data.table::as.IDate(lubridate::add_with_rollback(
    visit_date, 
    lubridate::period(1, units = "years"),
    roll_to_first = FALSE
  ))
  
  # Return end of eligibility
  data.table::fifelse(
    # Case 1: Next visit exists and occurs within the 1-year window
    !is.na(next_visit) & next_visit <= plus_1_yr,
    
    # If same day (e.g., multiple records), end date is today; otherwise, subtract 1 day
    data.table::fifelse(visit_date == next_visit, next_visit, next_visit - 1),
    
    # Case 2: No next visit within 1 year. Eligibility ends at the earliest of:
    # - 1 year after current visit
    # - One day before date of death (if applicable)
    # - The study end date
    pmin(
      plus_1_yr,
      date_of_death - 1,
      study_end,
      na.rm = TRUE
    )
  )
}

# enlist ICD codes
list_icds <- function() {
  # list the ICD patterns for comorbidities
  list(
    ihd       = "^I2[0-5]",
    pvd       = "^I7[023]",
    hf        = "^I110|^I130|^I50",
    dm        = "^E1[0-4]",
    scvd      = "^M3[0-2]",
    cancer    = paste0("^C(", paste(
      c(
        sprintf("0%d", 1:9),
        10:26,
        30:34,
        37:41,
        43,
        45:50,
        51:58,
        60:76,
        81:86,
        88,
        90:97
      ),
      collapse = "|"
    ), ")"),
    copd      = "^J44",
    cirr      = "^K74",
    psycho    = paste0("^F(", paste(c(
      sprintf("0%d", 0:9), 10:69, 99
    ), collapse = "|"), ")"),
    hiv       = "^B2[0-4]",
    dementia  = "^F0[0-4]",
    acs       = "^I200|^I2[1-2]",
    hyperten  = "^I1[0-5]",
    vhd       = "^I3[4-7]",
    cevd      = "^I6[5-9]|^G45[0-9]|^G46",
    af        = "^I48",
    arrh      = "^I4[4-7]|^I49",
    lung      = "^I27|^J(4[0-7]|6[0-9]|70|84|85|92|96|982|983)",
    thrombo   = "^I26|^I80[1-9]|^I81|^I820|^I82[2-9]",
    liver     = "^B18|^I85[09]|^I982|^K7[0-7]",
    fracture  = "^S0(2[0-9])|^S[1-9]2|^T0[2,8]|^T1[0,2]|^M48[45]|^M843",
    aki       = "^N17"
  )
}

# create dummy variables
create_dummies <- function(dt, var_names, patterns, id_name, col_name, date_name){
  # set name to retrieve dummies from
  setnames(dt, col_name, "VARIABLE")
  
  for (var_name in var_names) {
    # Extract pattern
    pattern <- patterns[[var_name]]
    
    # Create dummy for comorbidity
    dt[, (var_name) := as.numeric(stringi::stri_detect_regex(VARIABLE, pattern))]
  }

  # Remove unnecessary column
  dt[, VARIABLE := NULL]
  
  # Aggregate for each ID at each date
  dt <- dt[
    , lapply(.SD, max),
    by = c(id_name, date_name),
    .SDcols = var_names
  ]
  
  # Retrieve past info and propagate forwards
  for (var_name in var_names) {
    # retrieve info from past infinitely (except for cancer)
    max_roll_days <- ifelse(var_name == "cancer", 3 * 365.25, Inf)
    
    # Carry forward 1 within each patient, respecting max_roll_days
    dt[, (var_name) := {
      variable <- get(var_name)
      date <- get(date_name)
      
      # define the first one
      first_one <- match(1, variable)
      
      if (!is.na(first_one)) {
        # first date at which first one occurred
        first_one_date <- date[first_one]
        
        # define window
        within_window <- (date >= first_one_date) &
          ((date - first_one_date) <= max_roll_days)
        
        # set to 1 after first one
        variable <- pmax(variable, within_window * 1L, na.rm = TRUE)
      }
      
      # return variable
      variable
    }, by = id_name]
  }
  
  # rename date
  setnames(dt, date_name, "visit_date")
  
  return(dt)
}

# create dictionary of diagnoses
diagnoses.dictionary <- function(inpatient_dt,
                                 outpatient_dt,
                                 lopnr_obtain_diag,
                                 comorbidities,
                                 max_date_dict = NULL){
  ################################################################################
  ### Prepare inpatient and outpatient information
  ################################################################################
  for (dataset in c("inpatient", "outpatient")){
    # Load data
    if (dataset=="inpatient"){
      diagnoses_dt <- inpatient_dt
    } else{
      diagnoses_dt <- outpatient_dt
    }
    setDT(diagnoses_dt)
    
    # Select relevant outpatient columns
    nr_col <- ifelse(dataset=="inpatient", 19, 9)
    diagnosis_columns <- c("HDIA", paste0("DIA", 1:nr_col))
    diagnoses_dt <- diagnoses_dt[!is.na(INDATUMA) &
                                   LOPNR %in% lopnr_obtain_diag, 
                                 .SD, .SDcols = c("LOPNR", "INDATUMA", 
                                                  diagnosis_columns)]
    setkey(diagnoses_dt, LOPNR, INDATUMA)
    
    # Format date
    diagnoses_dt[, `:=`(INDATUMA = as.IDate(as.character(INDATUMA), format("%Y%m%d")))]
    
    # Remove NA's in INDATUMA because sometimes only the year is registered 
    diagnoses_dt <- diagnoses_dt[!is.na(INDATUMA)]
    
    # only look up diagnoses until max_date
    if (!is.null(max_date_dict)) {
      # add max_date information
      diagnoses_dt <- merge(diagnoses_dt,
                            max_date_dict,
                            by = "LOPNR",
                            all.x = TRUE)
      
      # only include those rows that do not have max_date or that are below max_date
      diagnoses_dt <- diagnoses_dt[is.na(max_date) | INDATUMA < max_date]
    }
    
    # Replace "[BLANKAD]" with NA_character_ 
    diagnoses_dt[, (diagnosis_columns) := lapply(.SD, function(x) fifelse(x == "[BLANKAD]", NA_character_, x)), 
                 .SDcols = diagnosis_columns]
    
    # Melt all diagnosis columns into one long column
    diagnoses_dt_long <- data.table::melt(
      diagnoses_dt,
      id.vars = c("LOPNR", "INDATUMA"),
      measure.vars = patterns(paste0("^HDIA$|^DIA[1-", nr_col, "]$")), # match HDIA and DIA1 ...
      value.name = "DIAGNOSIS",
      variable.name = "DIATYPE",
      na.rm = TRUE # drop NA
    )
    
    # Remove duplicate diagnoses codes on same date for ID for efficiency
    diagnoses_dt_long <- unique(diagnoses_dt_long, 
                                by = c("LOPNR", "INDATUMA", "DIAGNOSIS"))
    diagnoses_dt_long[, DIATYPE := NULL]
    setorder(diagnoses_dt_long, LOPNR, INDATUMA, DIAGNOSIS)
    
    # Save dt
    assign(paste0(dataset, "_dt"), diagnoses_dt_long)
  }
  
  # merge inpatient and outpatient
  diagnoses_dt <- rbindlist(list(inpatient_dt, outpatient_dt))
  
  # Keep only unique diagnoses
  diagnoses_dt <- unique(diagnoses_dt, by = c("LOPNR", "INDATUMA", "DIAGNOSIS"))
  setorder(diagnoses_dt, LOPNR, INDATUMA, DIAGNOSIS)
  
  # Create dummies for comorbidities
  diagnoses_dt <- create_dummies(dt = diagnoses_dt,
                                 var_names = comorbidities,
                                 patterns = list_icds(),
                                 id_name = "LOPNR", 
                                 col_name = "DIAGNOSIS",
                                 date_name = "INDATUMA")
  
  return(
    list(
      diagnoses_dt = diagnoses_dt,
      inpatient_dt = inpatient_dt,
      outpatient_dt = outpatient_dt
    )
  )
}

# Helper function to format a number
fmt <- function(x, digits = 1) {
  sprintf(paste0("%.", digits, "f"), round(x, digits))
}

# Helper function to format a value with CI
fmt_ci <- function(estimate, lower, upper, digits = 1) {
  paste0(fmt(estimate, digits),
         " (",
         fmt(lower, digits),
         ", ",
         fmt(upper, digits),
         ")")
}

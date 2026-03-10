################################################################################
### Decision for dialysis versus conservative care
### PART 2 - Coavariate and outcome derivation
################################################################################

# remove history
# rm(list = ls(all.names = TRUE))

# load data
setwd("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/")
load("Data/new_cohort.Rdata")
load("Data/merged_ckd.Rdata")
load("Data/cleaned/snr_inpatient.Rdata")
load("Data/cleaned/snr_outpatient.Rdata")
load("Data/cleaned/snr_lmed.Rdata")
load("Data/cleaned/snr_death.Rdata")

# load functions
pacman::p_load("dplyr", "tidyr", "readr", "lubridate", "stringr")
source("Code/utils/data_manipulation.R")
source("Code/utils/outcome_derivation.R")

# load libraries
library(data.table)

# convert inpatient
inpatient <- UT_R_PAR_SV_123160_2023
setDT(inpatient)
inpatient[, INDATUMA := as.IDate(as.character(INDATUMA), format="%Y%m%d")]

# convert inpatient
outpatient <- UT_R_PAR_OV_123160_2023
setDT(outpatient)

# convert lmed
lmed <- UT_R_LMED_123160_2023
setDT(lmed)
lmed[, EDATUM := as.IDate(EDATUM, format="%Y%m%d")]

# ID name
id_name <- "LOPNR"

################################################################################
### Add covariates to eligible cohort without treatment decision and analysis cohort
################################################################################
for (cohort_name in c("cohort", "elig_cohort")) {
  # evaluate for correct cohort
  cohort <- eval(parse(text = cohort_name))
  
  ################################################################################
  ### Add comorbidities acs, hyperten, vhd, cevd, af, arrh, lung, thrombo, liver, fracture, aki
  ################################################################################
  other_comorb <- c(
    "acs",
    "hyperten",
    "vhd",
    "cevd",
    "af",
    "arrh",
    "lung",
    "thrombo",
    "liver",
    "fracture",
    "aki"
  )
  
  # create dictionary of comorbidities from snr_inpatient and snr_outpatient
  in_out_dict <- diagnoses.dictionary(
    inpatient_dt = UT_R_PAR_SV_123160_2023,
    outpatient_dt = UT_R_PAR_OV_123160_2023,
    lopnr_obtain_diag = unique(cohort$LOPNR),
    comorbidities = other_comorb
  )
  
  # merge with cohort
  # TODO: do I really need merge or can I only use retrieve_past_info
  cohort_other_comorb <- merge(cohort,
                               in_out_dict$diagnoses_dt,
                               by = c(id_name, "visit_date"),
                               all.x = TRUE)
  
  # retrieve past info infinitely after merging with cohort
  for (comorbidity in other_comorb) {
    cohort_other_comorb <- retrieve_past_info(
      dt = cohort_other_comorb,
      dictionary = in_out_dict$diagnoses_dt,
      id_name = id_name,
      date_name = "visit_date",
      var_name = comorbidity,
      max_roll_days = Inf,
      setNA = TRUE # set remaining NA to zero
    )
  }
  
  ################################################################################
  ### Hospitalizations in past year
  ################################################################################
  # obtain inpatient info, only relevant patients
  dia_dt <- inpatient[LOPNR %in% cohort$LOPNR,
                      .SD[1], 
                      by = c(id_name,
                             "INDATUMA")][,
                                          c(id_name, "INDATUMA", "HDIA"),
                                          with = FALSE]
  
  # calculate the number of hospitalizations in past year (any + cardiovascular)
  hospital <- dia_dt[
    cohort[, c(id_name, "visit_date"), with = FALSE], 
    on = id_name, 
    allow.cartesian = TRUE
  ][
    !is.na(INDATUMA) &
      INDATUMA >= (visit_date - 365.25) & 
      INDATUMA <= visit_date,
    .(
      n_hospital = .N,
      n_cvd_hospital = sum(startsWith(HDIA, "I"), na.rm = TRUE)
    ),
    by = id_name
  ]
  
  # merge with cohort
  cohort_hosp <- merge(cohort_other_comorb,
                       hospital,
                       by = id_name,
                       all.x = TRUE)
  
  # replace NA by zero, i.e., if there was no hospitalization in the past year
  cohort_hosp[, `:=` (
    n_hospital = fifelse(is.na(n_hospital), 0L, n_hospital),
    n_cvd_hospital = fifelse(is.na(n_cvd_hospital), 0L, n_cvd_hospital)
  )]
  
  ################################################################################
  ### Medications and iron based on ATC code
  ################################################################################
  lmed_dict <- lmed[LOPNR %in% cohort$LOPNR,
                    c(id_name, "EDATUM", "ATC"),
                    with = FALSE][, unique(.SD)]
  
  # create medications data frames
  med_patterns <- list(
    bblock       = "^C07",
    calblock     = "^C08",
    diuretic     = "^C03",
    rasi         = "^C09A|^C09B|^C09C|^C09D",
    lipid        = "^C10",
    esa          = "^B03XA",
    phosbinder   = "^V03AE02|^V03AE03|^V03AE04|^V03AE05|^V03AE06|^V03AE07|^V03AE08",
    vitamind     = "^A11CC",
    digoxin      = "^C01AA05",
    vasodilator  = "^C01D",
    antiplatelet = "^B01AC",
    anticoag     = "^B01AA|^B01AE07|^B01AF|^B01AX05",
    # potasbinder  = "^V01AE01|^V03AE09|^V03AE10",
    iron_po      = "^B03AA|^B03AB|^B03AD|^B03AE",
    iron_iv      = "^B03AC"
  )

  # extract medications
  lmed_dt <- create_dummies(dt = lmed_dict,
                            var_names = names(med_patterns),
                            patterns = med_patterns,
                            id_name = id_name, 
                            col_name = "ATC",
                            date_name = "EDATUM")
  
  # append to cohort
  cohort_med <- merge(cohort_hosp,
                      lmed_dt,
                      by = c(id_name, "visit_date"),
                      all.x = TRUE)
  
  # retrieve past info infinitely after merging with cohort
  for (medication in names(med_patterns)) {
    cohort_med <- retrieve_past_info(
      dt = cohort_med,
      dictionary = lmed_dt,
      id_name = id_name,
      date_name = "visit_date",
      var_name = medication,
      max_roll_days = 365.25 / 2,
      setNA = TRUE # set remaining NA to zero
    )
  }
  
  ################################################################################
  ### Combine esa and iron from CKD and medications data
  ################################################################################
  # extract iron
  esa_iron_ckd_dt <- merged_ckd[LOPNR %in% cohort$LOPNR &
                                  (!is.na(iron_med) |
                                     !is.na(esa)),
                                c(id_name,
                                  "visit_date",
                                  "esa",
                                  "iron_med",
                                  "iron_type"),
                                with = FALSE][, unique(.SD)]
  
  # iron from medications dt
  esa_iron_med_dt <- cohort_med[!is.na(esa) | 
                                  !is.na(iron_iv) | 
                                  !is.na(iron_po), 
                                c(id_name, 
                                  "visit_date",
                                  "esa",
                                  "iron_iv",
                                  "iron_po"),
                                with = FALSE][, unique(.SD)]
  
  # combine iron_dt from CKD and medications data
  esa_iron_dt <- merge(esa_iron_ckd_dt,
                       esa_iron_med_dt,
                       by = c(id_name,
                              "visit_date"),
                       all.x = TRUE)
  esa_iron_dt[, `:=`
              (esa = fifelse((!is.na(esa.x) & esa.x == 1) | 
                               (!is.na(esa.y) & esa.y == 1), 1, 0),
                iron_iv = fifelse(
                  (!is.na(iron_type) & iron_type == "i.v.") |
                    (!is.na(iron_iv) & iron_iv == 1), 1, 0),
                iron_po = fifelse(
                  (!is.na(iron_type) & iron_type == "p.o.") |
                    (!is.na(iron_po) & iron_po == 1), 1, 0)
                )]
  esa_iron_dt <- esa_iron_dt[, c(id_name,
                                 "visit_date",
                                 "esa",
                                 "iron_iv",
                                 "iron_po"),
                             with = FALSE]
  
  # add esa and iron to cohort by one year look back for iron
  cohort_esa_iron <- copy(cohort_med)
  for (var_name in c("esa", "iron_iv", "iron_po")){
    cohort_esa_iron[, (var_name) := NA_real_]
    cohort_esa_iron <- retrieve_past_info(
      dt = cohort_esa_iron,
      dictionary = esa_iron_dt,
      id_name = id_name,
      date_name = "visit_date",
      var_name = var_name,
      max_roll_days = 365.25,
      setNA = TRUE # set remaining NA to zero
    )
  }
  
  ################################################################################
  ### Primary kidney disease from snr_ckd data
  ################################################################################
  # obtain information from snr_ckd
  # 0 = no primary kidney disease
  # 1 = Diabetesnefropati
  # 2 = Hyperoni
  # 3 = Other, i.e., Adult polycystisk njursjukdom, Glomerulonefrit, Pyelonefrit, Renovaskular, Uremi UNS
  prd_dt <- merged_ckd[LOPNR %in% cohort$LOPNR & !is.na(prd_cat),
                       c(id_name, "visit_date", "prd_cat"), 
                       with=FALSE][, unique(.SD)]
  prd_dt[, prd_cat := fifelse(prd_cat == "Diabetesnefropati",
                              1,
                              fifelse(prd_cat == "Hypertoni", 2, 3))]
  
  # add prd to cohort_iron, looking infinitely back
  cohort_prd <- prd_dt[cohort_esa_iron, on = c(id_name, "visit_date"), roll = TRUE]
  
  # define as factor
  cohort_prd[, prd_cat := factor(
    prd_cat,
    levels = c(1L, 2L, 3L),
    labels = c("Diabetesnefropati", "Hypertoni", "Other")
  )]
  
  ################################################################################
  ### Education category
  ################################################################################
  # extract education, only keep earliest education date for each patient
  edu_dt <- merged_ckd[
    LOPNR %in% cohort$LOPNR &
      !is.na(info_date1) &
      !is.na(info_type1), 
    c(id_name, "info_date1", "info_type1"),
    with = FALSE
    ][, `:=`
      (visit_date = as.IDate(info_date1, format = "%m/%d/%Y"),
      edu = 1)
      ][
        order(visit_date),
        .SD[1], 
        by = id_name
      ]
  
  # merge with cohort
  cohort_prd[, edu := NA_real_]
  cohort_edu <- retrieve_past_info(
    dt = cohort_prd,
    dictionary = edu_dt,
    id_name = id_name,
    date_name = "visit_date",
    var_name = "edu",
    max_roll_days = Inf,
    setNA = TRUE # set remaining NA to zero
  )
  
  ################################################################################
  ### Clinic level
  ################################################################################
  geo_dt <- merged_ckd[LOPNR %in% cohort$LOPNR &
                         (!is.na(clinic) |
                            !is.na(county)), # at least one is non-missing
                       c(id_name, 
                         "visit_date",
                         "clinic", 
                         "county"),
                       with = FALSE][, unique(.SD)]
  
  # Define mapping as a named vector
  county_to_region <- c(
    "Stockholm"       = "Stockholm",
    "Skane"           = "Sodra",
    "Blekinge"        = "Sodra",
    "Kronoberg"       = "Sodra",
    "Halland"         = "Sodra",
    "Orebro"          = "Orebro.Uppsala",
    "Uppsala"         = "Orebro.Uppsala",
    "Sodermanland"    = "Orebro.Uppsala",
    "Vastmanland"     = "Orebro.Uppsala",
    "Vastra Gotaland" = "Vastra",
    "Varmland"        = "Vastra",
    "Jonkoping"       = "Vastra",
    "Dalarna"         = "Vastra",
    "Gavleborg"       = "Vastra",
    "Kalmar"          = "Other regions", # "Sydostra",
    "Ostergotland"    = "Other regions", # "Sydostra",
    "Gotland"         = "Other regions", # "Sydostra",
    "Norrbotten"      = "Other regions", # "Norra",
    "Vasterbotten"    = "Other regions", # "Norra",
    "Vasternorrland"  = "Other regions", # "Norra",
    "Jamtland"        = "Other regions", # "Norra",
    "Ok\xe4nd"        = "Other regions", # meaning Unknown, these are all referred to other disciplines, merge with reference
    "Utrikes"         = "Other regions"  # meaning Emigrated, merge with reference
  )
  
  # Add region column based on mapping
  geo_dt[, region := county_to_region[county]]

  # create clinic levels
  # local clinics
  clinic_lev1 <- c(
    "Avesta",
    "Bollnas",
    "Eksjo",
    "Gallivare",
    "Hassleholm",
    "Karlshamn",
    "Karlskoga",
    "Koping",
    "Ljungby",
    # "Lyckesele",
    "Lycksele", # new
    "Mora",
    "Motala",
    "Nykoping",
    "Pitea",
    "Skelleftea",
    "Solleftea",
    "Varberg/Kungsbacka",
    "Varnamo",
    "Vastervik",
    "Ystad",
    "Angelholm",
    "Ornskoldsvik",
    "Visby",
    "Falkoping",
    "Gbg, Lundby",
    "Trelleborg",
    "Utrikes"         # emigration merged with local
  )
  
  # regional clinics
  clinic_lev2 <- c(
    "Boras",
    "Danderyd",
    "Eskilstuna",
    "Falun",
    # "Gävle",
    "Gavle", # new
    "Halmstad",
    "Helsingborg",
    "Jonkoping",
    "Kalmar",
    "Karlskrona",
    "Karlstad",
    "Kristianstad",
    "Norrkoping",
    "Skovde",
    "Sunderby",
    "Sundsvall",
    "Trollhattan, NAL",
    "Vasteras",
    "Vaxjo",
    "Ostersund",
    "Trollhattan",
    "Ej Njurmedicin"  # referred back to other discipline, merged with regional
  )
  
  # academic clinics
  clinic_lev3 <- c(
    # "Gbg SU/Ostra dialysmott",
    "Gbg SU/Ostra",       # new
    "Gbg, SU/Njurmed",
    "Gbg, SU/Trpl",       # new
    "Karolinska Njur med",
    "Linkoping",
    "Lund Njurmed",
    # "Malmo, njurmed",
    "Malmo, Heleneholms", # new
    "Molndal",
    "Uppsala, med",
    "Uppsala, Trpl",
    "Umea",
    "Orebro",
    "Gbg/Ostra",
    "Huddinge-K Njur med",
    # "Huddinge-K Njur med (Gammal)",
    "Huddinge-K, Trpl",   # new
    "Malmo",
    "Solna-K Njur med",
    # "Solna-K Njur med (Gammal)",
    "Solna, diaverum",    # new
    "Nacka",            # dialysis unit, academic 
    "Sodertalje",       # dialysis unit, academic 
    "J\xe4rf\xe4lla"    # dialysis unit, academic 
  )
  
  # merge and assign clinic levels, infinite look-back
  cohort_geo <- geo_dt[
    cohort_edu, on = .(LOPNR, visit_date),
    roll = TRUE                 
  ]
  
  # create factor for clinic level
  cohort_geo[, clinic_level := factor(
    fifelse(
      clinic %in% clinic_lev1, 1,
      fifelse(clinic %in% clinic_lev2, 2,
              fifelse(clinic %in% clinic_lev3, 3, NA))),
    levels = 1:3,
    labels = c("Local", "Regional", "Academic"))
  ]
  
  ################################################################################
  ### Nursing home
  ################################################################################
  # keep only nursing home, i.e., codes starting with 15
  nursing_dt <- outpatient[
    MVO == "020" | MVO == "243" | MVO == "246",
    c(id_name, "INDATUMA"), 
    with = FALSE
    ][, `:=`
      (visit_date = as.IDate(as.character(INDATUMA), format = "%Y%m%d"),
      nursing_home = 1)
      ][
        order(visit_date),
        .SD[1], 
        by = id_name
        ]
  
  # merge with cohort
  cohort_geo[, nursing_home := NA_real_]
  cohort_nursing <- retrieve_past_info(
    dt = cohort_geo,
    dictionary = nursing_dt,
    id_name = id_name,
    date_name = "visit_date",
    var_name = "nursing_home",
    max_roll_days = Inf,
    setNA = TRUE # set remaining NA to zero
  )
  
  ################################################################################
  ### Calendar year
  ################################################################################
  cohort_year <- cohort_nursing[, calendar_year := year(as.IDate(visit_date))]
  
  ################################################################################
  ### Primary outcome all-cause mortality (5-year, 6-month)
  ################################################################################
  # Prepare datasets
  # All-cause death
  death_dt <- cohort_year[, c(id_name, "DODSDAT"),
                          with = FALSE][, unique(.SD)]
  
  # CV death
  cv_death_dt <- merged_ckd[, c(id_name, "DODSDAT", "ULORSAK"),
                            with = FALSE][, unique(.SD)]
  
  # MI and Stroke
  mi_stroke_dt <- merge(inpatient, 
                        death_dt,
                        by = id_name, 
                        all.x = TRUE)
  
  # KRT date
  krt_dt <- merged_ckd[, c(id_name, "krt_startdate"),
                       with = FALSE][, unique(.SD)]
  
  # Define outcomes
  outcomes_list <- list(
    # All-cause death: uses DODSDAT from cohort
    death = list(
      # dataset with all-cause death dates
      dataset = death_dt,
      # column with date of death
      date_col = "DODSDAT",
      # no code column needed
      code_col = NULL,
      # no ICD filtering
      codes = NULL               
    ),
    
    # Cardiovascular death: uses DODSDAT and ULORSAK from merged CKD dataset
    cvdeath = list(
      # dataset with cause-specific death
      dataset = cv_death_dt,
      # date of death
      date_col = "DODSDAT",
      # ICD code column for cause of death
      code_col = "ULORSAK",
      # pattern to select cardiovascular deaths
      codes = "^I"               
    ),
    
    # Myocardial infarction (MI): ICD codes I21, I22, I23 from inpatient dataset
    mi = list(
      # inpatient dataset filtered for cohort
      dataset = mi_stroke_dt,
      # date of diagnoses registration
      date_col = "INDATUMA",
      # diagnosis code column
      code_col = "HDIA",
      # ICD codes for MI
      codes = "^I21|^I22|^I23"   
    ),
    
    # Stroke: ICD codes I60-I64 from inpatient dataset
    stroke = list(
      # inpatient dataset filtered for cohort
      dataset = mi_stroke_dt,
      # date of diagnoses registration
      date_col = "INDATUMA",
      # diagnosis code column
      code_col = "HDIA",
      # ICD codes for stroke
      codes = "^I60|^I61|^I62|^I63|^I64"
    ),
    
    # KRT: uses krt_startdate
    KRT = list(
      # dataset with all-cause death dates
      dataset = krt_dt,
      # column with date of death
      date_col = "krt_startdate",
      # no code column needed
      code_col = NULL,
      # no ICD filtering
      codes = NULL               
    )
  )
  
  # Define windows
  windows <- list(
    "1y" = 1,
    "2y" = 2,
    "inf" = 100
  )
  
  # add outcomes to cohort_dt
  cohort_outcomes <- add_multiple_outcomes(
    cohort_dt = cohort_year,
    windows = windows,
    id_name = id_name,
    end_follow_up = as.IDate("2024-12-31"),
    outcomes_list = outcomes_list
  )
  
  # add non-CV death and MACE outcomes
  for (window in names(windows)) {
    # Non-CV death: death event but not CV death
    cohort_outcomes[, paste0("noncvdeath_", window) :=
                      fifelse(get(paste0("event_death_", window)) == 1 &
                                get(paste0("event_cvdeath_", window)) == 0, 1L, 0L)]
    
    # MACE: composite of CV death, MI, Stroke
    # Event indicator columns
    mace_events <- c(
      paste0("event_cvdeath_", window),
      paste0("event_mi_", window),
      paste0("event_stroke_", window)
    )
    
    # Event date columns
    mace_dates <- c(
      paste0("event_dt_cvdeath_", window),
      paste0("event_dt_mi_", window),
      paste0("event_dt_stroke_", window)
    )
    
    # Compute MACE indicator: 1 if any of the components occurred
    cohort_outcomes[, paste0("event_mace_", window) := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = mace_events]
    
    # Compute MACE date: earliest date among components
    cohort_outcomes[, paste0("event_dt_mace_", window) := do.call(pmin, c(.SD, na.rm = TRUE)), .SDcols = mace_dates]
    
    # Compute time to MACE from visit date
    cohort_outcomes[, paste0("time2event_mace_", window) :=
                      as.numeric(get(paste0("event_dt_mace_", window)) - visit_date)]
  }
  
  ################################################################################
  ### Categorize some covariates
  ################################################################################
  # Categorize and set factors for all variables in one call
  cohort_final <- cohort_outcomes[, `:=`(
    age_cat = cut(
      age,
      breaks = c(0, 70, 75, 80, 10000),
      labels = c("65-69", "70-74", "75-79", ">=80"),
      right = FALSE
    ),
    
    Davies_score_cat = factor(cut(
      Davies_score,
      breaks = c(0, 2, 4, 10),
      labels = c("<2", "2-4", ">4"),
      right = FALSE
    )),
    
    egfr_cat = factor(cut(
      egfr2021,
      breaks = c(0, 10, 15, 20),
      labels = c("<10", "10-14", "15-20"),
      right = FALSE
    )),
    
    iron_cat = factor(
      fifelse(iron_iv==1, 1,
              fifelse(iron_po==1, 2, 0)),
      levels = c(0, 1, 2),
      labels = c("No iron", "IV iron", "PO iron")
    ),
    
    sbp_cat = factor(cut(
      sbp,
      breaks = c(0, 120, 140, 160, 1000),
      labels = c("<120", "120-139", "140-159", ">160"),
      right = FALSE
    )),
    
    dbp_cat = factor(cut(
      dbp,
      breaks = c(0, 80, 90, 100, 1000),
      labels = c("<80", "80-89", "90-99", ">100"),
      right = FALSE
    )),
    
    calendar_year_cat = factor(cut(
      calendar_year,
      breaks = c(0, 2013, 2018, 10000),
      labels = c("2007-2012", "2013-2017", "2018-2021"),
      right = FALSE
    )),
    
    female = factor(female),
    cancer = factor(cancer),
    ihd = factor(ihd),
    pvd = factor(pvd),
    hf = factor(hf),
    dm = factor(dm),
    scvd = factor(scvd),
    copd = factor(copd),
    cirr = factor(cirr),
    psycho = factor(psycho),
    acs = factor(acs),
    hyperten = factor(hyperten),
    vhd = factor(vhd),
    cevd = factor(cevd),
    af = factor(af),
    arrh = factor(arrh),
    lung = factor(lung),
    thrombo = factor(thrombo),
    liver = factor(liver),
    fracture = factor(fracture),
    aki = factor(aki),
    bblock = factor(bblock),
    calblock = factor(calblock),
    diuretic = factor(diuretic),
    rasi = factor(rasi),
    lipid = factor(lipid),
    phosbinder = factor(phosbinder),
    esa = factor(esa),
    vitamind = factor(vitamind),
    digoxin = factor(digoxin),
    vasodilator = factor(vasodilator),
    antiplatelet = factor(antiplatelet),
    anticoag = factor(anticoag),
    edu = factor(edu),
    clinic_level = factor(clinic_level),
    region = factor(region)
  )]
  
  # save cohort
  save(cohort_final, file = file.path(
    paste0(
      "Data/analysis_data_",
      cohort_name,
      "_new.Rdata"
    )
  ))
}


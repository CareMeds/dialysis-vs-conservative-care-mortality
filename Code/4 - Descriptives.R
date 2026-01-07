################################################################################
### Decision for dialysis versus conservative care
### PART 4 - Descriptive statistics
################################################################################

# remove history
rm(list = ls(all.names = TRUE))
knitr::opts_knit$set(root.dir = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/")
set.seed(1)

# set directory
setwd("P:/SCREAM2/SCREAM2_Research/Carolien Maas/")
results_path <- "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Results/"

# load libraries
library(patchwork) # combine figures
library(data.table)

# load functions
source("Code/utils/plots.R")
source("Code/utils/tables.R")
source("Code/utils/data_manipulation.R")

################################################################################
### Load data ##################################################################
################################################################################
load("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/cohort_with_weights.Rdata")
load("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/merged_ckd.Rdata")

################################################################################
### Time until dialysis ########################################################
################################################################################
id_name <- "LOPNR"
dia_dt <- merged_ckd[!is.na(krt_startdate) &
                       !is.na(krt_modality) &
                       krt_modality != "TX",
                     c(id_name, "krt_startdate", "krt_modality"), 
                     with = FALSE][, unique(.SD)]

# merge with cohort
time_to_dia_dt <- merge(baseline, 
                        dia_dt, 
                        by = id_name,
                        all.x = TRUE)
time_to_dia_dt[, `:=`(
  time_until_dialysis = krt_startdate - visit_date,
  time_until_HD = fifelse(krt_modality == "HD", krt_startdate - visit_date, NA),
  time_until_PD = fifelse(krt_modality == "PD", krt_startdate - visit_date, NA)
)]

# create time and event variable for dialysis
dialysis_df <- time_to_dia_dt[trt == 1]
dialysis_long <- dialysis_df[, .(LOPNR, krt_modality, 
                                 time_until_dialysis, 
                                 time_until_HD,
                                 time_until_PD)] |>
  tidyr::pivot_longer(
    cols = c(time_until_dialysis, time_until_HD, time_until_PD),
    names_to = "type",
    values_to = "time"
  ) |>
  dplyr::mutate(
    time_years = time / 365.25,
    type = dplyr::recode(type,
                         time_until_dialysis = "All Dialysis",
                         time_until_HD       = "Hemodialysis",
                         time_until_PD       = "Peritoneal Dialysis")
  ) |>
  dplyr::filter(is.finite((time_years))) # remove censored individuals (TODO: write)

time_hist <- ggplot2::ggplot(dialysis_long, ggplot2::aes(x = time_years, fill = type)) +
  ggplot2::geom_histogram(alpha = 0.6, 
                 binwidth = 1/12, 
                 boundary = 0,
                 color = "white",
                 linewidth = 0.1) +
  ggplot2::facet_wrap(~ type) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, max(dialysis_long$time_years, na.rm = TRUE), by = 1)
  ) +
  ggplot2::scale_fill_brewer(palette = "Set2") +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    legend.position = "none",
    panel.grid = ggplot2::element_blank(),
  ) +
  ggplot2::labs(
    x = "Time Until Dialysis (years)",
    y = "Count",
    title = ""
  )
ggplot2::ggsave(
  plot = time_hist,
  filename = paste0(results_path, "Supplemental/Figure_S2.png"),
  width = 10,
  height = 4,
  dpi = 300
)

################################################################################
### Baseline characteristics for full eligibility cohort and cohort after weighting
################################################################################
# Define pretty labels for the table
row_labels <- c(
  "Number of individuals",
  "Age (years) [median, IQR]",
  "Age in categories (%)",
  "65-69",
  "70-74",
  "75-79",
  ">=80",
  "Female (%)",
  "Davies comorbidity score [median, IQR]",
  "Davies comorbidity score category (%)",
  "<2",
  "2-4",
  ">4",
  "Region (%)",
  "Örebro/Uppsala", 
  "Other regions",
  "Södra",
  "Stockholm",  
  "Västra",
  "Clinic level (%)",
  "Local",
  "Regional",
  "Academic",
  "Calendar year in categories (%)",
  "2007-2012",
  "2013-2017",
  "2018-2021",
  "eGFR (ml/min/1.73 m2) [median, IQR]",
  "eGFR in categories (%)",
  "<10",
  "10-14",
  "15-20", 
  "SBP (mmHg) [mean, SD]",
  "SBP in categories (%)",
  "<120",
  "120-139",
  "140-159",
  ">160",
  "DBP (mmHg) [mean, SD]",
  "DBP in categories (%)",
  "<80",
  "80-89",
  "90-99",
  ">100",
  "Total calcium (mmol/L) [mean, SD]",
  "Phosphorus (mmol/L) [mean, SD]",
  "Albumin (g/L) [mean, SD]",
  "Haemoglobin (mmol/L) [mean, SD]",
  "Primary kidney disease (%)",
  "Diabetic nephropathy",
  "Hypertension",
  "Other",
  "Malignancy (%)",
  "Ischemic heart disease (%)",
  "Peripheral vascular disease (%)",
  "Heart failure (%)",
  "Diabetes mellitus (%)",
  "Systemic collagen vascular disease (%)",
  "COPD (%)",
  "Cirrhosis (%)",
  "Psychiatric illness (%)",
  "Acute coronary syndrome (%)",
  "Hypertension (%)",
  "Valvular heart disease (%)",
  "Other cerebrovascular disease (%)",
  "Atrial fibrillation (%)",
  "Other arrhythmias (%)",
  "Other lung disease (%)",
  "Venous thromboembolism (%)",
  "Liver disease (%)",
  "Fracture (%)",
  "Acute kidney injury (%)",
  "Beta blockers (%)",
  "Calcium channel blockers (%)",
  "Diuretics (%)",
  "Renin-angiotensin system inhibitors (%)",
  "Lipid lowering agents (%)",
  "Phosphate binder (%)",
  "Erythropoietin stimulating agents (%)",
  "Vitamin D (%)",
  "Digoxin (%)",
  "Vasodilator (%)",
  "Antiplatelet agents (%)",
  "Anticoagulants (%)",
  "Iron (%)",
  "No iron",
  "Intravenous",
  "Per os",
  "Hospitalizations in previous year [median, IQR]",
  "Cardiovascular hospitalizations in previous year [median, IQR]",
  "Education (%)"
)

# create descriptive statistics table for full eligible cohort
openxlsx::write.xlsx(
  create_baseline_table(
    data = elig_cohort,
    id_name = "LOPNR",
    weights = NULL,
    vars = listvar,
    categoricalVars = catvar,
    IQRVars = non_normal_vars,
    treatmentColumn = "S",
    # Column indicating if trt dec registered
    treatmentLabel = "Treatment registered",
    controlLabel = "Treatment not registered",
    tableCaption = paste("Baseline characteristics of all eligible patients")
  )$raw_table[, -c(3, 4)],
  rowNames = TRUE,
  file = paste0(
    results_path,
    "Supplemental/Descriptives_full_eligibility_cohort.xlsx"
  )
)
openxlsx::write.xlsx(
  create_baseline_table(
    data = elig_cohort,
    id_name = "LOPNR",
    weights = elig_cohort$sw_IPSW,
    vars = listvar,
    categoricalVars = catvar,
    IQRVars = non_normal_vars,
    treatmentColumn = "S",
    # Column indicating if trt dec registered
    treatmentLabel = "Treatment registered",
    controlLabel = "Treatment not registered",
    tableCaption = paste("Baseline characteristics of all eligible patients")
  )$raw_table[, -c(3, 4)],
  rowNames = TRUE,
  file = paste0(
    results_path,
    "Supplemental/Descriptives_full_eligibility_cohort_IPSW.xlsx"
  )
)

# Create baseline table without weighting and with the four methods of weighting
for (w_meth in w_meths) {
  if (w_meth == "") {
    weights_meth <- NULL
  } else {
    weights_meth <- baseline[[paste0("sw_", w_meth)]]
  }
  
  # create main descriptives table
  if (w_meth=="" | w_meth=="IPTW"){
    table_one <- create_baseline_table(
      data = baseline,
      id_name = "LOPNR",
      weights = weights_meth,
      vars = listvar_main,
      categoricalVars = catvar,
      IQRVars = non_normal_vars,
      treatmentColumn = "trt",
      # Column in dataframe with trt assignment
      treatmentLabel = treatment_label,
      controlLabel = control_label,
      tableCaption = paste(
        "Baseline characteristics of study patients stratified",
        " by the decision on dialysis vs. conservative care, ",
        "weighting method:",
        w_meth
      )
    )
    
    openxlsx::write.xlsx(
      table_one$raw_table,
      rowNames = TRUE,
      file = paste0(
        results_path,
        "Main/Descriptives_",
        ifelse(w_meth == "", "no_weighting", paste0("weighting_", w_meth)),
        ".xlsx"
      )
    )
  }
  
  # create descriptive statistics table
  if (w_meth != "IPSW" & w_meth != "IPSW_IPTW") {
    table_one <- create_baseline_table(
      data = baseline,
      id_name = "LOPNR",
      weights = weights_meth,
      vars = listvar,
      categoricalVars = catvar,
      IQRVars = non_normal_vars,
      treatmentColumn = "trt",
      # Column in dataframe with trt assignment
      treatmentLabel = treatment_label,
      controlLabel = control_label,
      tableCaption = paste(
        "Baseline characteristics of study patients stratified",
        " by the decision on dialysis vs. conservative care, ",
        "weighting method:",
        w_meth
      ),
      tableRowLabels = row_labels
    )
    
    # save table
    openxlsx::write.xlsx(
      table_one$raw_table,
      rowNames = TRUE,
      file = paste0(
        results_path,
        ifelse(w_meth == "" |
                 w_meth == "IPTW", "Supplemental", "Other"),
        "/Descriptives_",
        ifelse(w_meth == "", "no_weighting", paste0("weighting_", w_meth)),
        ".xlsx"
      )
    )
  }
  
  # save SMD values for love plot
  assign(paste0("SMD_", w_meth), table_one$smd_table)
}

################################################################################
### Create love plots of SMDs before and after weighting #######################
################################################################################
SMDs_dt <- cbind(SMD_, SMD_IPTW)
rownames(SMDs_dt) <- varnames
SMDs_dt <- SMDs_dt[order(SMDs_dt[, "SMD_IPTW"], decreasing = TRUE), ]
ggplot2::ggsave(
  plot = love_plot(
    SMDs_dt = SMDs_dt,
    SMD_names = c("SMD_", "SMD_IPTW"),
    plot_title = "Before and after IPTW",
    xmax = 1
  ),
  filename = paste0(results_path, "Supplemental/Figure_S5.png"),
  width = 15,
  height = 15,
  dpi = 300
)

################################################################################
### Kidney transplantation during follow-up
################################################################################
time_to_dia_dt$krt_startdate
table(time_to_dia_dt$krt_modality)
trans_dt <- merged_ckd[LOPNR %in% baseline$LOPNR &
             krt_modality=="TX" &
             !is.na(krt_startdate),
           .(LOPNR, krt_startdate, krt_modality)
           ][,
             unique(.SD)
             ]
baseline_TX <- merge(baseline,
                     trans_dt,
                     all.x = TRUE,
                     by = id_name)
baseline_TX[, 
            time_until_TX := krt_startdate - visit_date
]

# check
# baseline_TX[baseline_TX$time_until_TX<0,]
# baseline_TX[baseline_TX$time_until_TX>0 & baseline_TX$time_until_TX< 2 * 365.25, 
#             c(id_name, "trt", "visit_date", "krt_startdate", "krt_modality", "time_until_TX"),
#             with = FALSE]
# baseline_TX[baseline_TX$time_until_TX>0 & baseline_TX$time_until_TX< 4 * 365.25, 
#             c(id_name, "trt", "visit_date", "krt_startdate", "krt_modality", "time_until_TX"),
#             with = FALSE]

cat("Number of transplantations in 2 years after dialysis decision:",
    nrow(unique(baseline_TX[baseline_TX$time_until_TX>0 & 
                              baseline_TX$time_until_TX< 2 * 365.25, "LOPNR"])), "\n")

################################################################################
### Summarize number of decisions
################################################################################
# 1. select patients in baseline that have a second or third decision, but not a third decision
# 2. sort on decision date
# 3. keep only modality changes

# two decisions
two_decisions_dt <- merged_ckd[
  LOPNR %in% baseline$LOPNR &
    !is.na(decision_date1) & 
    !is.na(decision_date2) &
    !is.na(decision_modality1) &
    !is.na(decision_modality2),
  .(LOPNR, decision_date1, decision_date2,
    decision_modality1, decision_modality2)
][,
  # keep unique rows
  unique(.SD)
][
  # ensure ordered by date within patient
  order(LOPNR, decision_date1)
][
  # keep only if second decision is within 2 years of first
  decision_date2 <= decision_date1 + lubridate::years(2)
][, `:=`
  (decision_modality1_new = fifelse(decision_modality1=="Konservativ behandling", 0, 1),
    decision_modality2_new = fifelse(decision_modality2=="Konservativ behandling", 0, 1))
][
  # keep only if modality changes
  decision_modality1_new != decision_modality2_new
]
baseline[, n_decision_2 := fifelse(baseline$LOPNR %in% two_decisions_dt$LOPNR, 1, 0)]

# three decisions
three_decisions_dt <- merged_ckd[
  LOPNR %in% baseline$LOPNR &
    !is.na(decision_date1) & 
    !is.na(decision_date2) &
    !is.na(decision_date3) &
    !is.na(decision_modality1) &
    !is.na(decision_modality2) &
    !is.na(decision_modality3),
  .(LOPNR, decision_date1, decision_date2, decision_date3,
    decision_modality1, decision_modality2, decision_modality3)
][,
  # keep unique rows
  unique(.SD)
][
  # ensure ordered by date within patient
  order(LOPNR, decision_date1)
][
  # keep only if second decision is within 2 years of first
  decision_date2 <= decision_date1 + lubridate::years(2) &
  decision_date3 <= decision_date1 + lubridate::years(2)
][, `:=`
  (decision_modality1_new = fifelse(decision_modality1=="Konservativ behandling", 0, 1),
    decision_modality2_new = fifelse(decision_modality2=="Konservativ behandling", 0, 1),
    decision_modality3_new = fifelse(decision_modality3=="Konservativ behandling", 0, 1))
][
  # keep only if modality changes
  decision_modality1_new != decision_modality2_new & 
    decision_modality2_new != decision_modality3_new
]
baseline[, n_decision_3 := fifelse(baseline$LOPNR %in% three_decisions_dt$LOPNR, 1, 0)]
three_decisions_dt[LOPNR==176579584 | LOPNR==495890215,
                   .(LOPNR, decision_date1, decision_modality1,
                                       decision_date2, decision_modality2,
                                       decision_date3, decision_modality3)]

# Generate the summary table
ndec <- summarize_binary_vars(
  baseline,
  group_var = "trt",
  vars = c("n_decision_2", "n_decision_3"),
  labels = c(
    "Two treatment switches",
    "Three treatment switches"
  )
)

# Rename columns for readability
names(ndec) <- c("Variable", "Overall", "Dialysis", "Conservative Care")
print(ndec)

# save table
openxlsx::write.xlsx(ndec,
                     rowNames = FALSE,
                     file = paste0(results_path, "Main/Number_of_decisions.xlsx"))

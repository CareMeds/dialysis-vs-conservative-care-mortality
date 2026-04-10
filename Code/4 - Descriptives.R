################################################################################
### Decision for dialysis versus conservative care
### PART 4 - Descriptive statistics
################################################################################

# remove history
rm(list = ls(all.names = TRUE))
knitr::opts_knit$set(root.dir = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/")
set.seed(1)

# set directory
setwd(
  "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/"
)
results_path <- "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/Results/"

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
load("Data/cohort_with_weights.Rdata")
load("Data/merged_ckd.Rdata")

################################################################################
### Histogram for time until dialysis ##########################################
################################################################################
id_name <- "LOPNR"
dia_dt <- merged_ckd[!is.na(krt_startdate) &
                       !is.na(krt_modality) &
                       krt_modality != "TX", c(id_name, "krt_startdate", "krt_modality"), with = FALSE][, unique(.SD)]

# merge with cohort
time_to_dia_dt <- merge(baseline, dia_dt, by = id_name, all.x = TRUE)

# compute time to dialysis
time_to_dia_dt[, `:=`(
  time_until_dialysis = lubridate::time_length(lubridate::interval(visit_date, krt_startdate), "years"),
  time_until_dialysis_days = lubridate::time_length(lubridate::interval(visit_date, krt_startdate), "days"),
  time_until_HD = fifelse(
    krt_modality == "HD",
    lubridate::time_length(lubridate::interval(visit_date, krt_startdate), "years"),
    NA
  ),
  time_until_PD = fifelse(
    krt_modality == "PD",
    lubridate::time_length(lubridate::interval(visit_date, krt_startdate), "years"),
    NA
  )
)]

# create time and event variable for dialysis
dialysis_df <- time_to_dia_dt[trt == 1]
dialysis_long <- dialysis_df[, .(LOPNR,
                                 krt_modality,
                                 time_until_dialysis,
                                 time_until_HD,
                                 time_until_PD)] |>
  tidyr::pivot_longer(
    cols = c(time_until_dialysis, time_until_HD, time_until_PD),
    names_to = "type",
    values_to = "time"
  ) |>
  dplyr::mutate(
    time_years = time,
    type = dplyr::recode(
      type,
      time_until_dialysis = "All Dialysis",
      time_until_HD       = "Hemodialysis",
      time_until_PD       = "Peritoneal Dialysis"
    )
  ) |>
  dplyr::filter(is.finite((time_years)))

# histogram
time_hist <- ggplot2::ggplot(dialysis_long, ggplot2::aes(x = time_years, fill = type)) +
  ggplot2::geom_histogram(
    alpha = 0.6,
    binwidth = 1 / 12,
    boundary = 0,
    color = "white",
    linewidth = 0.1
  ) +
  ggplot2::facet_wrap( ~ type) +
  ggplot2::scale_x_continuous(breaks = seq(0, max(dialysis_long$time_years, na.rm = TRUE), by = 1)) +
  ggplot2::scale_fill_brewer(palette = "Set1") +
  ggplot2::theme_minimal(base_size = 18) +
  ggplot2::theme(
    legend.position = "none",
    panel.grid = ggplot2::element_blank(),
    plot.tag = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(hjust = 0.5)
  ) +
  ggplot2::labs(x = "Years", y = "Count", title = "Time until dialysis") +
  ggplot2::labs(tag = "A")
time_hist_zoomed <- time_hist +
  ggplot2::scale_x_continuous(
    limit = c(0, 2),
    breaks = seq(0, 2, 0.5),
    labels = c(0, 6, 12, 18, 24)
  ) +
  ggplot2::labs(tag = "B", title = "Time until dialysis, zoomed in on first two years", x = "Months") +
  ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(hjust = 0.5)
  )
ggplot2::ggsave(
  plot = time_hist / time_hist_zoomed,
  filename = paste0(results_path, "Supplemental/Figure_S4.png"),
  width = 10,
  height = 8,
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
table_full_elig <- create_baseline_table(
  data = elig_cohort,
  id_name = "LOPNR",
  weights = NULL,
  vars = listvar,
  categoricalVars = catvar,
  IQRVars = non_normal_vars,
  treatmentColumn = "S",
  treatmentLabel = "Treatment registered",
  controlLabel = "Treatment not registered",
  tableCaption = paste("Baseline characteristics of all eligible patients"),
  tableRowLabels = row_labels
)$raw_table

elig_cohort$sw_IPSW <- 1
elig_cohort[elig_cohort$S == 1, "sw_IPSW"] <- baseline$sw_IPSW
table_full_elig_IPSW <- create_baseline_table(
  data = elig_cohort,
  id_name = "LOPNR",
  weights = elig_cohort$sw_IPSW,
  vars = listvar,
  categoricalVars = catvar,
  IQRVars = non_normal_vars,
  treatmentColumn = "S",
  treatmentLabel = "Treatment registered",
  controlLabel = "Treatment not registered",
  tableCaption = paste("Baseline characteristics of all eligible patients, IPSW"),
  tableRowLabels = row_labels
)$raw_table

# Create baseline table without weighting and with the four methods of weighting
for (w_meth in w_meths) {
  if (w_meth == "unweighted") {
    weights_meth <- NULL
  } else {
    weights_meth <- baseline[[paste0("sw_", w_meth)]]
  }
  
  # create main descriptives table
  if (w_meth == "unweighted" | w_meth == "IPTW") {
    table_one <- create_baseline_table(
      data = baseline,
      id_name = "LOPNR",
      weights = weights_meth,
      vars = listvar_main,
      categoricalVars = catvar,
      IQRVars = non_normal_vars,
      treatmentColumn = trt_var,
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
        ifelse(w_meth == "unweighted", "no_weighting", paste0("weighting_", w_meth)),
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
      treatmentColumn = trt_var,
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
        ifelse(
          w_meth == "unweighted" |
            w_meth == "IPTW",
          "Supplemental",
          "Other"
        ),
        "/Descriptives_",
        ifelse(
          w_meth == "unweighted",
          "no_weighting",
          paste0("weighting_", w_meth)
        ),
        ".xlsx"
      )
    )
  }
  
  # save SMD values for love plot
  assign(paste0("SMD_", w_meth), table_one$smd_table)
}

################################################################################
### Create love plots of SMDs before and after IPTW ############################
################################################################################
# set colors
manual_colors <- RColorBrewer::brewer.pal(n = 4, name = "Set1")

# set names
varnames <- c(
  "Age",
  "Age category",
  "Female",
  "Davies score",
  "Davies score category",
  "Region",
  "Clinic level",
  "Calendar year category",
  "eGFR",
  "eGFR category",
  "SBP",
  "SBP category",
  "DBP",
  "DBP category",
  "Total calcium",
  "Phosphorus",
  "Albumin",
  "Haemoglobin",
  "Primary kidney disease",
  "Malignancy",
  "Ischemic heart disease",
  "Peripheral vascular disease",
  "Heart failure",
  "Diabetes mellitus",
  "Systemic collagen vascular disease",
  "COPD",
  "Cirrhosis",
  "Psychiatric illness",
  "Acute coronary syndrome",
  "Hypertension",
  "Valvular heart disease",
  "Other cerebrovascular disease",
  "Atrial fibrillation",
  "Other arrhythmias",
  "Other lung disease",
  "Venous thromboembolism",
  "Liver disease",
  "Fracture",
  "Acute kidney injury",
  "Beta blockers",
  "Calcium channel blockers",
  "Diuretic",
  "Renin-angiotensin system inhibitors",
  "Lipid lowering agents",
  "Phosphate binder",
  "Erythropoietin stimulating agents",
  "Vitamin D",
  "Digoxin",
  "Vasodilator",
  "Antiplatelet agents",
  "Anticoagulants",
  "Iron",
  "Hospitalizations in previous year",
  "Cardiovascular hospitalizations in previous year",
  "Education"
)

# combine in dt
SMDs_dt <- cbind(SMD_unweighted, SMD_IPTW)
rownames(SMDs_dt) <- varnames
SMDs_dt <- SMDs_dt[order(SMDs_dt[, "SMD_IPTW"], decreasing = TRUE), ]
ggplot2::ggsave(
  plot = love_plot(
    SMDs_dt = SMDs_dt,
    SMD_names = c("SMD_unweighted", "SMD_IPTW"),
    plotColors = manual_colors,
    xmax = 1
  ),
  filename = paste0(results_path, "Supplemental/Figure_S3.png"),
  width = 15,
  height = 15,
  dpi = 300
)

################################################################################
# calculate SMD for those with a treatment decision versus all eligible patients
################################################################################
# encode factors
elig_cohort_num <- encode_factors(
  dt = elig_cohort,
  catvar = catvar,
  contvar = contvar,
  expand = FALSE
)
means_all <- colMeans(elig_cohort_num[, ..listvar])
sd_all <- apply(elig_cohort_num[, ..listvar], 2, sd)

# encode factors
baseline_num <- encode_factors(
  dt = baseline,
  catvar = c(catvar, trt_var),
  contvar = contvar,
  expand = FALSE
)
means <- colMeans(baseline_num[, ..listvar])
means_treated <- colMeans(baseline_num[trt == 1, ..listvar])
means_untreated <- colMeans(baseline_num[trt == 0, ..listvar])
wmeans <- apply(baseline_num[, ..listvar], 2, weighted.mean, w = as.numeric(baseline[, "sw_IPSW"][[1]]))
wmeans_treated <- apply(baseline_num[trt == 1, ..listvar], 2, weighted.mean, w = as.numeric(baseline[trt ==
                                                                                                       1, "sw_IPSW"][[1]]))
wmeans_untreated <- apply(baseline_num[trt == 0, ..listvar], 2, weighted.mean, w = as.numeric(baseline[trt ==
                                                                                                         0, "sw_IPSW"][[1]]))

# balance before weighting
SMD_elig <- abs(means - means_all) / sd_all
SMD_elig_treated <- abs(means_treated - means_all) / sd_all
SMD_elig_untreated <- abs(means_untreated - means_all) / sd_all

# balance after weighting
SMD_elig_wt <- abs(wmeans - means_all) / sd_all
SMD_elig_treated_wt <- abs(wmeans_treated - means_all) / sd_all
SMD_elig_untreated_wt <- abs(wmeans_untreated - means_all) / sd_all

# save table
SMD_eligs_dt <- data.frame(
  SMD_elig = SMD_elig,
  untreated_unweighted = SMD_elig_untreated,
  treated_unweighted = SMD_elig_treated,
  SMD_elig_IPSW = SMD_elig_wt,
  untreated_weighted = SMD_elig_untreated_wt,
  treated_weighted = SMD_elig_treated_wt
)
row_names_without_levels <- row_labels[-c(1,
                                          4:7,
                                          11:13,
                                          15:19,
                                          21:23,
                                          25:27,
                                          30:32,
                                          35:38,
                                          41:44,
                                          50:52,
                                          86:88)]
rownames(SMD_eligs_dt) <- varnames

# save plot
SMD_eligs_dt <- SMD_eligs_dt[order(SMD_eligs_dt$SMD_elig_IPSW, decreasing = TRUE), ]
ggplot2::ggsave(
  plot = love_plot(
    SMDs_dt = SMD_eligs_dt,
    SMD_names = c("SMD_elig", "SMD_elig_IPSW"),
    plotColors = manual_colors[1:2],
    xlab_title = "Target absolute standardized mean difference",
    xmax = 1
  ),
  filename = paste0(results_path, "Supplemental/Figure_S9_IPSW.png"),
  width = 15,
  height = 15,
  dpi = 300
)

# descriptives full eligibility
table_full_elig[, "SMD"] <- ""
table_full_elig[rownames(table_full_elig) %in% row_names_without_levels, "SMD"] <- sprintf("%.3f", SMD_elig)
openxlsx::write.xlsx(
  table_full_elig,
  rowNames = TRUE,
  file = paste0(
    results_path,
    "Supplemental/Descriptives_full_eligibility_cohort.xlsx"
  )
)
table_full_elig_IPSW[, "SMD"] <- ""
table_full_elig_IPSW[rownames(table_full_elig_IPSW) %in% row_names_without_levels, "SMD"] <- sprintf("%.3f", SMD_elig_wt)
openxlsx::write.xlsx(
  table_full_elig_IPSW,
  rowNames = TRUE,
  file = paste0(
    results_path,
    "Supplemental/Descriptives_full_eligibility_cohort_IPSW.xlsx"
  )
)

################################################################################
### Kidney transplantation during follow-up
################################################################################
time_to_dia_dt$krt_startdate
table(time_to_dia_dt$krt_modality)
trans_dt <- merged_ckd[LOPNR %in% baseline$LOPNR &
                         krt_modality == "TX" &
                         !is.na(krt_startdate), .(LOPNR, krt_startdate, krt_modality)][, unique(.SD)]
baseline_TX <- merge(baseline, trans_dt, all.x = TRUE, by = id_name)
baseline_TX[, time_until_TX := krt_startdate - visit_date]

################################################################################
### Summarize number of decisions
################################################################################
# 1. select patients in baseline that have a second or third decision, but not a third decision
# 2. sort on decision date
# 3. keep only modality changes

# extract those with two decisions
# 1. extract relevant variables
two_decisions_dt <- merged_ckd[LOPNR %in% baseline$LOPNR &
                                 !is.na(decision_date1) &
                                 !is.na(decision_date2) &
                                 !is.na(decision_modality1) &
                                 !is.na(decision_modality2), .(LOPNR,
                                                               decision_date1,
                                                               decision_date2,
                                                               decision_modality1,
                                                               decision_modality2)]
# 2. keep unique rows and ensure ordered by date within patient
two_decisions_dt <- two_decisions_dt[, unique(.SD)][order(LOPNR, decision_date1)]
# 3. keep only if second decision is within 2 years of first
two_decisions_dt <- two_decisions_dt[decision_date2 <= decision_date1 + lubridate::years(2)]
# 4. keep only if modality changes
two_decisions_dt <- two_decisions_dt[, `:=`
                                     (
                                       decision_modality1_new = fifelse(decision_modality1 == "Konservativ behandling", 0, 1),
                                       decision_modality2_new = fifelse(decision_modality2 == "Konservativ behandling", 0, 1)
                                     )][decision_modality1_new != decision_modality2_new]
baseline[, n_decision_2 := fifelse(baseline$LOPNR %in% two_decisions_dt$LOPNR, 1, 0)]

# three decisions
# 1. extract relevant variables
three_decisions_dt <- merged_ckd[LOPNR %in% baseline$LOPNR &
                                   !is.na(decision_date1) &
                                   !is.na(decision_date2) &
                                   !is.na(decision_date3) &
                                   !is.na(decision_modality1) &
                                   !is.na(decision_modality2) &
                                   !is.na(decision_modality3), .(
                                     LOPNR,
                                     decision_date1,
                                     decision_date2,
                                     decision_date3,
                                     decision_modality1,
                                     decision_modality2,
                                     decision_modality3
                                   )]
# 2. keep unique rows and ensure ordered by date within patient
three_decisions_dt <- three_decisions_dt[, unique(.SD)][order(LOPNR, decision_date1)]
# 3. keep only if second decision is within 2 years of first
three_decisions_dt <- three_decisions_dt[decision_date2 <= decision_date1 + lubridate::years(2) &
                                           decision_date3 <= decision_date1 + lubridate::years(2)]
# 4. keep only if modality changes
three_decisions_dt <- three_decisions_dt[, `:=`
                                         (
                                           decision_modality1_new = fifelse(decision_modality1 == "Konservativ behandling", 0, 1),
                                           decision_modality2_new = fifelse(decision_modality2 == "Konservativ behandling", 0, 1),
                                           decision_modality3_new = fifelse(decision_modality3 == "Konservativ behandling", 0, 1)
                                         )][decision_modality1_new != decision_modality2_new &
                                              decision_modality2_new != decision_modality3_new]
baseline[, n_decision_3 := fifelse(baseline$LOPNR %in% three_decisions_dt$LOPNR, 1, 0)]
three_decisions_dt[LOPNR == 176579584 | LOPNR == 495890215, .(
  LOPNR,
  decision_date1,
  decision_modality1,
  decision_date2,
  decision_modality2,
  decision_date3,
  decision_modality3
)]

stats_dialysis <- data.frame(
  names = c("Maximum time until dialysis",
            "# patients choosing dialysis", 
            "Treatment switches for dialysis",
            "# treatment switches for dialysis",
            "Treatment switches for conservative care",
            "# treatment switches for conservative care",
            "# patients that chose dialysis starting dialysis within two years",
            "# patients that chose PD starting dialysis within two years",
            "# patients that chose HD starting dialysis within two years",
            "# patients that chose dialysis starting dialysis within 3 months",
            "# patients that chose dialysis starting dialysis between 3-6 months",
            "# patients that chose dialysis starting dialysis between 6-12 months",
            "# patients that chose dialysis starting dialysis between 12-24 months",
            "# patients that chose dialysis but died before starting dialysis within two years",
            "Number of transplantations in two years after dialysis decision:",
            "% transplantations in two years after dialysis decision:"
  ),
  counts = c(max(dialysis_df$time_until_dialysis, na.rm=TRUE),
             nrow(dialysis_df),
             sum(baseline$n_decision_2==1 & baseline$trt==1),
             sum(baseline$n_decision_2==1 & baseline$trt==1) / sum(baseline$trt==1) * 100,
             sum(baseline$n_decision_2==1 & baseline$trt==0),
             sum(baseline$n_decision_2==1 & baseline$trt==0) / sum(baseline$trt==0) * 100,
             dialysis_df[!is.na(time_until_dialysis) &
                           time_until_dialysis <= 2, .N],
             dialysis_df[!is.na(time_until_PD) &
                           time_until_PD <= 2, .N],
             dialysis_df[!is.na(time_until_HD) &
                           time_until_HD <= 2, .N],
             dialysis_df[!is.na(time_until_dialysis) &
                           time_until_dialysis <= 3 / 12, .N],
             dialysis_df[!is.na(time_until_dialysis) &
                           time_until_dialysis > 3 / 12 &
                           time_until_dialysis <= 0.5, .N],
             dialysis_df[!is.na(time_until_dialysis) &
                           time_until_dialysis > 0.5 &
                           time_until_dialysis <= 1, .N],
             dialysis_df[!is.na(time_until_dialysis) &
                           time_until_dialysis > 1 &
                           time_until_dialysis <= 2, .N],
             dialysis_df[time2event_death_2y < time_until_dialysis_days, .N],
             nrow(unique(baseline_TX[baseline_TX$time_until_TX > 0 &
                                       baseline_TX$time_until_TX < 2 * 365, "LOPNR"])),
             nrow(unique(baseline_TX[baseline_TX$time_until_TX > 0 &
                                       baseline_TX$time_until_TX < 2 * 365, "LOPNR"])) / sum(baseline$trt==1) * 100
  )
)
openxlsx::write.xlsx(
  stats_dialysis,
  rowNames = FALSE,
  colNames = FALSE,
  file = paste0(results_path, "Main/Statistics_text.xlsx")
)

# save cohort
save(
  id_name,
  listvar,
  listvar_main,
  catvar,
  contvar,
  non_normal_vars,
  treatment_label,
  control_label,
  baseline,
  model_PS,
  model_S,
  coef_PS_overall,
  elig_cohort,
  w_meths,
  trt_var,
  manual_colors,
  file = file.path("Data/cohort_with_weights.Rdata")
)

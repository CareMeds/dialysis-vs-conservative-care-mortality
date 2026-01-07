################################################################################
### Decision for dialysis versus conservative care
### PART 3 - Analysis
################################################################################

# remove history
rm(list = ls(all.names = TRUE))
knitr::opts_knit$set(root.dir = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/")
set.seed(1)

# set directory
setwd("P:/SCREAM2/SCREAM2_Research/Carolien Maas/")
results_path <- "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Results/"

# load libraries
library(data.table)
library(patchwork) # combine figures

# load functions
source("Code/utils/weighting.R")
source("Code/utils/plots.R")
source("Code/utils/tables.R")
source("Code/utils/data_manipulation.R")
source("Code/utils/outcomes_absolute_risks.R")

################################################################################
### Load data ##################################################################
################################################################################
load("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/analysis_data_elig_cohort_new.Rdata")
elig_cohort <- cohort_final
load("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/analysis_data_cohort_new.Rdata")
baseline <- cohort_final

# Set variables to include in the baseline table
listvar_main <- c(
  "age",
  "female",
  "Davies_score",
  "egfr2021",
  "egfr_cat",
  "sbp",
  "dbp",
  "calcium_total",
  "phosphate",
  "albumin",
  "hb",
  "cancer",
  "ihd",
  "pvd",
  "hf",
  "dm",
  "psycho",
  "hyperten",
  "aki",
  "bblock",
  "calblock",
  "diuretic",
  "rasi",
  "lipid",
  "phosbinder",
  "esa",
  "antiplatelet",
  "anticoag",
  "iron_cat",
  "n_hospital",
  "edu"
)

listvar <- c(
  "age",
  "age_cat",
  "female",
  "Davies_score",
  "Davies_score_cat",
  "region",
  "clinic_level",
  "calendar_year_cat",
  "egfr2021",
  "egfr_cat",
  "sbp",
  "sbp_cat",
  "dbp",
  "dbp_cat",
  "calcium_total",
  "phosphate",
  "albumin",
  "hb",
  "prd_cat",
  "cancer",
  "ihd",
  "pvd",
  "hf",
  "dm",
  "scvd",
  "copd",
  "cirr",
  "psycho",
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
  "aki",
  "bblock",
  "calblock",
  "diuretic",
  "rasi",
  "lipid",
  "phosbinder",
  "esa",
  "vitamind",
  "digoxin",
  "vasodilator",
  "antiplatelet",
  "anticoag",
  "iron_cat",
  "n_hospital",
  "n_cvd_hospital",
  "edu"
)

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

# define which variables are continuous and which are categorical
contvar <- c(
  "age",
  "egfr2021",
  "Davies_score",
  "sbp",
  "dbp",
  "calcium_total",
  "phosphate",
  "albumin",
  "hb",
  "n_hospital",
  "n_cvd_hospital"
  )
catvar <- listvar[!listvar %in% contvar]

# define labels
treatment_label <- "Dialysis"
control_label <- "Conservative care"
non_normal_vars <- c("age", 
                     "Davies_score",
                     "egfr2021",
                     "n_hospital",
                     "n_cvd_hospital") # Nonnormally distributed

# Define the RHS of the formula once
rhs_formula <- ~ rms::pol(age, 2) + age_cat +
  female +
  rms::pol(egfr2021, 2) + egfr_cat +
  Davies_score_cat +
  calendar_year_cat +
  sbp + sbp_cat +
  dbp + dbp_cat +
  calcium_total +
  phosphate +
  albumin +
  hb +
  prd_cat +
  cancer +
  ihd +
  pvd +
  hf +
  dm +
  scvd +
  copd +
  cirr +
  psycho +
  acs +
  hyperten +
  vhd +
  cevd +
  af +
  arrh +
  lung +
  thrombo +
  liver +
  fracture +
  aki +
  bblock +
  calblock +
  diuretic +
  rasi +
  lipid +
  phosbinder +
  esa +
  vitamind +
  digoxin +
  vasodilator +
  antiplatelet +
  anticoag +
  n_hospital +
  n_cvd_hospital +
  edu +
  iron_cat +
  clinic_level +
  region 

################################################################################
### Propensity score model #####################################################
################################################################################
# Define the propensity score model
model_PS <- update(rhs_formula, trt ~ . +
                     vasodilator:rms::pol(age, 2))  # TODO: write about it

# Create IPTW weights on full cohort
out_weights <- create_weights(
  data = baseline,
  model_PS = model_PS,
  w_meth = "IPTW",
  verbose = FALSE
)

# add weights to data
baseline <- out_weights$data

# check if SMDs after weighting below 0.1
id_name <- "LOPNR"
table_one <- create_baseline_table(
  data = baseline,
  id_name = id_name,
  weights = baseline$w,
  vars = listvar,
  categoricalVars = catvar,
  IQRVars = non_normal_vars,
  treatmentColumn = "trt",
  treatmentLabel = treatment_label,
  controlLabel = control_label,
  tableCaption = ""
)
print(table_one$smd_table[which(table_one$smd_table>0.1)])

# save coefficients of PS model
coef_PS_overall <- out_weights$coef_ps
PS_df <- data.frame(fmt_ci(
  as.numeric(out_weights$coef_ps),
  as.numeric(out_weights$lower_ps),
  as.numeric(out_weights$upper_ps),
  digits = 2
),
fmt_ci(exp(as.numeric(
  out_weights$coef_ps
)), exp(as.numeric(
  out_weights$lower_ps
)), exp(as.numeric(
  out_weights$upper_ps
))))
rownames(PS_df) <- c(
  "Intercept",
  "Age",
  "Age^2",
  "Age category 70-74 versus 65-69",
  "Age category 75-79 versus 65-69",
  "Age category >=80 versus 65-69",
  "Female",
  "eGFR",
  "eGFR^2",
  "eGFR category 10-14 versus < 10",
  "eGFR category 15-20 versus < 10",
  "Davies score category 2-4 versus < 2",
  "Davies score category >4 versus < 2",
  "Calendar year 2013-2017 versus 2007-2012",
  "Calendar year 2018-2021 versus 2007-2012",
  "SBP",
  "SBP category 120-139 versus <120",
  "SBP category 140-159 versus <120",
  "SBP category >160 versus <120",
  "DBP",
  "DBP category 80-89 versus <80",
  "DBP category 90-99 versus <80",
  "DBP category >100 versus <80",
  "Total calcium, mmol/L",
  "Phosphorus, mmol/L",
  "Albumin, g/L",
  "Haemoglobin , mmol/L",
  "Primary kidney disease hypertension versus diabetic nephropathy",
  "Primary kidney disease other versus diabetic nephropathy",
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
  "Diuretics",
  "Renin-angiotensin system inhibitors",
  "Lipid lowering agents",
  "Phosphate binder",
  "Erythropoietin stimulating agents",
  "Vitamin D",
  "Digoxin",
  "Vasodilator",
  "Antiplatelet agents",
  "Anticoagulants",
  "Number of hospitalizations in previous year",
  "Number of cardiovascular hospitalizations in previous year",
  "Education",
  "Iron intravenous versus no iron",
  "Iron per os versus no iron",
  "Clinic level regional versus local",
  "Clinic level academic versus local",
  "Other versus Örebro/Uppsala",
  "Södra versus Örebro/Uppsala",
  "Stockholm versus Örebro/Uppsala",
  "Västra versus Örebro/Uppsala",
  "Age * vasodilator",
  "Age^2 * vasodoliator"
)
colnames(PS_df) <- c("Coefficients (95% CI)", "HR (95% CI)")
# TODO: exploding vasodilator coefficient
openxlsx::write.xlsx(
  PS_df,
  rowNames = TRUE,
  file = paste0(results_path, "Other/Table_Propensity_Score_Model.xlsx")
)

# variance inflation factor
car::vif(glm(model_PS, family = "binomial", data = baseline), type = "terms")

# C-statistic of propensity score model is 0.86
pROC::auc(pROC::roc(
  baseline$trt,
  predict(glm(
    model_PS, family = "binomial", data = baseline
  )),
  levels = c(0, 1),
  direction = "<"
))

################################################################################
### Generalizability model #####################################################
################################################################################
# Define S = 1 if in analysis data set and 0 otherwise
elig_cohort$S <- ifelse(elig_cohort[, id_name, with = FALSE][[1]] %in% 
                          baseline[, id_name, with = FALSE][[1]], 1L, 0L)

################################################################################
### Compute weights using several techniques ###################################
################################################################################
# compute for all weighting methods
# 1) IPTW (ATE)
# 2) SMR (ATT)
# 3) Fine stratification weights (ATE)
# 4) Fine stratification weights (ATT)
# 5) Overlap weighting
# 6) IPTW and probability of belonging to eligibility cohort (generalizability)
w_meths <- c("", "IPTW", "overlap", "IPSW", "IPSW_IPTW", "SMR_ATT", "SMR_ATU")
for (w_meth in w_meths) {
  # create weights
  if (w_meth != "") {
    out_weights <- create_weights(
      data = baseline,
      elig_cohort = elig_cohort,
      model_PS = model_PS,
      w_meth = w_meth,
      catvar = catvar,
      contvar = contvar,
      verbose = FALSE
    )
    baseline[[paste0("sw_", w_meth)]] <- as.numeric(out_weights$data$w)
  }
  
  # set weights for each method
  if (w_meth == "") {
    weights_meth <- NULL
  } else {
    weights_meth <- baseline[[paste0("sw_", w_meth)]]
  }
}

# Summarize statistics on weights across methods
summ_weights <- data.table::rbindlist(lapply(w_meths[-1], function(w_meth) {
  # build a small data.table with trt and the corresponding weight
  weights_dt <- baseline[, .(trt, weights = get(paste0("sw_", w_meth)))]
  
  # summarize weights (data.table with Statistic, Overall, Control, Treated)
  summ_dt <- summarize.weights(weights_dt)
  
  # add Method column for this iteration
  summ_dt[, Method := w_meth]
  
  # reorder columns: Method first
  data.table::setcolorder(summ_dt,
                          c("Method", "Statistic", "Overall", "Control", "Treated"))
}), fill = TRUE)

# Write to Excel (one worksheet)
openxlsx::write.xlsx(
  x = summ_weights,
  file = file.path(results_path, "Supplemental/Table_S4.xlsx"),
  rowNames = FALSE,
  overwrite = TRUE
)

################################################################################
# calculate SMD for those with a treatment decision versus all eligible patients
################################################################################
# encode factors
elig_cohort_num <- encode_factors(dt = elig_cohort,
                                  catvar = catvar,
                                  contvar = contvar,
                                  expand = FALSE)
means_all <- colMeans(elig_cohort_num[, ..listvar])
sd_all <- apply(elig_cohort_num[, ..listvar], 2, sd)

# encode factors
baseline_num <- encode_factors(dt = baseline,
                               catvar = c(catvar, "trt"),
                               contvar = contvar,
                               expand = FALSE)
means <- colMeans(baseline_num[, ..listvar])
means_treated <- colMeans(baseline_num[trt==1, ..listvar])
means_untreated <- colMeans(baseline_num[trt==0, ..listvar])
wmeans <- apply(baseline_num[, ..listvar], 2, 
                weighted.mean, 
                w = as.numeric(baseline[, "sw_IPSW"][[1]]))
wmeans_treated <- apply(baseline_num[trt==1, ..listvar], 2, 
                        weighted.mean, 
                        w = as.numeric(baseline[trt==1, "sw_IPSW"][[1]]))
wmeans_untreated <- apply(baseline_num[trt==0, ..listvar], 2, 
                          weighted.mean, 
                          w = as.numeric(baseline[trt==0, "sw_IPSW"][[1]]))

# balance before weighting
TASMD <- abs(means - means_all) / sd_all
TASMD_treated <- abs(means_treated - means_all) / sd_all
TASMD_untreated <- abs(means_untreated - means_all) / sd_all

# balance after weighting
TASMD_wt <- abs(wmeans - means_all) / sd_all
TASMD_treated_wt <- abs(wmeans_treated - means_all) / sd_all
TASMD_untreated_wt <- abs(wmeans_untreated - means_all) / sd_all

# save table
TASMDs_dt <- data.frame(TASMD = TASMD, 
                       untreated_unweighted = TASMD_untreated,
                       treated_unweighted = TASMD_treated,
                       TASMD_IPSW = TASMD_wt, 
                       untreated_weighted = TASMD_untreated_wt, 
                       treated_weighted = TASMD_treated_wt)
rownames(TASMDs_dt) <- varnames
openxlsx::write.xlsx(
  TASMDs_dt,
  rowNames = TRUE,
  file = paste0(results_path, "Supplemental/TASMD.xlsx")
)

# save plot
TASMDs_dt <- TASMDs_dt[order(TASMDs_dt$TASMD_IPSW, decreasing = TRUE), ]
ggplot2::ggsave(
  plot = love_plot(
    SMDs_dt = TASMDs_dt,
    SMD_names = c("TASMD", "TASMD_IPSW"),
    plot_title = "Before and after IPSW",
    xlab_title = "Target absolute standardized mean difference",
    xmax = 1
  ),
  filename = paste0(results_path, "Supplemental/Figure_S6.png"),
  width = 15,
  height = 15,
  dpi = 300
)

# save cohort
save(
  id_name,
  listvar,
  listvar_main,
  varnames,
  catvar,
  contvar,
  non_normal_vars,
  treatment_label,
  control_label,
  baseline,
  model_PS,
  coef_PS_overall,
  elig_cohort,
  w_meths,
  file = file.path(
    "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/cohort_with_weights.Rdata"
  )
)

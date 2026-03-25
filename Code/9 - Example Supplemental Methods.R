# remove, set directory, and load data
rm(list = ls(all.names = TRUE))
setwd(
  "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Care/"
)
load("Data/cohort_with_models.Rdata")

################################################################################
### select patient
################################################################################
nr_patient <- 1
patients_risk <- baseline[, c("age",
                              "egfr2021",
                              "cancer",
                              "dm",
                              "ihd",
                              "vhd",
                              "pvd",
                              "female",
                              "albumin")]
patient_risk <- as.numeric(as.matrix(patients_risk[nr_patient, ]))

################################################################################
### fit risk model
################################################################################
elig_Surv <- survival::Surv(elig_cohort$time2event_death_2y, elig_cohort$event_death_2y)
dd <- rms::datadist(elig_cohort)
options(datadist = "dd")
risk_model <- rms::cph(
  elig_Surv ~ age + egfr2021 + cancer + dm + ihd +
    vhd + pvd + female + albumin,
  data = elig_cohort,
  method = c("breslow"),
  y = TRUE,
  x = TRUE
)

# extract coefficients
coef_risk <- coef(risk_model)

# extract means
centers_risk <- risk_model$means

# manual calculation linear predictor
manual_lp_risk <- sum(coef_risk * (patient_risk - centers_risk))

# formula
0.04740815 * (86 - 79.7733345) - 0.01213913 * (10.93256 - 14.7260768) +
  0.15658700 * (0 - 0.1282299) + 0.14536098 * (1 - 0.4978028) +
  0.14431054 * (0 - 0.4282827) + 0.29477188 * (0 - 0.1008086) +
  0.26784122 * (0 - 0.1784145) - 0.03105380 * (0 - 0.3711549) - 0.05511083 * (33 - 35.4782036)

# check using predict
predict(risk_model, newdata = baseline[nr_patient, ], type = "lp")

# baseline hazard
bh_risk <- survival::basehaz(risk_model)
h0_risk <- bh_risk$hazard[bh_risk$time == horizon]

# predicted probability
manual_prob_risk <- 1 - exp(-h0_risk * exp(manual_lp_risk))

# check 
1 - exp(-0.7692973 * exp(0.402959))

# check with function
PredictionTools::fun.event(h0 = h0_risk, lp = manual_lp_risk)

################################################################################
### ITE predicted risk model
################################################################################
baseline$lp_risk <- predict(risk_model, newdata = baseline)
patients_ITE <- baseline[, c("trt", "lp_risk")]
patient_ITE <- as.numeric(as.matrix(patients_ITE[nr_patient, ]))
patient_ITE <- c(patient_ITE, patient_ITE[1] * patient_ITE[2])

# ITE model
ITE_model <- rms::cph(
  survival::Surv(time2event_death_2y, event_death_2y) ~
    trt * lp_risk,
  data = baseline,
  weights = baseline$sw_IPTW,
  x = TRUE,
  y = TRUE
)
coef_ITE <- ITE_model$coefficients

# centers
centers_ITE <- ITE_model$means

# manual calculation linear predictor
manual_lp_ITE <- sum(coef_ITE * (patient_ITE - centers_ITE))

# check using predict
predict(ITE_model, newdata = baseline[nr_patient, ], type = "lp")

# baseline hazard
bh_ITE <- survival::basehaz(ITE_model)
h0_ITE <- bh_ITE$hazard[bh_ITE$time == horizon]

# predicted probability
manual_prob_ITE <- 1 - exp(-h0_ITE * exp(manual_lp_ITE))

# check with function
PredictionTools::fun.event(h0 = h0_ITE, lp = manual_lp_ITE)

# check with original code
data_1 <- data.frame(trt = 1, lp_risk = patient_ITE[2])
sf_1 <- survival::survfit(ITE_model, newdata = data_1)
risk_1 <- 1 - summary(sf_1, times = horizon)$surv

# counterfactual patient
patient_ITE_counterfactual <- data.frame(trt = 0,
                                         lp_risk = patient_ITE[2],
                                         `trt * lp_risk` = 0)

# counterfactual lp
manual_lp_ITE_counterfactual <- sum(coef_ITE * (patient_ITE_counterfactual - centers_ITE))

# check lp
predict(ITE_model, newdata = patient_ITE_counterfactual, type = "lp")

# counterfactual predicted probability
manual_prob_ITE_counterfactual <- 1 - exp(-h0_ITE * exp(manual_lp_ITE_counterfactual))
manual_prob_ITE_counterfactual

# check with function
PredictionTools::fun.event(h0 = h0_ITE, lp = manual_lp_ITE_counterfactual)

# check with original code
data_0 <- data.frame(trt = 0, lp_risk = patient_ITE[2])
sf_0 <- survival::survfit(ITE_model, newdata = data_0)
risk_0 <- 1 - summary(sf_0, times = horizon)$surv

# risk difference
RD <- manual_prob_ITE - manual_prob_ITE_counterfactual

# check
as.numeric(risk_1 - risk_0)

# dRMST
rmst_1 <- as.numeric(summary(sf_1, rmean = horizon)$table["rmean"]) / 30.5
rmst_0 <- as.numeric(summary(sf_0, rmean = horizon)$table["rmean"]) / 30.5

# hazard ratio
HR <- exp(manual_lp_ITE - manual_lp_ITE_counterfactual)

cat(
  "Patient is",
  baseline[nr_patient, ][["age"]],
  "years old \n",
  "eGFR                                 :",
  baseline[nr_patient, ][["egfr2021"]],
  "ml/min/m3 \n",
  "Albumin                              :",
  baseline[nr_patient, ][["albumin"]],
  "g/L \n",
  "Coefficients risk model              :",
  round(coef_risk, 4),
  "\n",
  "Centers risk model                   :",
  round(centers_risk, 4),
  "\n",
  "Baseline hazard at two years for risk:",
  round(h0_risk, 4),
  "\n",
  "Linear predictor value               :",
  round(manual_lp_risk, 4),
  "\n",
  "Probability risk                     :",
  round(manual_prob_risk, 4),
  "\n",
  "Baseline hazard at two years for ITE :",
  round(h0_ITE, 4),
  "\n",
  "Coefficients ITE model               :",
  round(coef_ITE, 4),
  "\n",
  "Centers ITE model                    :",
  round(centers_ITE, 4),
  "\n",
  "Linear predictor under dialysis      :",
  round(manual_lp_ITE, 4),
  "\n",
  "Linear predictor under CC            :",
  round(manual_lp_ITE_counterfactual, 4),
  "\n",
  "Probability risk under dialysis      :",
  round(manual_prob_ITE, 4),
  "\n",
  "Probability risk under CC            :",
  round(manual_prob_ITE_counterfactual, 4),
  "\n",
  "Risk difference                      :",
  round(RD, 4),
  "\n",
  "RMST under dialysis                  :",
  round(rmst_1, 1),
  "\n",
  "RMST under CC                        :",
  round(rmst_0, 1),
  "\n",
  "Difference in RMST                   :",
  round(rmst_1 - rmst_0, 1),
  "\n",
  "Hazard ratio                         :",
  round(HR, 3),
  "\n"
)

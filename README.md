# Effect of Choosing Dialysis versus Comprehensive Conservative Care on Survival in Older Adults with Advanced CKD: Nationwide Target Trial Emulation Study

#### Carolien C.H.M. Maas, PhD[1,2], Xuerui Zhang, MD[1], Juan-Jesus Carrero, PhD[2], Antoine Créon, MD, MSc[2], Roosa Lankinen, MD, PhD[2], Ilaria Prosepe, PhD[3], Willem Jan W. Bos, MD, PhD[1], Friedo W. Dekker, PhD[1], Marie Evans, MD, PhD[4], Edouard L. Fu, PhD[1,2]

1. Department of Clinical Epidemiology, Leiden University Medical Center, Leiden, The Netherlands
2. Department of Medical Epidemiology and Biostatistics, Karolinska Institute, Stockholm, Sweden
3. Department of Biomedical Data Sciences, Leiden University Medical Center, Leiden, The Netherlands
4. Department of Clinical Intervention and Technology, Karolinska Institute, Stockholm, Sweden

## 0 - Data preparation
This file combines the data files for CKD patients into one, ensuring the right encoding for each variable.

## 1 - Apply eligibility criteria
This file applies the eligibility criteria: 
1. Estimated glomerular filtration rate (eGFR) <20 mL/min/1.73m² 
2. Aged ≥65 with a Davies Comorbidity Score* ≥2 or aged ≥80
3. All lab measurements with maximum one-year look-back
4. No history of kidney transplantation or dialysis
5. No history of HIV or dementia
*Davies comorbidities are malignancies, ischemic heart disease, peripheral valvular disese, heart failure, diabetes mellitus, systemic collagen vascular disease, COPD, cirrhosis, psychiatric illness, and HIV.

First, define the time intervals in which patients are eligible. Second, expand the data table for each day within the time interval. Third, if the decision falls on an eligible date, include patient in cohort. 

## 2 - Covariate and outcome derivation
Obtain information on 
1. Other comorbidities from inpatient and outpatient files
2. Number of hospitalizations from inpatient file
3. Medications and iron based on ATC codes
4. ESA and iron from CKD file
5. Primary kidney disease from CKD file
6. Education on treatment choice from CKD file
7. Geographical clinic level from CKD file
8. Nursing home information from outpatient file
9. Outcomes (e.g., all-cause mortality) from death file

## 3 - Compute weights
Propensity score (PS) model
1. Fit PS model
2. Check if the coefficients of the PS model is not too extreme
4. Check if the AUC of the PS model is not too high

Generalizability model
1. Define S = 1 for eligible patients with a recorded treatment decision and S = 0 for eligible patients without a recorded treatment decision
2. Fit a model with S as the outcome and X the confounders

Compute weights
1. Inverse propensity treatment weights
2. Overlap weights
3. Selection weights
4. Selection and treatment weights
5. SMR weights for ATT
6. SMR weights for ATU

Describe weights
1. Check the minimum and maximum of the weights
2. IPTW: check if SMD below 0.1 for all confounders
3. Generalizibility: check if TASMD below 0.1, comparing those with recorded treatment decision versus all eligible patients

## 4 - Descriptives

## 5 - ATE

## 6 - SA positivity

## 7 - Subgroup analysis

## 8 - Continuous HTE

## 9 - Competing risk time-to-dialysis

# remove history
rm(list = ls(all.names = TRUE))
library(haven)

################################################################################
# Dialysis
################################################################################
SOS_DIALYSISDATA2024 <- readr::read_delim(
  "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/raw/SOS_DIALYSISDATA2024.txt",
  delim = "\t"
)
load(
  "P:/SCREAM2/SCREAM2_Research/Edouard Fu/Dialysis vs conservative care/data/cleaned/snr_hdpd.Rdata"
)
save(SOS_DIALYSISDATA2024, file = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/cleaned/snr_hdpd.Rdata")

################################################################################
# CKD
################################################################################
SOS_CKDFINAL2024 <- readr::read_delim(
  "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/raw/SOS_CKDFINAL2024.txt",
  delim = "\t"
)
load(
  "P:/SCREAM2/SCREAM2_Research/Edouard Fu/Dialysis vs conservative care/data/cleaned/snr_ckd.Rdata"
)
save(SOS_CKDFINAL2024, file = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/cleaned/snr_ckd.Rdata")

################################################################################
# RRT
################################################################################
SOS_KRTDATA2024 <- readr::read_delim("P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/raw/SOS_KRTDATA2024.txt",
                                     delim = "\t")
load(
  "P:/SCREAM2/SCREAM2_Research/Edouard Fu/Dialysis vs conservative care/data/cleaned/snr_rrt.Rdata"
)
save(SOS_KRTDATA2024, file = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/cleaned/snr_rrt.Rdata")

################################################################################
# Death
################################################################################
UT_R_DORS_123160_2023 <- readr::read_delim(
  "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/raw/UT_R_DORS_123160_2023.txt",
  delim = "\t"
)
load(
  "P:/SCREAM2/SCREAM2_Research/Edouard Fu/Dialysis vs conservative care/data/cleaned/snr_death.Rdata"
)
save(UT_R_DORS_123160_2023, file = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/cleaned/snr_death.Rdata")

################################################################################
# Medications
################################################################################
UT_R_LMED_123160_2023 <- readr::read_delim(
  "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/raw/UT_R_LMED_123160_2023.txt",
  delim = "\t"
)
load(
  "P:/SCREAM2/SCREAM2_Research/Edouard Fu/Dialysis vs conservative care/data/cleaned/snr_lmed.Rdata"
)
save(UT_R_LMED_123160_2023, file = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/cleaned/snr_lmed.Rdata")

################################################################################
# Outpatient comorbidities
################################################################################
UT_R_PAR_OV_123160_2023 <- readr::read_delim(
  "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/raw/UT_R_PAR_OV_123160_2023.txt",
  delim = "\t"
)
load(
  "P:/SCREAM2/SCREAM2_Research/Edouard Fu/Dialysis vs conservative care/data/cleaned/snr_outpatient.Rdata"
)
save(UT_R_PAR_OV_123160_2023, file = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/cleaned/snr_outpatient.Rdata")

################################################################################
# Inpatient comorbidities
################################################################################
UT_R_PAR_SV_123160_2023 <- readr::read_delim(
  "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/raw/UT_R_PAR_SV_123160_2023.txt",
  delim = "\t"
)
load(
  "P:/SCREAM2/SCREAM2_Research/Edouard Fu/Dialysis vs conservative care/data/cleaned/snr_inpatient.Rdata"
)
save(UT_R_PAR_SV_123160_2023, 
     file = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/cleaned/snr_inpatient.Rdata")

################################################################################
# Missing variables
################################################################################
# Create a named list of new datasets
new_datasets <- list(
  SOS_DIALYSISDATA2024 = colnames(SOS_DIALYSISDATA2024),
  SOS_CKDFINAL2024 = colnames(SOS_CKDFINAL2024),
  SOS_KRTDATA2024 = colnames(SOS_KRTDATA2024),
  UT_R_DORS_123160_2023 = colnames(UT_R_DORS_123160_2023),
  UT_R_LMED_123160_2023 = colnames(UT_R_LMED_123160_2023),
  UT_R_PAR_OV_123160_2023 = colnames(UT_R_PAR_OV_123160_2023),
  UT_R_PAR_SV_123160_2023 = colnames(UT_R_PAR_SV_123160_2023)
)

# Convert all column names to lowercase
new_datasets <- lapply(new_datasets, tolower)

# Combine old variable names
old_var <- sort(unique(tolower(
  c(
    colnames(snr_hdpd),
    colnames(snr_ckd),
    colnames(snr_rrt),
    colnames(snr_death),
    colnames(snr_lmed),
    colnames(snr_inpatient),
    colnames(snr_outpatient)
  )
)))

# Create a data frame mapping old_var to new datasets
mapping_df <- data.frame(
  old_var = old_var,
  new_datasets = sapply(old_var, function(var) {
    found_in <- names(Filter(function(cols)
      var %in% cols, new_datasets))
    if (length(found_in) == 0)
      return(NA)
    paste(found_in, collapse = ", ")
  }),
  stringsAsFactors = FALSE
)

# Write to Excel
# openxlsx::write.xlsx(mapping_df, file = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Data/raw/variables.xlsx", rowNames = FALSE)

# important covariate vectors for this emulation


# Table 1 covariates ------------------------------------------------------

# covariates and corresponding labels for table 1
table1_covariates <- tibble::tribble(
  ~covariate, ~label,
  "dem_age_index_cont", "Age at index date",
  "dem_sex_cont", "Sex",
  "dem_race", "Race",
  "c_smoking_history", "Smoking history",
  "c_ecog_cont", "ECOG",
  "c_year_index", "Index year"
  )

# covariates and corresponding labels for table 1
table1_labels <- list(
  "dem_age_index_cont" ~ "Age at index date",
  "dem_sex_cont" ~ "Sex",
  "dem_race" ~ "Race",
  "c_smoking_history" ~ "Smoking history",
  "c_ecog_cont" ~ "ECOG",
  "c_year_index" ~ "Index year"
  )


# global covariate vector for propensity score ----------------------------

# most comprehensive covariate vector for propensity score model
covariates_for_ps <- c(
  # all ropro covariates
  'dem_age_index_cont', 
  'dem_sex_cont', 
  'c_smoking_history', 
  'c_number_met_sites', 
  'c_hemoglobin_g_dl_cont', 
  'c_urea_nitrogen_mg_dl_cont', 
  'c_platelets_10_9_l_cont', 
  'c_calcium_mg_dl_cont', 
  'c_glucose_mg_dl_cont', 
  'c_lymphocyte_leukocyte_ratio_cont', 
  'c_alp_u_l_cont', 
  'c_protein_g_l_cont', 
  'c_alt_u_l_cont', 
  'c_albumin_g_l_cont', 
  'c_bilirubin_mg_dl_cont', 
  'c_chloride_mmol_l_cont', 
  'c_monocytes_10_9_l_cont', 
  'c_eosinophils_leukocytes_ratio_cont', 
  'c_ldh_u_l_cont', 
  'c_hr_cont', 
  'c_sbp_cont', 
  'c_oxygen_cont', 
  'c_ecog_cont', 
  'c_neutrophil_lymphocyte_ratio_cont', 
  'c_bmi_cont', 
  'c_ast_alt_ratio_cont', 
  'c_stage_initial_dx_cont',
  # additional covariates
  "dem_race",
  "dem_region",
  "dem_ses",
  "c_time_dx_to_index"
  )

# global covariate vector for imputation ----------------------------------

# most comprehensive covariate vector for imputation
# i.e., which covariates should be predictors in imputation
covariates_for_imputation <- c(
  # all ps covariates (including ropro covariates)
  'dem_age_index_cont', 
  'dem_sex_cont', 
  'c_smoking_history', 
  'c_number_met_sites', 
  'c_hemoglobin_g_dl_cont', 
  'c_urea_nitrogen_mg_dl_cont', 
  'c_platelets_10_9_l_cont', 
  'c_calcium_mg_dl_cont', 
  'c_glucose_mg_dl_cont', 
  'c_lymphocyte_leukocyte_ratio_cont', 
  'c_alp_u_l_cont', 
  'c_protein_g_l_cont', 
  'c_alt_u_l_cont', 
  'c_albumin_g_l_cont', 
  'c_bilirubin_mg_dl_cont', 
  'c_chloride_mmol_l_cont', 
  'c_monocytes_10_9_l_cont', 
  'c_eosinophils_leukocytes_ratio_cont', 
  'c_ldh_u_l_cont', 
  'c_hr_cont', 
  'c_sbp_cont', 
  'c_oxygen_cont', 
  'c_ecog_cont', 
  'c_neutrophil_lymphocyte_ratio_cont', 
  'c_bmi_cont', 
  'c_ast_alt_ratio_cont', 
  'c_stage_initial_dx_cont',
  # auxiliary covariates
  #"c_squamous",
  "c_de_novo_mets_dx",
  "c_height_cont",
  "c_weight_cont",
  "c_year_index",
  "c_dbp_cont",
  # exposure
  "treat",
  # outcome
  "death_itt",
  "fu_itt_months"
  )

# simulates a FLAURA cohort
# with covariate names and distributions
# approximated to match those of ENCORE

simulate_flaura <- function(n_total = 3500, 
                            treat_prevalence = .51, 
                            seed = 42,
                            include_id = TRUE,
                            imposeNA = TRUE 
                            ){
  
  # set the seed
  set.seed(seed)
  
  # Build the total cohort
  cohort <- tibble::tibble(
    # patient id
    patientid = paste0("Patient", seq(1, n_total, 1)),
    # treat
    treat = rbinom(n = n_total, size = 1, prob = treat_prevalence),
    # age
    dem_age_index_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 69, sd = 8),
      rnorm(n = n_total, mean = 70, sd = 9)
      ),
    # sex
    dem_sex_cont = ifelse(
      treat == 1,
      rbinom(n = n_total, size = 1, prob = .33),
      rbinom(n = n_total, size = 1, prob = .34)
      ),
    # smoking
    c_smoking_history = ifelse(
      treat == 1, 
      sample(c(TRUE, FALSE), size = n_total, replace = TRUE, prob = c(.44, .56)),
      sample(c(TRUE, FALSE), size = n_total, replace = TRUE, prob = c(.53, .47))
      ),
    # num met sites
    c_number_met_sites = sample(
      c(1, 2, 3, 4, 5), size = n_total, replace = TRUE, prob = c(.73, .22, .035, .006, .001)
      ),
    # c_ecog_cont
    c_ecog_cont = ifelse(
      treat == 1, 
      rbinom(n = n_total, size = 1, prob = .54),
      rbinom(n = n_total, size = 1, prob = .61)
      ),
    # c_stage_initial_dx_cont
    c_stage_initial_dx_cont = ifelse(
      treat == 1, 
      sample(c(4, 3, 2, 1), size = n_total, replace = TRUE, prob = c(.95, 0.04, 0.01, 0)),
      sample(c(4, 3, 2, 1), size = n_total, replace = TRUE, prob = c(.95, 0.02, 0.01, 0.1))
      ),
    # dem_race
    dem_race = factor(ifelse(
      treat == 1, 
      sample(c("White", "Asian", "Other"), size = n_total, replace = TRUE, prob = c(.57, .13, .3)),
      sample(c("White", "Asian", "Other"), size = n_total, replace = TRUE, prob = c(.60, .11, .29)))
      ),
    # dem_region
    dem_region = factor(ifelse(
      treat == 1, 
      sample(c("Midwest", "Northeast", "South", "West"), size = n_total, replace = TRUE, prob = c(.12, .31, .39, .18)),
      sample(c("Midwest", "Northeast", "South", "West"), size = n_total, replace = TRUE, prob = c(.1, .33, .37, .2)))
      ),
    # dem_ses
    dem_ses = ifelse(
      treat == 1, 
      sample(c(1, 2, 3, 4, 5), size = n_total, replace = TRUE, prob = c(.17, .14, .19, .22, .28)),
      sample(c(1, 2, 3, 4, 5), size = n_total, replace = TRUE, prob = c(.1, .14, .26, .22, .28))
      ),
    # c_hemoglobin_g_dl_cont
    c_hemoglobin_g_dl_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 12.80, sd = 1.3),
      rnorm(n = n_total, mean = 13, sd = 1.1)
      ),
    # c_urea_nitrogen_mg_dl_cont
    c_urea_nitrogen_mg_dl_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 2.77, sd = 0.3),
      rnorm(n = n_total, mean = 2.77, sd = 0.5)
      ),
    # c_platelets_10_9_l_cont
    c_platelets_10_9_l_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 266, sd = 53),
      rnorm(n = n_total, mean = 255, sd = 55)
      ),
    # c_calcium_mg_dl_cont
    c_calcium_mg_dl_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 2.24, sd = 0.04),
      rnorm(n = n_total, mean = 2.23, sd = 0.03)
      ),
    # c_glucose_mg_dl_cont
    c_glucose_mg_dl_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 4.65, sd = 0.11),
      rnorm(n = n_total, mean = 4.64, sd = 0.12)
      ),
    # c_lymphocyte_leukocyte_ratio_cont
    c_lymphocyte_leukocyte_ratio_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 2.93, sd = 0.17),
      rnorm(n = n_total, mean = 2.94, sd = 0.20)
      ),
    # c_alp_u_l_cont
    c_alp_u_l_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 4.51, sd = 0.13),
      rnorm(n = n_total, mean = 4.47, sd = 0.22)
      ),
    # c_protein_g_l_cont
    c_protein_g_l_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 69, sd = 4),
      rnorm(n = n_total, mean = 68, sd = 4)
      ),
    # c_alt_u_l_cont
    c_alt_u_l_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 2.89, sd = 0.33),
      rnorm(n = n_total, mean = 2.94, sd = 0.30)
      ),
    # c_albumin_g_l_cont
    c_albumin_g_l_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 40, sd = 3),
      rnorm(n = n_total, mean = 39, sd = 3)
      ),
    # c_bilirubin_mg_dl_cont
    c_bilirubin_mg_dl_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = -0.92, sd = 1.28),
      rnorm(n = n_total, mean = -0.69, sd = 1.21)
      ),
    # c_chloride_mmol_l_cont
    c_chloride_mmol_l_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 102.0, sd = 3),
      rnorm(n = n_total, mean = 102.0, sd = 3)
      ),
    # c_monocytes_10_9_l_cont
    c_monocytes_10_9_l_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = -0.51, sd = 0.25),
      rnorm(n = n_total, mean = -0.51, sd = 0.41)
      ),
    # c_eosinophils_leukocytes_ratio_cont
    c_eosinophils_leukocytes_ratio_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 0.69, sd = 0.59),
      rnorm(n = n_total, mean = 0.72, sd = 0.34)
      ),
    # c_ldh_u_l_cont
    c_ldh_u_l_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 1.69, sd = 0.05),
      rnorm(n = n_total, mean = 1.68, sd = 0.04)
      ),
    # c_hr_cont
    c_hr_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 4.43, sd = 0.04),
      rnorm(n = n_total, mean = 4.41, sd = 0.03)
      ),
    # c_sbp_cont
    c_sbp_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 4.85, sd = 0.09),
      rnorm(n = n_total, mean = 4.85, sd = 0.11)
      ),
    # c_oxygen_cont
    c_oxygen_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 97.00, sd = 0.01),
      rnorm(n = n_total, mean = 97.00, sd = 0.02)
      ),
    # c_neutrophil_lymphocyte_ratio_cont
    c_neutrophil_lymphocyte_ratio_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 1.29, sd = 0.38),
      rnorm(n = n_total, mean = 1.32, sd = 0.34)
      ),
    # c_bmi_cont
    c_bmi_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 3.23, sd = 0.13),
      rnorm(n = n_total, mean = 3.23, sd = 0.13)
      ),
    # c_ast_alt_ratio_cont
    c_ast_alt_ratio_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 0.12, sd = 0.29),
      rnorm(n = n_total, mean = 0.11, sd = 0.30)
      ),
    # c_time_dx_to_index
    c_time_dx_to_index = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 44, sd = 17),
      rnorm(n = n_total, mean = 74, sd = 43)
      ),
    # c_de_novo_mets_dx
    c_de_novo_mets_dx = ifelse(
      treat == 1, 
      sample(c(TRUE, FALSE), size = n_total, replace = TRUE, prob = c(.79, .21)),
      sample(c(TRUE, FALSE), size = n_total, replace = TRUE, prob = c(.69, .31))
      ),
    # c_height_cont
    c_height_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 1.65, sd = 0.08),
      rnorm(n = n_total, mean = 1.64, sd = 0.07)
      ),
    # c_weight_cont
    c_weight_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 69, sd = 11),
      rnorm(n = n_total, mean = 68, sd = 12)
      ),
    # c_dbp_cont
    c_dbp_cont = ifelse(
      treat == 1, 
      rnorm(n = n_total, mean = 76, sd = 6),
      rnorm(n = n_total, mean = 75, sd = 7)
      ),
    # c_year_index
    c_year_index = factor(ifelse(
      treat == 1, 
      sample(c("<2018", "2018+"), size = n_total, replace = TRUE, prob = c(.95, .05)),
      sample(c("<2018", "2018+"), size = n_total, replace = TRUE, prob = c(.05, .95)))
      )
    )
  
  # assign betas for hazard model
  betas_os <- c(
    treat = log(0.8)
    )
  
  set.seed(seed)
  cohort_outcome <- cohort |> 
    bind_cols(
      simsurv::simsurv(
        dist = "weibull",
        gammas = 1.2,
        lambdas = 0.015,
        betas = betas_os,
        x = cohort,
        maxt = 120
        )
      ) |> 
    dplyr::select(-id) |> 
    dplyr::rename(
      fu_itt_months = eventtime,
      death_itt = status
      )
  
  # summary(cohort_outcome$death_itt)
  # summary(cohort_outcome$fu_itt_months)
  # hist(cohort_outcome$fu_itt_months)
        
  if(!include_id) cohort_outcome <- cohort_outcome |> dplyr::select(-patientid)
  
  if(imposeNA){
    
    varComplete <-  c("patientid", "treat", "c_year_index", "fu_itt_months", "death_itt", "dem_age_index_cont", "dem_sex_cont")
    varNA <- setdiff(names(cohort_outcome), varComplete)
    
    # define column numbers for pattern and weight determination
    colsNA <- which(colnames(cohort_outcome) %in% varNA)
    
    # define missingness pattern
    default_pattern <- rep(1, ncol(cohort_outcome))
    pattern <- replace(default_pattern, colsNA, 0)
    
    set.seed(seed)
    cohort_outcome <- mice::ampute(
      data = cohort_outcome, 
      prop = .33,
      patterns = pattern,
      mech = "MCAR"
      )$amp
    
  }
  
  return(cohort_outcome)
  
}

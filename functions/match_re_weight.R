#' Custom function to match and re-weight a patient population
#'
#' @param i integer, the m_th imputed dataset (can be be used in conjunction with lapply)
#' @param mids_object mids object with m multiple imputed datasets
#' @param ps_formula formula of propensity score model fit
#' 
#' @importFrom mice complete
#' @importFrom MatchIt matchit match.data
#' @importFrom WeightIt weightit
#' 
#' @return matched and re-weighted dataset
#' 
#' @export
#'
#' @examples
#' data_miss <- simulate_flaura(
#'   n_total = 3500, 
#'   treat_prevalence = .51, 
#'   seed = 41, 
#'   include_id = FALSE, 
#'   imposeNA = TRUE
#'   )
#'   
#'  covariates <- data_miss |> 
#'    select(starts_with("c_"), starts_with("dem_")) |> 
#'    colnames()
#' 

match_re_weight <- function(i, 
                            mids_object = NULL, 
                            ps_formula = NULL,
                            ...
                            ){
  
  # perform matching on the mth imputed dataset -----------------------------
  
  # extract the mth imputed dataset
  imputed_data_m <- mice::complete(data = mids_object, action = i)
  
  # create a matchit object
  matchit_out <- MatchIt::matchit(
    formula = ps_formula,
    data = imputed_data_m,
    ...
    )
  
  # read data ---------------------------------------------------------------
  
  # get ALL patients with matching weights
  data_matchit_all <- MatchIt::match.data(
    object = matchit_out, 
    drop.unmatched = FALSE
    )
  
  # create a data frame with temporary patient/case ids
  data_matchit_all <- data_matchit_all |> 
    # create a case id
    mutate(caseid = seq(1, nrow(data_matchit_all), 1)) |> 
    # convert sex and ecog into a factor variable since anesrake doesn't accept 0/1 numeric encoding
    mutate(across(c(dem_sex_cont, c_ecog_cont), function(x) factor(as.character(x))))
  
  # create a sub-cohort with matching weights == 1 (= matched patients)
  # we will compute the new sampling/re-weighting weights only on matched patients
  data_matched <- data_matchit_all |> 
    dplyr::filter(weights == 1)
  
  # Define FLAURA distributions for key covariates --------------------------
  # order is as in Table 1

  ## sex ---------------------------------------------------------------------
  
  # male sex by exposure
  avg_prop_male <- sum(.36, .38)/2
  
  # female (0) to male (1) proportion:
  sex_target <- c(1-avg_prop_male, avg_prop_male) 
  names(sex_target) <- c("0", "1")

  ## race --------------------------------------------------------------------
  # white, asian, other
  race_target <- c(.36, .62, .02) # "other" is rounded up to 2% to sum up to 100%
  names(race_target) <- c("White", "Asian", "Other")

  ## smoking -----------------------------------------------------------------
  
  # never smoking by exposure
  avg_prop_never_smoker <- sum(.65, .63)/2
  
  # current/former smoker (TRUE) to never smoker (FALSE) proportion
  # note: logical variables in dataframe can be matched to a numeric vector of length 2 and ordered with the TRUE target as the first element and the FALSE target as the second element.
  smoker_target <- c(1-avg_prop_never_smoker, avg_prop_never_smoker)

  ## ecog --------------------------------------------------------------------
  
  # ecog 0 by exposure 
  avg_prop_ecog0 <- sum(.4, .42)/2
  
  # ecog 0 to ecog 1 proportion
  ecog_target <- c(avg_prop_ecog0, 1-avg_prop_ecog0)
  names(ecog_target) <- c("0", "1")
  

  # summarize target distributions in a named list vector --------------
  targets <- list(sex_target, race_target, smoker_target, ecog_target)
  names(targets) <- c("dem_sex_cont", "dem_race", "c_smoking_history", "c_ecog_cont")
  

  # apply anesrake function to re-weight -----------------------------------
  anesrake_out <- anesrake::anesrake(
    inputter = targets, 
    dataframe = data_matched, 
    caseid = data_matched$caseid,
    weightvec = data_matched$weights,
    verbose = FALSE
    )
  
  # create a temporary dataframe with id and sampling/re-weighting weights
  caseweights <- data.frame(caseid = anesrake_out$caseid, re_weights = anesrake_out$weightvec)
  
  # join back to matchit dataframe with ALL subjects
  # this is necessary to add the sampling weights to the final matchit object
  data_matchit_all <- data_matchit_all |> 
    dplyr::left_join(caseweights, by = "caseid") |> 
    # unmatched patients get a sampling/re-weighting weight of 0
    dplyr::mutate(re_weights = ifelse(is.na(re_weights), 0, re_weights)) |>
    # make sure we have the same original order of subject
    dplyr::arrange(caseid) |> 
    # now we can drop the temporary id
    dplyr::select(-caseid)
  
  matchit_out <- MatchIt::add_s.weights(
    m = matchit_out,
    s.weights = data_matchit_all$re_weights
    )
  
  return(matchit_out)
  
}


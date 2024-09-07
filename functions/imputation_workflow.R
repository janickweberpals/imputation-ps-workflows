#' Streamlines the imputation workflow
#'
#' @description
#' This function streamlines the imputation workflow from an eligible cohort
#' with missings to an imputed mids object with only a subset of the key covariates
#' and computed ropro
#'
#' @details
#' ...
#'
#' @param ard_eligible data frame with all eligible patients and covariates needed for imputation
#' @param database character, which database is used ("edb1", "edb2", "edb3" or "edb4")
#' @param cancer character, which cancer is investigated ("aNSCLC", "MetastaticBreast", "EarlyBreast", "MetastaticCRC", "MultipleMyeloma")
#' @param covars_for_imputation character vector with covariates for imputation (should include covariates for propensity score, treatment, outcome and optionally other auxiliary covariates)
#'
#' @return imputed mids object with ROPRO
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr select any_of all_of summarize pull
#' @importFrom glue glue
#' @importFrom mice futuremice quickpred complete as.mids 
#'
#' @export
#'
#' @examples
#'\dontrun{
#' if(require("smdi")){
#' 
#'   library(encore.io)
#'
#'   mids_data <- imputation_workflow(
#'     ard_eligible, 
#'     database = "edb1", 
#'     cancer = "aNSCLC", 
#'     covars_for_imputation = c("dem_age_index", "dem_sex")
#'     )
#'
#'  }
#'}
#'

imputation_workflow <- function(ard_eligible = NULL, 
                                database = NULL, 
                                cancer = NULL,
                                covars_for_imputation = NULL
                                ){

  # check if packages available
  install_on_demand("smdi")
  
  # input checks
  assertthat::assert_that(inherits(ard_eligible, c("data.frame", "tibble")), msg = "<ard_eligible> needs to be a data.frame or tibble")
  assertthat::assert_that(database %in% c("edb1", "edb2", "edb3", "edb4"), msg = "<database> needs to be correctly specified")
  assertthat::assert_that(cancer %in% c("aNSCLC", "MetastaticBreast", "EarlyBreast", "MetastaticCRC","MultipleMyeloma"), msg = "<cancer> needs to be correctly specified")
  assertthat::assert_that(is.vector(covars_for_imputation), msg = "<covars_for_imputation> needs to be a character vector")
  
  # vector of possible ids
  possible_ids <- c("patientid", "patient_id", "cpid", "chai_patient_id")
    
  # select key covariates ------------------------------------------------
  
  # select columns
  ard_eligible_select <- ard_eligible |> 
    # select only required variables
    dplyr::select(
      # patient id, so we can join variables later on
      dplyr::any_of(possible_ids),
      dplyr::all_of(covars_for_imputation)
      ) |>
    # drop unused levels of factor covariates
    droplevels()  

  # determine average proportion of missingness -------------------------------------
  prop_miss <- ard_eligible_select |> 
    smdi::smdi_summarize() |> 
    dplyr::summarize(mean(prop_miss)) |> 
    dplyr::pull() |> 
    ceiling()
  
  cat(glue::glue("The average proportion missing among all selected covariates is {prop_miss}% \n"))
  cat("\n\n")
  cat(glue::glue("=> Imputing m = {prop_miss} datasets \n"))
  
  
  # predictor matrix --------------------------------------------------------
  
  # create predictor matrix to not mistakenly use a predictor we don't want to use
  
  cat(glue::glue("Covariates included for imputation and as predictors: {paste(covars_for_imputation, collapse = ', ')} \n"))
  cat("\n\n")
  cat(glue::glue("Covariates NOT included as predictors: {paste(setdiff(names(ard_eligible_select), covars_for_imputation), collapse = ', ')} \n"))
  
  predictors_mice <- mice::quickpred(
    data = ard_eligible_select,
    include = covars_for_imputation,
    exclude = setdiff(names(ard_eligible_select), covars_for_imputation),
    mincor = 0,
    minpuc = 0
    )
  
  # actual imputation -------------------------------------------------------
  mids_data <- mice::futuremice(
    parallelseed = 42,
    n.core = prop_miss,
    data = ard_eligible_select,
    method = "rf",
    m = prop_miss,
    pred = predictors_mice,
    print = FALSE
    )

  # compute ropro -----------------------------------------------------------
  
  # convert imputed mids object into a long data object to compute ropro
  # here, action = "long" and include = TRUE are needed to be able to convert back
  # to a mids object
  data_long <- mice::complete(data = mids_data, action = "long", include = TRUE)
  
  # compute ropro and use formula according to what 
  # covariates are available
  if(database %in% c("edb1", "edb3", "edb4")){
    
    data_long <- data_long |> 
      encore.io::edb1_3_4_compute_ropro(cancer = cancer)
    
  }else{
    
    cancer_edb2 <- if(cancer %in% c("aNSCLC", "NSCLC")) "NSCLC" else "MM"
    
    data_long <- data_long |> 
      encore.io::edb2_compute_ropro(cancer = cancer_edb2)
    
  }

  # convert back to mids object
  mids_data <- mice::as.mids(long = data_long)
  
  # output any potential log
  if(!is.null(mids_data$loggedEvents)) print(mids_data$loggedEvents)
  
  # return mids dataset
  return(mids_data)

}

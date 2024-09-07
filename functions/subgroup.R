#' Custom function to perform a subgroup analysis within a
#' multiple imputation and propensity score workflow
#'
#' @details This function is ...
#' 
#' The function performs ...
#' 
#' @param x data.frame or mild object/list of data.frames if used in combination with lapply
#' @param subgroup character, name of the column with subgroups
#' @param matching_weighting character, one of "matching" or "weighting"
#' @param ... other arguments and specifications to propagate on to 
#' \code{\link[MatchIt]{matchit}} or \code{\link[WeightIt]{weighit}}, 
#' depending on chosen method (\code{matching_weighting}).
#' 
#' @importFrom mice complete
#' @importFrom MatchIt matchit match.data
#' @importFrom WeightIt weightit as.weightit
#' 
#' @return ...
#' 
#' 
#' @seealso 
#' \code{\link[MatchIt]{matchit}}
#' \code{\link[WeightIt]{weighit}}
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
#'  # impute data
#'  data_imp <- futuremice(
#'    parallelseed = 42,
#'    n.core = 7,
#'    data = data_miss,
#'    method = "rf",
#'    m = 10,
#'    print = FALSE
#'    )
#'    
#'  # define covariates for propensity score model
#'  covariates <- data_miss |> 
#'    select(starts_with("c_"), starts_with("dem_")) |> 
#'    colnames()
#'  
#'  # define propensity score model
#'  ps_form <- as.formula(paste("treat ~", paste(covariates, collapse = " + ")))
#'  
#'  # create a mild object containing lists of data.frames
#'  data_mild <- mice::complete(data = data_imp, action = "all", include = FALSE)
#'  
#' 

re_weight <- function(x,
                      subgroup = NULL,
                      matching_weighting = NULL,
                      ...
                      ){
  
  # checks
  assertthat::assert_that(inherits(x, c("data.frame", "list")), msg = "<x> needs to be a data frame or a list of data frames")
  assertthat::assert_that(inherits(subgroup, "character"), msg = "<subgroup> needs to be a character describing a column name in <x>")
  assertthat::assert_that(subgroup %in% names(x), msg = "<subgroup> needs to be a column in <x>")
  assertthat::assert_that(matching_weighting %in% c("matching", "weighting"), msg = "<matching_weighting> needs to be either matching or weighting")
  
  if(is.null(targets)) message("No target distributions specified, no re-weighting will be performed.")
  
  # perform matching on the imputed dataset -----------------------------
  if(matching_weighting == "matching"){
    
    # create a matchit object on x;
    # if no re-weighting is desired,
    # this is already it
    object_out <- MatchIt::matchit(
      data = x,
      ...
      )
    
    # if re-weighting is desired:
    # get ALL patients with corresponding
    # matching weights
    data_all <- MatchIt::match.data(
      object = object_out, 
      drop.unmatched = FALSE
      )
    
    # create temporary patient/case ids for raking procedure
    data_all <- data_all |> 
      mutate(caseid = seq(1, nrow(data_all), 1)) |> 
      # indicator for patients who re-weighting should be applied to
      # = matched patients (= 1; unmatched = 0)
      mutate(reweight = weights)
    
  }else if(matching_weighting == "weighting"){
    
    # create a weightit object on x;
    # if no additional re-weighting is desired,
    # this is already it
    object_out <- WeightIt::weightit(
      data = x,
      ...
      )
    
    # get ALL patients with corresonding weights
    data_all <- x |> 
      dplyr::mutate(weights = object_out$weights) |> 
      mutate(caseid = seq(1, nrow(x), 1)) |> 
      # indicator for patients who re-weighting should be applied to
      # = patients with weight > 0
      mutate(reweight = ifelse(weights == 0, 0, 1))
    
  }

  

  # if re-weighting by the anesrake procedure is desired: -------------------
  if(!is.null(targets)){
  
    # apply anesrake function to re-weight -----------------------------------
    anesrake_out <- anesrake::anesrake(
      # target distributions for raking procedure
      inputter = targets, 
      # data to be weighted
      dataframe = data_all, 
      # patient identifier
      caseid = data_all$caseid,
      # other weights that should be accounted for before re-weighting is conducted
      # for matching this is 0 (unmatched) and 1 (matched) and for weighting this is the actual weights
      weightvec = data_all$weights,
      # we only want to re-weight patients who were matched/have weights > 0
      filter = data_all$reweight,
      # no verbosity
      verbose = FALSE
      )
    
    # create a temporary dataframe with id and sampling/re-weighting weights
    caseweights <- data.frame(caseid = anesrake_out$caseid, re_weights = anesrake_out$weightvec)
    
    # join back to dataframe with ALL subjects
    # this is necessary to add the sampling weights to the final 
    # matchit/weightit object in the next step where all patients are required
    # unmatched patients get a sampling weight of 0
    data_all <- data_all |> 
      dplyr::left_join(caseweights, by = "caseid") |> 
      # unmatched patients or patients with weight = 0 get a re-weighting weight of 0
      dplyr::mutate(re_weights = ifelse(is.na(re_weights), 0, re_weights)) |>
      # make sure we have the same original order of subjects
      dplyr::arrange(caseid) |> 
      # now we can drop the temporary id
      dplyr::select(-caseid)
    
    if(matching_weighting == "matching"){
      
      # add re-weights to the matchit object so that 
      # they are incorporated into balance assessment 
      # and creation of the weights
      object_out <- MatchIt::add_s.weights(
        m = object_out,
        s.weights = data_all$re_weights
        )
      
    }else if(matching_weighting == "weighting"){
      
      object_out <- WeightIt::as.weightit(
        x = object_out$weights,
        covs = data_all,
        treat = data_all$treat,
        s.weights = data_all$re_weights
        )
      
      warning("<covariates> do not represent those that were used for weighting but all available covariates.")
      
      }
    
    }
  
  return(object_out)
  
}


---
title: "Subgroup analysis"
subtitle: "Workflow to establish a subgroup analysis after multiple imputation and propensity Score matching"
author: Janick Weberpals, RPh, PhD
date: last-modified
format: html
code-fold: false
toc: true
toc-depth: 3
code-tools: true
keep-md: true
editor: visual
embed-resources: true
bibliography: references.bib
---



This is a reproducible example on how to establish a workflow to run a subgroup analysis after multiple imputation and propensity Score matching.

Load packages:


::: {.cell}

```{.r .cell-code}
library(dplyr)
library(survival)
library(mice)
library(MatchThem)
library(survey)
library(here)
library(gtsummary)
library(parallelly)
library(ranger)
library(furrr)

source(here("functions", "source_encore.io_functions.R"))

# track time
runtime <- tictoc::tic()
```
:::


## About

This script is adapted from Noah Greifer's highly recommended blog post on "[Subgroup Analysis After Propensity Score Matching Using R](https://ngreifer.github.io/blog/subgroup-analysis-psm/)".

For a more formal manuscript on subgroup analysis with propensity scores, see Green and Stuart.[@Green2014]

## Data generation

We again use the `simulate_flaura()` function to simulate a realistic oncology comparative effectiveness cohort analytic dataset.


::: {.cell}

```{.r .cell-code}
# load example dataset with missing observations
data_miss <- simulate_flaura(
  n_total = 3500, 
  seed = 42, 
  include_id = FALSE, 
  imposeNA = TRUE,
  propNA = .13
  ) |> 
  # we simplify and assume a binary Asian/non-Asian race covariate where Asian is the reference group
  mutate(dem_race = ifelse(dem_race == "Asian", 0, 1))

covariates_for_ps <- setdiff(covariates_for_ps, "dem_race")
```
:::


## Moderator covariate

In this example, we assume heterogeneous treatment effect by sex and we aim to assess the average treatment effect among the treated for female and male patients separately. The effect size is time to all-cause mortality. In this dataset, sex is encoded with a binary covariate with 0 = female and 1 = male.


::: {.cell}

```{.r .cell-code}
table(data_miss$dem_race, useNA = "ifany")
```

::: {.cell-output .cell-output-stdout}

```

   0    1 <NA> 
1854 1166  480 
```


:::
:::


## Multiple imputation

Both the imputation and propensity score step Multiple imputation using `mice:`


::: {.cell}

```{.r .cell-code}
# impute data
data_imp <- futuremice(
  parallelseed = 42,
  n.core = parallel::detectCores()-1,
  data = data_miss,
  method = "rf",
  m = 10,
  print = FALSE
  )
```
:::


## Propensity score matching and weighting

Apply propensity score matching and weighting with replacement within in each imputed dataset.


::: {.cell}

```{.r .cell-code}
# apply propensity score matching on mids object
ps_form <- as.formula(paste("treat ~", paste(covariates_for_ps, collapse = " + ")))
ps_form
```

::: {.cell-output .cell-output-stdout}

```
treat ~ dem_age_index_cont + dem_sex_cont + c_smoking_history + 
    c_number_met_sites + c_hemoglobin_g_dl_cont + c_urea_nitrogen_mg_dl_cont + 
    c_platelets_10_9_l_cont + c_calcium_mg_dl_cont + c_glucose_mg_dl_cont + 
    c_lymphocyte_leukocyte_ratio_cont + c_alp_u_l_cont + c_protein_g_l_cont + 
    c_alt_u_l_cont + c_albumin_g_l_cont + c_bilirubin_mg_dl_cont + 
    c_chloride_mmol_l_cont + c_monocytes_10_9_l_cont + c_eosinophils_leukocytes_ratio_cont + 
    c_ldh_u_l_cont + c_hr_cont + c_sbp_cont + c_oxygen_cont + 
    c_ecog_cont + c_neutrophil_lymphocyte_ratio_cont + c_bmi_cont + 
    c_ast_alt_ratio_cont + c_stage_initial_dx_cont + dem_region + 
    dem_ses + c_time_dx_to_index
```


:::
:::


### Matching


::: {.cell}

```{.r .cell-code}
match_within_strata <- function(i, 
                                imputed_data = NULL, 
                                ps_formula = NULL, 
                                filter_expr = NULL
                                ){
  
  stratum <- mice::complete(imputed_data, i) |> 
    dplyr::filter(eval(filter_expr))
    
  matched <- MatchIt::matchit(
    formula = ps_formula, 
    data = stratum,
    method = "nearest",
    caliper = 0.01,
    ratio = 1,
    replace = F
    ) |> 
    MatchIt::match.data()
  
  return(matched)
  
}

asian_matched <- lapply(
  X = 1:data_imp$m, 
  FUN = match_within_strata, 
  imputed_data = data_imp,
  ps_formula = ps_form,
  filter_expr = expr(dem_race == 0)
  )

non_asian_matched <- lapply(
  X = 1:data_imp$m, 
  FUN = match_within_strata, 
  imputed_data = data_imp,
  ps_formula = ps_form,
  filter_expr = expr(dem_race == 1)
  )

# combine the mth imputed and matched datasets
combine_list <- function(i, data_0 = NULL, data_1 = NULL){
  
  data_combined <- rbind(data_0[[i]], data_1[[i]])
  
  return(data_combined)
  
}

matched_all <- lapply(
  X = 1:data_imp$m, 
  FUN = combine_list, 
  data_0 = non_asian_matched,
  data_1 = asian_matched
  )
```
:::


### Weighting


::: {.cell}

```{.r .cell-code}
weight_within_strata <- function(i, 
                                 imputed_data = NULL, 
                                 ps_formula = NULL,
                                 filter_expr = NULL
                                 ){
  
  stratum <- mice::complete(imputed_data, i) |> 
    dplyr::filter(eval(filter_expr))
  
  weighted <- WeightIt::weightit(
    formula = ps_formula, 
    data = stratum,
    method = "glm",
    estimand = "ATT"
    )
  
  # trim extreme weights
  weighted <- trim(
    x = weighted, 
    at = .95, 
    lower = TRUE
    )
  
  weighted_data <- mice::complete(imputed_data, i) |> 
    dplyr::filter(eval(filter_expr)) |> 
    mutate(weights = weighted$weights)
  
  return(weighted_data)
  
}

asian_weighted <- lapply(
  X = 1:data_imp$m, 
  FUN = weight_within_strata, 
  imputed_data = data_imp,
  ps_formula = ps_form,
  filter_expr = expr(dem_race == 0)
  )

non_asian_weighted <- lapply(
  X = 1:data_imp$m, 
  FUN = weight_within_strata, 
  imputed_data = data_imp,
  ps_formula = ps_form,
  filter_expr = expr(dem_race == 1)
  )

weighted_all <- lapply(
  X = 1:data_imp$m, 
  FUN = combine_list, 
  data_0 = asian_weighted,
  data_1 = non_asian_weighted
  )
```
:::


## Outcome model comparisons

### Matching


::: {.cell}

```{.r .cell-code}
cox_fit_matching <- function(i){
  
  survival_fit <- survival::coxph(
    data = i,
    formula = Surv(fu_itt_months, death_itt) ~ treat*dem_race, 
    weights = weights, 
    cluster = subclass,
    robust = TRUE
    )
  
}
```
:::

::: {.cell}

```{.r .cell-code}
matched_all |> 
  lapply(FUN = cox_fit_matching) |> 
  mice::pool() |> 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) |> 
  dplyr::select(term, estimate, std.error, conf.low, conf.high)
```

::: {.cell-output .cell-output-stdout}

```
            term  estimate  std.error  conf.low conf.high
1          treat 1.0164611 0.06285385 0.8964279 1.1525669
2       dem_race 1.0168610 0.06594461 0.8929551 1.1579600
3 treat:dem_race 0.5838534 0.10944305 0.4687215 0.7272652
```


:::
:::


### Weighting


::: {.cell}

```{.r .cell-code}
cox_fit_weighting <- function(i){
  
  survival_fit <- survival::coxph(
    data = i,
    formula = Surv(fu_itt_months, death_itt) ~ treat*dem_race, 
    weights = weights, 
    robust = TRUE
    )
  
}
```
:::

::: {.cell}

```{.r .cell-code}
weighted_all |> 
  lapply(FUN = cox_fit_weighting) |> 
  mice::pool() |> 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) |> 
  dplyr::select(term, estimate, std.error, conf.low, conf.high)
```

::: {.cell-output .cell-output-stdout}

```
            term  estimate  std.error  conf.low conf.high
1          treat 1.0129437 0.04540748 0.9265436  1.107401
2       dem_race 1.0140269 0.05675557 0.9070347  1.133640
3 treat:dem_race 0.5930234 0.07725232 0.5094262  0.690339
```


:::
:::


## References

::: {#refs}
:::

## Session info





Script runtime: 0.29 minutes.

::: panel-tabset
### Loaded packages


::: {.cell}

```{.r .cell-code}
pander::pander(subset(data.frame(sessioninfo::package_info()), attached==TRUE, c(package, loadedversion)))
```

::: {.cell-output-display}

---------------------------------------------
     &nbsp;        package     loadedversion 
---------------- ------------ ---------------
   **dplyr**        dplyr          1.1.4     

   **furrr**        furrr          0.3.1     

   **future**       future        1.34.0     

 **gtsummary**    gtsummary        2.0.1     

    **here**         here          1.0.1     

 **MatchThem**    MatchThem        1.2.1     

   **Matrix**       Matrix         1.7-0     

    **mice**         mice         3.16.0     

 **parallelly**   parallelly      1.38.0     

   **ranger**       ranger        0.16.0     

   **survey**       survey         4.4-2     

  **survival**     survival        3.5-8     
---------------------------------------------


:::
:::


### Session info


::: {.cell}

```{.r .cell-code}
pander::pander(sessionInfo())
```

::: {.cell-output-display}
**R version 4.4.0 (2024-04-24)**

**Platform:** aarch64-apple-darwin20 

**locale:**
en_US.UTF-8||en_US.UTF-8||en_US.UTF-8||C||en_US.UTF-8||en_US.UTF-8

**attached base packages:** 
_grid_, _stats_, _graphics_, _grDevices_, _datasets_, _utils_, _methods_ and _base_

**other attached packages:** 
_furrr(v.0.3.1)_, _future(v.1.34.0)_, _ranger(v.0.16.0)_, _parallelly(v.1.38.0)_, _gtsummary(v.2.0.1)_, _here(v.1.0.1)_, _survey(v.4.4-2)_, _Matrix(v.1.7-0)_, _MatchThem(v.1.2.1)_, _mice(v.3.16.0)_, _survival(v.3.5-8)_ and _dplyr(v.1.1.4)_

**loaded via a namespace (and not attached):** 
_gtable(v.0.3.5)_, _shape(v.1.4.6.1)_, _xfun(v.0.47)_, _ggplot2(v.3.5.1)_, _htmlwidgets(v.1.6.4)_, _MatchIt(v.4.5.5)_, _lattice(v.0.22-6)_, _simsurv(v.1.0.0)_, _vctrs(v.0.6.5)_, _tools(v.4.4.0)_, _generics(v.0.1.3)_, _parallel(v.4.4.0)_, _tibble(v.3.2.1)_, _fansi(v.1.0.6)_, _pan(v.1.9)_, _pkgconfig(v.2.0.3)_, _jomo(v.2.7-6)_, _assertthat(v.0.2.1)_, _lifecycle(v.1.0.4)_, _stringr(v.1.5.1)_, _compiler(v.4.4.0)_, _tictoc(v.1.2.1)_, _munsell(v.0.5.1)_, _mitools(v.2.4)_, _codetools(v.0.2-20)_, _htmltools(v.0.5.8.1)_, _yaml(v.2.3.10)_, _glmnet(v.4.1-8)_, _pillar(v.1.9.0)_, _nloptr(v.2.1.1)_, _crayon(v.1.5.3)_, _tidyr(v.1.3.1)_, _MASS(v.7.3-60.2)_, _sessioninfo(v.1.2.2)_, _iterators(v.1.0.14)_, _rpart(v.4.1.23)_, _boot(v.1.3-30)_, _foreach(v.1.5.2)_, _mitml(v.0.4-5)_, _nlme(v.3.1-164)_, _locfit(v.1.5-9.10)_, _WeightIt(v.1.3.0)_, _tidyselect(v.1.2.1)_, _digest(v.0.6.37)_, _stringi(v.1.8.4)_, _pander(v.0.6.5)_, _listenv(v.0.9.1)_, _purrr(v.1.0.2)_, _splines(v.4.4.0)_, _rprojroot(v.2.0.4)_, _fastmap(v.1.2.0)_, _colorspace(v.2.1-1)_, _cli(v.3.6.3)_, _magrittr(v.2.0.3)_, _utf8(v.1.2.4)_, _broom(v.1.0.6)_, _withr(v.3.0.1)_, _scales(v.1.3.0)_, _backports(v.1.5.0)_, _rmarkdown(v.2.28)_, _globals(v.0.16.3)_, _nnet(v.7.3-19)_, _lme4(v.1.1-35.5)_, _chk(v.0.9.2)_, _evaluate(v.0.24.0)_, _knitr(v.1.48)_, _rlang(v.1.1.4)_, _Rcpp(v.1.0.13)_, _glue(v.1.7.0)_, _DBI(v.1.2.3)_, _renv(v.1.0.7)_, _rstudioapi(v.0.16.0)_, _minqa(v.1.2.8)_, _jsonlite(v.1.8.8)_ and _R6(v.2.5.1)_
:::
:::


### Repositories


::: {.cell}

```{.r .cell-code}
pander::pander(options('repos'))
```

::: {.cell-output-display}


  * **repos**:

    ---------------------------------------------
                      REPO_NAME
    ---------------------------------------------
     https://packagemanager.posit.co/cran/latest
    ---------------------------------------------


<!-- end of list -->



:::
:::

:::

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
  treat_prevalence = .51, 
  seed = 42, 
  include_id = FALSE, 
  imposeNA = TRUE
  )

covariates <- data_miss |> 
  select(starts_with("c_"), starts_with("dem_")) |> 
  select(-dem_sex_cont) |> 
  colnames()

head(data_miss)
```

::: {.cell-output .cell-output-stdout}
```
  treat dem_age_index_cont dem_sex_cont c_smoking_history c_number_met_sites
1     0           45.72712            1              TRUE                  1
2     0           68.80948            0              TRUE                  1
3     1           69.08840            1              TRUE                  1
4     0           68.53324            0                NA                 NA
5     0           59.37128            1              TRUE                  2
6     0           75.85751            0              TRUE                  1
  c_ecog_cont c_stage_initial_dx_cont dem_race dem_region dem_ses
1           1                       4     TRUE      South       3
2           0                       4    FALSE  Northeast       3
3           0                       4    FALSE       West       2
4          NA                      NA       NA       <NA>      NA
5           0                       4    FALSE       West       3
6           0                       4     TRUE    Midwest       4
  c_hemoglobin_g_dl_cont c_urea_nitrogen_mg_dl_cont c_platelets_10_9_l_cont
1               12.66657                   2.522498                320.6275
2               13.19548                   2.998084                231.0663
3               13.94372                   2.878485                244.8208
4                     NA                         NA                      NA
5               13.20148                   2.068868                213.3328
6               10.69378                   2.826482                316.9183
  c_calcium_mg_dl_cont c_glucose_mg_dl_cont c_lymphocyte_leukocyte_ratio_cont
1             2.196983             4.621117                          3.024268
2             2.267144             4.654893                          2.702998
3             2.278628             4.644069                          3.032139
4                   NA                   NA                                NA
5             2.239781             4.645509                          2.950502
6             2.223204             4.673633                          2.627947
  c_alp_u_l_cont c_protein_g_l_cont c_alt_u_l_cont c_albumin_g_l_cont
1       4.316845           67.60541       2.484302           39.05167
2       4.300352           68.67764       2.638243           39.17591
3       4.340642           67.12073       2.813447           41.89578
4             NA                 NA             NA                 NA
5       4.050191           68.90887       3.271357           44.11435
6       4.836835           64.02802       2.732577           35.75704
  c_bilirubin_mg_dl_cont c_chloride_mmol_l_cont c_monocytes_10_9_l_cont
1          -1.4809641766               99.02275             -0.44092220
2           0.0007103677              105.58828             -0.40409105
3          -0.7584989810              100.31545             -0.31518310
4                     NA                     NA                      NA
5           0.2537484544              102.47050             -0.27677341
6          -2.2942350574              101.68754              0.06843379
  c_eosinophils_leukocytes_ratio_cont c_ldh_u_l_cont c_hr_cont c_sbp_cont
1                           1.1417198       1.722156  4.388467   4.773152
2                           0.7942357       1.728584  4.406478   4.814099
3                          -0.4574536       1.697640  4.479523   4.887520
4                                  NA             NA        NA         NA
5                           0.6337259       1.651701  4.393284   4.857790
6                           0.2562202       1.720717  4.350660   4.794559
  c_oxygen_cont c_neutrophil_lymphocyte_ratio_cont c_bmi_cont
1      97.01706                          1.8116984   3.329260
2      97.00790                          0.9970809   3.116247
3      96.99572                          1.4552387   3.307000
4            NA                                 NA         NA
5      96.99230                          1.1772781   3.138414
6      96.99236                          0.9956116   3.096611
  c_ast_alt_ratio_cont c_time_dx_to_index c_de_novo_mets_dx c_height_cont
1           0.25520166          22.333475              TRUE      1.645257
2           0.43717963          -8.733742              TRUE      1.613621
3          -0.27769998          32.036705              TRUE      1.534387
4                   NA                 NA                NA            NA
5           0.01972456          36.686128             FALSE      1.653337
6           0.33510661          37.093396              TRUE      1.513132
  c_weight_cont c_dbp_cont c_year_index fu_itt_months death_itt
1      65.79150   71.65043        2018+     45.673471         1
2      51.40890   85.69923        <2018     19.120776         1
3      57.89357   80.79833        <2018      5.492147         1
4            NA         NA        2018+     26.876450         1
5      85.48502   74.99556        <2018     23.348709         1
6      83.69665   82.97996        2018+      1.395827         1
```
:::
:::


## Moderator covariate

In this example, we assume heterogeneous treatment effect by sex and we aim to assess the average treatment effect among the treated for female and male patients separately. The effect size is time to all-cause mortality. In this dataset, sex is encoded with a binary covariate with 0 = female and 1 = male.


::: {.cell}

```{.r .cell-code}
table(data_miss$dem_sex_cont)
```

::: {.cell-output .cell-output-stdout}
```

   0    1 
2317 1183 
```
:::
:::


## Multiple imputation

Both the imputation and propensity score step Multiple imputation using `mice:`


::: {.cell}

```{.r .cell-code}
# impute data
female_imp <- futuremice(
  parallelseed = 42,
  n.core = parallel::detectCores()-1,
  data = data_miss |> filter(dem_sex_cont == 0),
  method = "rf",
  m = 10,
  print = FALSE
  )
```

::: {.cell-output .cell-output-stderr}
```
Warning: Number of logged events: 1
Warning: Number of logged events: 1
Warning: Number of logged events: 1
Warning: Number of logged events: 1
Warning: Number of logged events: 1
Warning: Number of logged events: 1
Warning: Number of logged events: 1
```
:::

```{.r .cell-code}
male_imp <- futuremice(
  parallelseed = 42,
  n.core = parallel::detectCores()-1,
  data = data_miss |> filter(dem_sex_cont == 1),
  method = "rf",
  m = 10,
  print = FALSE
  )
```

::: {.cell-output .cell-output-stderr}
```
Warning: Number of logged events: 1
Warning: Number of logged events: 1
Warning: Number of logged events: 1
Warning: Number of logged events: 1
Warning: Number of logged events: 1
Warning: Number of logged events: 1
Warning: Number of logged events: 1
```
:::
:::


## Propensity score matching and weighting

Apply propensity score matching and weighting with replacement within in each imputed dataset.


::: {.cell}

```{.r .cell-code}
# apply propensity score matching on mids object
ps_form <- as.formula(paste("treat ~", paste(covariates, collapse = " + ")))
ps_form
```

::: {.cell-output .cell-output-stdout}
```
treat ~ c_smoking_history + c_number_met_sites + c_ecog_cont + 
    c_stage_initial_dx_cont + c_hemoglobin_g_dl_cont + c_urea_nitrogen_mg_dl_cont + 
    c_platelets_10_9_l_cont + c_calcium_mg_dl_cont + c_glucose_mg_dl_cont + 
    c_lymphocyte_leukocyte_ratio_cont + c_alp_u_l_cont + c_protein_g_l_cont + 
    c_alt_u_l_cont + c_albumin_g_l_cont + c_bilirubin_mg_dl_cont + 
    c_chloride_mmol_l_cont + c_monocytes_10_9_l_cont + c_eosinophils_leukocytes_ratio_cont + 
    c_ldh_u_l_cont + c_hr_cont + c_sbp_cont + c_oxygen_cont + 
    c_neutrophil_lymphocyte_ratio_cont + c_bmi_cont + c_ast_alt_ratio_cont + 
    c_time_dx_to_index + c_de_novo_mets_dx + c_height_cont + 
    c_weight_cont + c_dbp_cont + c_year_index + dem_age_index_cont + 
    dem_race + dem_region + dem_ses
```
:::
:::


### Matching


::: {.cell}

```{.r .cell-code}
match_within_strata <- function(i, imputed_data = NULL, ps_formula = NULL){
  
  matched <- MatchIt::matchit(
    formula = ps_formula, 
    data = mice::complete(imputed_data, i),
    method = "nearest",
    caliper = 0.01,
    ratio = 1,
    replace = F
    ) |> 
    MatchIt::match.data()
  
  return(matched)
  
}

female_matched <- lapply(
  X = 1:female_imp$m, 
  FUN = match_within_strata, 
  imputed_data = female_imp,
  ps_formula = ps_form
  )

male_matched <- lapply(
  X = 1:male_imp$m, 
  FUN = match_within_strata, 
  imputed_data = male_imp,
  ps_formula = ps_form
  )

# combine the mth imputed and matched datasets
combine_list <- function(i, data_0 = NULL, data_1 = NULL){
  
  data_combined <- rbind(data_0[[i]], data_1[[i]])
  
  return(data_combined)
  
}

matched_all <- lapply(
  X = 1:female_imp$m, 
  FUN = combine_list, 
  data_0 = female_matched,
  data_1 = male_matched
  )
```
:::


### Weighting


::: {.cell}

```{.r .cell-code}
weight_within_strata <- function(i, imputed_data = NULL, ps_formula = NULL){
  
  weighted <- WeightIt::weightit(
    formula = ps_formula, 
    data = mice::complete(imputed_data, i),
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
    mutate(weights = weighted$weights)
  
  return(weighted_data)
  
}

female_weighted <- lapply(
  X = 1:female_imp$m, 
  FUN = weight_within_strata, 
  imputed_data = female_imp,
  ps_formula = ps_form
  )

male_weighted <- lapply(
  X = 1:male_imp$m, 
  FUN = weight_within_strata, 
  imputed_data = male_imp,
  ps_formula = ps_form
  )

weighted_all <- lapply(
  X = 1:female_imp$m, 
  FUN = combine_list, 
  data_0 = female_weighted,
  data_1 = male_weighted
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
    formula = Surv(fu_itt_months, death_itt) ~ treat*dem_sex_cont, 
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
                term  estimate std.error  conf.low conf.high
1              treat 0.6802007 0.1588437 0.4964419  0.931978
2       dem_sex_cont 0.8241669 0.2029580 0.5506235  1.233604
3 treat:dem_sex_cont 0.9520131 0.3065403 0.5139928  1.763310
```
:::
:::


### Weighting


::: {.cell}

```{.r .cell-code}
cox_fit_weighting <- function(i){
  
  survival_fit <- survival::coxph(
    data = i,
    formula = Surv(fu_itt_months, death_itt) ~ treat*dem_sex_cont, 
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
                term  estimate std.error  conf.low conf.high
1              treat 0.8034467 0.0742781 0.6942947 0.9297589
2       dem_sex_cont 0.9077231 0.1466264 0.6803678 1.2110527
3 treat:dem_sex_cont 0.9427430 0.1552176 0.6948934 1.2789939
```
:::
:::


## References

::: {#refs}
:::

## Session info





Script runtime: 0.43 minutes.

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
_gtable(v.0.3.5)_, _shape(v.1.4.6.1)_, _xfun(v.0.47)_, _ggplot2(v.3.5.1)_, _htmlwidgets(v.1.6.4)_, _MatchIt(v.4.5.5)_, _lattice(v.0.22-6)_, _simsurv(v.1.0.0)_, _vctrs(v.0.6.5)_, _tools(v.4.4.0)_, _generics(v.0.1.3)_, _parallel(v.4.4.0)_, _tibble(v.3.2.1)_, _fansi(v.1.0.6)_, _pan(v.1.9)_, _pkgconfig(v.2.0.3)_, _jomo(v.2.7-6)_, _assertthat(v.0.2.1)_, _lifecycle(v.1.0.4)_, _stringr(v.1.5.1)_, _compiler(v.4.4.0)_, _tictoc(v.1.2.1)_, _munsell(v.0.5.1)_, _mitools(v.2.4)_, _codetools(v.0.2-20)_, _htmltools(v.0.5.8.1)_, _yaml(v.2.3.10)_, _glmnet(v.4.1-8)_, _pillar(v.1.9.0)_, _nloptr(v.2.1.1)_, _crayon(v.1.5.3)_, _tidyr(v.1.3.1)_, _MASS(v.7.3-60.2)_, _sessioninfo(v.1.2.2)_, _iterators(v.1.0.14)_, _rpart(v.4.1.23)_, _boot(v.1.3-30)_, _foreach(v.1.5.2)_, _mitml(v.0.4-5)_, _nlme(v.3.1-164)_, _WeightIt(v.1.3.0)_, _tidyselect(v.1.2.1)_, _digest(v.0.6.37)_, _stringi(v.1.8.4)_, _pander(v.0.6.5)_, _listenv(v.0.9.1)_, _purrr(v.1.0.2)_, _splines(v.4.4.0)_, _rprojroot(v.2.0.4)_, _fastmap(v.1.2.0)_, _colorspace(v.2.1-1)_, _cli(v.3.6.3)_, _magrittr(v.2.0.3)_, _utf8(v.1.2.4)_, _broom(v.1.0.6)_, _withr(v.3.0.1)_, _scales(v.1.3.0)_, _backports(v.1.5.0)_, _rmarkdown(v.2.28)_, _globals(v.0.16.3)_, _nnet(v.7.3-19)_, _lme4(v.1.1-35.5)_, _chk(v.0.9.2)_, _evaluate(v.0.24.0)_, _knitr(v.1.48)_, _rlang(v.1.1.4)_, _Rcpp(v.1.0.13)_, _glue(v.1.7.0)_, _DBI(v.1.2.3)_, _renv(v.1.0.7)_, _rstudioapi(v.0.16.0)_, _minqa(v.1.2.8)_, _jsonlite(v.1.8.8)_ and _R6(v.2.5.1)_
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


---
title: "Re-weighting to a target population"
author: Janick Weberpals, RPh, PhD
date: last-modified
format: html
toc: true
toc-depth: 3
keep-md: true
embed-resources: true
bibliography: ../references.bib
csl: ../pharmacoepidemiology-and-drug-safety.csl
---

## About raking weights

This chapter is a reproducible example on how to incorporate raking weights to match distributions of a target population in multiple imputation \> matching/weighting \> balance assessment \> outcome analysis workflows.

In brief, raking is a procedure to compute weights for multiple variables of interest at the same time to re-weight a population to match pre-specified distributions in complex survey designs. The raking procedure is implemented in the [anesrake](https://surveyinsights.org/wp-content/uploads/2014/07/Full-anesrake-paper.pdf) package. Here, the anesrake function iteratively adjusts the weights to make the weighted sample percentages match the target population percentages for the selected variables. It does this by multiplying the current weight for each case by a factor based on the ratio of the target population proportion to the weighted sample proportion for a given category. This adjustment is performed sequentially for each category of each selected variable. Because adjusting for one variable can disrupt the match for previous variables, the process is repeated through all selected variables in cycles. This iterative process minimizes the Kullback-Leibler (KL) divergence and continues until the weighted sample proportions match the target population proportions for all categories ("full convergence"), or until no further change occurs.

Raking weights could be seen as an alternative to matching-adjusted indirect comparisons (MAIC) which are often used to adjust for differences in baseline characteristics between a trial and a real-world population. However, raking weights are not the same as MAIC, as they do not require a model to estimate the treatment effect in the target population. Instead, raking weights adjust the sample to match the target population without estimating a treatment effect.

## Objective

This chapter aims to introduce the [raking_weights()](https://janickweberpals.github.io/encore.analytics/reference/raking_weights.html) function, which performs the raking procedure in multiply imputed and matched (`mimids`) datasets. The function is a wrapper for the `anesrake()` function from the `anesrake` package and can be used in combination with multiple imputation and propensity score matching/weighting.

::: {.cell}

```{.r .cell-code}
library(here)
library(dplyr)
library(survival)
library(mice)
library(MatchThem)
library(MatchIt)
library(survey)
library(gtsummary)
library(encore.analytics)

source(here("functions", "covariate_vectors.R"))

# track time
runtime <- tictoc::tic()
```
:::

## Data generation

We use the `simulate_data()` function to simulate a realistic oncology comparative effectiveness cohort analytic dataset.



::: {.cell}

```{.r .cell-code}
# load example dataset with missing observations
data_miss <- simulate_data(
  n_total = 3500, 
  seed = 42, 
  include_id = FALSE, 
  imposeNA = TRUE,
  propNA = .33
  )
```
:::


If it is desired ro rake categorical variables, the function works best of those variables are encoded as factor variables. 

In this example, we aim to estimate raking weights for an EHR-derived real-world dataset that aims to mimick the distributions of the [FLAURA trial population, which is one of the prioritized clinical trials to be emulated in ENCORE](https://cdn.clinicaltrials.gov/large-docs/95/NCT06675695/Prot_000.pdf). 

## Defining target distributions

Before estimating raking weights, we need to define the target distributions of patient characteristics that we want to match from the clinical trial using the *raking* procedure. The following distributions are taken from Table 1 of the FLAURA trial.

::: {.cell}

```{.r .cell-code}
# Define FLAURA distributions for key covariates --------------------------
# order is as in Table 1

## age (taken from https://www.tagrissohcp.com/metastatic/flaura/efficacy.html)
# less than 65 years (54%, TRUE) to 65+ (46%, FALSE)
age_target <- c(.54, .46)
names(age_target) <- c("<65", "65+")

## sex ---------------------------------------------------------------------

# female (0) to male (1) proportion:
sex_target <- c(.63, .37) 
names(sex_target) <- c("Female", "Male")

## race --------------------------------------------------------------------
# asian, non-asian
# asian (TRUE) to non-asian (FALSE) proportion
# note: logical variables in dataframe can be matched to a numeric vector of length 2 and ordered with the TRUE target as the first element and the FALSE target as the second element.
race_target <- c(.62, .38)
names(race_target) <- c("Asian", "Non-Asian")

## smoking -----------------------------------------------------------------

# current/former smoker (TRUE) to never smoker (FALSE) proportion
# note: logical variables in dataframe can be matched to a numeric vector of length 2 and ordered with the TRUE target as the first element and the FALSE target as the second element.
smoker_target <- c(.35, .65)
names(smoker_target) <- c("Current/former", "Never")

## ecog --------------------------------------------------------------------

# ecog 0 to ecog 1 proportion
ecog_target <- c(.41, .59)
names(ecog_target) <- c("0", "1")

# summarize target distributions in a named list vector --------------
targets <- list(age_target, sex_target, race_target, smoker_target, ecog_target)
names(targets) <- c("dem_age_lt65", "dem_sex_cont", "dem_race", "c_smoking_history", "c_ecog_cont")

# print
targets
```

::: {.cell-output .cell-output-stdout}

```
$dem_age_lt65
 <65  65+ 
0.54 0.46 

$dem_sex_cont
Female   Male 
  0.63   0.37 

$dem_race
    Asian Non-Asian 
     0.62      0.38 

$c_smoking_history
Current/former          Never 
          0.35           0.65 

$c_ecog_cont
   0    1 
0.41 0.59 
```


:::
:::

Accordingly, we need to convert the categorical target variables in our simulated dataset into factors. The following code chunk shows how to do this for the `data_miss` dataset.

::: {.cell}

```{.r .cell-code}
data_miss <- data_miss |> 
 # anesrake works best with factor variables
  # create age category with age less than 65
  mutate(dem_age_lt65 = factor(ifelse(dem_age_index_cont < 65, "<65", "65+"))) |> 
  # convert dem_race into a binary Asian vs. non-Asian 
  mutate(dem_race = factor(ifelse(dem_race == "Asian", "Asian", "Non-Asian"))) |>
  # convert dem_sex_cont into a factor 
  mutate(dem_sex_cont = factor(ifelse(dem_sex_cont == "1", "Male", "Female"))) |> 
  # convert dem_sex_cont into a factor 
  mutate(c_smoking_history = factor(ifelse(c_smoking_history == TRUE, "Current/former", "Never"))) |> 
  # convert c_ecog_cont into a factor 
  mutate(across(c(c_ecog_cont), function(x) factor(as.character(x))))
```
:::

## Multiple imputation

Since we operate on a simulated dataset with missing observations, we first impute the missing data before performing the propensity score matching and raking. We use the `mice` package for this purpose.


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


## Propensity score matching/weighting

In the next step, we create `mimids` and `wimids` object which contain the imputed and 1:1 propensity score matched/SMR-weighted datasets that will serve as the input for the raking procedure. The propensity score model is specified as follows:


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
    c_ast_alt_ratio_cont + c_stage_initial_dx_cont + dem_race + 
    dem_region + dem_ses + c_time_dx_to_index
```


:::
:::



As already covered in the previous chapters, we can use the `matchthem()` function to perform the matching and weighting.

::: panel-tabset
### Propensity score matching


::: {.cell}

```{.r .cell-code}
# matching
mimids_data <- matchthem(
  formula = ps_form,
  datasets = data_imp,
  approach = 'within',
  method = 'nearest',
  distance = "glm",
  link = "logit",
  caliper = 0.01,
  ratio = 1,
  replace = F
  )

# print summary for matched dataset #1
mimids_data
```

::: {.cell-output .cell-output-stdout}

```
A `matchit` object
 - method: 1:1 nearest neighbor matching without replacement
 - distance: Propensity score [caliper]

             - estimated with logistic regression
 - caliper: <distance> (0.001)
 - number of obs.: 3500 (original), 2686 (matched)
 - target estimand: ATT
 - covariates: dem_age_index_cont, dem_sex_cont, c_smoking_history, c_number_met_sites, c_hemoglobin_g_dl_cont, c_urea_nitrogen_mg_dl_cont, c_platelets_10_9_l_cont, c_calcium_mg_dl_cont, c_glucose_mg_dl_cont, c_lymphocyte_leukocyte_ratio_cont, c_alp_u_l_cont, c_protein_g_l_cont, c_alt_u_l_cont, c_albumin_g_l_cont, c_bilirubin_mg_dl_cont, c_chloride_mmol_l_cont, c_monocytes_10_9_l_cont, c_eosinophils_leukocytes_ratio_cont, c_ldh_u_l_cont, c_hr_cont, c_sbp_cont, c_oxygen_cont, c_ecog_cont, c_neutrophil_lymphocyte_ratio_cont, c_bmi_cont, c_ast_alt_ratio_cont, c_stage_initial_dx_cont, dem_race, dem_region, dem_ses, c_time_dx_to_index
```


:::
:::


### Propensity score weighting


::: {.cell}

```{.r .cell-code}
# SMR weighting
wimids_data <- weightthem(
  formula = ps_form,
  datasets = data_imp,
  approach = 'within',
  method = "glm",
  estimand = "ATT"
  )

# trim extreme weights
wimids_data <- trim(
  x = wimids_data, 
  at = .95, 
  lower = TRUE
  )

wimids_data
```

::: {.cell-output .cell-output-stdout}

```
A weightit object
 - method: "glm" (propensity score weighting with GLM)
 - number of obs.: 3500
 - sampling weights: none
 - treatment: 2-category
 - estimand: ATT (focal: 1)
 - covariates: dem_age_index_cont, dem_sex_cont, c_smoking_history, c_number_met_sites, c_hemoglobin_g_dl_cont, c_urea_nitrogen_mg_dl_cont, c_platelets_10_9_l_cont, c_calcium_mg_dl_cont, c_glucose_mg_dl_cont, c_lymphocyte_leukocyte_ratio_cont, c_alp_u_l_cont, c_protein_g_l_cont, c_alt_u_l_cont, c_albumin_g_l_cont, c_bilirubin_mg_dl_cont, c_chloride_mmol_l_cont, c_monocytes_10_9_l_cont, c_eosinophils_leukocytes_ratio_cont, c_ldh_u_l_cont, c_hr_cont, c_sbp_cont, c_oxygen_cont, c_ecog_cont, c_neutrophil_lymphocyte_ratio_cont, c_bmi_cont, c_ast_alt_ratio_cont, c_stage_initial_dx_cont, dem_race, dem_region, dem_ses, c_time_dx_to_index
 - weights trimmed at 5% and 95%
```


:::
:::

:::

### Raking procedure

The `raking_weights()` function as implemented in encore.analytics takes the `mimids` and `wimids` objects and the target distributions as input and returns a list of imputed, matched/SMR-weighted, and raked datasets. The function uses the `anesrake()` function from the `anesrake` package to perform the raking procedure. The only two inputs it needs is the multiply imputed and matched/SMR-weighted datasets that results from the `matchthem()` function and the target distributions which are defined in a named list. The function will then iterate over the datasets and apply the raking procedure to each dataset. The raking procedure is performed by calling the `anesrake()` function for each dataset in the `mimids` and `wimids` object. 

After successful raking, the function returns list of imputed, matched and re-weighted datasets (which we refer to here as `mirwds` and `wirwds`) with an (updated) column `weights` which includes the raking weights.

::: {.cell}

```{.r .cell-code}
# raking on the imputed and 1:1 propensity score-matched datasets
mirwds <- raking_weights(
  x = mimids_data, 
  targets = targets
  )
```

::: {.cell-output .cell-output-stdout}

```
[1] "Raking converged in 11 iterations"
[1] "Raking converged in 10 iterations"
[1] "Raking converged in 13 iterations"
[1] "Raking converged in 12 iterations"
[1] "Raking converged in 10 iterations"
[1] "Raking converged in 10 iterations"
[1] "Raking converged in 14 iterations"
[1] "Raking converged in 11 iterations"
[1] "Raking converged in 10 iterations"
[1] "Raking converged in 13 iterations"
[1] "Raking converged in 8 iterations"
[1] "Raking converged in 10 iterations"
[1] "Raking converged in 11 iterations"
```


:::

```{.r .cell-code}
# raking on the imputed and SMR-weighted datasets
wirwds <- raking_weights(
  x = wimids_data, 
  targets = targets
  )
```

::: {.cell-output .cell-output-stdout}

```
[1] "Raking converged in 12 iterations"
[1] "Raking converged in 11 iterations"
[1] "Raking converged in 14 iterations"
[1] "Raking converged in 11 iterations"
[1] "Raking converged in 13 iterations"
[1] "Raking converged in 13 iterations"
[1] "Raking converged in 13 iterations"
[1] "Raking converged in 19 iterations"
[1] "Raking converged in 15 iterations"
[1] "Raking converged in 15 iterations"
[1] "Raking converged in 15 iterations"
```


:::
:::

## Comparing the final distributions in Table 1

To assess the balance of the covariates before and after raking, we can use the `tbl_summary()` and `tbl_svysummary` functions from the `gtsummary` package. To illustrate the distributions in the datasets before and after raking, we use the first imputed and matched dataset before and after raking as an example.

> Reminder : The target distributions look like this

::: {.cell}

```{.r .cell-code}
targets
```

::: {.cell-output .cell-output-stdout}

```
$dem_age_lt65
 <65  65+ 
0.54 0.46 

$dem_sex_cont
Female   Male 
  0.63   0.37 

$dem_race
    Asian Non-Asian 
     0.62      0.38 

$c_smoking_history
Current/former          Never 
          0.35           0.65 

$c_ecog_cont
   0    1 
0.41 0.59 
```


:::
:::

::: panel-tabset
### Propensity score-matched datasets

::: panel-tabset
#### Table 1 BEFORE raking

::: {#tbl-matched-before-raking .cell tbl-cap='Table 1 BEFORE raking (1:1 propensity score-matched data)'}

```{.r .cell-code}
# extract the first imputed and matched dataset BEFORE raking
first_dataset <- MatchThem::complete(mimids_data, action = 1, all = FALSE)

# print
first_dataset |>
  tbl_summary(
    by = treat,
    include = c(dem_age_index_cont, names(targets))
    ) |> 
  add_difference(test = dplyr::everything() ~ "smd") |>
  add_overall() |>
  modify_column_hide(columns = "conf.low") |> 
  modify_header(
    label ~ "**Patient characteristic**",
    stat_0 ~ "**Total** <br> N = {round(N, 2)}",
    stat_1 ~ "**{level}** <br> N = {round(n, 2)} <br> ({style_percent(p, digits=1)}%)",
    stat_2 ~ "**{level}** <br> N = {round(n, 2)} <br> ({style_percent(p, digits=1)}%)"
    ) |>
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Treatment received**")
```

::: {.cell-output-display}

```{=html}
<div id="ayjdrzqoxg" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#ayjdrzqoxg table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#ayjdrzqoxg thead, #ayjdrzqoxg tbody, #ayjdrzqoxg tfoot, #ayjdrzqoxg tr, #ayjdrzqoxg td, #ayjdrzqoxg th {
  border-style: none;
}

#ayjdrzqoxg p {
  margin: 0;
  padding: 0;
}

#ayjdrzqoxg .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#ayjdrzqoxg .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#ayjdrzqoxg .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#ayjdrzqoxg .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#ayjdrzqoxg .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#ayjdrzqoxg .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ayjdrzqoxg .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#ayjdrzqoxg .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#ayjdrzqoxg .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#ayjdrzqoxg .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#ayjdrzqoxg .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#ayjdrzqoxg .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#ayjdrzqoxg .gt_spanner_row {
  border-bottom-style: hidden;
}

#ayjdrzqoxg .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#ayjdrzqoxg .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#ayjdrzqoxg .gt_from_md > :first-child {
  margin-top: 0;
}

#ayjdrzqoxg .gt_from_md > :last-child {
  margin-bottom: 0;
}

#ayjdrzqoxg .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#ayjdrzqoxg .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#ayjdrzqoxg .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#ayjdrzqoxg .gt_row_group_first td {
  border-top-width: 2px;
}

#ayjdrzqoxg .gt_row_group_first th {
  border-top-width: 2px;
}

#ayjdrzqoxg .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ayjdrzqoxg .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#ayjdrzqoxg .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#ayjdrzqoxg .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ayjdrzqoxg .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ayjdrzqoxg .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#ayjdrzqoxg .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#ayjdrzqoxg .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#ayjdrzqoxg .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ayjdrzqoxg .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#ayjdrzqoxg .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#ayjdrzqoxg .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#ayjdrzqoxg .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#ayjdrzqoxg .gt_left {
  text-align: left;
}

#ayjdrzqoxg .gt_center {
  text-align: center;
}

#ayjdrzqoxg .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#ayjdrzqoxg .gt_font_normal {
  font-weight: normal;
}

#ayjdrzqoxg .gt_font_bold {
  font-weight: bold;
}

#ayjdrzqoxg .gt_font_italic {
  font-style: italic;
}

#ayjdrzqoxg .gt_super {
  font-size: 65%;
}

#ayjdrzqoxg .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#ayjdrzqoxg .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#ayjdrzqoxg .gt_indent_1 {
  text-indent: 5px;
}

#ayjdrzqoxg .gt_indent_2 {
  text-indent: 10px;
}

#ayjdrzqoxg .gt_indent_3 {
  text-indent: 15px;
}

#ayjdrzqoxg .gt_indent_4 {
  text-indent: 20px;
}

#ayjdrzqoxg .gt_indent_5 {
  text-indent: 25px;
}

#ayjdrzqoxg .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#ayjdrzqoxg div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings gt_spanner_row">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1" scope="col" id="label"><span data-qmd-base64="KipQYXRpZW50IGNoYXJhY3RlcmlzdGljKio="><span class='gt_from_md'><strong>Patient characteristic</strong></span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="2" colspan="1" scope="col" id="stat_0"><span data-qmd-base64="KipUb3RhbCoqIDxicj4gTiA9IDI2ODY="><span class='gt_from_md'><strong>Total</strong> <br> N = 2686</span></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2" scope="colgroup" id="level 1; stat_1">
        <div class="gt_column_spanner"><span data-qmd-base64="KipUcmVhdG1lbnQgcmVjZWl2ZWQqKg=="><span class='gt_from_md'><strong>Treatment received</strong></span></span></div>
      </th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="2" colspan="1" scope="col" id="estimate"><span data-qmd-base64="KipEaWZmZXJlbmNlKio="><span class='gt_from_md'><strong>Difference</strong></span></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>2</sup></span></th>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="stat_1"><span data-qmd-base64="KiowKiogPGJyPiBOID0gMTM0MyA8YnI+ICg1MC4wJSk="><span class='gt_from_md'><strong>0</strong> <br> N = 1343 <br> (50.0%)</span></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="stat_2"><span data-qmd-base64="KioxKiogPGJyPiBOID0gMTM0MyA8YnI+ICg1MC4wJSk="><span class='gt_from_md'><strong>1</strong> <br> N = 1343 <br> (50.0%)</span></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="label" class="gt_row gt_left">dem_age_index_cont</td>
<td headers="stat_0" class="gt_row gt_center">69 (64, 74)</td>
<td headers="stat_1" class="gt_row gt_center">69 (64, 74)</td>
<td headers="stat_2" class="gt_row gt_center">69 (63, 74)</td>
<td headers="estimate" class="gt_row gt_center">0.00</td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_age_lt65</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.01</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    &lt;65</td>
<td headers="stat_0" class="gt_row gt_center">840 (31%)</td>
<td headers="stat_1" class="gt_row gt_center">423 (31%)</td>
<td headers="stat_2" class="gt_row gt_center">417 (31%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    65+</td>
<td headers="stat_0" class="gt_row gt_center">1,846 (69%)</td>
<td headers="stat_1" class="gt_row gt_center">920 (69%)</td>
<td headers="stat_2" class="gt_row gt_center">926 (69%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_sex_cont</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.03</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Female</td>
<td headers="stat_0" class="gt_row gt_center">1,817 (68%)</td>
<td headers="stat_1" class="gt_row gt_center">898 (67%)</td>
<td headers="stat_2" class="gt_row gt_center">919 (68%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Male</td>
<td headers="stat_0" class="gt_row gt_center">869 (32%)</td>
<td headers="stat_1" class="gt_row gt_center">445 (33%)</td>
<td headers="stat_2" class="gt_row gt_center">424 (32%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_race</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.04</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Asian</td>
<td headers="stat_0" class="gt_row gt_center">990 (37%)</td>
<td headers="stat_1" class="gt_row gt_center">507 (38%)</td>
<td headers="stat_2" class="gt_row gt_center">483 (36%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Non-Asian</td>
<td headers="stat_0" class="gt_row gt_center">1,696 (63%)</td>
<td headers="stat_1" class="gt_row gt_center">836 (62%)</td>
<td headers="stat_2" class="gt_row gt_center">860 (64%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_smoking_history</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.01</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Current/former</td>
<td headers="stat_0" class="gt_row gt_center">1,255 (47%)</td>
<td headers="stat_1" class="gt_row gt_center">632 (47%)</td>
<td headers="stat_2" class="gt_row gt_center">623 (46%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Never</td>
<td headers="stat_0" class="gt_row gt_center">1,431 (53%)</td>
<td headers="stat_1" class="gt_row gt_center">711 (53%)</td>
<td headers="stat_2" class="gt_row gt_center">720 (54%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_ecog_cont</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.01</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    0</td>
<td headers="stat_0" class="gt_row gt_center">1,178 (44%)</td>
<td headers="stat_1" class="gt_row gt_center">587 (44%)</td>
<td headers="stat_2" class="gt_row gt_center">591 (44%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    1</td>
<td headers="stat_0" class="gt_row gt_center">1,508 (56%)</td>
<td headers="stat_1" class="gt_row gt_center">756 (56%)</td>
<td headers="stat_2" class="gt_row gt_center">752 (56%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
  </tbody>
  <tfoot class="gt_sourcenotes">
    <tr>
      <td class="gt_sourcenote" colspan="5"><span data-qmd-base64="QWJicmV2aWF0aW9uOiBDSSA9IENvbmZpZGVuY2UgSW50ZXJ2YWw="><span class='gt_from_md'>Abbreviation: CI = Confidence Interval</span></span></td>
    </tr>
  </tfoot>
  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="5"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span> <span data-qmd-base64="TWVkaWFuIChRMSwgUTMpOyBuICglKQ=="><span class='gt_from_md'>Median (Q1, Q3); n (%)</span></span></td>
    </tr>
    <tr>
      <td class="gt_footnote" colspan="5"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>2</sup></span> <span data-qmd-base64="U3RhbmRhcmRpemVkIE1lYW4gRGlmZmVyZW5jZQ=="><span class='gt_from_md'>Standardized Mean Difference</span></span></td>
    </tr>
  </tfoot>
</table>
</div>
```

:::
:::

#### Table 1 AFTER raking

::: {#tbl-matched-after-raking .cell tbl-cap='Table 1 AFTER raking (1:1 propensity score-matched data)'}

```{.r .cell-code}
# create survey object 
data_svy <- svydesign(ids = ~ 1, weights = ~ weights, data = mirwds[[1]])

# print
data_svy |>
  tbl_svysummary(
    by = treat,
    include = c(dem_age_index_cont, names(targets))
    ) |> 
  add_difference(test = dplyr::everything() ~ "smd") |>
  add_overall() |>
  modify_column_hide(columns = "conf.low") |> 
  modify_header(
    label ~ "**Patient characteristic**",
    stat_0 ~ "**Total** <br> N = {round(N, 2)}",
    stat_1 ~ "**{level}** <br> N = {round(n, 2)} <br> ({style_percent(p, digits=1)}%)",
    stat_2 ~ "**{level}** <br> N = {round(n, 2)} <br> ({style_percent(p, digits=1)}%)"
    ) |>
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Treatment received**")
```

::: {.cell-output-display}

```{=html}
<div id="bjwlziabxx" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#bjwlziabxx table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#bjwlziabxx thead, #bjwlziabxx tbody, #bjwlziabxx tfoot, #bjwlziabxx tr, #bjwlziabxx td, #bjwlziabxx th {
  border-style: none;
}

#bjwlziabxx p {
  margin: 0;
  padding: 0;
}

#bjwlziabxx .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#bjwlziabxx .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#bjwlziabxx .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#bjwlziabxx .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#bjwlziabxx .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#bjwlziabxx .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#bjwlziabxx .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#bjwlziabxx .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#bjwlziabxx .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#bjwlziabxx .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#bjwlziabxx .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#bjwlziabxx .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#bjwlziabxx .gt_spanner_row {
  border-bottom-style: hidden;
}

#bjwlziabxx .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#bjwlziabxx .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#bjwlziabxx .gt_from_md > :first-child {
  margin-top: 0;
}

#bjwlziabxx .gt_from_md > :last-child {
  margin-bottom: 0;
}

#bjwlziabxx .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#bjwlziabxx .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#bjwlziabxx .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#bjwlziabxx .gt_row_group_first td {
  border-top-width: 2px;
}

#bjwlziabxx .gt_row_group_first th {
  border-top-width: 2px;
}

#bjwlziabxx .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#bjwlziabxx .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#bjwlziabxx .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#bjwlziabxx .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#bjwlziabxx .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#bjwlziabxx .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#bjwlziabxx .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#bjwlziabxx .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#bjwlziabxx .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#bjwlziabxx .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#bjwlziabxx .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#bjwlziabxx .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#bjwlziabxx .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#bjwlziabxx .gt_left {
  text-align: left;
}

#bjwlziabxx .gt_center {
  text-align: center;
}

#bjwlziabxx .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#bjwlziabxx .gt_font_normal {
  font-weight: normal;
}

#bjwlziabxx .gt_font_bold {
  font-weight: bold;
}

#bjwlziabxx .gt_font_italic {
  font-style: italic;
}

#bjwlziabxx .gt_super {
  font-size: 65%;
}

#bjwlziabxx .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#bjwlziabxx .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#bjwlziabxx .gt_indent_1 {
  text-indent: 5px;
}

#bjwlziabxx .gt_indent_2 {
  text-indent: 10px;
}

#bjwlziabxx .gt_indent_3 {
  text-indent: 15px;
}

#bjwlziabxx .gt_indent_4 {
  text-indent: 20px;
}

#bjwlziabxx .gt_indent_5 {
  text-indent: 25px;
}

#bjwlziabxx .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#bjwlziabxx div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings gt_spanner_row">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1" scope="col" id="label"><span data-qmd-base64="KipQYXRpZW50IGNoYXJhY3RlcmlzdGljKio="><span class='gt_from_md'><strong>Patient characteristic</strong></span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="2" colspan="1" scope="col" id="stat_0"><span data-qmd-base64="KipUb3RhbCoqIDxicj4gTiA9IDI2ODY="><span class='gt_from_md'><strong>Total</strong> <br> N = 2686</span></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2" scope="colgroup" id="level 1; stat_1">
        <div class="gt_column_spanner"><span data-qmd-base64="KipUcmVhdG1lbnQgcmVjZWl2ZWQqKg=="><span class='gt_from_md'><strong>Treatment received</strong></span></span></div>
      </th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="2" colspan="1" scope="col" id="estimate"><span data-qmd-base64="KipEaWZmZXJlbmNlKio="><span class='gt_from_md'><strong>Difference</strong></span></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>2</sup></span></th>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="stat_1"><span data-qmd-base64="KiowKiogPGJyPiBOID0gMTM2My4yNCA8YnI+ICg1MC44JSk="><span class='gt_from_md'><strong>0</strong> <br> N = 1363.24 <br> (50.8%)</span></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="stat_2"><span data-qmd-base64="KioxKiogPGJyPiBOID0gMTMyMi43NiA8YnI+ICg0OS4yJSk="><span class='gt_from_md'><strong>1</strong> <br> N = 1322.76 <br> (49.2%)</span></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="label" class="gt_row gt_left">dem_age_index_cont</td>
<td headers="stat_0" class="gt_row gt_center">64 (60, 71)</td>
<td headers="stat_1" class="gt_row gt_center">64 (60, 71)</td>
<td headers="stat_2" class="gt_row gt_center">64 (60, 71)</td>
<td headers="estimate" class="gt_row gt_center">-0.02</td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_age_lt65</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.02</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    &lt;65</td>
<td headers="stat_0" class="gt_row gt_center">1,450 (54%)</td>
<td headers="stat_1" class="gt_row gt_center">744 (55%)</td>
<td headers="stat_2" class="gt_row gt_center">707 (53%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    65+</td>
<td headers="stat_0" class="gt_row gt_center">1,236 (46%)</td>
<td headers="stat_1" class="gt_row gt_center">619 (45%)</td>
<td headers="stat_2" class="gt_row gt_center">616 (47%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_sex_cont</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.11</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Female</td>
<td headers="stat_0" class="gt_row gt_center">1,692 (63%)</td>
<td headers="stat_1" class="gt_row gt_center">824 (60%)</td>
<td headers="stat_2" class="gt_row gt_center">868 (66%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Male</td>
<td headers="stat_0" class="gt_row gt_center">994 (37%)</td>
<td headers="stat_1" class="gt_row gt_center">539 (40%)</td>
<td headers="stat_2" class="gt_row gt_center">455 (34%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_race</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.05</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Asian</td>
<td headers="stat_0" class="gt_row gt_center">1,665 (62%)</td>
<td headers="stat_1" class="gt_row gt_center">862 (63%)</td>
<td headers="stat_2" class="gt_row gt_center">804 (61%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Non-Asian</td>
<td headers="stat_0" class="gt_row gt_center">1,021 (38%)</td>
<td headers="stat_1" class="gt_row gt_center">502 (37%)</td>
<td headers="stat_2" class="gt_row gt_center">519 (39%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_smoking_history</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.05</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Current/former</td>
<td headers="stat_0" class="gt_row gt_center">940 (35%)</td>
<td headers="stat_1" class="gt_row gt_center">493 (36%)</td>
<td headers="stat_2" class="gt_row gt_center">447 (34%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Never</td>
<td headers="stat_0" class="gt_row gt_center">1,746 (65%)</td>
<td headers="stat_1" class="gt_row gt_center">870 (64%)</td>
<td headers="stat_2" class="gt_row gt_center">875 (66%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_ecog_cont</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.05</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    0</td>
<td headers="stat_0" class="gt_row gt_center">1,101 (41%)</td>
<td headers="stat_1" class="gt_row gt_center">575 (42%)</td>
<td headers="stat_2" class="gt_row gt_center">526 (40%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    1</td>
<td headers="stat_0" class="gt_row gt_center">1,585 (59%)</td>
<td headers="stat_1" class="gt_row gt_center">788 (58%)</td>
<td headers="stat_2" class="gt_row gt_center">797 (60%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
  </tbody>
  <tfoot class="gt_sourcenotes">
    <tr>
      <td class="gt_sourcenote" colspan="5"><span data-qmd-base64="QWJicmV2aWF0aW9uOiBDSSA9IENvbmZpZGVuY2UgSW50ZXJ2YWw="><span class='gt_from_md'>Abbreviation: CI = Confidence Interval</span></span></td>
    </tr>
  </tfoot>
  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="5"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span> <span data-qmd-base64="TWVkaWFuIChRMSwgUTMpOyBuICglKQ=="><span class='gt_from_md'>Median (Q1, Q3); n (%)</span></span></td>
    </tr>
    <tr>
      <td class="gt_footnote" colspan="5"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>2</sup></span> <span data-qmd-base64="U3RhbmRhcmRpemVkIE1lYW4gRGlmZmVyZW5jZQ=="><span class='gt_from_md'>Standardized Mean Difference</span></span></td>
    </tr>
  </tfoot>
</table>
</div>
```

:::
:::
:::

### SMR-weighted datasets

::: panel-tabset
#### Table 1 BEFORE raking

::: {#tbl-smr-before-raking .cell tbl-cap='Table 1 BEFORE raking (SMR-weighted data)'}

```{.r .cell-code}
# extract the first imputed and matched dataset BEFORE raking
first_dataset <- MatchThem::complete(wimids_data, action = 1, all = FALSE)

# print
first_dataset |>
  tbl_summary(
    by = treat,
    include = c(dem_age_index_cont, names(targets))
    ) |> 
  add_difference(test = dplyr::everything() ~ "smd") |>
  add_overall() |>
  modify_column_hide(columns = "conf.low") |> 
  modify_header(
    label ~ "**Patient characteristic**",
    stat_0 ~ "**Total** <br> N = {round(N, 2)}",
    stat_1 ~ "**{level}** <br> N = {round(n, 2)} <br> ({style_percent(p, digits=1)}%)",
    stat_2 ~ "**{level}** <br> N = {round(n, 2)} <br> ({style_percent(p, digits=1)}%)"
    ) |>
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Treatment received**")
```

::: {.cell-output-display}

```{=html}
<div id="jvyvqulufo" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#jvyvqulufo table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#jvyvqulufo thead, #jvyvqulufo tbody, #jvyvqulufo tfoot, #jvyvqulufo tr, #jvyvqulufo td, #jvyvqulufo th {
  border-style: none;
}

#jvyvqulufo p {
  margin: 0;
  padding: 0;
}

#jvyvqulufo .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#jvyvqulufo .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#jvyvqulufo .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#jvyvqulufo .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#jvyvqulufo .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#jvyvqulufo .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#jvyvqulufo .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#jvyvqulufo .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#jvyvqulufo .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#jvyvqulufo .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#jvyvqulufo .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#jvyvqulufo .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#jvyvqulufo .gt_spanner_row {
  border-bottom-style: hidden;
}

#jvyvqulufo .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#jvyvqulufo .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#jvyvqulufo .gt_from_md > :first-child {
  margin-top: 0;
}

#jvyvqulufo .gt_from_md > :last-child {
  margin-bottom: 0;
}

#jvyvqulufo .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#jvyvqulufo .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#jvyvqulufo .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#jvyvqulufo .gt_row_group_first td {
  border-top-width: 2px;
}

#jvyvqulufo .gt_row_group_first th {
  border-top-width: 2px;
}

#jvyvqulufo .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#jvyvqulufo .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#jvyvqulufo .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#jvyvqulufo .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#jvyvqulufo .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#jvyvqulufo .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#jvyvqulufo .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#jvyvqulufo .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#jvyvqulufo .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#jvyvqulufo .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#jvyvqulufo .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#jvyvqulufo .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#jvyvqulufo .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#jvyvqulufo .gt_left {
  text-align: left;
}

#jvyvqulufo .gt_center {
  text-align: center;
}

#jvyvqulufo .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#jvyvqulufo .gt_font_normal {
  font-weight: normal;
}

#jvyvqulufo .gt_font_bold {
  font-weight: bold;
}

#jvyvqulufo .gt_font_italic {
  font-style: italic;
}

#jvyvqulufo .gt_super {
  font-size: 65%;
}

#jvyvqulufo .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#jvyvqulufo .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#jvyvqulufo .gt_indent_1 {
  text-indent: 5px;
}

#jvyvqulufo .gt_indent_2 {
  text-indent: 10px;
}

#jvyvqulufo .gt_indent_3 {
  text-indent: 15px;
}

#jvyvqulufo .gt_indent_4 {
  text-indent: 20px;
}

#jvyvqulufo .gt_indent_5 {
  text-indent: 25px;
}

#jvyvqulufo .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#jvyvqulufo div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings gt_spanner_row">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1" scope="col" id="label"><span data-qmd-base64="KipQYXRpZW50IGNoYXJhY3RlcmlzdGljKio="><span class='gt_from_md'><strong>Patient characteristic</strong></span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="2" colspan="1" scope="col" id="stat_0"><span data-qmd-base64="KipUb3RhbCoqIDxicj4gTiA9IDM1MDA="><span class='gt_from_md'><strong>Total</strong> <br> N = 3500</span></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2" scope="colgroup" id="level 1; stat_1">
        <div class="gt_column_spanner"><span data-qmd-base64="KipUcmVhdG1lbnQgcmVjZWl2ZWQqKg=="><span class='gt_from_md'><strong>Treatment received</strong></span></span></div>
      </th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="2" colspan="1" scope="col" id="estimate"><span data-qmd-base64="KipEaWZmZXJlbmNlKio="><span class='gt_from_md'><strong>Difference</strong></span></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>2</sup></span></th>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="stat_1"><span data-qmd-base64="KiowKiogPGJyPiBOID0gMTQ4NyA8YnI+ICg0Mi41JSk="><span class='gt_from_md'><strong>0</strong> <br> N = 1487 <br> (42.5%)</span></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="stat_2"><span data-qmd-base64="KioxKiogPGJyPiBOID0gMjAxMyA8YnI+ICg1Ny41JSk="><span class='gt_from_md'><strong>1</strong> <br> N = 2013 <br> (57.5%)</span></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="label" class="gt_row gt_left">dem_age_index_cont</td>
<td headers="stat_0" class="gt_row gt_center">69 (64, 74)</td>
<td headers="stat_1" class="gt_row gt_center">69 (64, 74)</td>
<td headers="stat_2" class="gt_row gt_center">69 (64, 74)</td>
<td headers="estimate" class="gt_row gt_center">-0.04</td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_age_lt65</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.03</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    &lt;65</td>
<td headers="stat_0" class="gt_row gt_center">1,085 (31%)</td>
<td headers="stat_1" class="gt_row gt_center">472 (32%)</td>
<td headers="stat_2" class="gt_row gt_center">613 (30%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    65+</td>
<td headers="stat_0" class="gt_row gt_center">2,415 (69%)</td>
<td headers="stat_1" class="gt_row gt_center">1,015 (68%)</td>
<td headers="stat_2" class="gt_row gt_center">1,400 (70%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_sex_cont</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.02</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Female</td>
<td headers="stat_0" class="gt_row gt_center">2,354 (67%)</td>
<td headers="stat_1" class="gt_row gt_center">993 (67%)</td>
<td headers="stat_2" class="gt_row gt_center">1,361 (68%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Male</td>
<td headers="stat_0" class="gt_row gt_center">1,146 (33%)</td>
<td headers="stat_1" class="gt_row gt_center">494 (33%)</td>
<td headers="stat_2" class="gt_row gt_center">652 (32%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_race</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.04</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Asian</td>
<td headers="stat_0" class="gt_row gt_center">1,330 (38%)</td>
<td headers="stat_1" class="gt_row gt_center">550 (37%)</td>
<td headers="stat_2" class="gt_row gt_center">780 (39%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Non-Asian</td>
<td headers="stat_0" class="gt_row gt_center">2,170 (62%)</td>
<td headers="stat_1" class="gt_row gt_center">937 (63%)</td>
<td headers="stat_2" class="gt_row gt_center">1,233 (61%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_smoking_history</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.10</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Current/former</td>
<td headers="stat_0" class="gt_row gt_center">1,616 (46%)</td>
<td headers="stat_1" class="gt_row gt_center">728 (49%)</td>
<td headers="stat_2" class="gt_row gt_center">888 (44%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Never</td>
<td headers="stat_0" class="gt_row gt_center">1,884 (54%)</td>
<td headers="stat_1" class="gt_row gt_center">759 (51%)</td>
<td headers="stat_2" class="gt_row gt_center">1,125 (56%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_ecog_cont</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.06</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    0</td>
<td headers="stat_0" class="gt_row gt_center">1,561 (45%)</td>
<td headers="stat_1" class="gt_row gt_center">636 (43%)</td>
<td headers="stat_2" class="gt_row gt_center">925 (46%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    1</td>
<td headers="stat_0" class="gt_row gt_center">1,939 (55%)</td>
<td headers="stat_1" class="gt_row gt_center">851 (57%)</td>
<td headers="stat_2" class="gt_row gt_center">1,088 (54%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
  </tbody>
  <tfoot class="gt_sourcenotes">
    <tr>
      <td class="gt_sourcenote" colspan="5"><span data-qmd-base64="QWJicmV2aWF0aW9uOiBDSSA9IENvbmZpZGVuY2UgSW50ZXJ2YWw="><span class='gt_from_md'>Abbreviation: CI = Confidence Interval</span></span></td>
    </tr>
  </tfoot>
  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="5"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span> <span data-qmd-base64="TWVkaWFuIChRMSwgUTMpOyBuICglKQ=="><span class='gt_from_md'>Median (Q1, Q3); n (%)</span></span></td>
    </tr>
    <tr>
      <td class="gt_footnote" colspan="5"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>2</sup></span> <span data-qmd-base64="U3RhbmRhcmRpemVkIE1lYW4gRGlmZmVyZW5jZQ=="><span class='gt_from_md'>Standardized Mean Difference</span></span></td>
    </tr>
  </tfoot>
</table>
</div>
```

:::
:::

#### Table 1 AFTER raking

::: {#tbl-smr-after-raking .cell tbl-cap='Table 1 AFTER raking (SMR-weighted data)'}

```{.r .cell-code}
# create survey object 
data_svy <- svydesign(ids = ~ 1, weights = ~ weights, data = wirwds[[1]])

# print
data_svy |>
  tbl_svysummary(
    by = treat,
    include = c(dem_age_index_cont, names(targets))
    ) |> 
  add_difference(test = dplyr::everything() ~ "smd") |>
  add_overall() |>
  modify_column_hide(columns = "conf.low") |> 
  modify_header(
    label ~ "**Patient characteristic**",
    stat_0 ~ "**Total** <br> N = {round(N, 2)}",
    stat_1 ~ "**{level}** <br> N = {round(n, 2)} <br> ({style_percent(p, digits=1)}%)",
    stat_2 ~ "**{level}** <br> N = {round(n, 2)} <br> ({style_percent(p, digits=1)}%)"
    ) |>
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Treatment received**")
```

::: {.cell-output-display}

```{=html}
<div id="tinncrjaht" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#tinncrjaht table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#tinncrjaht thead, #tinncrjaht tbody, #tinncrjaht tfoot, #tinncrjaht tr, #tinncrjaht td, #tinncrjaht th {
  border-style: none;
}

#tinncrjaht p {
  margin: 0;
  padding: 0;
}

#tinncrjaht .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#tinncrjaht .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#tinncrjaht .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#tinncrjaht .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#tinncrjaht .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#tinncrjaht .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#tinncrjaht .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#tinncrjaht .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#tinncrjaht .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#tinncrjaht .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#tinncrjaht .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#tinncrjaht .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#tinncrjaht .gt_spanner_row {
  border-bottom-style: hidden;
}

#tinncrjaht .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#tinncrjaht .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#tinncrjaht .gt_from_md > :first-child {
  margin-top: 0;
}

#tinncrjaht .gt_from_md > :last-child {
  margin-bottom: 0;
}

#tinncrjaht .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#tinncrjaht .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#tinncrjaht .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#tinncrjaht .gt_row_group_first td {
  border-top-width: 2px;
}

#tinncrjaht .gt_row_group_first th {
  border-top-width: 2px;
}

#tinncrjaht .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#tinncrjaht .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#tinncrjaht .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#tinncrjaht .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#tinncrjaht .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#tinncrjaht .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#tinncrjaht .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#tinncrjaht .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#tinncrjaht .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#tinncrjaht .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#tinncrjaht .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#tinncrjaht .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#tinncrjaht .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#tinncrjaht .gt_left {
  text-align: left;
}

#tinncrjaht .gt_center {
  text-align: center;
}

#tinncrjaht .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#tinncrjaht .gt_font_normal {
  font-weight: normal;
}

#tinncrjaht .gt_font_bold {
  font-weight: bold;
}

#tinncrjaht .gt_font_italic {
  font-style: italic;
}

#tinncrjaht .gt_super {
  font-size: 65%;
}

#tinncrjaht .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#tinncrjaht .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#tinncrjaht .gt_indent_1 {
  text-indent: 5px;
}

#tinncrjaht .gt_indent_2 {
  text-indent: 10px;
}

#tinncrjaht .gt_indent_3 {
  text-indent: 15px;
}

#tinncrjaht .gt_indent_4 {
  text-indent: 20px;
}

#tinncrjaht .gt_indent_5 {
  text-indent: 25px;
}

#tinncrjaht .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#tinncrjaht div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings gt_spanner_row">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1" scope="col" id="label"><span data-qmd-base64="KipQYXRpZW50IGNoYXJhY3RlcmlzdGljKio="><span class='gt_from_md'><strong>Patient characteristic</strong></span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="2" colspan="1" scope="col" id="stat_0"><span data-qmd-base64="KipUb3RhbCoqIDxicj4gTiA9IDM1MDA="><span class='gt_from_md'><strong>Total</strong> <br> N = 3500</span></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2" scope="colgroup" id="level 1; stat_1">
        <div class="gt_column_spanner"><span data-qmd-base64="KipUcmVhdG1lbnQgcmVjZWl2ZWQqKg=="><span class='gt_from_md'><strong>Treatment received</strong></span></span></div>
      </th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="2" colspan="1" scope="col" id="estimate"><span data-qmd-base64="KipEaWZmZXJlbmNlKio="><span class='gt_from_md'><strong>Difference</strong></span></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>2</sup></span></th>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="stat_1"><span data-qmd-base64="KiowKiogPGJyPiBOID0gMTczMC42OCA8YnI+ICg0OS40JSk="><span class='gt_from_md'><strong>0</strong> <br> N = 1730.68 <br> (49.4%)</span></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="stat_2"><span data-qmd-base64="KioxKiogPGJyPiBOID0gMTc2OS4zMiA8YnI+ICg1MC42JSk="><span class='gt_from_md'><strong>1</strong> <br> N = 1769.32 <br> (50.6%)</span></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="label" class="gt_row gt_left">dem_age_index_cont</td>
<td headers="stat_0" class="gt_row gt_center">64 (60, 71)</td>
<td headers="stat_1" class="gt_row gt_center">65 (61, 71)</td>
<td headers="stat_2" class="gt_row gt_center">64 (60, 72)</td>
<td headers="estimate" class="gt_row gt_center">0.00</td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_age_lt65</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.02</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    &lt;65</td>
<td headers="stat_0" class="gt_row gt_center">1,890 (54%)</td>
<td headers="stat_1" class="gt_row gt_center">928 (54%)</td>
<td headers="stat_2" class="gt_row gt_center">962 (54%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    65+</td>
<td headers="stat_0" class="gt_row gt_center">1,610 (46%)</td>
<td headers="stat_1" class="gt_row gt_center">803 (46%)</td>
<td headers="stat_2" class="gt_row gt_center">807 (46%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_sex_cont</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.05</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Female</td>
<td headers="stat_0" class="gt_row gt_center">2,205 (63%)</td>
<td headers="stat_1" class="gt_row gt_center">1,068 (62%)</td>
<td headers="stat_2" class="gt_row gt_center">1,137 (64%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Male</td>
<td headers="stat_0" class="gt_row gt_center">1,295 (37%)</td>
<td headers="stat_1" class="gt_row gt_center">663 (38%)</td>
<td headers="stat_2" class="gt_row gt_center">632 (36%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_race</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.01</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Asian</td>
<td headers="stat_0" class="gt_row gt_center">2,170 (62%)</td>
<td headers="stat_1" class="gt_row gt_center">1,069 (62%)</td>
<td headers="stat_2" class="gt_row gt_center">1,101 (62%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Non-Asian</td>
<td headers="stat_0" class="gt_row gt_center">1,330 (38%)</td>
<td headers="stat_1" class="gt_row gt_center">661 (38%)</td>
<td headers="stat_2" class="gt_row gt_center">669 (38%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_smoking_history</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.05</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Current/former</td>
<td headers="stat_0" class="gt_row gt_center">1,225 (35%)</td>
<td headers="stat_1" class="gt_row gt_center">626 (36%)</td>
<td headers="stat_2" class="gt_row gt_center">599 (34%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Never</td>
<td headers="stat_0" class="gt_row gt_center">2,275 (65%)</td>
<td headers="stat_1" class="gt_row gt_center">1,105 (64%)</td>
<td headers="stat_2" class="gt_row gt_center">1,170 (66%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_ecog_cont</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.03</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    0</td>
<td headers="stat_0" class="gt_row gt_center">1,435 (41%)</td>
<td headers="stat_1" class="gt_row gt_center">723 (42%)</td>
<td headers="stat_2" class="gt_row gt_center">712 (40%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    1</td>
<td headers="stat_0" class="gt_row gt_center">2,065 (59%)</td>
<td headers="stat_1" class="gt_row gt_center">1,007 (58%)</td>
<td headers="stat_2" class="gt_row gt_center">1,058 (60%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
  </tbody>
  <tfoot class="gt_sourcenotes">
    <tr>
      <td class="gt_sourcenote" colspan="5"><span data-qmd-base64="QWJicmV2aWF0aW9uOiBDSSA9IENvbmZpZGVuY2UgSW50ZXJ2YWw="><span class='gt_from_md'>Abbreviation: CI = Confidence Interval</span></span></td>
    </tr>
  </tfoot>
  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="5"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span> <span data-qmd-base64="TWVkaWFuIChRMSwgUTMpOyBuICglKQ=="><span class='gt_from_md'>Median (Q1, Q3); n (%)</span></span></td>
    </tr>
    <tr>
      <td class="gt_footnote" colspan="5"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>2</sup></span> <span data-qmd-base64="U3RhbmRhcmRpemVkIE1lYW4gRGlmZmVyZW5jZQ=="><span class='gt_from_md'>Standardized Mean Difference</span></span></td>
    </tr>
  </tfoot>
</table>
</div>
```

:::
:::
:::

:::

## Estimating final treatment effects

In the final step,  the `cox_pooling()` function that comes with encore.analytics is used to fit and pool the results of a Cox proportional hazards model on the imputed, matched, and raked datasets. The function takes list of imputed data frames with weights and cluster (matching) information (= our matched/smr-weighted and raked datasets) and a formula for the Cox proportional hazards model. The data frames must have column names weights and subclass (for matched datasets to indicate the cluster membership).

The `cox_pooling()` function is a convenience wrapper for fitting Cox models to multiple imputed datasets which do not come as a mimids or wimids object. This is useful when there are intermediate steps in the analysis pipeline, such as computing raking weights via `raking_weights()` which are so far not straightforward to implement in a mimids or wimids object.

The function follows the following logic:

1. Fit a Cox proportional hazards model to each imputed dataset. If a \code{subclass} column is present, it is used as a cluster variable for matched pairs.
2. Convert the list of fitted models into a \code{mira} object using \code{mice::as.mira}.
3. Pool the results using \code{mice::pool}, which applies Rubin's rules for combining estimates and variances across imputations.
4. Format the pooled results, including exponentiating the hazard ratios and calculating confidence intervals.

Our fitted Cox proportional hazards model is specified as follows:

::: {.cell}

```{.r .cell-code}
# fit a survival model
cox_fit <- as.formula(survival::Surv(fu_itt_months, death_itt) ~ treat)
cox_fit
```

::: {.cell-output .cell-output-stdout}

```
survival::Surv(fu_itt_months, death_itt) ~ treat
```


:::
:::

### Estimate treatment effects on matched and raked list of datasets

::: {.cell}

```{.r .cell-code}
# fit and pool Cox proportional hazards model results
cox_pooling(mirwds, surv_formula = cox_fit)
```

::: {.cell-output .cell-output-stdout}

```
   term  estimate  std.error statistic      p.value  conf.low conf.high
1 treat 0.7398062 0.05166569 -5.833021 1.060081e-08 0.6683716 0.8188757
             b      df dfcom       fmi    lambda  m       riv        ubar
1 0.0003134966 437.178  2654 0.1331443 0.1291877 10 0.1483531 0.002324497
```


:::
:::

### Estimate treatment effects on SMR-weighted and raked list of datasets

::: {.cell}

```{.r .cell-code}
# fit and pool Cox proportional hazards model results
cox_pooling(wirwds, surv_formula = cox_fit)
```

::: {.cell-output .cell-output-stdout}

```
   term  estimate std.error statistic      p.value  conf.low conf.high
1 treat 0.7438464 0.0454113 -6.516456 1.089725e-10 0.6804354 0.8131667
             b       df dfcom        fmi     lambda  m        riv        ubar
1 0.0001364739 1110.119  3459 0.07446312 0.07279717 10 0.07851267 0.001912065
```


:::
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

---------------------------------------------------------
        &nbsp;              package        loadedversion 
---------------------- ------------------ ---------------
      **dplyr**              dplyr             1.1.4     

 **encore.analytics**   encore.analytics       0.1.0     

    **gtsummary**          gtsummary           2.2.0     

       **here**               here             1.0.1     

     **MatchIt**            MatchIt            4.7.1     

    **MatchThem**          MatchThem           1.2.1     

      **Matrix**             Matrix            1.7-0     

       **mice**               mice            3.17.0     

      **survey**             survey            4.4-2     

     **survival**           survival           3.5-8     
---------------------------------------------------------


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
_encore.analytics(v.0.1.0)_, _gtsummary(v.2.2.0)_, _survey(v.4.4-2)_, _Matrix(v.1.7-0)_, _MatchIt(v.4.7.1)_, _MatchThem(v.1.2.1)_, _mice(v.3.17.0)_, _survival(v.3.5-8)_, _dplyr(v.1.1.4)_ and _here(v.1.0.1)_

**loaded via a namespace (and not attached):** 
_Rdpack(v.2.6.4)_, _DBI(v.1.2.3)_, _gridExtra(v.2.3)_, _sandwich(v.3.1-1)_, _rlang(v.1.1.6)_, _magrittr(v.2.0.3)_, _furrr(v.0.3.1)_, _compiler(v.4.4.0)_, _gdata(v.3.0.1)_, _vctrs(v.0.6.5)_, _stringr(v.1.5.1)_, _pkgconfig(v.2.0.3)_, _shape(v.1.4.6.1)_, _crayon(v.1.5.3)_, _fastmap(v.1.2.0)_, _backports(v.1.5.0)_, _pander(v.0.6.5)_, _rmarkdown(v.2.29)_, _sessioninfo(v.1.2.2)_, _markdown(v.2.0)_, _nloptr(v.2.2.1)_, _purrr(v.1.0.4)_, _xfun(v.0.52)_, _glmnet(v.4.1-8)_, _jomo(v.2.7-6)_, _litedown(v.0.7)_, _anesrake(v.0.80)_, _jsonlite(v.2.0.0)_, _tictoc(v.1.2.1)_, _chk(v.0.10.0)_, _pan(v.1.9)_, _broom(v.1.0.8)_, _parallel(v.4.4.0)_, _cluster(v.2.1.6)_, _R6(v.2.6.1)_, _simsurv(v.1.0.0)_, _stringi(v.1.8.7)_, _RColorBrewer(v.1.1-3)_, _parallelly(v.1.38.0)_, _boot(v.1.3-30)_, _rpart(v.4.1.23)_, _Rcpp(v.1.0.14)_, _assertthat(v.0.2.1)_, _iterators(v.1.0.14)_, _knitr(v.1.50)_, _zoo(v.1.8-14)_, _base64enc(v.0.1-3)_, _weights(v.1.0.4)_, _splines(v.4.4.0)_, _nnet(v.7.3-19)_, _tidyselect(v.1.2.1)_, _rstudioapi(v.0.17.1)_, _yaml(v.2.3.10)_, _codetools(v.0.2-20)_, _listenv(v.0.9.1)_, _lattice(v.0.22-6)_, _tibble(v.3.2.1)_, _withr(v.3.0.2)_, _evaluate(v.1.0.3)_, _foreign(v.0.8-86)_, _future(v.1.34.0)_, _xml2(v.1.3.8)_, _pillar(v.1.10.2)_, _WeightIt(v.1.4.0)_, _checkmate(v.2.3.2)_, _renv(v.1.0.7)_, _foreach(v.1.5.2)_, _reformulas(v.0.4.1)_, _generics(v.0.1.4)_, _rprojroot(v.2.0.4)_, _ggplot2(v.3.5.2)_, _commonmark(v.1.9.5)_, _scales(v.1.4.0)_, _minqa(v.1.2.8)_, _globals(v.0.16.3)_, _gtools(v.3.9.5)_, _glue(v.1.8.0)_, _Hmisc(v.5.2-3)_, _tools(v.4.4.0)_, _data.table(v.1.17.2)_, _lme4(v.1.1-37)_, _locfit(v.1.5-9.12)_, _tidyr(v.1.3.1)_, _mitools(v.2.4)_, _rbibutils(v.2.3)_, _cards(v.0.6.0)_, _colorspace(v.2.1-1)_, _nlme(v.3.1-164)_, _cardx(v.0.2.4)_, _htmlTable(v.2.4.3)_, _Formula(v.1.2-5)_, _cli(v.3.6.5)_, _smd(v.0.8.0)_, _gt(v.1.0.0)_, _gtable(v.0.3.6)_, _sass(v.0.4.10)_, _digest(v.0.6.37)_, _htmlwidgets(v.1.6.4)_, _farver(v.2.1.2)_, _htmltools(v.0.5.8.1)_, _lifecycle(v.1.0.4)_, _mitml(v.0.4-5)_ and _MASS(v.7.3-60.2)_
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


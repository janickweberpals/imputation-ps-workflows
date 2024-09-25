---
title: "Re-weighting to a target population"
author: Janick Weberpals, RPh, PhD
date: last-modified
format: html
toc: true
toc-depth: 3
code-tools: true
code-fold: true
keep-md: true
embed-resources: true
editor: visual
---



This is a reproducible example on how to incorporate population weights to match distributions of a target population in multiple imputation \> matching/weighting \> balance assessment \> outcome analysis workflows.

Load packages:


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
library(encore.io)

# calling functions
source(here::here("functions", "source_encore.io_functions.R"))

# track time
runtime <- tictoc::tic()
```
:::


## Data generation

We use the `simulate_flaura()` function to simulate a realistic oncology comparative effectiveness cohort analytic dataset.


::: {.cell}

```{.r .cell-code}
# load example dataset with missing observations
data_miss <- simulate_flaura(
  n_total = 3500, 
  seed = 41, 
  include_id = FALSE, 
  imposeNA = TRUE,
  propNA = .33
  ) |> 
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

Multiple imputation using `mice:`


::: {.cell}

```{.r .cell-code}
# impute data
data_imp <- futuremice(
  parallelseed = 42,
  n.core = 7,
  data = data_miss,
  method = "rf",
  m = 10,
  print = FALSE
  )
```
:::


## Defining target distributions

Before applying the re-weighting, we need to define the target distributions of patient characteristics that we want to match from the clinical trial using the *raking* procedure. The following distributions are taken from Table 1 of the FLAURA trial.

![FLAURA trial Table 1; in OS analysis race was simplified to Asian vs. non-Asian](/images/nejm_tbl1.png){#tbl-FLAURAtbl1 fig-align="center"}


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


## Propensity score matching and re-weighting

In this step, propensity score matching and re-weighting of key patient characteristics to match those of the original RCT is performed across all imputed datasets.

The propensity score model is specified as follows:


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


The matching and re-weighting is performed using the `re_weight()` function. This function is a wrapper for `matchit()` and `weightit()` in combination with the `anesrake()` function which performs the *raking* (= re-weighting) procedure.

We apply this function to each imputed dataset. Before doing so, the imputed datasets, which are currently stored as a `mids` object, needs to be converted to a list of dataframes:


::: {.cell}

```{.r .cell-code}
# create a mild object containing lists of data.frames
data_mild <- mice::complete(data = data_imp, action = "all", include = FALSE)

summary(data_mild)
```

::: {.cell-output .cell-output-stdout}

```
   Length Class      Mode
1  40     data.frame list
2  40     data.frame list
3  40     data.frame list
4  40     data.frame list
5  40     data.frame list
6  40     data.frame list
7  40     data.frame list
8  40     data.frame list
9  40     data.frame list
10 40     data.frame list
```


:::
:::


The lapply function loops the function through each dataframe and returns a list of `matchit` objects which contain imputed \> matched \> re-weighted datasets. To take advantage of the features that come with the `cobalt` and `matchthem` packages, the function stores the raking weights as sampling weights (`s.weights`).


::: {.cell}

```{.r .cell-code}
# call match re-weight
matchit_out_list <- lapply(
  # list of dataframes
  X = data_mild, 
  # call function
  FUN = re_weight,
  # target distributions
  targets = targets,
  # should matching or weighting be performed
  matching_weighting = "matching",
  # matching arguments passed on to matchit() function
  formula = ps_form,
  ratio = 1,
  method = "nearest",
  distance = "glm",
  link = "logit",
  caliper = 0.05,
  replace = F
  )
```

::: {.cell-output .cell-output-stdout}

```
[1] "Raking converged in 9 iterations"
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] "Raking converged in 12 iterations"
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] "Raking converged in 14 iterations"
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] "Raking converged in 9 iterations"
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] "Raking converged in 10 iterations"
[1] "Raking converged in 8 iterations"
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] "Raking converged in 9 iterations"
[1] "Raking converged in 10 iterations"
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] "Raking converged in 14 iterations"
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] "Raking converged in 10 iterations"
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] "Raking converged in 11 iterations"
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] "Raking converged in 9 iterations"
[1] "Raking converged in 8 iterations"
```


:::
:::


We can inspect the output of the first imputed \> matched \> re-weighted dataset.


::: {.cell}

```{.r .cell-code}
matchit_out_list[[1]]
```

::: {.cell-output .cell-output-stdout}

```
A matchit object
 - method: 1:1 nearest neighbor matching without replacement
 - distance: Propensity score [caliper]
             - estimated with logistic regression
             - sampling weights not included in estimation
 - caliper: <distance> (0.003)
 - number of obs.: 3500 (original), 2970 (matched)
 - sampling weights: present
 - target estimand: ATT
 - covariates: dem_age_index_cont, dem_sex_cont, c_smoking_history, c_number_met_sites, c_hemoglobin_g_dl_cont, c_urea_nitrogen_mg_dl_cont, c_platelets_10_9_l_cont, c_calcium_mg_dl_cont, c_glucose_mg_dl_cont, c_lymphocyte_leukocyte_ratio_cont, c_alp_u_l_cont, c_protein_g_l_cont, c_alt_u_l_cont, c_albumin_g_l_cont, c_bilirubin_mg_dl_cont, c_chloride_mmol_l_cont, c_monocytes_10_9_l_cont, c_eosinophils_leukocytes_ratio_cont, c_ldh_u_l_cont, c_hr_cont, c_sbp_cont, c_oxygen_cont, c_ecog_cont, c_neutrophil_lymphocyte_ratio_cont, c_bmi_cont, c_ast_alt_ratio_cont, c_stage_initial_dx_cont, dem_race, dem_region, dem_ses, c_time_dx_to_index
```


:::
:::


## Table 1

To check if the re-weighting process worked, we can extract the matched patients and compare a Table 1 that does not include the weights vs. a Table that considers the weights. For this example, we look at the first imputed \> matched \> re-weighted dataset.


::: {.cell}

```{.r .cell-code}
# extract the matched of
first_dataset <- get_matches(
  object = matchit_out_list[[1]]
  )
```
:::


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
### Unweighted Table 1


::: {#tbl-tbl1-unweighted .cell tbl-cap='Table 1 BEFORE re-weighting'}

```{.r .cell-code}
library(cardx)
library(smd)

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
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KipQYXRpZW50IGNoYXJhY3RlcmlzdGljKio=&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;Patient characteristic&lt;/strong&gt;&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;"><div data-qmd-base64="KipQYXRpZW50IGNoYXJhY3RlcmlzdGljKio="><div class='gt_from_md'><p><strong>Patient characteristic</strong></p>
</div></div></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="2" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KipUb3RhbCoqIDxicj4gTiA9IDI5NzA=&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;Total&lt;/strong&gt; &lt;br&gt; N = 2970&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><div data-qmd-base64="KipUb3RhbCoqIDxicj4gTiA9IDI5NzA="><div class='gt_from_md'><p><strong>Total</strong> <br> N = 2970</p>
</div></div><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span></th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2" scope="colgroup" id="&lt;div data-qmd-base64=&quot;KipUcmVhdG1lbnQgcmVjZWl2ZWQqKg==&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;Treatment received&lt;/strong&gt;&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;">
        <span class="gt_column_spanner"><div data-qmd-base64="KipUcmVhdG1lbnQgcmVjZWl2ZWQqKg=="><div class='gt_from_md'><p><strong>Treatment received</strong></p>
</div></div></span>
      </th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="2" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KipEaWZmZXJlbmNlKio=&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;Difference&lt;/strong&gt;&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;2&lt;/sup&gt;&lt;/span&gt;"><div data-qmd-base64="KipEaWZmZXJlbmNlKio="><div class='gt_from_md'><p><strong>Difference</strong></p>
</div></div><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>2</sup></span></th>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KiowKiogPGJyPiBOID0gMTQ4NSA8YnI+ICg1MC4wJSk=&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;0&lt;/strong&gt; &lt;br&gt; N = 1485 &lt;br&gt; (50.0%)&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><div data-qmd-base64="KiowKiogPGJyPiBOID0gMTQ4NSA8YnI+ICg1MC4wJSk="><div class='gt_from_md'><p><strong>0</strong> <br> N = 1485 <br> (50.0%)</p>
</div></div><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KioxKiogPGJyPiBOID0gMTQ4NSA8YnI+ICg1MC4wJSk=&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;1&lt;/strong&gt; &lt;br&gt; N = 1485 &lt;br&gt; (50.0%)&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><div data-qmd-base64="KioxKiogPGJyPiBOID0gMTQ4NSA8YnI+ICg1MC4wJSk="><div class='gt_from_md'><p><strong>1</strong> <br> N = 1485 <br> (50.0%)</p>
</div></div><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="label" class="gt_row gt_left">dem_age_index_cont</td>
<td headers="stat_0" class="gt_row gt_center">69 (63, 74)</td>
<td headers="stat_1" class="gt_row gt_center">69 (63, 74)</td>
<td headers="stat_2" class="gt_row gt_center">69 (64, 74)</td>
<td headers="estimate" class="gt_row gt_center">-0.03</td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_age_lt65</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.01</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    &lt;65</td>
<td headers="stat_0" class="gt_row gt_center">942 (32%)</td>
<td headers="stat_1" class="gt_row gt_center">468 (32%)</td>
<td headers="stat_2" class="gt_row gt_center">474 (32%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    65+</td>
<td headers="stat_0" class="gt_row gt_center">2,028 (68%)</td>
<td headers="stat_1" class="gt_row gt_center">1,017 (68%)</td>
<td headers="stat_2" class="gt_row gt_center">1,011 (68%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_sex_cont</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.00</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Female</td>
<td headers="stat_0" class="gt_row gt_center">1,957 (66%)</td>
<td headers="stat_1" class="gt_row gt_center">978 (66%)</td>
<td headers="stat_2" class="gt_row gt_center">979 (66%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Male</td>
<td headers="stat_0" class="gt_row gt_center">1,013 (34%)</td>
<td headers="stat_1" class="gt_row gt_center">507 (34%)</td>
<td headers="stat_2" class="gt_row gt_center">506 (34%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_race</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.01</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Asian</td>
<td headers="stat_0" class="gt_row gt_center">1,146 (39%)</td>
<td headers="stat_1" class="gt_row gt_center">568 (38%)</td>
<td headers="stat_2" class="gt_row gt_center">578 (39%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Non-Asian</td>
<td headers="stat_0" class="gt_row gt_center">1,824 (61%)</td>
<td headers="stat_1" class="gt_row gt_center">917 (62%)</td>
<td headers="stat_2" class="gt_row gt_center">907 (61%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_smoking_history</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.03</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Current/former</td>
<td headers="stat_0" class="gt_row gt_center">1,409 (47%)</td>
<td headers="stat_1" class="gt_row gt_center">714 (48%)</td>
<td headers="stat_2" class="gt_row gt_center">695 (47%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Never</td>
<td headers="stat_0" class="gt_row gt_center">1,561 (53%)</td>
<td headers="stat_1" class="gt_row gt_center">771 (52%)</td>
<td headers="stat_2" class="gt_row gt_center">790 (53%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_ecog_cont</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.02</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    0</td>
<td headers="stat_0" class="gt_row gt_center">1,308 (44%)</td>
<td headers="stat_1" class="gt_row gt_center">660 (44%)</td>
<td headers="stat_2" class="gt_row gt_center">648 (44%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    1</td>
<td headers="stat_0" class="gt_row gt_center">1,662 (56%)</td>
<td headers="stat_1" class="gt_row gt_center">825 (56%)</td>
<td headers="stat_2" class="gt_row gt_center">837 (56%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
  </tbody>
  
  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="5"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span> <div data-qmd-base64="TWVkaWFuIChRMSwgUTMpOyBuICglKQ=="><div class='gt_from_md'><p>Median (Q1, Q3); n (%)</p>
</div></div></td>
    </tr>
    <tr>
      <td class="gt_footnote" colspan="5"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>2</sup></span> <div data-qmd-base64="U3RhbmRhcmRpemVkIE1lYW4gRGlmZmVyZW5jZQ=="><div class='gt_from_md'><p>Standardized Mean Difference</p>
</div></div></td>
    </tr>
  </tfoot>
</table>
</div>
```


:::
:::


### Weighted Table 1


::: {#tbl-tbl1-weighted .cell tbl-cap='Table 1 AFTER re-weighting'}

```{.r .cell-code}
# create survey object 
data_svy <- svydesign(ids = ~ 1, weights = ~ weights, data = first_dataset)

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
<div id="jheaqgdmji" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#jheaqgdmji table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#jheaqgdmji thead, #jheaqgdmji tbody, #jheaqgdmji tfoot, #jheaqgdmji tr, #jheaqgdmji td, #jheaqgdmji th {
  border-style: none;
}

#jheaqgdmji p {
  margin: 0;
  padding: 0;
}

#jheaqgdmji .gt_table {
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

#jheaqgdmji .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#jheaqgdmji .gt_title {
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

#jheaqgdmji .gt_subtitle {
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

#jheaqgdmji .gt_heading {
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

#jheaqgdmji .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#jheaqgdmji .gt_col_headings {
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

#jheaqgdmji .gt_col_heading {
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

#jheaqgdmji .gt_column_spanner_outer {
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

#jheaqgdmji .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#jheaqgdmji .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#jheaqgdmji .gt_column_spanner {
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

#jheaqgdmji .gt_spanner_row {
  border-bottom-style: hidden;
}

#jheaqgdmji .gt_group_heading {
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

#jheaqgdmji .gt_empty_group_heading {
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

#jheaqgdmji .gt_from_md > :first-child {
  margin-top: 0;
}

#jheaqgdmji .gt_from_md > :last-child {
  margin-bottom: 0;
}

#jheaqgdmji .gt_row {
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

#jheaqgdmji .gt_stub {
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

#jheaqgdmji .gt_stub_row_group {
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

#jheaqgdmji .gt_row_group_first td {
  border-top-width: 2px;
}

#jheaqgdmji .gt_row_group_first th {
  border-top-width: 2px;
}

#jheaqgdmji .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#jheaqgdmji .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#jheaqgdmji .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#jheaqgdmji .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#jheaqgdmji .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#jheaqgdmji .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#jheaqgdmji .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#jheaqgdmji .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#jheaqgdmji .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#jheaqgdmji .gt_footnotes {
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

#jheaqgdmji .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#jheaqgdmji .gt_sourcenotes {
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

#jheaqgdmji .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#jheaqgdmji .gt_left {
  text-align: left;
}

#jheaqgdmji .gt_center {
  text-align: center;
}

#jheaqgdmji .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#jheaqgdmji .gt_font_normal {
  font-weight: normal;
}

#jheaqgdmji .gt_font_bold {
  font-weight: bold;
}

#jheaqgdmji .gt_font_italic {
  font-style: italic;
}

#jheaqgdmji .gt_super {
  font-size: 65%;
}

#jheaqgdmji .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#jheaqgdmji .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#jheaqgdmji .gt_indent_1 {
  text-indent: 5px;
}

#jheaqgdmji .gt_indent_2 {
  text-indent: 10px;
}

#jheaqgdmji .gt_indent_3 {
  text-indent: 15px;
}

#jheaqgdmji .gt_indent_4 {
  text-indent: 20px;
}

#jheaqgdmji .gt_indent_5 {
  text-indent: 25px;
}

#jheaqgdmji .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#jheaqgdmji div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings gt_spanner_row">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KipQYXRpZW50IGNoYXJhY3RlcmlzdGljKio=&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;Patient characteristic&lt;/strong&gt;&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;"><div data-qmd-base64="KipQYXRpZW50IGNoYXJhY3RlcmlzdGljKio="><div class='gt_from_md'><p><strong>Patient characteristic</strong></p>
</div></div></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="2" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KipUb3RhbCoqIDxicj4gTiA9IDI5NzA=&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;Total&lt;/strong&gt; &lt;br&gt; N = 2970&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><div data-qmd-base64="KipUb3RhbCoqIDxicj4gTiA9IDI5NzA="><div class='gt_from_md'><p><strong>Total</strong> <br> N = 2970</p>
</div></div><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span></th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2" scope="colgroup" id="&lt;div data-qmd-base64=&quot;KipUcmVhdG1lbnQgcmVjZWl2ZWQqKg==&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;Treatment received&lt;/strong&gt;&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;">
        <span class="gt_column_spanner"><div data-qmd-base64="KipUcmVhdG1lbnQgcmVjZWl2ZWQqKg=="><div class='gt_from_md'><p><strong>Treatment received</strong></p>
</div></div></span>
      </th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="2" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KipEaWZmZXJlbmNlKio=&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;Difference&lt;/strong&gt;&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;2&lt;/sup&gt;&lt;/span&gt;"><div data-qmd-base64="KipEaWZmZXJlbmNlKio="><div class='gt_from_md'><p><strong>Difference</strong></p>
</div></div><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>2</sup></span></th>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KiowKiogPGJyPiBOID0gMTQ3MS41MyA8YnI+ICg0OS41JSk=&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;0&lt;/strong&gt; &lt;br&gt; N = 1471.53 &lt;br&gt; (49.5%)&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><div data-qmd-base64="KiowKiogPGJyPiBOID0gMTQ3MS41MyA8YnI+ICg0OS41JSk="><div class='gt_from_md'><p><strong>0</strong> <br> N = 1471.53 <br> (49.5%)</p>
</div></div><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KioxKiogPGJyPiBOID0gMTQ5OC40NyA8YnI+ICg1MC41JSk=&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;1&lt;/strong&gt; &lt;br&gt; N = 1498.47 &lt;br&gt; (50.5%)&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><div data-qmd-base64="KioxKiogPGJyPiBOID0gMTQ5OC40NyA8YnI+ICg1MC41JSk="><div class='gt_from_md'><p><strong>1</strong> <br> N = 1498.47 <br> (50.5%)</p>
</div></div><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="label" class="gt_row gt_left">dem_age_index_cont</td>
<td headers="stat_0" class="gt_row gt_center">64 (61, 72)</td>
<td headers="stat_1" class="gt_row gt_center">64 (61, 71)</td>
<td headers="stat_2" class="gt_row gt_center">64 (61, 72)</td>
<td headers="estimate" class="gt_row gt_center">-0.02</td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_age_lt65</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.02</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    &lt;65</td>
<td headers="stat_0" class="gt_row gt_center">1,604 (54%)</td>
<td headers="stat_1" class="gt_row gt_center">787 (53%)</td>
<td headers="stat_2" class="gt_row gt_center">817 (55%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    65+</td>
<td headers="stat_0" class="gt_row gt_center">1,366 (46%)</td>
<td headers="stat_1" class="gt_row gt_center">684 (47%)</td>
<td headers="stat_2" class="gt_row gt_center">682 (45%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_sex_cont</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.04</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Female</td>
<td headers="stat_0" class="gt_row gt_center">1,871 (63%)</td>
<td headers="stat_1" class="gt_row gt_center">940 (64%)</td>
<td headers="stat_2" class="gt_row gt_center">931 (62%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Male</td>
<td headers="stat_0" class="gt_row gt_center">1,099 (37%)</td>
<td headers="stat_1" class="gt_row gt_center">531 (36%)</td>
<td headers="stat_2" class="gt_row gt_center">567 (38%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_race</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.01</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Asian</td>
<td headers="stat_0" class="gt_row gt_center">1,841 (62%)</td>
<td headers="stat_1" class="gt_row gt_center">909 (62%)</td>
<td headers="stat_2" class="gt_row gt_center">932 (62%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Non-Asian</td>
<td headers="stat_0" class="gt_row gt_center">1,129 (38%)</td>
<td headers="stat_1" class="gt_row gt_center">563 (38%)</td>
<td headers="stat_2" class="gt_row gt_center">566 (38%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_smoking_history</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.01</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Current/former</td>
<td headers="stat_0" class="gt_row gt_center">1,040 (35%)</td>
<td headers="stat_1" class="gt_row gt_center">518 (35%)</td>
<td headers="stat_2" class="gt_row gt_center">522 (35%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Never</td>
<td headers="stat_0" class="gt_row gt_center">1,931 (65%)</td>
<td headers="stat_1" class="gt_row gt_center">954 (65%)</td>
<td headers="stat_2" class="gt_row gt_center">977 (65%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_ecog_cont</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.02</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    0</td>
<td headers="stat_0" class="gt_row gt_center">1,218 (41%)</td>
<td headers="stat_1" class="gt_row gt_center">596 (41%)</td>
<td headers="stat_2" class="gt_row gt_center">621 (41%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    1</td>
<td headers="stat_0" class="gt_row gt_center">1,752 (59%)</td>
<td headers="stat_1" class="gt_row gt_center">875 (59%)</td>
<td headers="stat_2" class="gt_row gt_center">877 (59%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
  </tbody>
  
  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="5"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span> <div data-qmd-base64="TWVkaWFuIChRMSwgUTMpOyBuICglKQ=="><div class='gt_from_md'><p>Median (Q1, Q3); n (%)</p>
</div></div></td>
    </tr>
    <tr>
      <td class="gt_footnote" colspan="5"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>2</sup></span> <div data-qmd-base64="U3RhbmRhcmRpemVkIE1lYW4gRGlmZmVyZW5jZQ=="><div class='gt_from_md'><p>Standardized Mean Difference</p>
</div></div></td>
    </tr>
  </tfoot>
</table>
</div>
```


:::
:::

:::

### Comparison of custom function vs. `matchthem()`

Lastly, we want to make sure that our custom function results in the same matched datasets as the `matchthem()` function which we use in the main analysis - not considering the re-weighting.

For this demonstration, we use the same matching parameters, but without re-weighting after matching in our custom function.

::: panel.tabset
### Custom function

We run again our custom function but with `targets = NULL` to not re-weight any of the included covariates. To convert the returned output of a list of matchit objects into an object of type `mimids` we use the `MatchThem::as.mimids()` function.


::: {.cell}

```{.r .cell-code}
# call match re-weight
set.seed(42)
matchit_out_list <- lapply(
  X = data_mild, 
  FUN = re_weight,
  targets = NULL,
  matching_weighting = "matching",
  formula = ps_form,
  ratio = 1,
  method = "nearest",
  distance = "glm",
  link = "logit",
  caliper = 0.05,
  replace = F
  )
```

::: {.cell-output .cell-output-stderr}

```
No target distributions specified, no re-weighting will be performed.
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```
No target distributions specified, no re-weighting will be performed.
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```
No target distributions specified, no re-weighting will be performed.
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```
No target distributions specified, no re-weighting will be performed.
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```
No target distributions specified, no re-weighting will be performed.
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```
No target distributions specified, no re-weighting will be performed.
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```
No target distributions specified, no re-weighting will be performed.
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```
No target distributions specified, no re-weighting will be performed.
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```
No target distributions specified, no re-weighting will be performed.
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```
No target distributions specified, no re-weighting will be performed.
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

```{.r .cell-code}
# convert the output into a mimids object
mimids_data_from_function <- MatchThem::as.mimids(
  x = matchit_out_list, 
  datasets = data_imp
  )

mimids_data_from_function
```

::: {.cell-output .cell-output-stderr}

```
Printing               | dataset: #1
```


:::

::: {.cell-output .cell-output-stdout}

```
A matchit object
 - method: 1:1 nearest neighbor matching without replacement
 - distance: Propensity score [caliper]
             - estimated with logistic regression
 - caliper: <distance> (0.003)
 - number of obs.: 3500 (original), 2970 (matched)
 - target estimand: ATT
 - covariates: dem_age_index_cont, dem_sex_cont, c_smoking_history, c_number_met_sites, c_hemoglobin_g_dl_cont, c_urea_nitrogen_mg_dl_cont, c_platelets_10_9_l_cont, c_calcium_mg_dl_cont, c_glucose_mg_dl_cont, c_lymphocyte_leukocyte_ratio_cont, c_alp_u_l_cont, c_protein_g_l_cont, c_alt_u_l_cont, c_albumin_g_l_cont, c_bilirubin_mg_dl_cont, c_chloride_mmol_l_cont, c_monocytes_10_9_l_cont, c_eosinophils_leukocytes_ratio_cont, c_ldh_u_l_cont, c_hr_cont, c_sbp_cont, c_oxygen_cont, c_ecog_cont, c_neutrophil_lymphocyte_ratio_cont, c_bmi_cont, c_ast_alt_ratio_cont, c_stage_initial_dx_cont, dem_race, dem_region, dem_ses, c_time_dx_to_index
```


:::
:::


### `matchthem()` function

The following code resembles the code we would use in the main analysis by implementing the generic `matchthem()` function.


::: {.cell}

```{.r .cell-code}
# matching
set.seed(42)
mimids_data <- matchthem(
  datasets = data_imp,
  formula = ps_form,
  ratio = 1,
  method = "nearest",
  distance = "glm",
  link = "logit",
  caliper = 0.05,
  replace = F
  )
```

::: {.cell-output .cell-output-stderr}

```

Matching Observations  | dataset: #1
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```
#2
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```
#3
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```
#4
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```
#5
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```
#6
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```
#7
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```
#8
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```
#9
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```
#10
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Fewer control units than treated units; not all treated units will get
a match.
```


:::

::: {.cell-output .cell-output-stderr}

```

```


:::

```{.r .cell-code}
mimids_data
```

::: {.cell-output .cell-output-stderr}

```
Printing               | dataset: #1
```


:::

::: {.cell-output .cell-output-stdout}

```
A matchit object
 - method: 1:1 nearest neighbor matching without replacement
 - distance: Propensity score [caliper]
             - estimated with logistic regression
 - caliper: <distance> (0.003)
 - number of obs.: 3500 (original), 2970 (matched)
 - target estimand: ATT
 - covariates: dem_age_index_cont, dem_sex_cont, c_smoking_history, c_number_met_sites, c_hemoglobin_g_dl_cont, c_urea_nitrogen_mg_dl_cont, c_platelets_10_9_l_cont, c_calcium_mg_dl_cont, c_glucose_mg_dl_cont, c_lymphocyte_leukocyte_ratio_cont, c_alp_u_l_cont, c_protein_g_l_cont, c_alt_u_l_cont, c_albumin_g_l_cont, c_bilirubin_mg_dl_cont, c_chloride_mmol_l_cont, c_monocytes_10_9_l_cont, c_eosinophils_leukocytes_ratio_cont, c_ldh_u_l_cont, c_hr_cont, c_sbp_cont, c_oxygen_cont, c_ecog_cont, c_neutrophil_lymphocyte_ratio_cont, c_bmi_cont, c_ast_alt_ratio_cont, c_stage_initial_dx_cont, dem_race, dem_region, dem_ses, c_time_dx_to_index
```


:::
:::


### Comparison of stacked datasets

We can now stack the datasets (= vertically append them) and compare the resulting 10 x 10 datasets for any differences:


::: {.cell}

```{.r .cell-code}
#| echo: true
waldo::compare(
  MatchThem::complete(mimids_data_from_function), 
  MatchThem::complete(mimids_data)
  )
```

::: {.cell-output .cell-output-stdout}

```
✔ No differences
```


:::
:::

:::

## Session info





Script runtime: 0.46 minutes.

::: panel-tabset
### Loaded packages


::: {.cell}

```{.r .cell-code}
pander::pander(subset(data.frame(sessioninfo::package_info()), attached==TRUE, c(package, loadedversion)))
```

::: {.cell-output-display}

-------------------------------------------
    &nbsp;        package    loadedversion 
--------------- ----------- ---------------
   **cardx**       cardx         0.2.1     

   **dplyr**       dplyr         1.1.4     

 **encore.io**   encore.io       0.2.0     

 **gtsummary**   gtsummary       2.0.2     

   **here**        here          1.0.1     

  **MatchIt**     MatchIt        4.5.5     

 **MatchThem**   MatchThem       1.2.1     

  **Matrix**      Matrix         1.7-0     

   **mice**        mice         3.16.0     

    **smd**         smd          0.7.0     

  **survey**      survey         4.4-2     

 **survival**    survival        3.5-8     
-------------------------------------------


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
_smd(v.0.7.0)_, _cardx(v.0.2.1)_, _encore.io(v.0.2.0)_, _gtsummary(v.2.0.2)_, _survey(v.4.4-2)_, _Matrix(v.1.7-0)_, _MatchIt(v.4.5.5)_, _MatchThem(v.1.2.1)_, _mice(v.3.16.0)_, _survival(v.3.5-8)_, _dplyr(v.1.1.4)_ and _here(v.1.0.1)_

**loaded via a namespace (and not attached):** 
_DBI(v.1.2.3)_, _pROC(v.1.18.5)_, _gridExtra(v.2.3)_, _rlang(v.1.1.4)_, _magrittr(v.2.0.3)_, _furrr(v.0.3.1)_, _compiler(v.4.4.0)_, _gdata(v.3.0.0)_, _vctrs(v.0.6.5)_, _stringr(v.1.5.1)_, _pkgconfig(v.2.0.3)_, _shape(v.1.4.6.1)_, _crayon(v.1.5.3)_, _fastmap(v.1.2.0)_, _backports(v.1.5.0)_, _pander(v.0.6.5)_, _utf8(v.1.2.4)_, _rmarkdown(v.2.28)_, _sessioninfo(v.1.2.2)_, _markdown(v.1.13)_, _nloptr(v.2.1.1)_, _waldo(v.0.5.3)_, _purrr(v.1.0.2)_, _bit(v.4.5.0)_, _xfun(v.0.47)_, _glmnet(v.4.1-8)_, _jomo(v.2.7-6)_, _anesrake(v.0.80)_, _jsonlite(v.1.8.9)_, _tictoc(v.1.2.1)_, _chk(v.0.9.2)_, _pan(v.1.9)_, _broom(v.1.0.6)_, _parallel(v.4.4.0)_, _cluster(v.2.1.6)_, _R6(v.2.5.1)_, _simsurv(v.1.0.0)_, _stringi(v.1.8.4)_, _parallelly(v.1.38.0)_, _boot(v.1.3-30)_, _rpart(v.4.1.23)_, _lubridate(v.1.9.3)_, _Rcpp(v.1.0.13)_, _assertthat(v.0.2.1)_, _iterators(v.1.0.14)_, _knitr(v.1.48)_, _base64enc(v.0.1-3)_, _weights(v.1.0.4)_, _splines(v.4.4.0)_, _nnet(v.7.3-19)_, _timechange(v.0.3.0)_, _tidyselect(v.1.2.1)_, _rstudioapi(v.0.16.0)_, _yaml(v.2.3.10)_, _codetools(v.0.2-20)_, _listenv(v.0.9.1)_, _lattice(v.0.22-6)_, _tibble(v.3.2.1)_, _plyr(v.1.8.9)_, _withr(v.3.0.1)_, _evaluate(v.1.0.0)_, _foreign(v.0.8-86)_, _future(v.1.34.0)_, _xml2(v.1.3.6)_, _pillar(v.1.9.0)_, _WeightIt(v.1.3.0)_, _checkmate(v.2.3.2)_, _renv(v.1.0.7)_, _foreach(v.1.5.2)_, _generics(v.0.1.3)_, _rprojroot(v.2.0.4)_, _ggplot2(v.3.5.1)_, _commonmark(v.1.9.1)_, _munsell(v.0.5.1)_, _scales(v.1.3.0)_, _minqa(v.1.2.8)_, _globals(v.0.16.3)_, _gtools(v.3.9.5)_, _glue(v.1.7.0)_, _Hmisc(v.5.1-3)_, _tools(v.4.4.0)_, _data.table(v.1.16.0)_, _lme4(v.1.1-35.5)_, _locfit(v.1.5-9.10)_, _forcats(v.1.0.0)_, _tidyr(v.1.3.1)_, _mitools(v.2.4)_, _cards(v.0.2.2)_, _colorspace(v.2.1-1)_, _nlme(v.3.1-164)_, _htmlTable(v.2.4.3)_, _Formula(v.1.2-5)_, _cli(v.3.6.3)_, _fansi(v.1.0.6)_, _gt(v.0.11.0)_, _arrow(v.17.0.0.1)_, _gtable(v.0.3.5)_, _sass(v.0.4.9)_, _digest(v.0.6.37)_, _htmlwidgets(v.1.6.4)_, _htmltools(v.0.5.8.1)_, _lifecycle(v.1.0.4)_, _mitml(v.0.4-5)_, _bit64(v.4.5.2)_ and _MASS(v.7.3-60.2)_
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

## 

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
  treat_prevalence = .51, 
  seed = 41, 
  include_id = FALSE, 
  imposeNA = TRUE
  ) |> 
  # we have to convert sex and ecog into a factor variable 
  # since anesrake doesn't accept 0/1 numeric encoding for 
  # binary variables
  mutate(across(c(dem_sex_cont, c_ecog_cont), function(x) factor(as.character(x))))

covariates <- data_miss |> 
  select(starts_with("dem_"), starts_with("c_")) |> 
  colnames()
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

## sex ---------------------------------------------------------------------

# female (0) to male (1) proportion:
sex_target <- c(.63, .37) 
names(sex_target) <- c("0", "1")

## race --------------------------------------------------------------------
# asian, non-asian
# asian (TRUE) to non-asian (FALSE) proportion
# note: logical variables in dataframe can be matched to a numeric vector of length 2 and ordered with the TRUE target as the first element and the FALSE target as the second element.
race_target <- c(.62, .38)

## smoking -----------------------------------------------------------------

# current/former smoker (TRUE) to never smoker (FALSE) proportion
# note: logical variables in dataframe can be matched to a numeric vector of length 2 and ordered with the TRUE target as the first element and the FALSE target as the second element.
smoker_target <- c(.35, .65)

## ecog --------------------------------------------------------------------

# ecog 0 by exposure 
avg_prop_ecog0 <- .41

# ecog 0 to ecog 1 proportion
ecog_target <- c(.41, .59)
names(ecog_target) <- c("0", "1")


# summarize target distributions in a named list vector --------------
targets <- list(sex_target, race_target, smoker_target, ecog_target)
names(targets) <- c("dem_sex_cont", "dem_race", "c_smoking_history", "c_ecog_cont")

# print
targets
```

::: {.cell-output .cell-output-stdout}
```
$dem_sex_cont
   0    1 
0.63 0.37 

$dem_race
[1] 0.62 0.38

$c_smoking_history
[1] 0.35 0.65

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
ps_form <- as.formula(paste("treat ~", paste(covariates, collapse = " + ")))
ps_form
```

::: {.cell-output .cell-output-stdout}
```
treat ~ dem_age_index_cont + dem_sex_cont + dem_race + dem_region + 
    dem_ses + c_smoking_history + c_number_met_sites + c_ecog_cont + 
    c_stage_initial_dx_cont + c_hemoglobin_g_dl_cont + c_urea_nitrogen_mg_dl_cont + 
    c_platelets_10_9_l_cont + c_calcium_mg_dl_cont + c_glucose_mg_dl_cont + 
    c_lymphocyte_leukocyte_ratio_cont + c_alp_u_l_cont + c_protein_g_l_cont + 
    c_alt_u_l_cont + c_albumin_g_l_cont + c_bilirubin_mg_dl_cont + 
    c_chloride_mmol_l_cont + c_monocytes_10_9_l_cont + c_eosinophils_leukocytes_ratio_cont + 
    c_ldh_u_l_cont + c_hr_cont + c_sbp_cont + c_oxygen_cont + 
    c_neutrophil_lymphocyte_ratio_cont + c_bmi_cont + c_ast_alt_ratio_cont + 
    c_time_dx_to_index + c_de_novo_mets_dx + c_height_cont + 
    c_weight_cont + c_dbp_cont + c_year_index
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
1  39     data.frame list
2  39     data.frame list
3  39     data.frame list
4  39     data.frame list
5  39     data.frame list
6  39     data.frame list
7  39     data.frame list
8  39     data.frame list
9  39     data.frame list
10 39     data.frame list
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
  link = "linear.logit",
  caliper = 0.05,
  replace = F
  )
```

::: {.cell-output .cell-output-stdout}
```
[1] "Raking converged in 15 iterations"
[1] "Raking converged in 9 iterations"
[1] "Raking converged in 9 iterations"
[1] "Raking converged in 9 iterations"
[1] "Raking converged in 12 iterations"
[1] "Raking converged in 14 iterations"
[1] "Raking converged in 11 iterations"
[1] "Raking converged in 10 iterations"
[1] "Raking converged in 11 iterations"
[1] "Raking converged in 8 iterations"
[1] "Raking converged in 16 iterations"
[1] "Raking converged in 12 iterations"
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
             - estimated with logistic regression and linearized
             - sampling weights not included in estimation
 - caliper: <distance> (0.179)
 - number of obs.: 3500 (original), 366 (matched)
 - sampling weights: present
 - target estimand: ATT
 - covariates: dem_age_index_cont, dem_sex_cont, dem_race, dem_region, dem_ses, c_smoking_history, c_number_met_sites, c_ecog_cont, c_stage_initial_dx_cont, c_hemoglobin_g_dl_cont, c_urea_nitrogen_mg_dl_cont, c_platelets_10_9_l_cont, c_calcium_mg_dl_cont, c_glucose_mg_dl_cont, c_lymphocyte_leukocyte_ratio_cont, c_alp_u_l_cont, c_protein_g_l_cont, c_alt_u_l_cont, c_albumin_g_l_cont, c_bilirubin_mg_dl_cont, c_chloride_mmol_l_cont, c_monocytes_10_9_l_cont, c_eosinophils_leukocytes_ratio_cont, c_ldh_u_l_cont, c_hr_cont, c_sbp_cont, c_oxygen_cont, c_neutrophil_lymphocyte_ratio_cont, c_bmi_cont, c_ast_alt_ratio_cont, c_time_dx_to_index, c_de_novo_mets_dx, c_height_cont, c_weight_cont, c_dbp_cont, c_year_index
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
$dem_sex_cont
   0    1 
0.63 0.37 

$dem_race
[1] 0.62 0.38

$c_smoking_history
[1] 0.35 0.65

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
    include = names(targets)
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
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="2" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KipUb3RhbCoqIDxicj4gTiA9IDM2Ng==&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;Total&lt;/strong&gt; &lt;br&gt; N = 366&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><div data-qmd-base64="KipUb3RhbCoqIDxicj4gTiA9IDM2Ng=="><div class='gt_from_md'><p><strong>Total</strong> <br> N = 366</p>
</div></div><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span></th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2" scope="colgroup" id="&lt;div data-qmd-base64=&quot;KipUcmVhdG1lbnQgcmVjZWl2ZWQqKg==&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;Treatment received&lt;/strong&gt;&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;">
        <span class="gt_column_spanner"><div data-qmd-base64="KipUcmVhdG1lbnQgcmVjZWl2ZWQqKg=="><div class='gt_from_md'><p><strong>Treatment received</strong></p>
</div></div></span>
      </th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="2" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KipEaWZmZXJlbmNlKio=&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;Difference&lt;/strong&gt;&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;2&lt;/sup&gt;&lt;/span&gt;"><div data-qmd-base64="KipEaWZmZXJlbmNlKio="><div class='gt_from_md'><p><strong>Difference</strong></p>
</div></div><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>2</sup></span></th>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KiowKiogPGJyPiBOID0gMTgzIDxicj4gKDUwLjAlKQ==&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;0&lt;/strong&gt; &lt;br&gt; N = 183 &lt;br&gt; (50.0%)&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><div data-qmd-base64="KiowKiogPGJyPiBOID0gMTgzIDxicj4gKDUwLjAlKQ=="><div class='gt_from_md'><p><strong>0</strong> <br> N = 183 <br> (50.0%)</p>
</div></div><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KioxKiogPGJyPiBOID0gMTgzIDxicj4gKDUwLjAlKQ==&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;1&lt;/strong&gt; &lt;br&gt; N = 183 &lt;br&gt; (50.0%)&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><div data-qmd-base64="KioxKiogPGJyPiBOID0gMTgzIDxicj4gKDUwLjAlKQ=="><div class='gt_from_md'><p><strong>1</strong> <br> N = 183 <br> (50.0%)</p>
</div></div><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="label" class="gt_row gt_left">dem_sex_cont</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.06</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    0</td>
<td headers="stat_0" class="gt_row gt_center">247 (67%)</td>
<td headers="stat_1" class="gt_row gt_center">121 (66%)</td>
<td headers="stat_2" class="gt_row gt_center">126 (69%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    1</td>
<td headers="stat_0" class="gt_row gt_center">119 (33%)</td>
<td headers="stat_1" class="gt_row gt_center">62 (34%)</td>
<td headers="stat_2" class="gt_row gt_center">57 (31%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_race</td>
<td headers="stat_0" class="gt_row gt_center">217 (59%)</td>
<td headers="stat_1" class="gt_row gt_center">105 (57%)</td>
<td headers="stat_2" class="gt_row gt_center">112 (61%)</td>
<td headers="estimate" class="gt_row gt_center">-0.08</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_smoking_history</td>
<td headers="stat_0" class="gt_row gt_center">179 (49%)</td>
<td headers="stat_1" class="gt_row gt_center">91 (50%)</td>
<td headers="stat_2" class="gt_row gt_center">88 (48%)</td>
<td headers="estimate" class="gt_row gt_center">0.03</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_ecog_cont</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.15</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    0</td>
<td headers="stat_0" class="gt_row gt_center">160 (44%)</td>
<td headers="stat_1" class="gt_row gt_center">73 (40%)</td>
<td headers="stat_2" class="gt_row gt_center">87 (48%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    1</td>
<td headers="stat_0" class="gt_row gt_center">206 (56%)</td>
<td headers="stat_1" class="gt_row gt_center">110 (60%)</td>
<td headers="stat_2" class="gt_row gt_center">96 (52%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
  </tbody>
  
  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="5"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span> <div data-qmd-base64="biAoJSk="><div class='gt_from_md'><p>n (%)</p>
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
    include = names(targets)
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
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="2" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KipUb3RhbCoqIDxicj4gTiA9IDM2Ng==&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;Total&lt;/strong&gt; &lt;br&gt; N = 366&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><div data-qmd-base64="KipUb3RhbCoqIDxicj4gTiA9IDM2Ng=="><div class='gt_from_md'><p><strong>Total</strong> <br> N = 366</p>
</div></div><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span></th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2" scope="colgroup" id="&lt;div data-qmd-base64=&quot;KipUcmVhdG1lbnQgcmVjZWl2ZWQqKg==&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;Treatment received&lt;/strong&gt;&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;">
        <span class="gt_column_spanner"><div data-qmd-base64="KipUcmVhdG1lbnQgcmVjZWl2ZWQqKg=="><div class='gt_from_md'><p><strong>Treatment received</strong></p>
</div></div></span>
      </th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="2" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KipEaWZmZXJlbmNlKio=&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;Difference&lt;/strong&gt;&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;2&lt;/sup&gt;&lt;/span&gt;"><div data-qmd-base64="KipEaWZmZXJlbmNlKio="><div class='gt_from_md'><p><strong>Difference</strong></p>
</div></div><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>2</sup></span></th>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KiowKiogPGJyPiBOID0gMTgyLjE4IDxicj4gKDQ5LjglKQ==&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;0&lt;/strong&gt; &lt;br&gt; N = 182.18 &lt;br&gt; (49.8%)&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><div data-qmd-base64="KiowKiogPGJyPiBOID0gMTgyLjE4IDxicj4gKDQ5LjglKQ=="><div class='gt_from_md'><p><strong>0</strong> <br> N = 182.18 <br> (49.8%)</p>
</div></div><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KioxKiogPGJyPiBOID0gMTgzLjgyIDxicj4gKDUwLjIlKQ==&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;1&lt;/strong&gt; &lt;br&gt; N = 183.82 &lt;br&gt; (50.2%)&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><div data-qmd-base64="KioxKiogPGJyPiBOID0gMTgzLjgyIDxicj4gKDUwLjIlKQ=="><div class='gt_from_md'><p><strong>1</strong> <br> N = 183.82 <br> (50.2%)</p>
</div></div><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="label" class="gt_row gt_left">dem_sex_cont</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.04</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    0</td>
<td headers="stat_0" class="gt_row gt_center">231 (63%)</td>
<td headers="stat_1" class="gt_row gt_center">113 (62%)</td>
<td headers="stat_2" class="gt_row gt_center">117 (64%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    1</td>
<td headers="stat_0" class="gt_row gt_center">135 (37%)</td>
<td headers="stat_1" class="gt_row gt_center">69 (38%)</td>
<td headers="stat_2" class="gt_row gt_center">66 (36%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_race</td>
<td headers="stat_0" class="gt_row gt_center">227 (62%)</td>
<td headers="stat_1" class="gt_row gt_center">110 (61%)</td>
<td headers="stat_2" class="gt_row gt_center">117 (63%)</td>
<td headers="estimate" class="gt_row gt_center">-0.06</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_smoking_history</td>
<td headers="stat_0" class="gt_row gt_center">128 (35%)</td>
<td headers="stat_1" class="gt_row gt_center">66 (36%)</td>
<td headers="stat_2" class="gt_row gt_center">62 (34%)</td>
<td headers="estimate" class="gt_row gt_center">0.05</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_ecog_cont</td>
<td headers="stat_0" class="gt_row gt_center"><br /></td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td>
<td headers="estimate" class="gt_row gt_center">0.07</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    0</td>
<td headers="stat_0" class="gt_row gt_center">150 (41%)</td>
<td headers="stat_1" class="gt_row gt_center">72 (39%)</td>
<td headers="stat_2" class="gt_row gt_center">78 (43%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    1</td>
<td headers="stat_0" class="gt_row gt_center">216 (59%)</td>
<td headers="stat_1" class="gt_row gt_center">111 (61%)</td>
<td headers="stat_2" class="gt_row gt_center">105 (57%)</td>
<td headers="estimate" class="gt_row gt_center"><br /></td></tr>
  </tbody>
  
  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="5"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span> <div data-qmd-base64="biAoJSk="><div class='gt_from_md'><p>n (%)</p>
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
No target distributions specified, no re-weighting will be performed.
No target distributions specified, no re-weighting will be performed.
No target distributions specified, no re-weighting will be performed.
No target distributions specified, no re-weighting will be performed.
No target distributions specified, no re-weighting will be performed.
No target distributions specified, no re-weighting will be performed.
No target distributions specified, no re-weighting will be performed.
No target distributions specified, no re-weighting will be performed.
No target distributions specified, no re-weighting will be performed.
```
:::

```{.r .cell-code}
# convert the output into a mimids object
data_mimids_from_function <- MatchThem::as.mimids(
  x = matchit_out_list, 
  datasets = data_imp
  )

data_mimids_from_function
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
 - caliper: <distance> (0.023)
 - number of obs.: 3500 (original), 364 (matched)
 - target estimand: ATT
 - covariates: dem_age_index_cont, dem_sex_cont, dem_race, dem_region, dem_ses, c_smoking_history, c_number_met_sites, c_ecog_cont, c_stage_initial_dx_cont, c_hemoglobin_g_dl_cont, c_urea_nitrogen_mg_dl_cont, c_platelets_10_9_l_cont, c_calcium_mg_dl_cont, c_glucose_mg_dl_cont, c_lymphocyte_leukocyte_ratio_cont, c_alp_u_l_cont, c_protein_g_l_cont, c_alt_u_l_cont, c_albumin_g_l_cont, c_bilirubin_mg_dl_cont, c_chloride_mmol_l_cont, c_monocytes_10_9_l_cont, c_eosinophils_leukocytes_ratio_cont, c_ldh_u_l_cont, c_hr_cont, c_sbp_cont, c_oxygen_cont, c_neutrophil_lymphocyte_ratio_cont, c_bmi_cont, c_ast_alt_ratio_cont, c_time_dx_to_index, c_de_novo_mets_dx, c_height_cont, c_weight_cont, c_dbp_cont, c_year_index
```
:::
:::


### `matchthem()` function

The following code resembles the code we would use in the main analysis by implementing the generic `matchthem()` function.


::: {.cell}

```{.r .cell-code}
# matching
set.seed(42)
data_mimids <- matchthem(
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

Matching Observations  | dataset: #1 #2 #3 #4 #5 #6 #7 #8 #9 #10
```
:::

```{.r .cell-code}
data_mimids
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
 - caliper: <distance> (0.023)
 - number of obs.: 3500 (original), 364 (matched)
 - target estimand: ATT
 - covariates: dem_age_index_cont, dem_sex_cont, dem_race, dem_region, dem_ses, c_smoking_history, c_number_met_sites, c_ecog_cont, c_stage_initial_dx_cont, c_hemoglobin_g_dl_cont, c_urea_nitrogen_mg_dl_cont, c_platelets_10_9_l_cont, c_calcium_mg_dl_cont, c_glucose_mg_dl_cont, c_lymphocyte_leukocyte_ratio_cont, c_alp_u_l_cont, c_protein_g_l_cont, c_alt_u_l_cont, c_albumin_g_l_cont, c_bilirubin_mg_dl_cont, c_chloride_mmol_l_cont, c_monocytes_10_9_l_cont, c_eosinophils_leukocytes_ratio_cont, c_ldh_u_l_cont, c_hr_cont, c_sbp_cont, c_oxygen_cont, c_neutrophil_lymphocyte_ratio_cont, c_bmi_cont, c_ast_alt_ratio_cont, c_time_dx_to_index, c_de_novo_mets_dx, c_height_cont, c_weight_cont, c_dbp_cont, c_year_index
```
:::
:::


### Comparison of stacked datasets

We can now stack the datasets (= vertically append them) and compare the resulting 10 x 10 datasets for any differences:


::: {.cell}

```{.r .cell-code}
waldo::compare(
  MatchThem::complete(data_mimids_from_function), 
  MatchThem::complete(data_mimids)
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





Script runtime: 0.41 minutes.

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
   **cardx**       cardx         0.2.0     

   **dplyr**       dplyr         1.1.4     

 **gtsummary**   gtsummary       2.0.1     

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
_smd(v.0.7.0)_, _cardx(v.0.2.0)_, _gtsummary(v.2.0.1)_, _survey(v.4.4-2)_, _Matrix(v.1.7-0)_, _MatchIt(v.4.5.5)_, _MatchThem(v.1.2.1)_, _mice(v.3.16.0)_, _survival(v.3.5-8)_, _dplyr(v.1.1.4)_ and _here(v.1.0.1)_

**loaded via a namespace (and not attached):** 
_tidyselect(v.1.2.1)_, _fastmap(v.1.2.0)_, _digest(v.0.6.37)_, _rpart(v.4.1.23)_, _lifecycle(v.1.0.4)_, _cluster(v.2.1.6)_, _waldo(v.0.5.3)_, _gdata(v.3.0.0)_, _magrittr(v.2.0.3)_, _compiler(v.4.4.0)_, _sass(v.0.4.9)_, _rlang(v.1.1.4)_, _Hmisc(v.5.1-3)_, _tools(v.4.4.0)_, _gt(v.0.11.0)_, _utf8(v.1.2.4)_, _yaml(v.2.3.10)_, _data.table(v.1.16.0)_, _knitr(v.1.48)_, _htmlwidgets(v.1.6.4)_, _xml2(v.1.3.6)_, _withr(v.3.0.1)_, _foreign(v.0.8-86)_, _purrr(v.1.0.2)_, _nnet(v.7.3-19)_, _fansi(v.1.0.6)_, _jomo(v.2.7-6)_, _colorspace(v.2.1-1)_, _future(v.1.34.0)_, _ggplot2(v.3.5.1)_, _gtools(v.3.9.5)_, _globals(v.0.16.3)_, _scales(v.1.3.0)_, _iterators(v.1.0.14)_, _MASS(v.7.3-60.2)_, _cli(v.3.6.3)_, _rmarkdown(v.2.28)_, _crayon(v.1.5.3)_, _generics(v.0.1.3)_, _rstudioapi(v.0.16.0)_, _sessioninfo(v.1.2.2)_, _commonmark(v.1.9.1)_, _minqa(v.1.2.8)_, _DBI(v.1.2.3)_, _pander(v.0.6.5)_, _stringr(v.1.5.1)_, _splines(v.4.4.0)_, _assertthat(v.0.2.1)_, _parallel(v.4.4.0)_, _base64enc(v.0.1-3)_, _mitools(v.2.4)_, _vctrs(v.0.6.5)_, _WeightIt(v.1.3.0)_, _boot(v.1.3-30)_, _glmnet(v.4.1-8)_, _jsonlite(v.1.8.8)_, _mitml(v.0.4-5)_, _Formula(v.1.2-5)_, _htmlTable(v.2.4.3)_, _listenv(v.0.9.1)_, _weights(v.1.0.4)_, _foreach(v.1.5.2)_, _tidyr(v.1.3.1)_, _glue(v.1.7.0)_, _parallelly(v.1.38.0)_, _nloptr(v.2.1.1)_, _pan(v.1.9)_, _chk(v.0.9.2)_, _codetools(v.0.2-20)_, _stringi(v.1.8.4)_, _shape(v.1.4.6.1)_, _gtable(v.0.3.5)_, _lme4(v.1.1-35.5)_, _munsell(v.0.5.1)_, _tibble(v.3.2.1)_, _anesrake(v.0.80)_, _pillar(v.1.9.0)_, _furrr(v.0.3.1)_, _htmltools(v.0.5.8.1)_, _R6(v.2.5.1)_, _rprojroot(v.2.0.4)_, _evaluate(v.0.24.0)_, _lattice(v.0.22-6)_, _markdown(v.1.13)_, _cards(v.0.2.1)_, _backports(v.1.5.0)_, _tictoc(v.1.2.1)_, _broom(v.1.0.6)_, _renv(v.1.0.7)_, _simsurv(v.1.0.0)_, _Rcpp(v.1.0.13)_, _gridExtra(v.2.3)_, _nlme(v.3.1-164)_, _checkmate(v.2.3.2)_, _xfun(v.0.47)_ and _pkgconfig(v.2.0.3)_
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

---
title: "`coxph` versus `svycoxph`"
subtitle: "Comparison of `coxph` versus `svycoxph` after multiple imputation and propensity score matching"
author: Janick Weberpals, RPh, PhD
date: last-modified
format: html
code-fold: false
toc: true
toc-depth: 3
code-tools: true
keep-md: true
editor: visual
---

This is a reproducible example on how to use `coxph` and `svycoxph` in combination with multiple imputation and propensity score matching using a `mimids` object from the MatchThem package.

Load packages:

::: cell
``` {.r .cell-code}
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
```
:::

## Data generation

We use the `simulate_flaura()` function to simulate a realistic oncology comparative effectiveness cohort analytic dataset.

::: cell
``` {.r .cell-code}
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
  colnames()

data_miss |> 
  tbl_summary(
    by = "treat", 
    include = covariates
    )
```

::: cell-output-display
```{=html}
<div id="mhvesddiol" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#mhvesddiol table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#mhvesddiol thead, #mhvesddiol tbody, #mhvesddiol tfoot, #mhvesddiol tr, #mhvesddiol td, #mhvesddiol th {
  border-style: none;
}

#mhvesddiol p {
  margin: 0;
  padding: 0;
}

#mhvesddiol .gt_table {
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

#mhvesddiol .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#mhvesddiol .gt_title {
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

#mhvesddiol .gt_subtitle {
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

#mhvesddiol .gt_heading {
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

#mhvesddiol .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#mhvesddiol .gt_col_headings {
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

#mhvesddiol .gt_col_heading {
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

#mhvesddiol .gt_column_spanner_outer {
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

#mhvesddiol .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#mhvesddiol .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#mhvesddiol .gt_column_spanner {
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

#mhvesddiol .gt_spanner_row {
  border-bottom-style: hidden;
}

#mhvesddiol .gt_group_heading {
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

#mhvesddiol .gt_empty_group_heading {
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

#mhvesddiol .gt_from_md > :first-child {
  margin-top: 0;
}

#mhvesddiol .gt_from_md > :last-child {
  margin-bottom: 0;
}

#mhvesddiol .gt_row {
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

#mhvesddiol .gt_stub {
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

#mhvesddiol .gt_stub_row_group {
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

#mhvesddiol .gt_row_group_first td {
  border-top-width: 2px;
}

#mhvesddiol .gt_row_group_first th {
  border-top-width: 2px;
}

#mhvesddiol .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#mhvesddiol .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#mhvesddiol .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#mhvesddiol .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#mhvesddiol .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#mhvesddiol .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#mhvesddiol .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#mhvesddiol .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#mhvesddiol .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#mhvesddiol .gt_footnotes {
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

#mhvesddiol .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#mhvesddiol .gt_sourcenotes {
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

#mhvesddiol .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#mhvesddiol .gt_left {
  text-align: left;
}

#mhvesddiol .gt_center {
  text-align: center;
}

#mhvesddiol .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#mhvesddiol .gt_font_normal {
  font-weight: normal;
}

#mhvesddiol .gt_font_bold {
  font-weight: bold;
}

#mhvesddiol .gt_font_italic {
  font-style: italic;
}

#mhvesddiol .gt_super {
  font-size: 65%;
}

#mhvesddiol .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#mhvesddiol .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#mhvesddiol .gt_indent_1 {
  text-indent: 5px;
}

#mhvesddiol .gt_indent_2 {
  text-indent: 10px;
}

#mhvesddiol .gt_indent_3 {
  text-indent: 15px;
}

#mhvesddiol .gt_indent_4 {
  text-indent: 20px;
}

#mhvesddiol .gt_indent_5 {
  text-indent: 25px;
}

#mhvesddiol .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#mhvesddiol div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KipDaGFyYWN0ZXJpc3RpYyoq&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;Characteristic&lt;/strong&gt;&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;"><div data-qmd-base64="KipDaGFyYWN0ZXJpc3RpYyoq"><div class='gt_from_md'><p><strong>Characteristic</strong></p>
</div></div></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KiowKiogIApOID0gMSw3MTI=&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;0&lt;/strong&gt;&lt;br /&gt;&#10;N = 1,712&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><div data-qmd-base64="KiowKiogIApOID0gMSw3MTI="><div class='gt_from_md'><p><strong>0</strong><br />
N = 1,712</p>
</div></div><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;div data-qmd-base64=&quot;KioxKiogIApOID0gMSw3ODg=&quot;&gt;&lt;div class='gt_from_md'&gt;&lt;p&gt;&lt;strong&gt;1&lt;/strong&gt;&lt;br /&gt;&#10;N = 1,788&lt;/p&gt;&#10;&lt;/div&gt;&lt;/div&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><div data-qmd-base64="KioxKiogIApOID0gMSw3ODg="><div class='gt_from_md'><p><strong>1</strong><br />
N = 1,788</p>
</div></div><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="label" class="gt_row gt_left">c_smoking_history</td>
<td headers="stat_1" class="gt_row gt_center">579 (51%)</td>
<td headers="stat_2" class="gt_row gt_center">520 (43%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_number_met_sites</td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    1</td>
<td headers="stat_1" class="gt_row gt_center">837 (74%)</td>
<td headers="stat_2" class="gt_row gt_center">899 (75%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    2</td>
<td headers="stat_1" class="gt_row gt_center">249 (22%)</td>
<td headers="stat_2" class="gt_row gt_center">255 (21%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    3</td>
<td headers="stat_1" class="gt_row gt_center">41 (3.6%)</td>
<td headers="stat_2" class="gt_row gt_center">42 (3.5%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    4</td>
<td headers="stat_1" class="gt_row gt_center">7 (0.6%)</td>
<td headers="stat_2" class="gt_row gt_center">8 (0.7%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_ecog_cont</td>
<td headers="stat_1" class="gt_row gt_center">714 (63%)</td>
<td headers="stat_2" class="gt_row gt_center">637 (53%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_stage_initial_dx_cont</td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    1</td>
<td headers="stat_1" class="gt_row gt_center">101 (8.9%)</td>
<td headers="stat_2" class="gt_row gt_center">0 (0%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    2</td>
<td headers="stat_1" class="gt_row gt_center">17 (1.5%)</td>
<td headers="stat_2" class="gt_row gt_center">12 (1.0%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    3</td>
<td headers="stat_1" class="gt_row gt_center">15 (1.3%)</td>
<td headers="stat_2" class="gt_row gt_center">38 (3.2%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    4</td>
<td headers="stat_1" class="gt_row gt_center">1,001 (88%)</td>
<td headers="stat_2" class="gt_row gt_center">1,154 (96%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_hemoglobin_g_dl_cont</td>
<td headers="stat_1" class="gt_row gt_center">12.97 (12.16, 13.72)</td>
<td headers="stat_2" class="gt_row gt_center">12.89 (11.98, 13.75)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_urea_nitrogen_mg_dl_cont</td>
<td headers="stat_1" class="gt_row gt_center">2.76 (2.42, 3.08)</td>
<td headers="stat_2" class="gt_row gt_center">2.77 (2.58, 2.97)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_platelets_10_9_l_cont</td>
<td headers="stat_1" class="gt_row gt_center">255 (218, 290)</td>
<td headers="stat_2" class="gt_row gt_center">266 (230, 302)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_calcium_mg_dl_cont</td>
<td headers="stat_1" class="gt_row gt_center">2.23 (2.21, 2.25)</td>
<td headers="stat_2" class="gt_row gt_center">2.24 (2.21, 2.27)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_glucose_mg_dl_cont</td>
<td headers="stat_1" class="gt_row gt_center">4.64 (4.55, 4.72)</td>
<td headers="stat_2" class="gt_row gt_center">4.65 (4.58, 4.73)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_lymphocyte_leukocyte_ratio_cont</td>
<td headers="stat_1" class="gt_row gt_center">2.94 (2.81, 3.08)</td>
<td headers="stat_2" class="gt_row gt_center">2.93 (2.82, 3.04)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_alp_u_l_cont</td>
<td headers="stat_1" class="gt_row gt_center">4.47 (4.33, 4.61)</td>
<td headers="stat_2" class="gt_row gt_center">4.51 (4.42, 4.59)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_protein_g_l_cont</td>
<td headers="stat_1" class="gt_row gt_center">67.8 (65.1, 70.6)</td>
<td headers="stat_2" class="gt_row gt_center">69.0 (66.0, 72.1)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_alt_u_l_cont</td>
<td headers="stat_1" class="gt_row gt_center">2.93 (2.72, 3.14)</td>
<td headers="stat_2" class="gt_row gt_center">2.86 (2.64, 3.10)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_albumin_g_l_cont</td>
<td headers="stat_1" class="gt_row gt_center">39.01 (36.94, 41.01)</td>
<td headers="stat_2" class="gt_row gt_center">39.98 (37.76, 42.07)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_bilirubin_mg_dl_cont</td>
<td headers="stat_1" class="gt_row gt_center">-0.71 (-1.45, 0.07)</td>
<td headers="stat_2" class="gt_row gt_center">-0.85 (-1.74, -0.03)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_chloride_mmol_l_cont</td>
<td headers="stat_1" class="gt_row gt_center">102.08 (100.08, 104.20)</td>
<td headers="stat_2" class="gt_row gt_center">102.05 (100.20, 104.12)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_monocytes_10_9_l_cont</td>
<td headers="stat_1" class="gt_row gt_center">-0.53 (-0.78, -0.24)</td>
<td headers="stat_2" class="gt_row gt_center">-0.51 (-0.68, -0.35)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_eosinophils_leukocytes_ratio_cont</td>
<td headers="stat_1" class="gt_row gt_center">0.72 (0.50, 0.97)</td>
<td headers="stat_2" class="gt_row gt_center">0.69 (0.28, 1.07)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_ldh_u_l_cont</td>
<td headers="stat_1" class="gt_row gt_center">1.68 (1.65, 1.71)</td>
<td headers="stat_2" class="gt_row gt_center">1.69 (1.65, 1.72)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_hr_cont</td>
<td headers="stat_1" class="gt_row gt_center">4.41 (4.39, 4.43)</td>
<td headers="stat_2" class="gt_row gt_center">4.43 (4.40, 4.46)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_sbp_cont</td>
<td headers="stat_1" class="gt_row gt_center">4.85 (4.77, 4.92)</td>
<td headers="stat_2" class="gt_row gt_center">4.85 (4.79, 4.92)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_oxygen_cont</td>
<td headers="stat_1" class="gt_row gt_center">97.000 (96.986, 97.014)</td>
<td headers="stat_2" class="gt_row gt_center">97.001 (96.993, 97.007)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_neutrophil_lymphocyte_ratio_cont</td>
<td headers="stat_1" class="gt_row gt_center">1.33 (1.11, 1.56)</td>
<td headers="stat_2" class="gt_row gt_center">1.28 (1.02, 1.52)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_bmi_cont</td>
<td headers="stat_1" class="gt_row gt_center">3.23 (3.14, 3.31)</td>
<td headers="stat_2" class="gt_row gt_center">3.23 (3.14, 3.31)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_ast_alt_ratio_cont</td>
<td headers="stat_1" class="gt_row gt_center">0.09 (-0.10, 0.30)</td>
<td headers="stat_2" class="gt_row gt_center">0.12 (-0.07, 0.32)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_time_dx_to_index</td>
<td headers="stat_1" class="gt_row gt_center">73 (43, 103)</td>
<td headers="stat_2" class="gt_row gt_center">43 (32, 55)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_de_novo_mets_dx</td>
<td headers="stat_1" class="gt_row gt_center">774 (68%)</td>
<td headers="stat_2" class="gt_row gt_center">945 (78%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_height_cont</td>
<td headers="stat_1" class="gt_row gt_center">1.64 (1.59, 1.69)</td>
<td headers="stat_2" class="gt_row gt_center">1.65 (1.59, 1.70)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_weight_cont</td>
<td headers="stat_1" class="gt_row gt_center">68 (60, 76)</td>
<td headers="stat_2" class="gt_row gt_center">69 (61, 76)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_dbp_cont</td>
<td headers="stat_1" class="gt_row gt_center">75 (70, 79)</td>
<td headers="stat_2" class="gt_row gt_center">76 (72, 80)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">c_year_index</td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    &lt;2018</td>
<td headers="stat_1" class="gt_row gt_center">87 (5.1%)</td>
<td headers="stat_2" class="gt_row gt_center">1,703 (95%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    2018+</td>
<td headers="stat_1" class="gt_row gt_center">1,625 (95%)</td>
<td headers="stat_2" class="gt_row gt_center">85 (4.8%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_age_index_cont</td>
<td headers="stat_1" class="gt_row gt_center">70 (63, 76)</td>
<td headers="stat_2" class="gt_row gt_center">69 (64, 74)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_sex_cont</td>
<td headers="stat_1" class="gt_row gt_center">614 (36%)</td>
<td headers="stat_2" class="gt_row gt_center">569 (32%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_race</td>
<td headers="stat_1" class="gt_row gt_center">664 (59%)</td>
<td headers="stat_2" class="gt_row gt_center">753 (63%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_region</td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Midwest</td>
<td headers="stat_1" class="gt_row gt_center">123 (11%)</td>
<td headers="stat_2" class="gt_row gt_center">137 (11%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Northeast</td>
<td headers="stat_1" class="gt_row gt_center">382 (34%)</td>
<td headers="stat_2" class="gt_row gt_center">390 (32%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    South</td>
<td headers="stat_1" class="gt_row gt_center">409 (36%)</td>
<td headers="stat_2" class="gt_row gt_center">456 (38%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    West</td>
<td headers="stat_1" class="gt_row gt_center">220 (19%)</td>
<td headers="stat_2" class="gt_row gt_center">221 (18%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
    <tr><td headers="label" class="gt_row gt_left">dem_ses</td>
<td headers="stat_1" class="gt_row gt_center"><br /></td>
<td headers="stat_2" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    1</td>
<td headers="stat_1" class="gt_row gt_center">102 (9.0%)</td>
<td headers="stat_2" class="gt_row gt_center">217 (18%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    2</td>
<td headers="stat_1" class="gt_row gt_center">152 (13%)</td>
<td headers="stat_2" class="gt_row gt_center">157 (13%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    3</td>
<td headers="stat_1" class="gt_row gt_center">302 (27%)</td>
<td headers="stat_2" class="gt_row gt_center">210 (17%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    4</td>
<td headers="stat_1" class="gt_row gt_center">271 (24%)</td>
<td headers="stat_2" class="gt_row gt_center">291 (24%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    5</td>
<td headers="stat_1" class="gt_row gt_center">307 (27%)</td>
<td headers="stat_2" class="gt_row gt_center">329 (27%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">578</td>
<td headers="stat_2" class="gt_row gt_center">584</td></tr>
  </tbody>
  
  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="3"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span> <div data-qmd-base64="biAoJSk7IE1lZGlhbiAoUTEsIFEzKQ=="><div class='gt_from_md'><p>n (%); Median (Q1, Q3)</p>
</div></div></td>
    </tr>
  </tfoot>
</table>
</div>
```
:::
:::

## Multiple imputation

Multiple imputation using `mice:`

::: cell
``` {.r .cell-code}
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

Apply propensity score matching and weighting with replacement within in each imputed dataset:

::: cell
``` {.r .cell-code}
# apply propensity score matching on mids object
ps_form <- as.formula(paste("treat ~", paste(covariates, collapse = " + ")))

# matching
data_mimids <- matchthem(
  formula = ps_form,
  datasets = data_imp,
  approach = 'within',
  method = 'nearest',
  caliper = 0.01,
  ratio = 1,
  replace = F
  )

# SMR weighting
data_wimids <- weightthem(
  formula = ps_form,
  datasets = data_imp,
  approach = 'within',
  method = "glm",
  estimand = "ATT"
  )

# trim extreme weights
data_wimids <- trim(
  x = data_wimids, 
  at = .95, 
  lower = TRUE
  )
```
:::

## Outcome model comparisons

Next, we compare the marginal treatment effect estimates coming from a Cox proportional hazards model after propensity score matching and weighting as implemented in the `coxph()` and in the `svycoxph()` functions.

From the `MatchThem` documentation:

::: callout-important
-   `with()` applies the supplied model in `expr` to the (matched or weighted) multiply imputed datasets, automatically incorporating the (matching) weights when possible. The argument to `expr` should be of the form `glm(y ~ z, family = quasibinomial)`, for example, excluding the data or weights argument, which are automatically supplied.

-   Functions from the **survey** package, such as `svyglm()`, are treated a bit differently. No `svydesign` object needs to be supplied because `with()` automatically constructs and supplies it with the imputed dataset and estimated weights. When `cluster = TRUE` (or `with()` detects that pairs should be clustered; see the `cluster` argument above), pair membership is supplied to the `ids` argument of `svydesign()`.

-   After weighting using `weightthem()`, `glm_weightit()` should be used as the modeling function to fit generalized linear models. It correctly produces robust standard errors that account for estimation of the weights, if possible. See [`WeightIt::glm_weightit()`](http://127.0.0.1:31281/help/library/WeightIt/help/glm_weightit) for details. Otherwise, `svyglm()` should be used rather than `glm()` in order to correctly compute standard errors.

-   **For Cox models, `coxph()` will produce approximately correct standard errors when used with weighting, but `svycoxph()` will produce more accurate standard errors when matching is used.**
:::

::: panel-tabset
### Matching

We now want to compare treatment effect estimates for `treat` when computed (a) using `coxph` (survival package) and (b) `svycoxph` (survey package). More information on estimating treatment effects after matching is provided in <https://kosukeimai.github.io/MatchIt/articles/estimating-effects.html#survival-outcomes>

#### `coxph`

::: cell
``` {.r .cell-code}
# coxph result
coxph_results <- with(
  data = data_mimids,
  expr = coxph(formula = Surv(fu_itt_months, death_itt) ~ treat, 
               weights = weights, 
               cluster = subclass,
               robust = TRUE
               )
  ) |> 
  pool() |> 
  tidy(exponentiate = TRUE, conf.int = TRUE) |> 
  mutate(package = "survival") |> 
  select(package, term, estimate, std.error, conf.low, conf.high) 

coxph_results
```
:::

#### `svycoxph`

::: cell
``` {.r .cell-code}
# svycoxph result
svycoxph_results <- with(
  data = data_mimids,
  expr = svycoxph(formula = Surv(fu_itt_months, death_itt) ~ treat),
  cluster = TRUE
  ) |> 
  pool() |> 
  tidy(exponentiate = TRUE, conf.int = TRUE) |> 
  mutate(package = "survey") |> 
  select(package, term, estimate, std.error, conf.low, conf.high)

svycoxph_results
```
:::

#### Summary

::: cell
``` {.r .cell-code}
rbind(coxph_results, svycoxph_results)
```

::: {.cell-output .cell-output-stdout}
```         
   package  term  estimate std.error  conf.low conf.high
1 survival treat 0.6807406 0.1584256 0.4922441 0.9414186
2   survey treat 0.6807406 0.1586700 0.4934002 0.9392128
```
:::
:::

### Weighting

#### `coxph`

::: cell
``` {.r .cell-code}
# coxph result
coxph_results <- with(
  data = data_wimids,
  expr = coxph(formula = Surv(fu_itt_months, death_itt) ~ treat,
               weights = weights, 
               robust = TRUE
               )
  ) |> 
  pool() |> 
  tidy(exponentiate = TRUE, conf.int = TRUE) |> 
  mutate(package = "survival") |> 
  select(package, term, estimate, std.error, conf.low, conf.high) 

coxph_results
```
:::

#### `svycoxph`

::: cell
``` {.r .cell-code}
# svycoxph result
svycoxph_results <- with(
  data = data_wimids,
  expr = svycoxph(formula = Surv(fu_itt_months, death_itt) ~ treat),
  cluster = TRUE
  ) |> 
  pool() |> 
  tidy(exponentiate = TRUE, conf.int = TRUE) |> 
  mutate(package = "survey") |> 
  select(package, term, estimate, std.error, conf.low, conf.high) 

svycoxph_results
```
:::

#### Summary

::: cell
``` {.r .cell-code}
rbind(coxph_results, svycoxph_results)
```

::: {.cell-output .cell-output-stdout}
```         
   package  term  estimate  std.error  conf.low conf.high
1 survival treat 0.7761114 0.07023749 0.6758481 0.8912490
2   survey treat 0.7761114 0.07024572 0.6758771 0.8912107
```
:::
:::
:::

---
subtitle: "Agreement metrics for comparing RCT and RWE results"
author: Janick Weberpals, RPh, PhD
date: last-modified
format: 
  html:
    css: styles.css
    page-layout: full
toc: true
toc-depth: 3
keep-md: true
embed-resources: true
bibliography: ../references.bib
csl: ../pharmacoepidemiology-and-drug-safety.csl
---

# Agreement Metrics {#sec-agreement-metrics}

In this vignette, we demonstrate how to use the `agreement_metrics()` and `smd_agreement()` functions from the `encore.analytics` package to quantitatively assess the agreement between results from randomized controlled trials (RCTs) and their real-world evidence (RWE) emulations. These metrics are particularly important in the context of trial emulation studies, where we aim to understand how well RWE analyses can reproduce or complement RCT findings.

## Background

The use of real-world evidence for regulatory and clinical decision-making has gained significant attention, particularly with the FDA's Real World Evidence Program [@fda2018]. However, a key question remains: How well do RWE studies align with RCT results when attempting to answer similar clinical questions? Recent research by @wang2023 in the RCT-DUPLICATE initiative has shown that when RWE studies can closely emulate trial design elements, they can achieve high concordance with RCT results, with correlation coefficients as high as 0.93 (95% CI: 0.79-0.97) in well-emulated studies.

## Overview

The `agreement_metrics()` function provides a comprehensive framework to evaluate the concordance between RCT and RWE results using three complementary metrics that have been validated in large-scale emulation studies [@wang2023]:

1.  Statistical significance agreement
2.  Estimate agreement
3.  Standardized mean difference (SMD) agreement

## Methodology

### Types of Agreement

The package implements three complementary approaches to assess agreement:

1.  **Statistical Significance Agreement**: Evaluates whether the RCT and RWE results align in terms of both direction and statistical significance.

2.  **Estimate Agreement**: Determines whether the RWE point estimate falls within the confidence interval of the RCT result.

3.  **SMD Agreement**: Calculates a standardized mean difference that accounts for both the magnitude of difference between estimates and their uncertainty.

### Mathematical Details

The SMD calculation, which is particularly important for quantifying agreement, follows this methodology:

For estimates $\theta_{RCT}$ and $\theta_{RWE}$ with their respective variances, the SMD is calculated as:

$$SMD = \frac{\theta_{RCT} - \theta_{RWE}}{\sqrt{Var(\theta_{RCT}) + Var(\theta_{RWE})}}$$

where variances are derived from confidence intervals assuming normal distribution:

$$Var(\theta) = \left(\frac{upper - lower}{2 \times 1.96}\right)^2$$

The default threshold for SMD agreement is ±1.96, corresponding to α=0.05.

## Example Application

Let's walk through examples of using these functions:

::: {.cell}

```{.r .cell-code}
library(dplyr)
library(tibble)
library(encore.analytics)
```
:::

### Simple Comparison

First, let's look at a simple comparison between one RCT and RWE result:

::: {.cell}

```{.r .cell-code}
# Create example data
x <- tribble(
  ~Analysis, ~rct_estimate, ~rct_lower, ~rct_upper, ~rwe_estimate, ~rwe_lower, ~rwe_upper,
  "Main analysis", 0.87, 0.78, 0.97, 0.82, 0.76, 0.87
  )

# Calculate agreement metrics
agreement_metrics(x, analysis_col = "Analysis")
```

::: {.cell-output-display}

```{=html}
<div id="xgjtdraski" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#xgjtdraski table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#xgjtdraski thead, #xgjtdraski tbody, #xgjtdraski tfoot, #xgjtdraski tr, #xgjtdraski td, #xgjtdraski th {
  border-style: none;
}

#xgjtdraski p {
  margin: 0;
  padding: 0;
}

#xgjtdraski .gt_table {
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

#xgjtdraski .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#xgjtdraski .gt_title {
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

#xgjtdraski .gt_subtitle {
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

#xgjtdraski .gt_heading {
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

#xgjtdraski .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#xgjtdraski .gt_col_headings {
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

#xgjtdraski .gt_col_heading {
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

#xgjtdraski .gt_column_spanner_outer {
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

#xgjtdraski .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#xgjtdraski .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#xgjtdraski .gt_column_spanner {
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

#xgjtdraski .gt_spanner_row {
  border-bottom-style: hidden;
}

#xgjtdraski .gt_group_heading {
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

#xgjtdraski .gt_empty_group_heading {
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

#xgjtdraski .gt_from_md > :first-child {
  margin-top: 0;
}

#xgjtdraski .gt_from_md > :last-child {
  margin-bottom: 0;
}

#xgjtdraski .gt_row {
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

#xgjtdraski .gt_stub {
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

#xgjtdraski .gt_stub_row_group {
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

#xgjtdraski .gt_row_group_first td {
  border-top-width: 2px;
}

#xgjtdraski .gt_row_group_first th {
  border-top-width: 2px;
}

#xgjtdraski .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#xgjtdraski .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#xgjtdraski .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#xgjtdraski .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#xgjtdraski .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#xgjtdraski .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#xgjtdraski .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#xgjtdraski .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#xgjtdraski .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#xgjtdraski .gt_footnotes {
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

#xgjtdraski .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#xgjtdraski .gt_sourcenotes {
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

#xgjtdraski .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#xgjtdraski .gt_left {
  text-align: left;
}

#xgjtdraski .gt_center {
  text-align: center;
}

#xgjtdraski .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#xgjtdraski .gt_font_normal {
  font-weight: normal;
}

#xgjtdraski .gt_font_bold {
  font-weight: bold;
}

#xgjtdraski .gt_font_italic {
  font-style: italic;
}

#xgjtdraski .gt_super {
  font-size: 65%;
}

#xgjtdraski .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#xgjtdraski .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#xgjtdraski .gt_indent_1 {
  text-indent: 5px;
}

#xgjtdraski .gt_indent_2 {
  text-indent: 10px;
}

#xgjtdraski .gt_indent_3 {
  text-indent: 15px;
}

#xgjtdraski .gt_indent_4 {
  text-indent: 20px;
}

#xgjtdraski .gt_indent_5 {
  text-indent: 25px;
}

#xgjtdraski .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#xgjtdraski div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings gt_spanner_row">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1" style="font-weight: bold;" scope="col" id="Analysis">Analysis</th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2" style="font-weight: bold;" scope="colgroup" id="HR (95% CI)">
        <div class="gt_column_spanner">HR (95% CI)</div>
      </th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1" style="font-weight: bold;" scope="col" id="significance_agreement"><span data-qmd-base64="U3RhdGlzdGljYWwgPGJyPiBzaWduaWZpY2FuY2UgPGJyPiBhZ3JlZW1lbnQ="><span class='gt_from_md'>Statistical <br> significance <br> agreement</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1" style="font-weight: bold;" scope="col" id="estimate_agreement"><span data-qmd-base64="RXN0aW1hdGUgPGJyPiBhZ3JlZW1lbnQ="><span class='gt_from_md'>Estimate <br> agreement</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="2" colspan="1" style="font-weight: bold;" scope="col" id="smd_agreement">SMD</th>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" style="font-weight: bold;" scope="col" id="RCT">RCT</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" style="font-weight: bold;" scope="col" id="RWE">RWE</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Analysis" class="gt_row gt_left">Main analysis</td>
<td headers="RCT" class="gt_row gt_center">0.87 (0.78 - 0.97)</td>
<td headers="RWE" class="gt_row gt_center">0.82 (0.76 - 0.87)</td>
<td headers="significance_agreement" class="gt_row gt_left" style="color: #006400;">Yes</td>
<td headers="estimate_agreement" class="gt_row gt_left" style="color: #006400;">Yes</td>
<td headers="smd_agreement" class="gt_row gt_center" style="color: #006400;">Yes (0.90)</td></tr>
  </tbody>
  <tfoot>
    <tr class="gt_footnotes">
      <td class="gt_footnote" colspan="6"> Abbreviations: CI = Confidence interval, HR = Hazard ratio, RCT = Randomized controlled trial, RWE = Real-world evidence, SMD = standardized mean difference (based on log hazard ratios)</td>
    </tr>
  </tfoot>
</table>
</div>
```

:::
:::

### Multi-Database Comparison

Now let's examine agreement across multiple databases:

::: {.cell}

```{.r .cell-code}
# Create multi-database example
x_multi <- tribble(
  ~Analysis, ~Database, ~rct_estimate, ~rct_lower, ~rct_upper, ~rwe_estimate, ~rwe_lower, ~rwe_upper,
  "Main analysis", "Database 1", 0.87, 0.78, 0.97, 0.82, 0.76, 0.87,
  "Main analysis", "Database 2", 0.50, 0.40, 0.60, 2.00, 1.80, 2.20,
  "Main analysis", "Database 3", 0.80, 0.70, 0.90, 1.50, 1.40, 1.60
  )

# Calculate agreement metrics with grouping
agreement_metrics(x_multi, 
                 analysis_col = "Analysis", 
                 group_col = "Database")
```

::: {.cell-output-display}

```{=html}
<div id="wmymirmzbq" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#wmymirmzbq table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#wmymirmzbq thead, #wmymirmzbq tbody, #wmymirmzbq tfoot, #wmymirmzbq tr, #wmymirmzbq td, #wmymirmzbq th {
  border-style: none;
}

#wmymirmzbq p {
  margin: 0;
  padding: 0;
}

#wmymirmzbq .gt_table {
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

#wmymirmzbq .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#wmymirmzbq .gt_title {
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

#wmymirmzbq .gt_subtitle {
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

#wmymirmzbq .gt_heading {
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

#wmymirmzbq .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#wmymirmzbq .gt_col_headings {
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

#wmymirmzbq .gt_col_heading {
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

#wmymirmzbq .gt_column_spanner_outer {
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

#wmymirmzbq .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#wmymirmzbq .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#wmymirmzbq .gt_column_spanner {
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

#wmymirmzbq .gt_spanner_row {
  border-bottom-style: hidden;
}

#wmymirmzbq .gt_group_heading {
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

#wmymirmzbq .gt_empty_group_heading {
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

#wmymirmzbq .gt_from_md > :first-child {
  margin-top: 0;
}

#wmymirmzbq .gt_from_md > :last-child {
  margin-bottom: 0;
}

#wmymirmzbq .gt_row {
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

#wmymirmzbq .gt_stub {
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

#wmymirmzbq .gt_stub_row_group {
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

#wmymirmzbq .gt_row_group_first td {
  border-top-width: 2px;
}

#wmymirmzbq .gt_row_group_first th {
  border-top-width: 2px;
}

#wmymirmzbq .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#wmymirmzbq .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#wmymirmzbq .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#wmymirmzbq .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#wmymirmzbq .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#wmymirmzbq .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#wmymirmzbq .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#wmymirmzbq .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#wmymirmzbq .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#wmymirmzbq .gt_footnotes {
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

#wmymirmzbq .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#wmymirmzbq .gt_sourcenotes {
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

#wmymirmzbq .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#wmymirmzbq .gt_left {
  text-align: left;
}

#wmymirmzbq .gt_center {
  text-align: center;
}

#wmymirmzbq .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#wmymirmzbq .gt_font_normal {
  font-weight: normal;
}

#wmymirmzbq .gt_font_bold {
  font-weight: bold;
}

#wmymirmzbq .gt_font_italic {
  font-style: italic;
}

#wmymirmzbq .gt_super {
  font-size: 65%;
}

#wmymirmzbq .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#wmymirmzbq .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#wmymirmzbq .gt_indent_1 {
  text-indent: 5px;
}

#wmymirmzbq .gt_indent_2 {
  text-indent: 10px;
}

#wmymirmzbq .gt_indent_3 {
  text-indent: 15px;
}

#wmymirmzbq .gt_indent_4 {
  text-indent: 20px;
}

#wmymirmzbq .gt_indent_5 {
  text-indent: 25px;
}

#wmymirmzbq .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#wmymirmzbq div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings gt_spanner_row">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1" style="font-weight: bold;" scope="col" id="Analysis">Analysis</th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2" style="font-weight: bold;" scope="colgroup" id="HR (95% CI)">
        <div class="gt_column_spanner">HR (95% CI)</div>
      </th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1" style="font-weight: bold;" scope="col" id="significance_agreement"><span data-qmd-base64="U3RhdGlzdGljYWwgPGJyPiBzaWduaWZpY2FuY2UgPGJyPiBhZ3JlZW1lbnQ="><span class='gt_from_md'>Statistical <br> significance <br> agreement</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1" style="font-weight: bold;" scope="col" id="estimate_agreement"><span data-qmd-base64="RXN0aW1hdGUgPGJyPiBhZ3JlZW1lbnQ="><span class='gt_from_md'>Estimate <br> agreement</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="2" colspan="1" style="font-weight: bold;" scope="col" id="smd_agreement">SMD</th>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" style="font-weight: bold;" scope="col" id="RCT">RCT</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" style="font-weight: bold;" scope="col" id="RWE">RWE</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr class="gt_group_heading_row">
      <th colspan="6" class="gt_group_heading" style="font-weight: bold;" scope="colgroup" id="Database 1">Database 1</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="Database 1  Analysis" class="gt_row gt_left">Main analysis</td>
<td headers="Database 1  RCT" class="gt_row gt_center">0.87 (0.78 - 0.97)</td>
<td headers="Database 1  RWE" class="gt_row gt_center">0.82 (0.76 - 0.87)</td>
<td headers="Database 1  significance_agreement" class="gt_row gt_left" style="color: #006400;">Yes</td>
<td headers="Database 1  estimate_agreement" class="gt_row gt_left" style="color: #006400;">Yes</td>
<td headers="Database 1  smd_agreement" class="gt_row gt_center" style="color: #006400;">Yes (  0.90)</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="6" class="gt_group_heading" style="font-weight: bold;" scope="colgroup" id="Database 2">Database 2</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="Database 2  Analysis" class="gt_row gt_left">Main analysis</td>
<td headers="Database 2  RCT" class="gt_row gt_center">0.50 (0.40 - 0.60)</td>
<td headers="Database 2  RWE" class="gt_row gt_center">2.00 (1.80 - 2.20)</td>
<td headers="Database 2  significance_agreement" class="gt_row gt_left" style="color: #8B0000;">No</td>
<td headers="Database 2  estimate_agreement" class="gt_row gt_left" style="color: #8B0000;">No</td>
<td headers="Database 2  smd_agreement" class="gt_row gt_center" style="color: #8B0000;">No (-12.01)</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="6" class="gt_group_heading" style="font-weight: bold;" scope="colgroup" id="Database 3">Database 3</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="Database 3  Analysis" class="gt_row gt_left">Main analysis</td>
<td headers="Database 3  RCT" class="gt_row gt_center">0.80 (0.70 - 0.90)</td>
<td headers="Database 3  RWE" class="gt_row gt_center">1.50 (1.40 - 1.60)</td>
<td headers="Database 3  significance_agreement" class="gt_row gt_left" style="color: #8B0000;">No</td>
<td headers="Database 3  estimate_agreement" class="gt_row gt_left" style="color: #8B0000;">No</td>
<td headers="Database 3  smd_agreement" class="gt_row gt_center" style="color: #8B0000;">No ( -8.66)</td></tr>
  </tbody>
  <tfoot>
    <tr class="gt_footnotes">
      <td class="gt_footnote" colspan="6"> Abbreviations: CI = Confidence interval, HR = Hazard ratio, RCT = Randomized controlled trial, RWE = Real-world evidence, SMD = standardized mean difference (based on log hazard ratios)</td>
    </tr>
  </tfoot>
</table>
</div>
```

:::
:::

### Detailed SMD Calculation

To understand the SMD calculation in detail:

::: {.cell}

```{.r .cell-code}
# Calculate SMD for one comparison
smd <- smd_agreement(
  rct_estimate = log(0.87),
  rct_lower = log(0.78),
  rct_upper = log(0.97),
  rwe_estimate = log(0.82),
  rwe_lower = log(0.76),
  rwe_upper = log(0.87)
  )

print(paste("SMD value:", round(smd, 2)))
```

::: {.cell-output .cell-output-stdout}

```
[1] "SMD value: 0.9"
```


:::
:::

## Interpretation

The output table from `agreement_metrics()` provides a comprehensive view of agreement:

-   **RCT and RWE Estimates**: Shows point estimates and confidence intervals
-   **Statistical Agreement**: Indicates if results agree in direction and significance
-   **Estimate Agreement**: Shows if RWE estimate falls within RCT confidence interval
-   **SMD Agreement**: Provides standardized difference with threshold-based assessment

Results should be interpreted considering:

1.  Clinical relevance of differences
2.  Quality and characteristics of both RCT and RWE studies
3.  Context-specific tolerance for disagreement

## Technical Details

Important considerations when using these functions:

-   All estimates should be positive (e.g., hazard ratios, odds ratios)
-   Estimates are log-transformed for SMD calculation
-   Confidence intervals are assumed to be 95% intervals
-   The default SMD threshold of 1.96 corresponds to α=0.05

## Interpretation Guidelines

When interpreting agreement metrics, @wang2023 suggest several important considerations:

1.  **Context Matters**: Agreement thresholds may vary depending on the clinical context and the feasibility of emulating specific trial design elements.

2.  **Design Emulation**: Higher agreement is typically observed when RWE studies can closely emulate key trial design elements (population, intervention, comparator, outcome, and timing).

3.  **Multiple Metrics**: Using multiple agreement metrics provides a more complete picture than any single metric alone. For example:

    -   Statistical significance agreement captures directional alignment
    -   Estimate agreement ensures magnitude compatibility
    -   SMD provides a standardized measure of difference accounting for uncertainty

4.  **Limitations**: Discrepancies between RCT and RWE results may arise from:

    -   Residual confounding in RWE studies
    -   Different patient populations
    -   Measurement challenges in real-world data
    -   Random variation

## References

::: {#refs}
:::

## Session info



Script runtime:  minutes.

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

 **encore.analytics**   encore.analytics       0.2.1     

      **tibble**             tibble            3.3.0     
---------------------------------------------------------


:::
:::

### Session info

::: {.cell}

```{.r .cell-code}
pander::pander(sessionInfo())
```

::: {.cell-output-display}
**R version 4.5.1 (2025-06-13)**

**Platform:** aarch64-apple-darwin20 

**locale:**
en_US.UTF-8||en_US.UTF-8||en_US.UTF-8||C||en_US.UTF-8||en_US.UTF-8

**attached base packages:** 
_stats_, _graphics_, _grDevices_, _utils_, _datasets_, _methods_ and _base_

**other attached packages:** 
_encore.analytics(v.0.2.1)_, _tibble(v.3.3.0)_ and _dplyr(v.1.1.4)_

**loaded via a namespace (and not attached):** 
_jsonlite(v.2.0.0)_, _compiler(v.4.5.1)_, _Rcpp(v.1.1.0)_, _tidyselect(v.1.2.1)_, _xml2(v.1.5.1)_, _stringr(v.1.6.0)_, _assertthat(v.0.2.1)_, _tidyr(v.1.3.2)_, _yaml(v.2.3.12)_, _fastmap(v.1.2.0)_, _R6(v.2.6.1)_, _commonmark(v.2.0.0)_, _generics(v.0.1.4)_, _knitr(v.1.51)_, _htmlwidgets(v.1.6.4)_, _MASS(v.7.3-65)_, _pander(v.0.6.6)_, _pillar(v.1.11.1)_, _rlang(v.1.1.6)_, _stringi(v.1.8.7)_, _litedown(v.0.9)_, _xfun(v.0.55)_, _fs(v.1.6.6)_, _sass(v.0.4.10)_, _smd(v.0.8.0)_, _cli(v.3.6.5)_, _withr(v.3.0.2)_, _magrittr(v.2.0.4)_, _tictoc(v.1.2.1)_, _digest(v.0.6.39)_, _rstudioapi(v.0.17.1)_, _markdown(v.2.0)_, _base64enc(v.0.1-3)_, _lifecycle(v.1.0.4)_, _vctrs(v.0.6.5)_, _evaluate(v.1.0.5)_, _glue(v.1.8.0)_, _sessioninfo(v.1.2.3)_, _gt(v.1.2.0)_, _rmarkdown(v.2.30)_, _purrr(v.1.2.0)_, _tools(v.4.5.1)_, _pkgconfig(v.2.0.3)_ and _htmltools(v.0.5.9)_
:::
:::

### Repositories

::: {.cell}

```{.r .cell-code}
pander::pander(options('repos'))
```

::: {.cell-output-display}


  * **repos**:

    --------
      CRAN
    --------
     @CRAN@
    --------


<!-- end of list -->



:::
:::
:::


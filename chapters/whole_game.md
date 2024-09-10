---
title: "The whole game"
subtitle: "Multiple imputation workflows with propensity score and survival analysis"
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





## Background

-   In 2022, nearly every third drug approval was granted in the field of oncology (tendency ↑)[@mullard2022]

-   Decision-makers increasingly rely on real-world evidence (RWE) generated from routine-care health data such as electronic health records (EHR) to evaluate the comparative safety and effectiveness of novel cancer therapies

**ENCORE**

-   The [ENCORE project](https://www.fda.gov/about-fda/oncology-center-excellence/calibrating-real-world-evidence-studies-oncology-against-randomized-trials-encore) is an [RCT DUPLICATE](https://www.rctduplicate.org/) expansion to oncology which is going to emulate 12 randomized clinical trials using multiple EHR data sources. The process includes an emphasis on transparency with documented assessment of data fitness of the RWD source for each trial and conducting extensive sensitivity analyses to assess robustness of findings and trial eligibility criteria.

-   Partially observed covariates/confounders are a common and pervasive challenge

-   To date, most oncology studies utilizing RWD have relied on complete case analysis although assumptions for a complete case analysis (missing completely at random \[MCAR\]) are even stronger than those (missing at random \[MAR\]) for multiple imputation (MI). Besides this, MI has additional advantages:

    -   All patients are retained

    -   Flexible modeling (parametric, non-parametric)

    -   Can incorporate additional information (auxiliary covariates) to make the MAR assumption more likely

    -   Realistic variance estimation (Rubin's rule)

-   However:

    -   Not much is known about how to use multiple imputation in combination with propensity score-based approaches

    -   Computational implementation can be complex

## Objective

::: callout-tip
To establish a computationally reproducible workflow that streamlines multiple imputation \> propensity score matching/weighting \> survival analysis workflows in a transparent fashion
:::

![Streamlined workflow to approach partially observed covariate data in oncology trial emulations.](/images/workflow.png){fig-align="center"}

## Leyrat et al. simulation study

One of the most comprehensive and influental simulation studies that addressed the question on how to combine multiple imputation with propensity scores (IPTW weighting) was published in 2019 by Leyrat et al.[@leyrat2019]. In this study, the authors looked at three different potential ways:

-   **MIte**: MI \> PS estimation \> Outcome model for each PS model \> Pooling of results

-   **MIps**: MI \> PS estimation \> PS pooling across datasets \> single outcome model

-   **MIpar**: MI \> Pooling of covariate parameters \> single PS model \> single outcome model

Additional questions that were also addressed:

-   Should outcome be included in imputation model?

-   How to estimate variance of IPTW estimator in context of MIte or MIps or MIpar?

![Illustration of potential approaches that could be considered after multiple imputation (MI) of the partially observed covariates are missing values on the original dataset.](/images/mi_ps.png){fig-align="center"}

### Simulation study results

-   **MIte** performed best in terms of bias, standardized differences/balancing, coverage rate and variance estimation

    -   MI \> PS estimation \> Outcome model for each PS model \> Pooling of results

-   Standard IPTW variance estimation is valid for MIte

-   Outcome must be included in imputation model

![Leyrat et al. simulation study results.](/images/leyrat_results.png){fig-align="center"}

## References

::: {#refs}
:::

## Session info





Script runtime: 0.00 minutes.

::: panel-tabset
### Loaded packages


::: {.cell}

```{.r .cell-code}
pander::pander(subset(data.frame(sessioninfo::package_info()), attached==TRUE, c(package, loadedversion)))
```

::: {.cell-output-display}
------------- -------------------
 **package**   **loadedversion** 

------------- -------------------
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
_stats_, _graphics_, _grDevices_, _datasets_, _utils_, _methods_ and _base_

**loaded via a namespace (and not attached):** 
_digest(v.0.6.37)_, _fastmap(v.1.2.0)_, _xfun(v.0.47)_, _tictoc(v.1.2.1)_, _knitr(v.1.48)_, _htmltools(v.0.5.8.1)_, _rmarkdown(v.2.28)_, _cli(v.3.6.3)_, _pander(v.0.6.5)_, _sessioninfo(v.1.2.2)_, _renv(v.1.0.7)_, _compiler(v.4.4.0)_, _rstudioapi(v.0.16.0)_, _tools(v.4.4.0)_, _evaluate(v.0.24.0)_, _Rcpp(v.1.0.13)_, _yaml(v.2.3.10)_, _rlang(v.1.1.4)_, _jsonlite(v.1.8.8)_ and _htmlwidgets(v.1.6.4)_
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


#!/bin/sh
set -e

# Bootstrap: pak is installed via install.packages, everything else via pak
R -e "install.packages('pak')"

# Core / data wrangling
R -e "pak::pak('tidyverse')"
R -e "pak::pak('tibble')"
R -e "pak::pak('here')"
R -e "pak::pak('yaml')"

# Survival analysis
R -e "pak::pak('survival')"

# Imputation, matching & weighting
R -e "pak::pak('mice')"
R -e "pak::pak('MatchThem')"
R -e "pak::pak('cobalt')"
R -e "pak::pak('ranger')"
R -e "pak::pak('survey')"
R -e "pak::pak('marginaleffects')"

# Trial design
R -e "pak::pak('gsDesign')"

# Tables & reporting
R -e "pak::pak('gt')"
R -e "pak::pak('gtsummary')"
R -e "pak::pak('cardx')"
R -e "pak::pak('pander')"
R -e "pak::pak('rmarkdown')"

# Parallelization & timing
R -e "pak::pak('parallelly')"
R -e "pak::pak('furrr')"
R -e "pak::pak('tictoc')"

# Reproducibility
R -e "pak::pak('sessioninfo')"
R -e "pak::pak('remotes')"

# GitHub packages
R -e "pak::pak('janickweberpals/encore.analytics@dev')"

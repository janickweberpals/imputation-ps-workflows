#!/usr/bin/env Rscript

# Initialize renv if not already initialized
if (!dir.exists("renv")) {
  renv::init()
}

# Install all required packages
packages <- c(
  'pak',
  'tidyverse',
  'survival',
  'mice',
  'MatchThem',
  'survey',
  'here',
  'cardx',
  'gtsummary',
  'parallelly',
  'ranger',
  'furrr',
  'cobalt',
  'gsDesign',
  'yaml',
  'gt',
  'marginaleffects',
  'tictoc',
  'tibble',
  'rmarkdown',
  'pander',
  'sessioninfo',
  'remotes'
)

# Install packages
for (pkg in packages) {
  renv::install(pkg)
}

# Install GitHub packages
remotes::install_github('janickweberpals/encore.analytics', ref = 'dev')

# Lock dependencies
renv::snapshot()

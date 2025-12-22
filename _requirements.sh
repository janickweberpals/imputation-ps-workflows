#!/bin/sh

R -e "install.packages(c(
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
))"

R -e "remotes::install_github('janickweberpals/encore.analytics', ref = 'dev')"

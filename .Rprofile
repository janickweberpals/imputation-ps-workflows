source("renv/activate.R")

# For reproducibility, freeze package versions

# Mac OS
if(Sys.info()['sysname'][[1]]=="Darwin") options(repos = c(REPO_NAME = "https://packagemanager.posit.co/cran/2024-04-25"))

# Linux RHEL 9
if(Sys.info()['sysname'][[1]]=="Linux")options(repos = c(REPO_NAME = "https://packagemanager.posit.co/cran/__linux__/rhel9/2024-04-25"))

# encore.io functions
source(here::here("functions", "install_on_demand.R"))
source(here::here("functions", "simulate_flaura.R"))
source(here::here("functions", "match_re_weight.R"))

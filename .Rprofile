
# For reproducibility, freeze package versions

# Mac OS
if(Sys.info()['sysname'][[1]]=="Darwin") options(repos = c(REPO_NAME = "https://packagemanager.posit.co/cran/latest"))

# Linux Ubuntu
if(Sys.info()['sysname'][[1]]=="Linux")options(repos = c(REPO_NAME = "https://packagemanager.posit.co/cran/__linux__/noble/latest"))

# The following code creates a local directory for your projects packages and
# removes the users home package directory to stop issues with using packages
# from other projects. It is only executed if neither packrat nor renv are present
# as this approach is not compatible with these packages
if(!(dir.exists("./renv") | dir.exists("./packrat"))){

  if(!dir.exists(".rpkg")){

    dir.create(".rpkg")

    .libPaths(c(".rpkg", .libPaths()[!grepl("/home/", .libPaths())]))

  }

  # set lib paths
  .libPaths(c(".rpkg", .libPaths()[!grepl("/home/", .libPaths())]))
  Sys.setenv(R_LIBS_SITE = .libPaths()[1])

  # notify user upon start up
  cat(paste("Using libraries:", paste(.libPaths(), collapse = ", ")))
}


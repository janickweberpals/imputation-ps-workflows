# R packages --------------------------------------------------------------

# on local machine
.rpkg_dir <- '.rpkg'

# if neither packrat nor renv are present
# as this approach is not compatible with these packages
if (!(dir.exists("./renv") | dir.exists("./packrat"))) {
  if (!dir.exists(.rpkg_dir)) {
    dir.create(.rpkg_dir)

    .libPaths(c(.rpkg_dir, .libPaths()))
  }

  # add .rpkg to lib paths
  .libPaths(c(.rpkg_dir, .libPaths()))
  Sys.setenv(R_LIBS_SITE = .libPaths()[1])
}

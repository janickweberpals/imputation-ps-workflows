#' @importFrom utils askYesNo install.packages
#' 
install_on_demand <- function(pkg, quiet = FALSE, ...) {
  # internal function that checks whether package pkg is
  # in the library. If not found, it asks the user permission
  # to install.
  if (requireNamespace(pkg, quietly = TRUE)) {
    return()
  }
  if (interactive()) {
    answer <- askYesNo(paste("Package", pkg, "needed. Install from CRAN?"))
    if (answer) install.packages(pkg, quiet = quiet)
  }
}

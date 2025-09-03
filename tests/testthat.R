# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)

# Try to load the package from source if not installed
tryCatch({
  library(GeometricMorphometricsMix)
}, error = function(e) {
  # If package is not installed, try loading from source
  source_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
  for (file in source_files) {
    source(file)
  }
})

test_check("GeometricMorphometricsMix")

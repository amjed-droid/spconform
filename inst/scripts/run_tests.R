library(spconform)
cat("Running tests on installed spconform package (version:", as.character(packageVersion("spconform")), ")...\n")
testthat::test_package("spconform")
cat("\nAll installed package tests completed!\n")

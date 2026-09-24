library(spconform)
cat("Running complete testthat directory test suite...\n")
res <- testthat::test_dir("tests/testthat", stop_on_failure = TRUE)
cat("\n=== ALL TESTS PASSED SUCCESSFULLY! ===\n")

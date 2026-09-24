rtools_path <- "D:\\rtools45\\ucrt64\\bin;D:\\rtools45\\usr\\bin;"
Sys.setenv(PATH = paste0(rtools_path, Sys.getenv("PATH")))

cat("Running roxygen2::roxygenise()...\n")
roxygen2::roxygenise()
cat("roxygen2 completed successfully!\n")

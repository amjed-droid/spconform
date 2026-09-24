library(spconform)

cat("===============================================================\n")
cat("   SPEED BENCHMARK: PURE R ENGINE vs COMPILED C++ ENGINE\n")
cat("===============================================================\n\n")

# Benchmark Function
benchmark_run <- function(n_cal, m_test) {
  cat(sprintf("--- Testing with N_cal = %s points, M_test = %s queries ---\n",
              format(n_cal, big.mark=","), format(m_test, big.mark=",")))
  cat(sprintf("    Total Pairwise Interactions: %s\n", format(as.numeric(n_cal) * m_test, big.mark=",")))
  
  set.seed(42)
  s_train <- cbind(runif(n_cal), runif(n_cal))
  y_train <- sin(3 * s_train[, 1]) + cos(3 * s_train[, 2]) + rnorm(n_cal, sd = 0.2)
  s0 <- cbind(runif(m_test), runif(m_test))
  
  pfun <- function(s_tr, y_tr, s_new) rep(mean(y_tr), nrow(s_new))
  
  # 1. Pure R Engine
  cat("    Running Pure R Engine...")
  t_r <- system.time({
    res_r <- scp_geostatistical(s_train, y_train, s0, pfun, split = 0.5, engine = "r", seed = 123)
  })
  cat(sprintf("  Done in %.3f seconds.\n", t_r[3]))
  
  # 2. Compiled C++ Engine
  cat("    Running Compiled C++ Engine...")
  t_cpp <- system.time({
    res_cpp <- scp_geostatistical(s_train, y_train, s0, pfun, split = 0.5, engine = "cpp", seed = 123)
  })
  cat(sprintf(" Done in %.3f seconds.\n", t_cpp[3]))
  
  speedup <- t_r[3] / t_cpp[3]
  cat(sprintf("    >>> SPEEDUP FACTOR: %.1fx FASTER WITH C++! <<<\n\n", speedup))
  
  return(data.frame(
    Pairs = paste0(round((as.numeric(n_cal) * m_test) / 1e6, 1), "M"),
    Time_R = round(t_r[3], 3),
    Time_CPP = round(t_cpp[3], 3),
    Speedup = paste0(round(speedup, 1), "x")
  ))
}

# Run 3 test tiers
tier1 <- benchmark_run(n_cal = 4000, m_test = 1000)   # 4 Million Pairs
tier2 <- benchmark_run(n_cal = 10000, m_test = 2000)  # 20 Million Pairs
tier3 <- benchmark_run(n_cal = 20000, m_test = 5000)  # 100 Million Pairs

cat("===============================================================\n")
cat("                     BENCHMARK SUMMARY TABLE\n")
cat("===============================================================\n")
summary_df <- rbind(tier1, tier2, tier3)
print(summary_df)
cat("===============================================================\n")

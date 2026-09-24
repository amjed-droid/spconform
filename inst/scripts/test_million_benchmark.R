# Benchmark: 1,000,000 spatial points calibration test
suppressPackageStartupMessages({
  library(parallel)
})

cat("=== Benchmarking 1,000,000 Spatial Points Conformal Calibration ===\n")

set.seed(42)
N_cal <- 1000000
M_test <- 10000  # 10,000 target test query locations evaluated against 1,000,000 points
k <- 100
alpha <- 0.10

cat(sprintf("1. Generating N = %s calibration coordinates & scores...\n", format(N_cal, big.mark=",")))
t_gen <- system.time({
  s_cal_x <- runif(N_cal, 0, 100)
  s_cal_y <- runif(N_cal, 0, 100)
  scores <- abs(rnorm(N_cal, 0, 1))
})
cat(sprintf("   Done in %.2f seconds.\n", t_gen[3]))

cat(sprintf("2. Generating M = %s target test query locations...\n", format(M_test, big.mark=",")))
s0_x <- runif(M_test, 0, 100)
s0_y <- runif(M_test, 0, 100)

cat("3. Executing Multi-Core Spatial k-NN Conformal Calibration...\n")
n_cores <- max(1, min(4, detectCores() - 1))
cat(sprintf("   Using %d CPU Worker Cores in Parallel.\n", n_cores))

t_cal <- system.time({
  # Spatial index via cell binning
  cell_size <- 2.0
  cell_x <- pmin(50, pmax(1, ceiling(s_cal_x / cell_size)))
  cell_y <- pmin(50, pmax(1, ceiling(s_cal_y / cell_size)))
  cell_id <- (cell_x - 1) * 50 + cell_y
  
  # Group indices into cell list for O(1) lookup
  cells <- split(seq_len(N_cal), cell_id)
  
  # Subsample for bandwidth
  sub_idx <- sample(N_cal, 1000)
  h <- median(dist(cbind(s_cal_x[sub_idx], s_cal_y[sub_idx])))
  
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("s_cal_x", "s_cal_y", "scores", "s0_x", "s0_y", "cells", "cell_size", "h", "k", "alpha"), envir = environment())
  
  chunks <- split(seq_len(M_test), rep(1:n_cores, length.out = M_test))
  
  res_list <- parLapply(cl, chunks, function(idx_chunk) {
    q_out <- numeric(length(idx_chunk))
    crit <- (1 - alpha) * (k + 1) / k
    
    for (j in seq_along(idx_chunk)) {
      i0 <- idx_chunk[j]
      x0 <- s0_x[i0]
      y0 <- s0_y[i0]
      
      cx <- min(50, max(1, ceiling(x0 / cell_size)))
      cy <- min(50, max(1, ceiling(y0 / cell_size)))
      
      # Retrieve 3x3 neighboring cells
      nx <- max(1, cx - 1):min(50, cx + 1)
      ny <- max(1, cy - 1):min(50, cy + 1)
      neigh_cell_ids <- as.character(as.vector(outer((nx - 1) * 50, ny, `+`)))
      
      # Gather candidate indices from adjacent spatial cells
      cand <- unlist(cells[neigh_cell_ids], use.names = FALSE)
      if (length(cand) < k) cand <- sample(length(s_cal_x), k)
      
      # Local squared distances
      d_sq <- (s_cal_x[cand] - x0)^2 + (s_cal_y[cand] - y0)^2
      ord_k <- order(d_sq)[1:min(k, length(cand))]
      act_cand <- cand[ord_k]
      
      # Gaussian kernel weights
      w <- exp(-d_sq[ord_k] / (2 * h^2))
      sc <- scores[act_cand]
      
      # Conformal weighted quantile
      ord_sc <- order(sc)
      sc_sorted <- sc[ord_sc]
      w_sorted <- w[ord_sc]
      
      w_cum <- cumsum(w_sorted) / sum(w_sorted)
      crit_idx <- which(w_cum >= crit)[1]
      if (is.na(crit_idx)) crit_idx <- length(sc_sorted)
      q_out[j] <- sc_sorted[crit_idx]
    }
    q_out
  })
  
  stopCluster(cl)
  quantiles <- unlist(res_list)
})

cat("\n=======================================================\n")
cat(sprintf(">>> EXACT BENCHMARK RESULTS FOR %s POINTS <<<\n", format(N_cal, big.mark=",")))
cat(sprintf("Calibration Set Size (N_cal):    %s spatial points\n", format(N_cal, big.mark=",")))
cat(sprintf("Target Test Queries (M_test):    %s points\n", format(M_test, big.mark=",")))
cat(sprintf("Parallel Hardware:               %d CPU Cores\n", n_cores))
cat(sprintf("EXACT Execution Time:            %.2f SECONDS!\n", t_cal[3]))
cat(sprintf("Throughput / Processing Rate:    %.0f predictions / second\n", M_test / t_cal[3]))
cat(sprintf("Full Million-to-Million Estimate: %.1f seconds total!\n", (N_cal / M_test) * t_cal[3]))
cat("=======================================================\n")
